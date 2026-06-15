"""
ligand_parameters.amber
Run the post-Gaussian AMBER parameterization pipeline for a ligand:

  espgen      – extract ESP grid from Gaussian HF log
  resp        – fit RESP charges to ESP
  antechamber – assign GAFF atom types and read RESP charges -> .ac
  remap names – restore original PDB atom naming in the .ac file
  antechamber – convert .ac -> .mol2 (preserving GAFF types + RESP charges)
  parmchk2    – generate missing GAFF parameters -> .frcmod
  tleap       – build AMBER library file -> .lib

All tools must be available in $PATH (AmberTools installed).

Output files per ligand (written to workdir):
  {lig}.mol2   – GAFF atom types + RESP charges
  {lig}.frcmod – missing GAFF parameters
  {lig}.lib    – AMBER library file
"""

import re
import subprocess
import sys
from pathlib import Path


_GAFF2_TO_AMBER = {
    'c':  'C',   'c1': 'CX', 'c2': 'CA', 'c3': 'CT', 'ca': 'CA',
    'cb': 'CB',  'cc': 'CA', 'cd': 'CA', 'ce': 'CA', 'cf': 'CA',
    'cp': 'CA',  'cq': 'CA',
    'n':  'N',   'n1': 'N1', 'n2': 'N2', 'n3': 'N3', 'n4': 'N3',
    'na': 'NA',  'nb': 'NA', 'nh': 'N2',
    'o':  'O',   'oh': 'OH', 'os': 'OS', 'op': 'OS', 'oq': 'OS',
    'sh': 'SH',  'ss': 'S',  's2': 'S',  's4': 'S',  's6': 'S',
    'sx': 'S',   'sy': 'S',
    'p2': 'P',   'p3': 'P',  'p4': 'P',  'p5': 'P',
    'hc': 'HC',  'ha': 'HA', 'hn': 'H',  'ho': 'HO', 'hs': 'HS',
    'hp': 'HP',  'h1': 'H1', 'h2': 'H2', 'h3': 'H3', 'h4': 'H4', 'h5': 'H5',
    'f':  'F',   'cl': 'Cl', 'br': 'Br', 'i':  'I',
}


def _run(cmd, cwd):
    """Run a shell command, raise on failure, print output."""
    print(f"  $ {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout.rstrip())
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed (exit {result.returncode}):\n"
            f"  {' '.join(cmd)}\n"
            f"{result.stderr}"
        )


def _postprocess_frcmod(frcmod_path):
    """
    Post-process a parmchk2 frcmod in place:
      - Remove 'ATTN, need revision' lines (zero placeholder parameters).
      - Warn on parameters with penalty score > 100.
    """
    path  = Path(frcmod_path)
    lines = path.read_text().splitlines(keepends=True)

    clean        = []
    high_penalty = []

    for line in lines:
        if 'ATTN' in line:
            continue
        m = re.search(r'penalty score=\s*([\d.]+)', line)
        if m and float(m.group(1)) > 100:
            high_penalty.append(line.rstrip())
        clean.append(line)

    path.write_text(''.join(clean))

    if high_penalty:
        print(f"\n  Warning: parameters with penalty score > 100 in {path.name}:")
        for ln in high_penalty:
            print(f"    {ln.strip()}")


def remap_ac_atom_names(ac_path, ref_pdb_path):
    """
    Rename atoms in an antechamber .ac file to match the original ligand PDB.

    Builds a positional mapping (ac_name[i] -> pdb_name[i]) from ATOM record
    order and applies it to all ATOM and BOND lines.  Overwrites ac_path in place.

    Both ATOM and HETATM records are read from the PDB (ligand PDB files
    typically use HETATM rather than ATOM).

    Parameters
    ----------
    ac_path      : str | Path -- antechamber .ac file to update
    ref_pdb_path : str | Path -- original ligand PDB whose atom names are the target
    """
    ac_path      = Path(ac_path)
    ref_pdb_path = Path(ref_pdb_path)

    ac_names = [
        ln.split()[2]
        for ln in ac_path.read_text().splitlines()
        if ln.startswith('ATOM')
    ]
    # Read both ATOM and HETATM — ligand PDBs use HETATM for non-standard residues
    pdb_names = [
        ln.split()[2]
        for ln in ref_pdb_path.read_text().splitlines()
        if ln.startswith(('ATOM', 'HETATM'))
    ]

    if len(ac_names) != len(pdb_names):
        raise ValueError(
            f"Atom count mismatch: {ac_path.name} has {len(ac_names)} atoms, "
            f"{ref_pdb_path.name} has {len(pdb_names)} atoms."
        )

    name_map = dict(zip(ac_names, pdb_names))

    out_lines = []
    for line in ac_path.read_text().splitlines(keepends=True):
        nl = '\n' if line.endswith('\n') else ''
        s  = line.rstrip('\n')

        if s.startswith('ATOM'):
            m = re.match(r'^(ATOM\s+\d+)\s+(\S+)\s+(\S+\s+.*)$', s)
            if m:
                new_name = name_map.get(m.group(2), m.group(2))
                s = f"{m.group(1)}  {new_name:<4}{m.group(3)}"

        elif s.startswith('BOND'):
            m = re.match(r'^(BOND\s+\d+\s+\d+\s+\d+\s+\d+\s+)(\S+)(\s+)(\S+)(.*)$', s)
            if m:
                new1 = name_map.get(m.group(2), m.group(2))
                new2 = name_map.get(m.group(4), m.group(4))
                s = m.group(1) + new1 + m.group(3) + new2 + m.group(5)

        out_lines.append(s + nl)

    ac_path.write_text(''.join(out_lines))
    print(f"  Remapped {len(name_map)} atom names in {ac_path.name}")
    return name_map


def _read_ac_types(ac_path):
    """Return {atom_name: atom_type} for every ATOM record in an .ac file."""
    types = {}
    for line in Path(ac_path).read_text().splitlines():
        if line.startswith('ATOM'):
            parts = line.split()
            if len(parts) >= 10:
                types[parts[2]] = parts[-1]
    return types


def _patch_ac_types(ac_path, corrections):
    """Replace atom types in-place for the specified atom names."""
    ac_path = Path(ac_path)
    out = []
    for line in ac_path.read_text().splitlines(keepends=True):
        stripped = line.rstrip()
        if stripped.startswith('ATOM'):
            parts = stripped.split()
            if len(parts) >= 10 and parts[2] in corrections:
                trailing = line[len(stripped):]
                stripped  = re.sub(r'\S+$', corrections[parts[2]], stripped)
                line      = stripped + trailing
        out.append(line)
    ac_path.write_text(''.join(out))


def resolve_du_atom_types(ac_path, hf_log, wd, resname, total_charge,
                          atom_type, user_overrides=None):
    """
    Scan ac_path for DU atom types, suggest corrections via a gaff2 probe run,
    merge with user_overrides, and patch the .ac file in place.
    Hard-stops if any DU atom cannot be resolved.
    """
    user_overrides = user_overrides or {}
    ac_path        = Path(ac_path)

    ac_types = _read_ac_types(ac_path)
    du_names = [name for name, t in ac_types.items() if t == 'DU']
    if not du_names:
        return

    print("\n-- resolving DU atom types -------------------------------------------")
    print(f"  antechamber assigned DU to: {', '.join(du_names)}")

    gaff2_types = {}
    suggestions  = {}

    if atom_type != 'gaff2':
        print("  re-running with -at gaff2 to generate type suggestions ...")
        probe_name = '_du_probe.ac'
        probe_path = Path(wd) / probe_name
        try:
            _run([
                'antechamber',
                '-fi', 'gout', '-i', hf_log,
                '-bk', resname,
                '-fo', 'ac',   '-o', probe_name,
                '-c',  'rc',   '-cf', 'resp.chg',
                '-at', 'gaff2',
                '-nc', str(total_charge),
            ], cwd=wd)
            gaff2_types = _read_ac_types(probe_path)
        except RuntimeError:
            print("  Warning: gaff2 probe failed — auto-suggestions unavailable.")
        finally:
            probe_path.unlink(missing_ok=True)

        for name in du_names:
            g2 = gaff2_types.get(name, 'DU')
            if g2 == 'DU':
                continue
            if atom_type == 'amber':
                amber_t = _GAFF2_TO_AMBER.get(g2.lower())
                if amber_t:
                    suggestions[name] = amber_t
            else:
                suggestions[name] = g2
    else:
        print("  atom_type is gaff2 — probe would be identical; skipping auto-suggest.")

    w = max((len(n) for n in du_names), default=4) + 2
    print()
    if gaff2_types:
        print(f"  {'Atom':<{w}} {'GAFF2':<8} {'Applied ({})'.format(atom_type):<22} Source")
        print(f"  {'-'*w} {'-'*8} {'-'*22} {'-'*28}")
        for name in du_names:
            g2  = gaff2_types.get(name, '--')
            if name in user_overrides:
                apl = user_overrides[name]
                src = f"config override  (auto: {suggestions.get(name, '--')})"
            elif name in suggestions:
                apl = suggestions[name]
                src = "auto-suggestion"
            else:
                apl = "--"
                src = "UNRESOLVED"
            print(f"  {name:<{w}} {g2:<8} {apl:<22} {src}")
    else:
        print(f"  {'Atom':<{w}} {'Applied ({})'.format(atom_type):<22} Source")
        print(f"  {'-'*w} {'-'*22} {'-'*28}")
        for name in du_names:
            if name in user_overrides:
                print(f"  {name:<{w}} {user_overrides[name]:<22} config override")
            else:
                print(f"  {name:<{w}} {'--':<22} UNRESOLVED")

    corrections = {}
    unresolved  = []
    for name in du_names:
        if name in user_overrides:
            corrections[name] = user_overrides[name]
        elif name in suggestions:
            corrections[name] = suggestions[name]
        else:
            unresolved.append(name)

    if unresolved:
        print()
        print(f"  Error: unresolved DU type(s): {', '.join(unresolved)}")
        print(f"  Add atom_type_overrides to your config:")
        print(f"    amber:")
        print(f"      atom_type_overrides:")
        for name in unresolved:
            g2_hint = gaff2_types.get(name, '?')
            print(f"        {name}: <type>  # gaff2 assigned: {g2_hint}")
        sys.exit(1)

    auto_applied = {n: t for n, t in corrections.items() if n not in user_overrides}
    if auto_applied:
        print()
        print("  Auto-suggestions applied. To override, add to config:")
        print("    amber:")
        print("      atom_type_overrides:")
        for name, atype in auto_applied.items():
            print(f"        {name}: {atype}")

    _patch_ac_types(ac_path, corrections)
    applied = ',  '.join(f"{n}: DU → {t}" for n, t in corrections.items())
    print(f"\n  Patched {ac_path.name}: {applied}")


def _write_tleap(workdir, resname, mol2_name, frcmod_name, lib_name, forcefield):
    """Write a tleap input file for loading mol2 + frcmod and saving a lib."""
    (Path(workdir) / 'tleap_resp.in').write_text(
        f"source {forcefield}\n"
        f"{resname} = loadmol2 {mol2_name}\n"
        f"check {resname}\n"
        f"loadamberparams {frcmod_name}\n"
        f"saveoff {resname} {lib_name}\n"
        "quit\n"
    )


def run_amber_pipeline(hf_log, resname, charge, ligand_pdb,
                       workdir='.', atom_type='gaff2',
                       forcefield='leaprc.gaff2',
                       atom_type_overrides=None):
    """
    Run the full post-Gaussian AMBER parameterization pipeline for a ligand.

    Parameters
    ----------
    hf_log       : str | Path -- Gaussian HF/ESP log ({lig}_hf.log)
    resname      : str        -- residue name read from the ligand PDB
    charge       : int        -- net molecular charge
    ligand_pdb   : str | Path -- original ligand PDB, used to restore atom names
    workdir      : str | Path -- directory where intermediate and output files live
                                 (Phase 1 files must already be present here)
    atom_type    : str        -- antechamber -at flag: 'gaff2' (default), 'gaff', 'amber'
    forcefield   : str        -- tleap source line: 'leaprc.gaff2', 'leaprc.gaff', etc.
    atom_type_overrides : dict | None -- {atom_name: type} to force-fix DU types
    """
    lig      = Path(ligand_pdb).stem
    hf_log   = str(Path(hf_log).resolve())
    wd       = str(Path(workdir).resolve())

    ac_file   = f"{lig}.ac"
    mol2_out  = f"{lig}.mol2"
    frcmod    = f"{lig}.frcmod"
    lib       = f"{lig}.lib"

    print("\n-- espgen ------------------------------------------------------------")
    _run(['espgen', '-i', hf_log, '-o', 'esp.dat'], cwd=wd)

    print("\n-- resp --------------------------------------------------------------")
    _run([
        'resp', '-O',
        '-i', 'resp.in',
        '-o', 'resp.out',
        '-p', 'resp.pch',
        '-t', 'resp.chg',
        '-q', 'resp.qin',
        '-e', 'esp.dat',
    ], cwd=wd)

    print("\n-- antechamber (gout → ac) -------------------------------------------")
    _run([
        'antechamber',
        '-fi', 'gout',
        '-i',  hf_log,
        '-bk', resname,
        '-fo', 'ac',
        '-o',  ac_file,
        '-c',  'rc',
        '-cf', 'resp.chg',
        '-at', atom_type,
        '-nc', str(charge),
    ], cwd=wd)

    print("\n-- remap atom names --------------------------------------------------")
    remap_ac_atom_names(Path(wd) / ac_file, ligand_pdb)

    resolve_du_atom_types(
        Path(wd) / ac_file, hf_log, wd,
        resname, charge, atom_type, atom_type_overrides,
    )

    print("\n-- antechamber (ac → mol2) -------------------------------------------")
    _run([
        'antechamber',
        '-fi', 'ac',
        '-i',  ac_file,
        '-fo', 'mol2',
        '-o',  mol2_out,
        '-rn', resname,
    ], cwd=wd)

    print("\n-- parmchk2 ----------------------------------------------------------")
    _run([
        'parmchk2',
        '-i', mol2_out,
        '-f', 'mol2',
        '-o', frcmod,
    ], cwd=wd)
    _postprocess_frcmod(Path(wd) / frcmod)

    print("\n-- tleap -------------------------------------------------------------")
    _write_tleap(wd, resname, mol2_out, frcmod, lib, forcefield)
    _run(['tleap', '-f', 'tleap_resp.in'], cwd=wd)

    output_files = [mol2_out, frcmod, lib]
    print(f"\nDone. Output files in {wd}:")
    for f in output_files:
        p = Path(wd) / f
        status = "ok" if p.exists() else "MISSING"
        print(f"  [{status}]  {f}")
