import glob
from pathlib import Path

import numpy as np

# Residue names added by tleap that do not appear in the input dimer PDB
_SOLVENT_RES = frozenset({
    "WAT", "HOH", "Na+", "Cl-", "K+", "Mg2+", "Ca2+",
    "Na", "Cl", "MG", "CA", "ZN", "FE", "SOD", "CLA", "IP",
})


# ---------------------------------------------------------------------------
# Chain map: build (chain_id, pdb_resnum) -> 1-based AMBER residue number
# ---------------------------------------------------------------------------

def _parse_pdb_residues(pdb_path: Path) -> list:
    """Return ordered (chain, resnum, resname) for each unique non-solvent residue in the PDB."""
    seen = set()
    residues = []
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            chain   = line[21]
            resnum  = int(line[22:26])
            resname = line[17:20].strip()
            if resname in _SOLVENT_RES:
                continue
            key = (chain, resnum)
            if key not in seen:
                seen.add(key)
                residues.append((chain, resnum, resname))
    return residues


def _build_chain_map(pdb_path: Path, top) -> dict:
    """
    Map (chain, pdb_resnum) -> 1-based AMBER residue number by positional
    matching of the dimer PDB against the parm7 topology.

    PDB side: water and ions are excluded via _SOLVENT_RES (residue name filter,
    not record type -- both HOH and ligands are HETATM records).

    Topology side: instead of a predefined solvent list we keep only residues
    whose name appears in the filtered PDB set.  This excludes whatever name
    tleap assigns to water/ions (WAT, HOH, TP3, Na+, ...) without needing to
    enumerate them.

    Note: when a HETATM residue shares a (chain, resnum) with an ATOM residue
    (e.g. a ligand given the same resnum as a protein residue), only the first
    record encountered is kept.  Such ligands will be absent from chain_map;
    their AMBER resids must be declared explicitly via substrate_sequence.
    """
    pdb_res      = _parse_pdb_residues(pdb_path)
    pdb_resnames = {rname for _, _, rname in pdb_res}  # protein + ligand names only

    amber_res = [(r.index, r.name) for r in top.residues if r.name in pdb_resnames]

    if len(pdb_res) != len(amber_res):
        raise ValueError(
            f"Residue count mismatch: PDB has {len(pdb_res)} non-solvent residues, "
            f"topology has {len(amber_res)} matching residues. "
            f"Ensure {pdb_path.name} is the exact file fed to tleap."
        )

    # Store 1-based AMBER residue numbers to match PDB and AMBER convention
    return {
        (chain, resnum): amber_idx + 1
        for (chain, resnum, _), (amber_idx, _) in zip(pdb_res, amber_res)
    }


# ---------------------------------------------------------------------------
# Substrate chain map: geometric proximity
# ---------------------------------------------------------------------------

def _build_substrate_chain_map(
    substrate_A_set:     set,
    all_substrate_amber: set,
    sequence_specs:      list,
    chain_map:           dict,
    protein_chains:      list,
    top_atoms:           list,
    top_path:            Path,
    sim_base:            Path,
) -> dict:
    """
    For every chain, assign the correct substrate copy by geometric proximity.

    Proximity is used for ALL chains so that PDB chain ordering in the dimer
    does not affect which substrate copy is used per chain.  Copies are assigned
    greedily in chain order; each copy is used exactly once.

    substrate_A_set:     user-declared 1-based AMBER resids (keys in the returned map)
    all_substrate_amber: all substrate copies found in the topology (by residue name)
    sequence_specs:      [{"pdb_resnum": N, "name": aname}, ...] -- protein reference atoms
    """
    import pytraj as pt

    nc_files = sorted(glob.glob(str(sim_base / "replica_1" / "04_NVT" / "structure_NVT_*.nc")))
    if not nc_files:
        raise FileNotFoundError(
            f"No NVT trajectories found in {sim_base / 'replica_1' / '04_NVT'}; "
            "cannot perform geometric substrate assignment"
        )
    coords = pt.load(nc_files[0], top=str(top_path), frame_indices=[0]).xyz[0]

    # Pre-compute atom coordinates for each substrate copy
    cand_coords: dict = {
        ar: coords[[a.index for a in top_atoms if a.resid + 1 == ar]]
        for ar in all_substrate_amber
    }

    result:    dict = {}
    remaining: set  = set(all_substrate_amber)

    for ch in protein_chains:
        # Reference: protein atoms in this chain resolved via PDB resnum
        ref_list = []
        for spec in sequence_specs:
            amber_ch = chain_map.get((ch, spec["pdb_resnum"]))
            if amber_ch is None:
                continue
            for atom in top_atoms:
                if atom.resid + 1 == amber_ch and atom.name == spec["name"]:
                    ref_list.append(coords[atom.index])
                    break

        if not ref_list:
            raise ValueError(
                f"No reference sequence atoms could be resolved for chain {ch!r}. "
                f"Specs tried: {sequence_specs}. "
                f"Check that the PDB residue numbers match entries in the dimer PDB."
            )
        ref_arr = np.array(ref_list)

        # Pick the closest unassigned substrate copy
        best_amber, best_dist = None, float("inf")
        for ar in remaining:
            cand_c = cand_coords[ar]
            if len(cand_c) == 0:
                continue
            diffs = cand_c[:, np.newaxis, :] - ref_arr[np.newaxis, :, :]
            min_d = float(np.sqrt((diffs ** 2).sum(axis=2)).min())
            if min_d < best_dist:
                best_dist  = min_d
                best_amber = ar

        if best_amber is None:
            raise ValueError(f"No substrate candidate found for chain {ch!r}")

        remaining.discard(best_amber)
        # Map every user-declared substrate_sequence key to this chain's copy
        result[ch] = {amber_A: best_amber for amber_A in substrate_A_set}

    return result


# ---------------------------------------------------------------------------
# Atom index resolution
# ---------------------------------------------------------------------------

def _resolve_atom_index(
    chain_map:           dict,
    top_atoms:           list,
    spec:                dict,
    chain:               str,
    substrate_chain_map: dict,
) -> int:
    """
    Return 1-based AMBER atom index for a spec in the given chain.

    spec fields:
      sequence           -- PDB residue number (same in all monomers; protein residues)
      substrate_sequence -- 1-based AMBER residue number declared in chains[0]; the actual
                            resid used per chain is looked up via substrate_chain_map,
                            which is assigned by geometric proximity for every chain.
    """
    aname = spec["name"]

    if "substrate_sequence" in spec:
        amber_A = spec["substrate_sequence"]
        if chain not in substrate_chain_map:
            raise ValueError(
                f"No substrate mapping for chain {chain!r}. "
                "Ensure substrate_sequence specs are present and NVT frames exist."
            )
        amber_resid = substrate_chain_map[chain][amber_A]
    else:
        pdb_resnum  = spec["sequence"]
        amber_resid = chain_map.get((chain, pdb_resnum))
        if amber_resid is None:
            raise ValueError(
                f"PDB residue {pdb_resnum} not found for chain {chain!r} in the dimer PDB. "
                "Check that the sequence number matches the PDB file."
            )

    for atom in top_atoms:
        if atom.resid + 1 == amber_resid and atom.name == aname:
            return atom.index + 1  # 1-based for AMBER masks

    raise ValueError(
        f"Atom {aname!r} not found at AMBER residue {amber_resid} (chain {chain!r})"
    )


# ---------------------------------------------------------------------------
# Topology and PDB path helpers
# ---------------------------------------------------------------------------

def _find_topology(rep_dir: Path) -> Path:
    for name in ("structure_HMR.parm7", "structure.parm7"):
        p = rep_dir / "00_prep" / name
        if p.exists():
            return p
    raise FileNotFoundError(f"No topology found in {rep_dir / '00_prep'}")


def _find_pdb(rep1_dir: Path, protein_dir: Path, inh: str, mut: str) -> Path:
    # Prefer the copy in 00_prep -- that is exactly what tleap received
    p = rep1_dir / "00_prep" / f"{mut}_{inh}_dimer.pdb"
    if p.exists():
        return p
    p = protein_dir / f"{mut}_{inh}_dimer.pdb"
    if p.exists():
        return p
    raise FileNotFoundError(
        f"Dimer PDB not found in {rep1_dir / '00_prep'} or {protein_dir}"
    )


# ---------------------------------------------------------------------------
# Distribution analysis
# ---------------------------------------------------------------------------

def _analyse_column(values: np.ndarray):
    counts, edges = np.histogram(values, bins=100)
    bin_centers   = 0.5 * (edges[:-1] + edges[1:])
    mean = float(np.mean(values))
    std  = float(np.std(values))

    try:
        from scipy.ndimage import uniform_filter1d
        from scipy.signal import find_peaks
        smoothed        = uniform_filter1d(counts.astype(float), size=5)
        peak_idx, props = find_peaks(smoothed, height=smoothed.max() * 0.1, distance=5)
        multimodal      = False
        if len(peak_idx) > 1:
            heights  = props["peak_heights"]
            sorted_h = np.sort(heights)[::-1]
            if sorted_h[1] > sorted_h[0] * 0.30:
                multimodal = True
        mode = float(bin_centers[np.argmax(smoothed)])
    except ImportError:
        multimodal = False
        mode       = float(bin_centers[np.argmax(counts)])

    return mode, mean, std, counts, bin_centers, multimodal


# ---------------------------------------------------------------------------
# Frame selection
# ---------------------------------------------------------------------------

def _select_frame(pooled: np.ndarray, modes: list, stds: list) -> int:
    safe_stds  = np.array([max(s, 1e-10) for s in stds])
    deviations = np.abs(pooled - np.array(modes)) / safe_stds
    scores     = deviations.sum(axis=1)

    for sigma in (1, 2, 3):
        mask = np.all(deviations <= sigma, axis=1)
        if mask.any():
            if sigma > 1:
                print(f"  WARNING: no frame within {sigma - 1}-sigma; relaxing to {sigma}-sigma")
            candidates = np.where(mask)[0]
            return int(candidates[np.argmin(scores[candidates])])

    return int(np.argmin(scores))


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _plot(out_path, dist_specs, pooled, modes, means, stds, multimodal_flags, selected_row):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n     = pooled.shape[1]
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4))
    if n == 1:
        axes = [axes]

    for i, ax in enumerate(axes):
        values  = pooled[:, i]
        label   = dist_specs[i].get("label", f"dist_{i}")
        mode, mean, std = modes[i], means[i], stds[i]
        mm      = multimodal_flags[i]
        sel_val = pooled[selected_row, i]

        ax.hist(values, bins=100, density=True, color="steelblue", alpha=0.5, label="density")
        ax.axvline(mode,    color="red",    linestyle="--",     label=f"mode {mode:.2f}")
        ax.axvline(mean,    color="orange", linestyle="dotted", label=f"mean {mean:.2f}")
        ax.axvspan(mode - std, mode + std,  color="red", alpha=0.1, label="mode +/- 1sigma")
        ax.axvline(sel_val, color="green",  linestyle="-",      label=f"selected {sel_val:.2f}")

        title = label + (" [MULTIMODAL]" if mm else "")
        if mm:
            ax.set_facecolor("lightyellow")
        ax.set_title(title)
        ax.set_xlabel("Distance (A)")
        ax.set_ylabel("Density")
        ax.legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(str(out_path), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Per-chain analysis
# ---------------------------------------------------------------------------

def _analyse_chain(
    chain:               str,
    chain_specs:         list,
    n_replicas:          int,
    sim_base:            Path,
    inh:                 str,
    mut:                 str,
    top_path:            Path,
    chain_map:           dict,
    substrate_chain_map: dict,
    plots_dir:           Path,
    data_dir:            Path,
) -> None:
    import pytraj as pt

    top          = pt.load_topology(str(top_path))
    top_atoms    = list(top.atoms)
    top_residues = list(top.residues)

    labels     = []
    atom_pairs = []
    for spec in chain_specs:
        atoms = spec["atoms"]
        if len(atoms) != 2:
            raise ValueError(f"Each distance must have exactly 2 atoms, got {len(atoms)}")
        idx1 = _resolve_atom_index(chain_map, top_atoms, atoms[0], chain, substrate_chain_map)
        idx2 = _resolve_atom_index(chain_map, top_atoms, atoms[1], chain, substrate_chain_map)
        atom_pairs.append((idx1, idx2))
        label = spec.get("label", f"{atoms[0]['name']}-{atoms[1]['name']}")
        labels.append(label)

        a1 = top_atoms[idx1 - 1]
        a2 = top_atoms[idx2 - 1]
        r1 = top_residues[a1.resid]
        r2 = top_residues[a2.resid]
        print(f"  [DEBUG] {label}:")
        print(f"          atom1 @{idx1}  {a1.name} res {r1.index + 1} ({r1.name})")
        print(f"          atom2 @{idx2}  {a2.name} res {r2.index + 1} ({r2.name})")

    # Pool distances across all replicas
    all_dist_cols = [[] for _ in atom_pairs]
    frame_lookup  = []

    for rep_n in range(1, n_replicas + 1):
        rep_dir  = sim_base / f"replica_{rep_n}"
        nc_files = sorted(glob.glob(str(rep_dir / "04_NVT" / "structure_NVT_*.nc")))
        if not nc_files:
            print(f"    replica_{rep_n}: no NVT trajectories found, skipping")
            continue

        traj     = pt.load(nc_files, top=str(top_path))
        pt.autoimage(traj)
        n_frames = len(traj)
        print(f"    replica_{rep_n}: {n_frames} frames")

        for col_i, (idx1, idx2) in enumerate(atom_pairs):
            dists = pt.distance(traj, f"@{idx1} @{idx2}")
            all_dist_cols[col_i].extend(dists.tolist())

        for local_i in range(1, n_frames + 1):
            frame_lookup.append((rep_n, local_i))

    if not frame_lookup:
        print(f"  No frames collected for chain {chain}, skipping")
        return

    pooled = np.column_stack([np.array(col) for col in all_dist_cols])

    # Distribution analysis
    modes, means, stds, multimodal_flags = [], [], [], []
    print(f"\n  {'Label':<24} {'Mode':>8} {'Mean':>8} {'Std':>8}  Multimodal")
    print(f"  {'-'*24} {'-'*8} {'-'*8} {'-'*8}  ----------")
    for i, label in enumerate(labels):
        mode, mean, std, _, _, mm = _analyse_column(pooled[:, i])
        modes.append(mode)
        means.append(mean)
        stds.append(std)
        multimodal_flags.append(mm)
        print(f"  {label + ('*' if mm else ''):<24} {mode:>8.3f} {mean:>8.3f} {std:>8.3f}")
        if mm:
            print(f"    WARNING: multimodal distribution for {label!r} -- using highest-probability peak")

    # Frame selection
    best_row           = _select_frame(pooled, modes, stds)
    best_rep, best_loc = frame_lookup[best_row]
    print(f"\n  Selected: replica {best_rep}, frame {best_loc}")
    print(f"\n  {'Label':<24} {'Selected':>10}")
    print(f"  {'-'*24} {'-'*10}")
    for i, label in enumerate(labels):
        print(f"  {label:<24} {pooled[best_row, i]:>10.3f}")

    # Write rst7
    nc_files_sel = sorted(glob.glob(str(sim_base / f"replica_{best_rep}" / "04_NVT" / "structure_NVT_*.nc")))
    traj_sel     = pt.load(nc_files_sel, top=str(top_path))
    rst7_path    = sim_base / f"{inh}_{mut}_chain{chain}_representative.rst7"
    pt.write_traj(str(rst7_path), traj_sel[best_loc - 1:best_loc], format="rst7", overwrite=True)
    print(f"\n  Restart written: {rst7_path}")

    # Save CSV data
    stem     = f"{inh}_{mut}_chain{chain}_distances"
    csv_path = data_dir / f"{stem}.csv"
    header   = "replica,frame," + ",".join(labels)
    rows     = [
        f"{rep},{frm}," + ",".join(f"{pooled[i, j]:.6f}" for j in range(len(labels)))
        for i, (rep, frm) in enumerate(frame_lookup)
    ]
    csv_path.write_text(header + "\n" + "\n".join(rows) + "\n")
    print(f"  Data written:    {csv_path}")

    # Plot
    try:
        import matplotlib  # noqa: F401
        plot_path = plots_dir / f"{stem}.png"
        _plot(plot_path, chain_specs, pooled, modes, means, stds, multimodal_flags, best_row)
        print(f"  Plot written:    {plot_path}")
    except ImportError:
        print("  NOTE: matplotlib not installed, skipping plot")


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def analyse(
    inh:         str,
    mut:         str,
    dist_specs:  list,
    n_replicas:  int,
    chains:      list,
    protein_dir: Path,
    cwd:         Path,
) -> None:
    try:
        import pytraj as pt
    except ImportError:
        raise ImportError("pytraj is required for analysis: pip install pytraj")

    sim_base = cwd / "simulations" / inh / mut
    rep1_dir = sim_base / "replica_1"

    top_path     = _find_topology(rep1_dir)
    pdb_path     = _find_pdb(rep1_dir, protein_dir, inh, mut)
    top          = pt.load_topology(str(top_path))
    top_atoms    = list(top.atoms)
    top_residues = list(top.residues)
    chain_map    = _build_chain_map(pdb_path, top)

    # Collect unique sequence specs (PDB resnum + atom name) and substrate resids
    seen_seq:        set  = set()
    sequence_specs:  list = []
    substrate_A_set: set  = set()

    for ds in dist_specs:
        for spec in ds["atoms"]:
            if "substrate_sequence" in spec:
                substrate_A_set.add(spec["substrate_sequence"])
            elif "sequence" in spec:
                key = (spec["sequence"], spec["name"])
                if key not in seen_seq:
                    seen_seq.add(key)
                    sequence_specs.append({"pdb_resnum": spec["sequence"], "name": spec["name"]})

    substrate_chain_map: dict = {}
    if substrate_A_set:
        # Find ALL copies of the substrate in the topology by residue name.
        # This is chain-order-independent: we don't assume substrate_sequence
        # belongs to chains[0]; proximity assigns the right copy to each chain.
        sub_name            = top_residues[next(iter(substrate_A_set)) - 1].name
        all_substrate_amber = {r.index + 1 for r in top_residues if r.name == sub_name}
        substrate_chain_map = _build_substrate_chain_map(
            substrate_A_set, all_substrate_amber, sequence_specs, chain_map,
            chains, top_atoms, top_path, sim_base,
        )

    print(f"  Topology : {top_path.name}")
    print(f"  PDB      : {pdb_path.name}")
    print(f"  Chains   : {chains}")
    if substrate_chain_map:
        for ch, mapping in substrate_chain_map.items():
            pairs = ", ".join(f"{a}->{b}" for a, b in sorted(mapping.items()))
            print(f"  Substrate chain {ch}: {pairs}")

    plots_dir = cwd / "Analysis" / "plots" / "distances"
    data_dir  = cwd / "Analysis" / "data"  / "distances"
    plots_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    for chain in chains:
        print(f"\n  --- Chain {chain} ---")
        _analyse_chain(
            chain, dist_specs, n_replicas, sim_base, inh, mut,
            top_path, chain_map, substrate_chain_map,
            plots_dir, data_dir,
        )
