"""
senda.integration.string

PMF/free-energy integration of stage 07's (adaptive string method)
sampling output. Port of two originally standalone scripts
(sampling_to_wham.py, mbar_10chunks.py) into parameterized functions.

Reads simulations/{inh}/{mut}/07_QMMM_string/results/*_final.dat and writes,
inside that same results/ directory:
  wham/N.dat, wham/meta            -- WHAM-format per-window data + meta
  mbar/meta                        -- ndfes-format meta (also used by vfep/)
  vfep/meta                        -- same as mbar/meta
  mbar/chunks/chunk_NN/             -- per-chunk data/meta/MBAR output
  mbar/mbar_PMF_<n_chunks>.PMF      -- averaged PMF with SE (5 columns)
  mbar/mbar_PMF_<n_chunks>.png      -- PMF plot

Requires the external ndfes/ndfes-PrintFES.py tools in $PATH, and
matplotlib (pip install -e ".[analysis]") for the PMF plot.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

from ..qmmm.common.atoms import resolve_stage_cfg


def convert_sampling_to_wham(results_dir: Path, temp: float) -> None:
    """
    Convert results_dir/*_final.dat files to WHAM/MBAR/vFEP input files.

    *_final.dat layout:
      line 1 (header): force_constant  rc_center  ...
      line 2+:         rc_value  ...  (many columns; only col 0 is needed)
    """
    wham_dir = results_dir / "wham"
    mbar_dir = results_dir / "mbar"
    vfep_dir = results_dir / "vfep"
    for d in (wham_dir, mbar_dir, vfep_dir):
        d.mkdir(exist_ok=True)

    sampling_files = sorted(
        (p for p in results_dir.glob("*_final.dat") if not p.name.startswith(".")),
        key=lambda p: int(p.stem.split("_")[0]),
    )
    if not sampling_files:
        raise FileNotFoundError(f"No *_final.dat files found in {results_dir}")

    wham_meta_rows  = []
    ndfes_meta_rows = []

    for src in sampling_files:
        idx   = src.stem.split("_")[0]
        lines = src.read_text().splitlines()
        header    = lines[0].split()
        force_k   = float(header[0])
        rc_center = float(header[1])

        dat_name = f"{idx}.dat"
        dat_path = wham_dir / dat_name
        with open(dat_path, "w") as f:
            for step, line in enumerate(lines[1:], start=1):
                rc_val = line.split()[0]
                f.write(f"{step} {rc_val}\n")

        wham_meta_rows.append(f"{dat_name} {rc_center:.5E} {force_k:.5E}")
        ndfes_meta_rows.append(
            f"0 {temp:.2f} ../wham/{dat_name} {rc_center:.5E} {force_k / 2:.4f}"
        )
        print(f"    {src.name} -> {dat_path.name}  rc={rc_center:.4f}  k={force_k:.4f}")

    (wham_dir / "meta").write_text("\n".join(wham_meta_rows) + "\n")
    (mbar_dir / "meta").write_text("\n".join(ndfes_meta_rows) + "\n")
    (vfep_dir / "meta").write_text("\n".join(ndfes_meta_rows) + "\n")

    print(f"    Wrote wham/meta ({len(wham_meta_rows)} windows)")
    print(f"    Wrote mbar/meta and vfep/meta ({len(ndfes_meta_rows)} windows)")


def run_mbar_chunks(results_dir: Path, n_chunks: int, temp: float) -> Path:
    """
    Split each umbrella window's sampling data into n_chunks consecutive
    chunks, run MBAR (via the external ndfes tools) on each chunk, average
    the resulting PMFs, and write a 5-column PMF file (rc, fe, se, weight,
    n_chunks). Returns the path to that file.
    """
    import numpy as np

    wham_dir   = results_dir / "wham"
    mbar_dir   = results_dir / "mbar"
    chunks_dir = mbar_dir / "chunks"

    meta_lines = (mbar_dir / "meta").read_text().splitlines()
    windows = []
    for line in meta_lines:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        win_name = Path(parts[2]).name   # "N.dat"
        windows.append({
            "name": win_name,
            "rc":   parts[3],
            "k":    parts[4],
            "dat":  wham_dir / win_name,
        })

    print(f"    Found {len(windows)} umbrella windows")

    print("    Splitting data into chunks...")
    chunk_lines = {c: {} for c in range(1, n_chunks + 1)}
    for win in windows:
        all_lines = win["dat"].read_text().splitlines()
        n  = len(all_lines)
        cs = n // n_chunks
        for c in range(n_chunks):
            start = c * cs
            end   = (c + 1) * cs if c < n_chunks - 1 else n
            chunk_lines[c + 1][win["name"]] = all_lines[start:end]

    pmfs = []
    for c in range(1, n_chunks + 1):
        cdir = chunks_dir / f"chunk_{c:02d}"
        cdir.mkdir(parents=True, exist_ok=True)

        meta_rows = []
        for win in windows:
            (cdir / win["name"]).write_text(
                "\n".join(chunk_lines[c][win["name"]]) + "\n"
            )
            meta_rows.append(f"0 {temp:.2f} ./{win['name']} {win['rc']} {win['k']}")
        (cdir / "meta").write_text("\n".join(meta_rows) + "\n")

        print(f"    [{c}/{n_chunks}] Running MBAR in {cdir.name}...")
        subprocess.run(
            f"ndfes --temp={temp:.0f} --mbar -w 0.4 --nboot=0 -c mbar.xml meta",
            shell=True, check=True, cwd=cdir,
        )
        with open(cdir / "mbar.PMF", "w") as fout:
            subprocess.run(
                "ndfes-PrintFES.py --ci mbar.xml",
                shell=True, check=True, cwd=cdir, stdout=fout,
            )

        data = np.loadtxt(cdir / "mbar.PMF")
        pmfs.append(data)
        print(f"      -> {len(data)} RC points")

    # Find common RC grid across all chunks (floating-point safe)
    def rc_key(r):
        return round(float(r), 6)

    common_rc = sorted(
        set.intersection(*[{rc_key(row[0]) for row in pmf} for pmf in pmfs])
    )
    print(f"    Common RC points across all chunks: {len(common_rc)}")

    fe_matrix = np.zeros((n_chunks, len(common_rc)))
    for i, pmf in enumerate(pmfs):
        rc_map = {rc_key(row[0]): row[1] for row in pmf}
        fe_matrix[i] = [rc_map[rc] for rc in common_rc]

    mean_fe = fe_matrix.mean(axis=0)

    # Zero energy to the first local minimum on the reactant (low-RC) side.
    # A short moving average suppresses point-to-point MBAR noise (sub-kcal/mol
    # wiggles) so the genuine reactant basin is picked rather than a noise dip.
    def find_reactant_min_idx(fe, smooth_window=5):
        kernel   = np.ones(smooth_window) / smooth_window
        pad      = smooth_window // 2
        smoothed = np.convolve(np.pad(fe, pad, mode="edge"), kernel, mode="valid")
        for i in range(1, len(smoothed) - 1):
            if smoothed[i] < smoothed[i - 1] and smoothed[i] < smoothed[i + 1]:
                lo, hi = max(0, i - pad), min(len(fe), i + pad + 1)
                return lo + int(np.argmin(fe[lo:hi]))
        return int(np.argmin(fe))  # fallback: no interior minimum found

    ref_idx = find_reactant_min_idx(mean_fe)
    ref_val = mean_fe[ref_idx]

    # Re-reference each chunk's PMF to zero at the reactant minimum before
    # computing the mean and SE. Every chunk then agrees exactly (=0) at
    # ref_idx, so the SE there is exactly 0 by construction, and the error
    # propagates backward/forward from that RC point based on each chunk's
    # deviation relative to its own anchor there.
    fe_matrix -= fe_matrix[:, ref_idx:ref_idx + 1]
    mean_fe = fe_matrix.mean(axis=0)
    se_fe   = fe_matrix.std(axis=0, ddof=1) / np.sqrt(n_chunks)

    print(f"    Zero reference: first reactant-side minimum at RC = {common_rc[ref_idx]:.4f}, "
          f"was {ref_val:.4f} kcal/mol; SE at reference = {se_fe[ref_idx]:.3e} kcal/mol")

    out_path = mbar_dir / f"mbar_PMF_{n_chunks}.PMF"
    with open(out_path, "w") as f:
        for rc, fe, err in zip(common_rc, mean_fe, se_fe):
            f.write(f"  {rc:14.8f}  {fe:14.6e}  {err:10.3e}  {1.000:7.3f}  {n_chunks:6d}\n")

    print(f"    Wrote {out_path}")
    print(f"    RC range : {common_rc[0]:.4f} to {common_rc[-1]:.4f}")
    print(f"    FE range : {mean_fe.min():.3f} to {mean_fe.max():.3f} kcal/mol")
    print(f"    SE range : {se_fe.min():.4f} to {se_fe.max():.4f} kcal/mol")

    return out_path


def plot_pmf(pmf_path: Path, out_path: Path, label: str) -> None:
    """Plot the averaged PMF (mean +/- SE band) from an mbar_PMF_*.PMF file."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    data = np.loadtxt(pmf_path)
    rc, fe, se = data[:, 0], data[:, 1], data[:, 2]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.fill_between(rc, fe - se, fe + se, color="steelblue", alpha=0.2, label="SE")
    ax.plot(rc, fe, color="steelblue", lw=1.5, label="mean PMF")
    ax.set_xlabel("Reaction coordinate")
    ax.set_ylabel(r"Free energy (kcal$\bullet$mol$^{-1}$)")
    ax.set_title(f"PMF: {label}")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(str(out_path), dpi=150)
    plt.close(fig)


def setup(inh: str, mut: str, inh_cfg: dict, cfg: dict, cwd: Path) -> None:
    """
    Run PMF integration for one inhibitor/mutant pair's stage 07 sampling
    output: convert to WHAM/MBAR/vFEP input, run chunked MBAR, and plot.
    """
    sim_base    = cwd / "simulations" / inh / mut
    results_dir = sim_base / "07_QMMM_string" / "results"
    if not results_dir.is_dir():
        raise FileNotFoundError(
            f"{results_dir} not found -- run 'senda-qmmm string' and let the "
            "ASM simulation produce sampling output first."
        )

    string_cfg_top = (cfg.get("qmmm") or {}).get("string") or {}
    string_cfg = resolve_stage_cfg(inh_cfg, string_cfg_top, "string")
    equil_cfg  = resolve_stage_cfg(inh_cfg, string_cfg_top, "equil")
    # Same temperature the string simulation itself actually ran at, not a
    # separately-declared value that could drift out of sync with it.
    temp = float(string_cfg.get("temp", equil_cfg.get("temp", 300.0)))

    integration_cfg = resolve_stage_cfg(inh_cfg, string_cfg_top, "integration")
    n_chunks = int(integration_cfg.get("n_chunks", 10))

    print(f"  Converting sampling data to WHAM/MBAR/vFEP input (T={temp:.2f} K)...")
    convert_sampling_to_wham(results_dir, temp)

    print(f"  Running chunked MBAR ({n_chunks} chunks)...")
    pmf_path = run_mbar_chunks(results_dir, n_chunks, temp)

    plot_path = pmf_path.with_suffix(".png")
    plot_pmf(pmf_path, plot_path, label=f"{inh}/{mut}")
    print(f"  Wrote {plot_path}")
