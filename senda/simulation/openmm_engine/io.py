"""
senda.simulation.openmm_engine.io
Coordinate/state I/O and reporters that keep the AMBER file layout.

  write_rst7 / load_start   -- AMBER ASCII restart (coords, velocities, box)
  save_checkpoint           -- portable OpenMM state XML + step sidecar
  AmberNetCDFWriter         -- streaming AMBER NetCDF trajectory (CDF-2), no
                               netCDF4/scipy dependency; readable by cpptraj,
                               pytraj and OpenMM
  NetCDFReporter            -- OpenMM reporter wrapping AmberNetCDFWriter
                               (wraps molecules into the box, like iwrap=1)
  write_done / is_done      -- JSON completion markers (AMBER's "Total wall time")
"""

from __future__ import annotations

import json
import math
import struct
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import app, unit

VELSCALE = 20.455  # AMBER velocity unit: A per (1/20.455) ps


# ---------------------------------------------------------------------------
# Box helpers
# ---------------------------------------------------------------------------

def box_lengths_angles(box_vectors) -> tuple[list[float], list[float]]:
    """(a, b, c) in A and (alpha, beta, gamma) in degrees from OpenMM box vectors."""
    v = np.array([vec.value_in_unit(unit.angstrom) for vec in box_vectors])
    lengths = [float(np.linalg.norm(x)) for x in v]

    def ang(x, y):
        c = float(np.dot(x, y) / (np.linalg.norm(x) * np.linalg.norm(y)))
        return math.degrees(math.acos(max(-1.0, min(1.0, c))))

    return lengths, [ang(v[1], v[2]), ang(v[0], v[2]), ang(v[0], v[1])]


# ---------------------------------------------------------------------------
# AMBER ASCII restart
# ---------------------------------------------------------------------------

def _six_per_line(values: np.ndarray) -> str:
    flat = values.reshape(-1)
    lines = []
    for i in range(0, len(flat), 6):
        lines.append("".join(f"{x:12.7f}" for x in flat[i:i + 6]))
    return "\n".join(lines) + "\n"


def write_rst7(path: Path, state: mm.State, title: str = "senda OpenMM restart",
               velocities: bool = True) -> None:
    """Write an AMBER ASCII .rst7 (coordinates [A], velocities [AMBER units], box)."""
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    natom = len(pos)
    time_ps = state.getTime().value_in_unit(unit.picosecond)
    head = f"{natom:6d}" if natom > 99999 else f"{natom:5d}"
    out = [f"{title}\n", f"{head}{time_ps:15.7e}\n", _six_per_line(pos)]
    if velocities:
        vel = state.getVelocities(asNumpy=True).value_in_unit(unit.angstrom / unit.picosecond)
        out.append(_six_per_line(vel / VELSCALE))
    lengths, angles = box_lengths_angles(state.getPeriodicBoxVectors())
    out.append("".join(f"{x:12.7f}" for x in lengths + angles) + "\n")
    Path(path).write_text("".join(out))


def load_start(context: mm.Context, path: Path) -> None:
    """Load positions/box (+ velocities, time) from an OpenMM state .xml or AMBER .rst7.

    The box always comes from the file, never from the prmtop, so NPT box
    changes propagate between stages exactly as with AMBER restarts.
    """
    path = Path(path)
    if path.suffix == ".xml":
        state = mm.XmlSerializer.deserialize(path.read_text())
        context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
        context.setPositions(state.getPositions())
        context.setVelocities(state.getVelocities())
        context.setTime(state.getTime())
        return
    crd = app.AmberInpcrdFile(str(path))
    if crd.boxVectors is None:
        raise ValueError(f"{path} has no box -- the OpenMM engine needs a periodic system")
    context.setPeriodicBoxVectors(*crd.boxVectors)
    context.setPositions(crd.positions)
    if crd.velocities is not None:
        context.setVelocities(crd.velocities)


def full_state(context: mm.Context, enforce_box: bool = True) -> mm.State:
    return context.getState(getPositions=True, getVelocities=True, getEnergy=True,
                            enforcePeriodicBox=enforce_box)


def save_checkpoint(context: mm.Context, path: Path, step: int) -> None:
    """Portable state XML (+ ``.step`` sidecar). Written atomically."""
    path = Path(path)
    state = context.getState(getPositions=True, getVelocities=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(mm.XmlSerializer.serialize(state))
    tmp.replace(path)
    path.with_name(path.name + ".step").write_text(f"{step}\n")


def checkpoint_step(path: Path) -> int:
    return int(Path(str(path) + ".step").read_text().strip())


# ---------------------------------------------------------------------------
# Completion markers
# ---------------------------------------------------------------------------

def done_path(stem: Path) -> Path:
    return Path(f"{stem}.done")


def is_done(stem: Path) -> bool:
    return done_path(stem).exists()


def write_done(stem: Path, **info) -> None:
    done_path(stem).write_text(json.dumps(info, indent=2) + "\n")


# ---------------------------------------------------------------------------
# AMBER NetCDF trajectory writer (NetCDF classic, 64-bit offset)
# ---------------------------------------------------------------------------

_NC_DIMENSION, _NC_VARIABLE, _NC_ATTRIBUTE = 0x0A, 0x0B, 0x0C
_NC_CHAR, _NC_FLOAT, _NC_DOUBLE = 2, 5, 6
_TYPE_SIZE = {_NC_CHAR: 1, _NC_FLOAT: 4, _NC_DOUBLE: 8}


def _pad4(b: bytes) -> bytes:
    return b + b"\0" * (-len(b) % 4)


def _name(s: str) -> bytes:
    raw = s.encode()
    return struct.pack(">i", len(raw)) + _pad4(raw)


def _attrs(attrs: dict[str, str]) -> bytes:
    if not attrs:
        return struct.pack(">ii", 0, 0)
    out = struct.pack(">ii", _NC_ATTRIBUTE, len(attrs))
    for k, v in attrs.items():
        raw = v.encode()
        out += _name(k) + struct.pack(">ii", _NC_CHAR, len(raw)) + _pad4(raw)
    return out


class AmberNetCDFWriter:
    """Append-only AMBER-convention NetCDF trajectory, flushed after every frame.

    ``append=True`` reopens an existing file written by this class (resume
    after a crash/walltime) and truncates it to *keep_frames* frames so no
    frame is duplicated past the checkpoint that is being resumed.
    """

    _DIMS = [("frame", 0), ("spatial", 3), ("atom", None), ("cell_spatial", 3),
             ("label", 5), ("cell_angular", 3)]

    def __init__(self, path: Path, natoms: int, append: bool = False,
                 keep_frames: int | None = None, title: str = "senda OpenMM"):
        self.path, self.natoms = Path(path), natoms
        dims = [(n, natoms if size is None else size) for n, size in self._DIMS]
        dim_id = {n: i for i, (n, _) in enumerate(dims)}
        # name, dims, type, attributes, static data (None = record variable)
        variables = [
            ("spatial",      ["spatial"],                 _NC_CHAR,   {}, b"xyz"),
            ("cell_spatial", ["cell_spatial"],            _NC_CHAR,   {}, b"abc"),
            ("cell_angular", ["cell_angular", "label"],   _NC_CHAR,   {}, b"alphabeta gamma"),
            ("time",         ["frame"],                   _NC_FLOAT,  {"units": "picosecond"}, None),
            ("coordinates",  ["frame", "atom", "spatial"], _NC_FLOAT, {"units": "angstrom"}, None),
            ("cell_lengths", ["frame", "cell_spatial"],   _NC_DOUBLE, {"units": "angstrom"}, None),
            ("cell_angles",  ["frame", "cell_angular"],   _NC_DOUBLE, {"units": "degree"}, None),
        ]
        dim_size = dict(dims)

        def vsize(v):
            n = _TYPE_SIZE[v[2]]
            for d in v[1]:
                if dim_size[d]:
                    n *= dim_size[d]
            return n + (-n % 4)

        def header(begins):
            out = b"CDF\x02" + struct.pack(">i", 0)
            out += struct.pack(">ii", _NC_DIMENSION, len(dims))
            for n, size in dims:
                out += _name(n) + struct.pack(">i", size)
            out += _attrs({"title": title, "application": "AMBER", "program": "senda",
                           "programVersion": "openmm_engine", "Conventions": "AMBER",
                           "ConventionVersion": "1.0"})
            out += struct.pack(">ii", _NC_VARIABLE, len(variables))
            for v, b in zip(variables, begins):
                out += _name(v[0]) + struct.pack(">i", len(v[1]))
                out += b"".join(struct.pack(">i", dim_id[d]) for d in v[1])
                out += _attrs(v[3]) + struct.pack(">iiq", v[2], vsize(v), b)
            return out

        hdr_len = len(header([0] * len(variables)))
        begins, off = [], hdr_len
        for v in variables:            # static variables first
            if v[4] is not None:
                begins.append(off)
                off += vsize(v)
            else:
                begins.append(None)
        self._rec_begin = off
        for i, v in enumerate(variables):
            if v[4] is None:
                begins[i] = off
                off += vsize(v)
        self._rec_size = off - self._rec_begin
        hdr = header(begins)

        if append and self.path.exists():
            self._fh = self.path.open("r+b")
            if self._fh.read(len(hdr))[8:] != hdr[8:]:
                raise ValueError(f"{self.path}: not a compatible senda NetCDF trajectory")
            self._fh.seek(4)
            n = struct.unpack(">i", self._fh.read(4))[0]
            self.nframes = n if keep_frames is None else min(n, keep_frames)
            self._fh.truncate(self._rec_begin + self.nframes * self._rec_size)
            self._write_numrecs()
        else:
            self._fh = self.path.open("wb")
            self._fh.write(hdr)
            for v in variables:
                if v[4] is not None:
                    self._fh.write(_pad4(v[4]))
            self.nframes = 0
        self._fh.flush()

    def _write_numrecs(self) -> None:
        self._fh.seek(4)
        self._fh.write(struct.pack(">i", self.nframes))

    def write_frame(self, coords_A: np.ndarray, time_ps: float,
                    lengths: list[float], angles: list[float]) -> None:
        self._fh.seek(self._rec_begin + self.nframes * self._rec_size)
        self._fh.write(struct.pack(">f", time_ps))
        self._fh.write(np.asarray(coords_A, dtype=">f4").tobytes())
        self._fh.write(struct.pack(">3d", *lengths))
        self._fh.write(struct.pack(">3d", *angles))
        self.nframes += 1
        self._write_numrecs()
        self._fh.flush()

    def close(self) -> None:
        if not self._fh.closed:
            self._fh.close()


class NetCDFReporter:
    """OpenMM reporter writing wrapped coordinates every *interval* steps."""

    def __init__(self, path: Path, interval: int, natoms: int,
                 append: bool = False, keep_frames: int | None = None):
        self._interval = interval
        self._writer = AmberNetCDFWriter(path, natoms, append=append, keep_frames=keep_frames)

    def describeNextReport(self, simulation):
        steps = self._interval - simulation.currentStep % self._interval
        # steps, positions, velocities, forces, energies, enforcePeriodicBox
        return (steps, True, False, False, False, True)

    def report(self, simulation, state):
        lengths, angles = box_lengths_angles(state.getPeriodicBoxVectors())
        self._writer.write_frame(
            state.getPositions(asNumpy=True).value_in_unit(unit.angstrom),
            state.getTime().value_in_unit(unit.picosecond), lengths, angles,
        )

    def close(self):
        self._writer.close()

    def __del__(self):
        self._writer.close()


class CheckpointReporter:
    """Writes a portable state XML every *interval* steps (AMBER ntwr)."""

    def __init__(self, path: Path, interval: int):
        self._path, self._interval = Path(path), interval

    def describeNextReport(self, simulation):
        steps = self._interval - simulation.currentStep % self._interval
        return (steps, False, False, False, False, False)

    def report(self, simulation, state):
        save_checkpoint(simulation.context, self._path, simulation.currentStep)


def trim_log(path: Path, max_step: int) -> None:
    """Drop StateDataReporter lines past *max_step* before resuming a run."""
    path = Path(path)
    if not path.exists():
        return
    kept = []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            kept.append(line)
            continue
        try:
            if int(line.split(",")[0]) <= max_step:
                kept.append(line)
        except ValueError:
            pass
    path.write_text("\n".join(kept) + ("\n" if kept else ""))
