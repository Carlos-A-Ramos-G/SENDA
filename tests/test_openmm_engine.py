"""Unit tests for senda.simulation.openmm_engine (need openmm; skipped otherwise)."""

import math

import numpy as np
import pytest

mm = pytest.importorskip("openmm")
from openmm import app, unit  # noqa: E402

from senda.simulation.openmm_engine.io import (  # noqa: E402
    AmberNetCDFWriter, trim_log, write_rst7,
)
from senda.simulation.openmm_engine.restraints import (  # noqa: E402
    KCAL_TO_KJ, nmr_forces, parse_disang, positional_force, read_mask_atoms,
)


def amber_flat_bottom(r, r1, r2, r3, r4, rk2, rk3):
    """AMBER nmropt flat-bottomed restraint energy (manual, &rst iat/r1-r4/rk2/rk3)."""
    if r < r1:
        return rk2 * (r1 - r2) ** 2 + 2 * rk2 * (r1 - r2) * (r - r1)
    if r < r2:
        return rk2 * (r - r2) ** 2
    if r <= r3:
        return 0.0
    if r <= r4:
        return rk3 * (r - r3) ** 2
    return rk3 * (r4 - r3) ** 2 + 2 * rk3 * (r4 - r3) * (r - r4)


def _energy_kcal(positions_A, forces, box_nm=5.0):
    system = mm.System()
    for _ in positions_A:
        system.addParticle(12.0)
    system.setDefaultPeriodicBoxVectors(mm.Vec3(box_nm, 0, 0), mm.Vec3(0, box_nm, 0),
                                        mm.Vec3(0, 0, box_nm))
    for f in forces:
        system.addForce(f)
    ctx = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(np.asarray(positions_A) * 0.1)
    return ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)


def _dihedral_deg(p):
    b0, b1, b2 = p[0] - p[1], p[2] - p[1], p[3] - p[2]
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return math.degrees(math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


WELL = dict(r1=2.0, r2=2.5, r3=3.0, r4=3.5, rk2=20.0, rk3=30.0)


@pytest.mark.parametrize("d", [1.2, 2.0, 2.3, 2.7, 3.0, 3.2, 3.5, 4.4])
def test_distance_restraint_matches_amber(d):
    pos = [[0, 0, 0], [d, 0, 0]]
    r = {"atoms": [0, 1], **WELL}
    e = _energy_kcal(pos, nmr_forces([r]))
    assert e == pytest.approx(amber_flat_bottom(d, **WELL), abs=1e-4)


def test_distance_restraint_uses_minimum_image():
    # 4.6 nm apart in a 5 nm box -> 0.4 nm = 4 A through the boundary
    pos = [[0, 0, 0], [46.0, 0, 0]]
    e = _energy_kcal(pos, nmr_forces([{"atoms": [0, 1], **WELL}]))
    assert e == pytest.approx(amber_flat_bottom(4.0, **WELL), abs=1e-4)


@pytest.mark.parametrize("theta", [60.0, 95.0, 109.5, 120.0, 150.0])
def test_angle_restraint_matches_amber(theta):
    t = math.radians(theta)
    pos = [[1.5, 0, 0], [0, 0, 0], [1.5 * math.cos(t), 1.5 * math.sin(t), 0]]
    well = dict(r1=90.0, r2=105.0, r3=112.0, r4=130.0, rk2=50.0, rk3=40.0)
    e = _energy_kcal(pos, nmr_forces([{"atoms": [0, 1, 2], **well}]))
    # AMBER: r in degrees in the input, energy evaluated in radians
    expected = amber_flat_bottom(t, *map(math.radians, (well["r1"], well["r2"], well["r3"], well["r4"])),
                                 well["rk2"], well["rk3"])
    assert e == pytest.approx(expected, abs=1e-4)


def _dihedral_positions(phi_deg):
    phi = math.radians(phi_deg)
    return np.array([[0.0, 1.5, 0.0], [0, 0, 0], [1.5, 0, 0],
                     [1.5, 1.5 * math.cos(phi), 1.5 * math.sin(phi)]])


@pytest.mark.parametrize("phi", [-170.0, -60.0, -8.0, 0.0, 10.0, 17.0, 45.0, 175.0])
def test_dihedral_restraint_matches_amber(phi):
    well = dict(r1=-10.0, r2=-5.0, r3=15.0, r4=20.0, rk2=500.0, rk3=500.0)  # example config
    pos = _dihedral_positions(phi)
    measured = _dihedral_deg(pos)
    e = _energy_kcal(pos, nmr_forces([{"atoms": [0, 1, 2, 3], **well}]))
    mid = 0.5 * (well["r2"] + well["r3"])
    shifted = measured - 360.0 * math.floor((measured - mid + 180.0) / 360.0)
    expected = amber_flat_bottom(math.radians(shifted),
                                 *map(math.radians, (well["r1"], well["r2"], well["r3"], well["r4"])),
                                 well["rk2"], well["rk3"])
    assert e == pytest.approx(expected, rel=1e-5, abs=1e-4)


def test_dihedral_window_across_180():
    # window centred on 180: -175 deg must be 5 deg from the well, not 355
    well = dict(r1=160.0, r2=170.0, r3=178.0, r4=190.0, rk2=10.0, rk3=10.0)
    pos = _dihedral_positions(-175.0)
    e = _energy_kcal(pos, nmr_forces([{"atoms": [0, 1, 2, 3], **well}]))
    expected = 10.0 * math.radians(185.0 - 178.0) ** 2
    assert e == pytest.approx(expected, abs=1e-4)


def test_positional_restraint_units():
    # AMBER: E = wt * |dr|^2 (no 1/2); 0.5 A displacement, wt = 10 -> 2.5 kcal/mol
    ref = np.array([[1.0, 1.0, 1.0]])           # nm
    f = positional_force([0], ref, 10.0)
    assert f.getGlobalParameterDefaultValue(0) == pytest.approx(10.0 * KCAL_TO_KJ * 100)
    e = _energy_kcal([[10.0 + 0.5, 10.0, 10.0]], [f])
    assert e == pytest.approx(2.5, abs=1e-6)


def test_parse_disang_restrainer_format(tmp_path):
    p = tmp_path / "restr"
    p.write_text(
        "&rst\niat=5, 7, 9, 0,\nr1=1.0, r2=2.0, r3=3.0, r4=4.0\nrk2=10, rk3=20\n/\n\n"
        "# comment\n"
        "&rst\niat=1, 2, 3, 4,\nr1=-10.0, r2=-5.0, r3=15.0, r4=20.0\nrk2=500, rk3=500\n/\n"
    )
    rs = parse_disang(p)
    assert [r["atoms"] for r in rs] == [[4, 6, 8], [0, 1, 2, 3]]
    assert rs[1]["r1"] == -10.0 and rs[0]["rk3"] == 20.0


def test_parse_disang_rejects_unsupported_keys(tmp_path):
    p = tmp_path / "restr"
    p.write_text("&rst iat=1,2, r1=1, r2=2, r3=3, r4=4, rk2=1, rk3=1, igr1=3 /\n")
    with pytest.raises(ValueError, match="igr1"):
        parse_disang(p)


def test_read_mask_atoms(tmp_path):
    p = tmp_path / "restraint_atoms.dat"
    p.write_text("#Frame    AtomNum Atom   ResNum  Res   MolNum\n"
                 "       1        5 C           1 ACE         1\n"
                 "       1        7 N           2 ALA         1\n")
    assert read_mask_atoms(p) == [4, 6]


def test_netcdf_writer_roundtrip_and_resume(tmp_path):
    netcdf = pytest.importorskip("scipy.io")
    path = tmp_path / "t.nc"
    rng = np.random.default_rng(0)
    frames = rng.uniform(-50, 50, size=(5, 7, 3)).astype(np.float32)
    w = AmberNetCDFWriter(path, 7)
    for i, f in enumerate(frames[:4]):
        w.write_frame(f, float(i), [40.0, 41.0, 42.0], [90.0, 90.0, 90.0])
    w.close()
    # resume as if from a checkpoint after frame 2 -> frames 3,4 dropped
    w = AmberNetCDFWriter(path, 7, append=True, keep_frames=2)
    w.write_frame(frames[4], 99.0, [40.0, 41.0, 42.0], [90.0, 90.0, 90.0])
    w.close()

    nc = netcdf.netcdf_file(str(path), "r", mmap=False)
    assert nc.Conventions == b"AMBER"
    coords = nc.variables["coordinates"][:]
    assert coords.shape == (3, 7, 3)
    np.testing.assert_allclose(coords[:2], frames[:2])
    np.testing.assert_allclose(coords[2], frames[4])
    np.testing.assert_allclose(nc.variables["time"][:], [0.0, 1.0, 99.0])
    np.testing.assert_allclose(nc.variables["cell_lengths"][1], [40.0, 41.0, 42.0])
    nc.close()


def test_rst7_roundtrip(tmp_path):
    n = 5
    system = mm.System()
    for _ in range(n):
        system.addParticle(1.0)
    system.setDefaultPeriodicBoxVectors(mm.Vec3(3, 0, 0), mm.Vec3(0, 3.1, 0), mm.Vec3(0, 0, 3.2))
    ctx = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("Reference"))
    rng = np.random.default_rng(1)
    pos = rng.uniform(0, 3, size=(n, 3))
    vel = rng.normal(0, 0.5, size=(n, 3))
    ctx.setPositions(pos)
    ctx.setVelocities(vel)
    ctx.setTime(12.5)
    write_rst7(tmp_path / "x.rst7",
               ctx.getState(getPositions=True, getVelocities=True))
    crd = app.AmberInpcrdFile(str(tmp_path / "x.rst7"))
    np.testing.assert_allclose(crd.getPositions(asNumpy=True).value_in_unit(unit.nanometer), pos, atol=1e-7)
    np.testing.assert_allclose(crd.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond),
                               vel, atol=1e-6)
    np.testing.assert_allclose([crd.boxVectors[i][i].value_in_unit(unit.nanometer) for i in range(3)],
                               [3, 3.1, 3.2], atol=1e-7)


def test_trim_log(tmp_path):
    p = tmp_path / "x.log"
    p.write_text('#"Step","Time"\n500,1.0\n1000,2.0\n1500,3.0\n')
    trim_log(p, 1000)
    assert p.read_text() == '#"Step","Time"\n500,1.0\n1000,2.0\n'
