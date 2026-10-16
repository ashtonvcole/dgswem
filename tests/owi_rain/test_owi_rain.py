#!/usr/bin/env python3
"""Tests of the OWI_RAIN module, through the owi_rain_driver executable.

Each test writes a tiny set of OWI files into a temporary directory, runs the
driver in a fresh process (the module holds global state), and compares the
interpolated rainfall against hand-computed values.

Run from the tests directory with:
    pytest owi_rain --binpath ../build
"""
import os
import subprocess

import numpy as np
import pytest

MMHR = 3.6e6  # (mm/hr) -> (m/s) divisor


def write_owi(path, records, swlon=0.0, swlat=0.0, dx=1.0, dy=1.0, bad=False):
    """Write an OWI file. Each record is an array of shape (ilat, ilon)."""
    with open(path, "w") as f:
        f.write("Test OWI precipitation".ljust(55)
                + f"{2026010100:10d}" + " " * 5 + f"{2026010123:10d}\n")
        for k, vals in enumerate(records):
            ilat, ilon = vals.shape
            ilat_str = "  ab" if bad else f"{ilat:4d}"
            f.write("iLat=" + ilat_str + "iLong=" + f"{ilon:4d}"
                    + "DX=" + f"{dx:6.4f}" + "DY=" + f"{dy:6.4f}"
                    + "SWLat=" + f"{swlat:8.4f}" + "SWLon=" + f"{swlon:8.4f}"
                    + "DT=" + f"{202601010000 + 100 * k:12d}\n")
            flat = vals.flatten()  # Longitude varies fastest
            for i in range(0, len(flat), 8):
                f.write("".join(f"{v:10.4f}" for v in flat[i:i + 8]) + "\n")


def uniform(value, ilat=2, ilon=2):
    """A spatially constant field."""
    return np.full((ilat, ilon), float(value))


def run_driver(binpath, workdir, points, times, num_sets=1, offset=0,
               t0=0.0, timinc=3600.0, mask=0, mmhr=1):
    """Run the driver and return an array of shape (len(times), len(points))."""
    exe = os.path.join(binpath, "owi_rain_driver")
    if not os.path.exists(exe):
        pytest.skip("owi_rain_driver executable not found. Skipping...")
    with open(os.path.join(workdir, "fort.27"), "w") as f:
        f.write(f"{num_sets}\n{offset}\n")
    with open(os.path.join(workdir, "points.txt"), "w") as f:
        f.write(f"{len(points)}\n")
        f.writelines(f"{lon} {lat}\n" for lon, lat in points)
    with open(os.path.join(workdir, "times.txt"), "w") as f:
        f.write(f"{len(times)}\n")
        f.writelines(f"{t}\n" for t in times)
    args = [exe, str(workdir), str(t0), str(timinc), str(mask), str(mmhr)]
    return subprocess.run(args, cwd=workdir, check=False,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          text=True)


def results(proc):
    """Parse the RESULT lines of the driver output into an array."""
    rows = [list(map(float, line.split()[2:]))
            for line in proc.stdout.splitlines() if line.startswith("RESULT")]
    return np.array(rows)


def test_linear_in_time(binpath, tmp_path):
    """Rain is linear between snapshots, and the start time is honored.

    Three uniform snapshots of 36, 72, and 108 mm/hr (1e-5, 2e-5, 3e-5 m/s)
    at t0, t0 + 3600 s, and t0 + 7200 s, with a nonzero t0 to exercise STATIM.
    Querying across a snapshot boundary also exercises the snapshot cycling.
    """
    t0 = 86400.0
    write_owi(tmp_path / "fort.425", [uniform(v) for v in (36, 72, 108)])
    times = t0 + np.array([0, 1800, 3600, 5400, 7200])
    proc = run_driver(binpath, tmp_path, [(0.5, 0.5)], times, t0=t0)
    assert proc.returncode == 0, proc.stderr
    expected = np.array([1.0, 1.5, 2.0, 2.5, 3.0]) * 1e-5
    np.testing.assert_allclose(results(proc)[:, 0], expected, rtol=1e-9)


def test_bilinear_in_space(binpath, tmp_path):
    """Bilinear interpolation reproduces a linear field exactly.

    The field is f = 1 + 2 (lon - swlon) + 3 (lat - swlat), which bilinear
    weights represent exactly. Points outside the grid get zero.
    """
    swlon, swlat = -95.0, 28.0
    lons = swlon + np.arange(4)
    lats = swlat + np.arange(4)
    field = 1 + 2 * (lons[None, :] - swlon) + 3 * (lats[:, None] - swlat)
    write_owi(tmp_path / "fort.425", [field], swlon=swlon, swlat=swlat)
    points = [(-94.5, 28.5), (-92.75, 29.75), (10.0, 10.0)]
    proc = run_driver(binpath, tmp_path, points, [0.0], mmhr=0)
    assert proc.returncode == 0, proc.stderr
    f = lambda lon, lat: 1 + 2 * (lon - swlon) + 3 * (lat - swlat)
    expected = [f(*points[0]), f(*points[1]), 0.0]
    np.testing.assert_allclose(results(proc)[0], expected, rtol=1e-9)


@pytest.mark.parametrize("mask, expected", [(1, 3.0), (0, -2496.75)])
def test_negative_masking(binpath, tmp_path, mask, expected):
    """A negative land flag is masked to zero only when requested.

    A 2 x 2 cell with a -9999 flag at the southwest corner and 4 elsewhere,
    sampled at the cell center (all weights 0.25). Masked: 0.25 * (0 + 4 + 4 +
    4) = 3. Unmasked: 0.25 * (-9999 + 12) = -2496.75.
    """
    field = np.array([[-9999.0, 4.0], [4.0, 4.0]])
    write_owi(tmp_path / "fort.425", [field])
    proc = run_driver(binpath, tmp_path, [(0.5, 0.5)], [0.0], mask=mask,
                      mmhr=0)
    assert proc.returncode == 0, proc.stderr
    np.testing.assert_allclose(results(proc)[0, 0], expected, rtol=1e-9)


def test_eof_decays_to_zero(binpath, tmp_path):
    """The end of the file inserts blank snapshots instead of failing.

    Two records of 36 and 72 mm/hr, queried past the end. Rain falls linearly
    from 2e-5 to zero over the interval after the last record, then stays 0.
    """
    write_owi(tmp_path / "fort.425", [uniform(36), uniform(72)])
    times = [3600.0, 5400.0, 7200.0, 10800.0]
    proc = run_driver(binpath, tmp_path, [(0.5, 0.5)], times)
    assert proc.returncode == 0, proc.stderr
    np.testing.assert_allclose(results(proc)[:, 0],
                               [2e-5, 1e-5, 0.0, 0.0], atol=1e-15)


@pytest.mark.parametrize("offset, values, times, expected", [
    # One blank snapshot is inserted before the first record
    (1, (36, 72), [0, 3600, 7200], [0.0, 1e-5, 2e-5]),
    # The first record is skipped
    (-1, (36, 72, 108), [0, 3600], [2e-5, 3e-5]),
])
def test_blank_and_skip_offset(binpath, tmp_path, offset, values, times,
                               expected):
    """The fort.27 offset inserts blank snapshots or skips records."""
    write_owi(tmp_path / "fort.425", [uniform(v) for v in values])
    proc = run_driver(binpath, tmp_path, [(0.5, 0.5)], times, offset=offset)
    assert proc.returncode == 0, proc.stderr
    np.testing.assert_allclose(results(proc)[:, 0], expected, atol=1e-15)


def test_region_overrides_basin(binpath, tmp_path):
    """The region grid overwrites the basin grid where the two overlap.

    A 4 x 4 basin at 36 mm/hr and a 2 x 2 region over [1, 2] x [1, 2] at 108
    mm/hr. Points in the region get the region value, others the basin value.
    """
    write_owi(tmp_path / "fort.425", [uniform(36, 4, 4)])
    write_owi(tmp_path / "fort.427", [uniform(108)], swlon=1.0, swlat=1.0)
    points = [(1.5, 1.5), (0.5, 0.5), (2.5, 2.5)]
    proc = run_driver(binpath, tmp_path, points, [0.0], num_sets=2)
    assert proc.returncode == 0, proc.stderr
    np.testing.assert_allclose(results(proc)[0], [3e-5, 1e-5, 1e-5],
                               rtol=1e-9)


def test_format_error_is_fatal(binpath, tmp_path):
    """A malformed record header ends the run with a nonzero exit status."""
    write_owi(tmp_path / "fort.425", [uniform(36)], bad=True)
    proc = run_driver(binpath, tmp_path, [(0.5, 0.5)], [0.0])
    assert proc.returncode != 0
    assert "format error" in proc.stderr
