#!/usr/bin/env python3
"""Cap a peat map at the smallest natural fraction LUH3 ever has in each cell.

LPJ-GUESS takes the fixed peat area out of natural land first. Where the map
holds more peat than LUH3 has natural land, natural is zero and the rest comes
out of crops and pasture. The corrected map is

    natural = min(primf + primn + secdf + secdn, 1 - crops - pasture - urban)
    peat    = min(peat, min over all years of natural)

so peat always fits inside natural. The second term is the natural land
LPJ-GUESS is left with after it closes the separately rounded LUH3 classes to
one. It lowers the cap by at most 6.6e-7, which keeps peat from borrowing even
a rounding error from managed land.

This is the rule of a270270's correct_peat_lund_method.py, which made the
TCO95 and TCO319 maps. That script needs a LUH3 file with a coordinate per
cell in the order of the peat map. This one also takes the reduced Gaussian
and the native layout, and reproduces both of those maps.

Cells LUH3 does not cover (fill value in every year) are set to zero, as in
that file. Pass "keep" as the last argument to leave them as they are.
The result is rounded down, never up, so the written peat cannot exceed natural.

Coordinates are written as they were read. The regridded LUH3 files carry
theirs rounded to four decimals, so they are no source for better ones.

The LUH3 file can be a regridded one, on the grid of the peat map, or the
native 0.25 degree one. The regrids are nearest neighbour, and a minimum over
time commutes with that, so sampling the native minimum at the nearest cell
gives what a regrid would. Where a point lies on the edge between two native
cells the smaller of the two is taken. It otherwise differs at the coast only: cdo's remapnn skips
missing source cells, this does not, so a few coastal cells a regrid would
fill are uncovered here, and get zero peat. That errs on the safe side.

Several LUH3 files can be given. The cap is then the smallest of them, so the
map fits whichever of those versions the model is later run with.

Usage:
  correct_peat_map_luh3.py <peat_frac.txt> <out.txt> <luh3_states.nc> [...] [keep]
"""
import sys

import numpy as np
import netCDF4 as nc
from scipy.spatial import cKDTree

NATURAL = ("primf", "primn", "secdf", "secdn")
MANAGED = ("c3ann", "c4ann", "c3per", "c4per", "c3nfx", "pastr", "range", "urban")
# Degrees. The peat map and a regridded LUH3 file are on the same grid, so a
# match is exact to the precision the coordinates were written with, which is
# two decimals in the older maps. Far below the spacing of any grid used here.
MATCH_TOL = 0.01
# In cells of the native grid. Closer than this to an edge counts as on it.
TIE = 1e-6
# Rows of LUH3 read at a time, in bytes of float64.
BLOCK_BYTES = 2.0e9
DECIMALS = 8


def natural_minimum(ds):
    """Smallest natural fraction over all years, inf where LUH3 has no data."""
    shape = ds[NATURAL[0]].shape
    ntime, ncell = shape[0], int(np.prod(shape[1:]))
    block = max(1, int(BLOCK_BYTES / (8.0 * ncell)))
    nat_min = np.full(shape[1:], np.inf)
    for t0 in range(0, ntime, block):
        t1 = min(t0 + block, ntime)
        nat = np.zeros((t1 - t0,) + shape[1:])
        managed = np.zeros_like(nat)
        valid = np.ones(nat.shape, dtype=bool)
        for name in NATURAL + MANAGED:
            x = np.asarray(ds[name][t0:t1], dtype=np.float64)
            # A fill value in any class makes the sums meaningless.
            ok = np.isfinite(x) & (x > -1e-4) & (x < 1.0 + 1e-4)
            valid &= ok
            x = np.where(ok, x, 0.0)
            if name in NATURAL:
                nat += x
            else:
                managed += x
        nat = np.clip(nat, 0.0, 1.0)
        # What LPJ-GUESS keeps as natural once it has closed LUH3 to one.
        closed = np.clip(np.where(managed > 1.0, 0.0, 1.0 - managed), 0.0, 1.0)
        nat = np.where(valid, np.minimum(nat, closed), np.inf)
        nat_min = np.minimum(nat_min, nat.min(axis=0))
        print(f"  years {t0}-{t1 - 1} of {ntime}", flush=True)
    return nat_min


def cell_coordinates(ds):
    """Longitude and latitude of every cell of a regridded file."""
    if "lon" in ds.variables:
        lon = np.asarray(ds["lon"][:], dtype=np.float64) % 360.0
        lat = np.asarray(ds["lat"][:], dtype=np.float64)
        return lon, lat
    # cdo's reduced Gaussian layout: a latitude per row and the number of
    # points in it, the points starting at the Greenwich meridian.
    npts = np.asarray(ds["reduced_points"][:], dtype=np.int64)
    row_lat = np.asarray(ds["lat"][:], dtype=np.float64)
    if not (np.all(np.diff(row_lat) < 0) and abs(row_lat[0]) <= 90.0):
        # The regrids in this layout carry an unwritten latitude array.
        # The rows are the Gaussian latitudes, north to south.
        nodes = np.polynomial.legendre.leggauss(npts.size)[0]
        row_lat = np.degrees(np.arcsin(nodes[::-1]))
    lat = np.repeat(row_lat, npts)
    lon = np.concatenate([360.0 * np.arange(n) / n for n in npts])
    return lon, lat


def cap_from(states_path, plon, plat):
    """Natural minimum at each peat row, and whether the row is on the grid."""
    ds = nc.Dataset(states_path)
    ds.set_auto_mask(False)
    nat_min = natural_minimum(ds)
    if nat_min.ndim == 2:
        # Native regular grid: the nearest cell is index arithmetic.
        glat = np.asarray(ds["lat"][:], dtype=np.float64)
        glon = np.asarray(ds["lon"][:], dtype=np.float64)
        dlat, dlon = glat[1] - glat[0], glon[1] - glon[0]
        y = (plat - glat[0]) / dlat
        x = ((plon - glon[0]) % 360.0) / dlon
        # Gaussian longitudes often fall exactly on the edge between two
        # 0.25 degree cells. Either is a nearest neighbour and a regrid may
        # have taken either, so take the smaller of the two.
        cap = np.full(plon.size, np.inf)
        for ey in (-TIE, TIE):
            j = np.clip(np.floor(y + 0.5 + ey).astype(np.int64), 0, glat.size - 1)
            for ex in (-TIE, TIE):
                i = np.floor(x + 0.5 + ex).astype(np.int64) % glon.size
                cap = np.minimum(cap, nat_min[j, i])
        matched = np.ones(plon.size, dtype=bool)
        kind = f"native {glat.size}x{glon.size}"
    else:
        lon, lat = cell_coordinates(ds)

        # Longitude wraps, so match on the sphere's chord rather than on degrees.
        def xyz(lo, la):
            lo, la = np.radians(lo), np.radians(la)
            return np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]

        dist, idx = cKDTree(xyz(lon, lat)).query(xyz(plon, plat))
        matched = np.degrees(dist) < MATCH_TOL
        cap = np.where(matched, nat_min[idx], np.inf)
        kind = f"regridded, {lon.size} cells"
    ds.close()
    print(f"  {states_path}: {kind}")
    return cap, matched


def main():
    args = sys.argv[1:]
    keep_uncovered = bool(args) and args[-1] == "keep"
    if keep_uncovered:
        args = args[:-1]
    if len(args) < 3:
        raise SystemExit(__doc__)
    peat_path, out_path, states_paths = args[0], args[1], args[2:]

    with open(peat_path) as f:
        header = f.readline().rstrip("\n")
        coords = [line.split()[:2] for line in f if line.strip()]
    pfrac = np.loadtxt(peat_path, skiprows=1, usecols=2)
    plon = np.array([float(c[0]) for c in coords]) % 360.0
    plat = np.array([float(c[1]) for c in coords])

    cap = np.full(pfrac.size, np.inf)
    matched = np.ones(pfrac.size, dtype=bool)
    for path in states_paths:
        c, m = cap_from(path, plon, plat)
        cap, matched = np.minimum(cap, c), matched & m
    uncovered = matched & ~np.isfinite(cap)
    if not keep_uncovered:
        cap[uncovered] = 0.0

    scale = 10.0 ** DECIMALS
    capped = np.floor(np.minimum(pfrac, cap) * scale) / scale
    new = np.where(cap < pfrac, capped, pfrac)

    with open(out_path, "w") as f:
        f.write(header + "\n")
        f.writelines(f"{lo}\t{la}\t{v:.8f}\n" for (lo, la), v in zip(coords, new))

    changed = new < pfrac
    print(f"{peat_path} -> {out_path}")
    print(f"  peat rows {pfrac.size}, not on the LUH3 grid {int((~matched).sum())}, "
          f"not covered by LUH3 {int(uncovered.sum())} "
          f"(with peat {int((uncovered & (pfrac > 0)).sum())}, "
          f"{'kept' if keep_uncovered else 'set to zero'})")
    print(f"  cells with peat {int((pfrac > 0).sum())} -> {int((new > 0).sum())}, "
          f"lowered {int(changed.sum())}, set to zero {int(((new == 0) & (pfrac > 0)).sum())}")
    print(f"  sum of peat fractions {pfrac.sum():.3f} -> {new.sum():.3f} "
          f"({100 * (new.sum() / pfrac.sum() - 1):+.1f}%)")


if __name__ == "__main__":
    main()
