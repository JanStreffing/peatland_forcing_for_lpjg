#!/usr/bin/env python3
"""Write the full reduced Gaussian grid of a regridded LUH3 file, at full precision.

generate_peat_map.py needs a target grid with a longitude and latitude per
cell. The regridded LUH3 files cannot serve as they are: the ncells ones carry
coordinates rounded to four decimals, and the TL511 and TCO639 ones carry a
latitude array that was never written. LPJ-GUESS finds a cell in the peat map
to within 5e-5 degrees, so four decimals are not enough.

Only the row lengths are taken from the LUH3 file. The latitudes are the
Gaussian ones, computed, and each row starts at the Greenwich meridian.

Usage:
  make_gaussian_target_grid.py <luh3_states.nc> <out_grid.nc>
"""
import sys

import numpy as np
import netCDF4 as nc


def row_lengths(ds):
    if "reduced_points" in ds.variables:
        return np.asarray(ds["reduced_points"][:], dtype=np.int64)
    lat = np.asarray(ds["lat"][:], dtype=np.float64)
    # Cells are stored row by row, so a row is a run of equal latitude.
    edges = np.flatnonzero(np.diff(lat) != 0.0) + 1
    return np.diff(np.concatenate([[0], edges, [lat.size]]))


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    src, out = sys.argv[1:3]
    ds = nc.Dataset(src)
    ds.set_auto_mask(False)
    npts = row_lengths(ds)
    if npts.size % 2 or not np.array_equal(npts, npts[::-1]):
        raise SystemExit(f"{src}: row lengths {npts.size} rows are not a "
                         "symmetric Gaussian grid")
    nodes = np.polynomial.legendre.leggauss(npts.size)[0]
    row_lat = np.degrees(np.arcsin(nodes[::-1]))
    lat = np.repeat(row_lat, npts)
    lon = np.concatenate([360.0 * np.arange(n) / n for n in npts])

    if "lon" in ds.variables:
        # The rounded coordinates in the file must agree with the computed ones.
        off = max(np.abs(np.asarray(ds["lat"][:]) - lat).max(),
                  np.abs((np.asarray(ds["lon"][:]) - lon + 180.0) % 360.0 - 180.0).max())
        if off > 1e-3:
            raise SystemExit(f"{src}: computed grid is {off} degrees off the file's own")
        print(f"  agrees with the file's rounded coordinates to {off:.1e} degrees")
    ds.close()

    with nc.Dataset(out, "w") as tgt:
        tgt.createDimension("ncells", lat.size)
        for name, vals, unit in (("lon", lon, "degrees_east"), ("lat", lat, "degrees_north")):
            v = tgt.createVariable(name, "f8", ("ncells",))
            v.units = unit
            v[:] = vals
    print(f"{out}: {npts.size} rows, {lat.size} cells, row lengths {npts.min()}-{npts.max()}")


if __name__ == "__main__":
    main()
