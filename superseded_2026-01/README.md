# Superseded, do not use

Peat maps made in January 2026 by `create_peat_files.py`, and their plots.

- Coordinates were written with two decimals. LPJ-GUESS looks a cell up in the map to
  within 5e-5 degrees, so these do not load under the fixed-peat code.
- The TL511 and TCO639 maps were laid out on the latitude array of the regridded LUH3
  files, which was never written for those two grids. TL511 has latitude 0.00 on every
  row, TCO639 has half its rows on two latitudes. Their peat values sit in the wrong places.

Kept only so that a run made with them can be traced. The maps one level up replace them.
