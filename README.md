# peatland_forcing_for_lpjg

This repository contains a small utility to generate LPJ-GUESS peatland fraction forcing files on various model grids.

## Contents

- `TL255_peat_frac.txt`: source peat fraction file (Lon/Lat/Peat_Frac).
- `create_peat_files.py`: interpolates the TL255 peat fractions onto target grids using nearest-neighbor interpolation.

## Usage

Run:

```bash
python3 create_peat_files.py
```

The script reads the source peat file from the local folder and derives target grids from available regridded land-use NetCDFs under:

- `../land_use/regridded/multiple-transitions_input4MIPs_landState_CMIP_UofMD-landState-3-1-1_gn_0850-2023_<RES>.nc`

If those files are not present, it will also attempt to interpolate to any `../gridlists/gridlist_*_CORE2.txt` grids.

## Outputs

The script writes `*_peat_frac.txt` files:

- into this directory, and
- into `generated/`

It also writes quick-look plots into `plots/` (skipping plots for very large grids).
