#!/usr/bin/env python3
"""
Create peat fraction files for different atmospheric resolutions.

Source: TL255_peat_frac.txt (TL255 resolution, 88838 cells)
Targets: TCO95, TCO319 (interpolated via nearest-neighbor)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.interpolate import griddata
import netCDF4
import os

INPUT_DIR = "/work/ab0246/a270092/input/lpj-guess"
PEAT_DIR = os.path.join(INPUT_DIR, "peat")
GRIDLIST_DIR = os.path.join(INPUT_DIR, "gridlists")
LAND_USE_REGRID_DIR = os.path.join(INPUT_DIR, "land_use", "regridded")
OUTPUT_DIR = os.path.join(PEAT_DIR, "generated")
PLOT_DIR = os.path.join(PEAT_DIR, "plots")

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)


def read_peat_file(filepath):
    return pd.read_csv(filepath, sep=r'\s+')


def read_gridlist(filepath):
    return pd.read_csv(filepath, sep=r'\s+', header=None, names=['Lat', 'Lon'])


def read_grid_from_netcdf(filepath):
    ds = netCDF4.Dataset(filepath, mode="r")
    try:
        if "lon" in ds.variables and "lat" in ds.variables:
            lon = np.asarray(ds.variables["lon"][:])
            lat = np.asarray(ds.variables["lat"][:])
        elif "lat" in ds.variables and "reduced_points" in ds.variables:
            lats = np.asarray(ds.variables["lat"][:])
            reduced_points = np.asarray(ds.variables["reduced_points"][:]).astype(int)
            lon_list = []
            lat_list = []
            for lat_val, nlon in zip(lats, reduced_points):
                if nlon <= 0:
                    continue
                lon_row = (np.arange(nlon, dtype=float) * (360.0 / float(nlon)))
                lon_list.append(lon_row)
                lat_list.append(np.full(nlon, float(lat_val)))
            if not lon_list:
                raise ValueError(f"Could not build reduced grid lon/lat from file: {filepath}")
            lon = np.concatenate(lon_list)
            lat = np.concatenate(lat_list)
        else:
            raise KeyError(f"Could not find recognizable coordinate variables in: {filepath}")
    finally:
        ds.close()

    return pd.DataFrame({'Lat': lat, 'Lon': lon})


def interpolate_peat_to_grid(source_peat, target_grid):
    """Interpolate peat fractions using nearest neighbor."""
    source_lon = source_peat['Lon'].values
    source_lat = source_peat['Lat'].values
    source_peat_frac = source_peat['Peat_Frac'].values
    source_mask = np.isfinite(source_lon) & np.isfinite(source_lat) & np.isfinite(source_peat_frac)
    source_coords = np.column_stack([source_lon[source_mask], source_lat[source_mask]])
    source_vals = source_peat_frac[source_mask]

    target_lon = target_grid['Lon'].values
    target_lat = target_grid['Lat'].values
    target_mask = np.isfinite(target_lon) & np.isfinite(target_lat)
    target_coords = np.column_stack([target_lon[target_mask], target_lat[target_mask]])

    peat_interp_valid = griddata(source_coords, source_vals,
                                 target_coords, method='nearest')
    peat_interp_valid = np.nan_to_num(peat_interp_valid, nan=0.0)

    peat_interp = np.zeros(len(target_grid), dtype=float)
    peat_interp[target_mask] = peat_interp_valid
    return peat_interp


def write_peat_file(grid, peat_values, filepath):
    """Write peat fraction file in LPJ-GUESS format."""
    with open(filepath, 'w') as f:
        f.write("Lon\tLat\tPeat_Frac\n")
        for i in range(len(grid)):
            f.write(f"{grid['Lon'].iloc[i]:.2f}\t{grid['Lat'].iloc[i]:.2f}\t{peat_values[i]:.8f}\n")


def plot_peat_map(peat_df, title, output_path):
    """Plot peat fraction on a global map."""
    fig, ax = plt.subplots(figsize=(14, 8), subplot_kw={'projection': ccrs.Robinson()})
    ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor='lightgray', alpha=0.3)
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue', alpha=0.3)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, alpha=0.5)

    lons, lats, peat = peat_df['Lon'].values, peat_df['Lat'].values, peat_df['Peat_Frac'].values
    mask = peat > 0.001
    
    if mask.sum() > 0:
        scatter = ax.scatter(lons[mask], lats[mask], c=peat[mask],
                             s=2 if len(peat_df) > 50000 else 5,
                             cmap='YlOrBr', vmin=0, vmax=1,
                             transform=ccrs.PlateCarree(), alpha=0.8)
        cbar = plt.colorbar(scatter, ax=ax, orientation='horizontal', pad=0.05, shrink=0.7)
        cbar.set_label('Peat Fraction', fontsize=12)

    ax.set_title(title, fontsize=14, fontweight='bold')
    stats = f"Total: {len(peat_df):,} | With peat: {mask.sum():,} | Max: {peat.max():.3f}"
    ax.text(0.02, 0.02, stats, transform=ax.transAxes, fontsize=9,
            verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

def main():
    # Read source TL255 peat data
    tl255_peat_file = os.path.join(PEAT_DIR, "TL255_peat_frac.txt")
    tl255_peat = read_peat_file(tl255_peat_file)

    # Plot TL255 peat map
    plot_peat_map(tl255_peat, "Peatland Fraction - TL255",
                  os.path.join(PLOT_DIR, "peat_frac_TL255_map.png"))

    if os.path.isdir(LAND_USE_REGRID_DIR):
        prefix = "multiple-transitions_input4MIPs_landState_CMIP_UofMD-landState-3-1-1_gn_0850-2023_"
        for filename in sorted(os.listdir(LAND_USE_REGRID_DIR)):
            if not (filename.startswith(prefix) and filename.endswith(".nc")):
                continue

            grid_id = filename[len(prefix):-len(".nc")]
            if grid_id == "TL255":
                continue

            target_grid = read_grid_from_netcdf(os.path.join(LAND_USE_REGRID_DIR, filename))
            peat_values = interpolate_peat_to_grid(tl255_peat, target_grid)
            peat_df = target_grid.copy()
            peat_df['Peat_Frac'] = peat_values

            write_peat_file(target_grid, peat_values, os.path.join(OUTPUT_DIR, f"{grid_id}_peat_frac.txt"))
            write_peat_file(target_grid, peat_values, os.path.join(PEAT_DIR, f"{grid_id}_peat_frac.txt"))
            if len(peat_df) <= 300000:
                plot_peat_map(peat_df, f"Peatland Fraction - {grid_id}",
                              os.path.join(PLOT_DIR, f"peat_frac_{grid_id}_map.png"))

    for filename in sorted(os.listdir(GRIDLIST_DIR)):
        if not (filename.startswith("gridlist_") and filename.endswith("_CORE2.txt")):
            continue

        grid_id = filename[len("gridlist_"):-len("_CORE2.txt")]
        if grid_id == "TL255":
            continue

        target_gridlist = read_gridlist(os.path.join(GRIDLIST_DIR, filename))
        peat_values = interpolate_peat_to_grid(tl255_peat, target_gridlist)
        peat_df = target_gridlist.copy()
        peat_df['Peat_Frac'] = peat_values

        write_peat_file(target_gridlist, peat_values, os.path.join(OUTPUT_DIR, f"{grid_id}_peat_frac.txt"))
        write_peat_file(target_gridlist, peat_values, os.path.join(PEAT_DIR, f"{grid_id}_peat_frac.txt"))
        plot_peat_map(peat_df, f"Peatland Fraction - {grid_id}",
                      os.path.join(PLOT_DIR, f"peat_frac_{grid_id}_map.png"))


if __name__ == "__main__":
    main()
