import os
import json
import warnings
from datetime import datetime, timezone

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import metpy.calc as mpcalc
from metpy.units import units
from herbie import Herbie
from cartopy import crs as ccrs, feature as cfeature
from scipy.ndimage import gaussian_filter

# silence the harmless cfgrib/xarray FutureWarning
warnings.filterwarnings(
    "ignore",
    message=".*default value for compat will change.*",
    category=FutureWarning,
)

output_dir = 'assets/maps/q_vector'
os.makedirs(output_dir, exist_ok=True)


# subsets the GFS data so you don't process/plot the whole globe
def subset(ds, lon_min=-130, lon_max=-60, lat_min=15, lat_max=58):
    lon_min_convert = lon_min % 360
    lon_max_convert = lon_max % 360
    return ds.sel(
        latitude=slice(lat_max, lat_min),  # GFS latitude is descending
        longitude=slice(lon_min_convert, lon_max_convert),
    )


# loads one pressure level, cropped to CONUS and pulled into memory
def load_level(H, level):
    ds = H.xarray(f":(?:HGT|UGRD|VGRD|TMP):{level} mb:", remove_grib=True)
    return subset(ds).load()


# computes the q-vector and its divergence for one level
def q_for_level(ds, level, dx, dy):
    u = ds['u']
    v = ds['v']
    t = ds['t']
    uq, vq = mpcalc.q_vector(u, v, t, level * units.hPa, dx=dx, dy=dy)
    return uq, vq


# function that builds the images from the Herbie object
def extract_build_data(H, fxx):
    ds_700 = load_level(H, 700)
    ds_600 = load_level(H, 600)
    ds_500 = load_level(H, 500)
    ds_400 = load_level(H, 400)

    lat = ds_700['latitude']
    lon = ds_700['longitude']
    gh_500 = ds_500['gh'] / 10

    # strip time (same valid time across levels)
    time_valid = pd.Timestamp(ds_700['time'].values)
    time_str = time_valid.strftime('%HZ')
    valid_time = pd.Timestamp(ds_700['valid_time'].values).strftime('%HZ %a %b %d %Y')

    # grid spacing computed once; all levels share the same cropped grid
    dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat)

    uq_700, vq_700 = q_for_level(ds_700, 700, dx, dy)
    uq_600, vq_600 = q_for_level(ds_600, 600, dx, dy)
    uq_500, vq_500 = q_for_level(ds_500, 500, dx, dy)
    uq_400, vq_400 = q_for_level(ds_400, 400, dx, dy)

    #### Averaging layer ######
    uqvect_mean = (uq_700 + uq_600 + uq_500 + uq_400) / 4
    vqvect_mean = (vq_700 + vq_600 + vq_500 + vq_400) / 4
    q_div_mean = 2 * mpcalc.divergence(uqvect_mean, vqvect_mean, dx=dx, dy=dy)
    q_div_smooth = gaussian_filter(q_div_mean * 1e18, 4)
    gh_smooth_500 = gaussian_filter(gh_500.values, 2)

    # plot settings
    proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35,
                                 standard_parallels=(30, 60))
    fig = plt.figure(figsize=(20, 12))
    ax = plt.axes(projection=proj)
    ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())

    # q-vector convergence/divergence
    clevs_qdiv = list(range(-35, -4, 5)) + list(range(5, 36, 5))
    norm = mcolors.BoundaryNorm(clevs_qdiv, ncolors=256, extend='both')
    cf = ax.contourf(lon, lat, q_div_smooth, clevs_qdiv, cmap='bwr_r', norm=norm,
                     extend='both', transform=ccrs.PlateCarree())
    line1 = ax.contour(lon, lat, q_div_smooth, colors='black', linewidths=1.5,
                       transform=ccrs.PlateCarree())
    ax.clabel(line1, inline=True, colors='black', fontsize=14, fmt='%d')

    cb = plt.colorbar(cf, ax=ax, orientation='horizontal', shrink=.4, pad=.02, aspect=25)
    cb.set_label('m^2kg^-1s^-1', fontsize=16)
    cb.ax.tick_params(labelsize=16)

    # 500 mb GH contours
    line = ax.contour(lon, lat, gh_smooth_500, levels=list(range(500, 600, 6)),
                      colors='black', linewidths=3, transform=ccrs.PlateCarree())
    ax.clabel(line, inline=True, colors='black', fontsize=14, fmt='%d')

    # vectors: lon/lat are 1D, so thin each axis separately
    step = 5
    ax.quiver(lon.values[::step], lat.values[::step],
              uqvect_mean.values[::step, ::step],
              vqvect_mean.values[::step, ::step],
              pivot='mid', color='black',
              scale=5e-12, scale_units='inches',
              transform=ccrs.PlateCarree())

    # additional plot settings
    ax.coastlines('50m')
    ax.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax.add_feature(cfeature.STATES.with_scale('50m'))
    ax.set_title(f'Model: {time_str} {H.model.upper()} | F{H.fxx:03d}', fontsize=15, loc='left')
    ax.set_title('400-700mb Q-Vectors\n (m^2kg^-1s^-1)', fontsize=20, loc='center')
    ax.set_title(f'\nValid: {valid_time}', fontsize=15, loc='right')

    out_path = os.path.join(output_dir, f'{fxx:02d}.png')
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)

    return out_path, valid_time


# pull from model
target_time = (pd.Timestamp.now(tz='UTC') - pd.Timedelta(hours=4)).tz_localize(None).floor('6h')
print('looking for run', target_time)

fxx_list = list(range(0, 25, 1)) + list(range(27, 121, 3))
frames = []
for fxx in fxx_list:
    try:
        H = Herbie(target_time, model="gfs", product="pgrb2.0p25", fxx=fxx)
        out_path, valid_time = extract_build_data(H, fxx)
        frames.append({"fxx": fxx, "file": os.path.basename(out_path), "valid_time": valid_time})
        print(f"Saved frame f{fxx:02d}")
    except Exception as e:
        print(f"Skipping fxx={fxx}: {e}")
        continue

# write document dictionary so the website JS knows what frames exist
document = {
    "run_time": target_time.strftime("%Y-%m-%d %H:%M"),
    "generated_at": datetime.now(timezone.utc).isoformat().replace('+00:00', 'Z'),
    "frames": frames,
}
with open(os.path.join(output_dir, "document.json"), 'w') as f:
    json.dump(document, f, indent=2)
print(f"Done {len(frames)} frames saved")