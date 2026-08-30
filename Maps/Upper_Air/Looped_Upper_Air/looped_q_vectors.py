import os

import matplotlib
# matplotlib.use("Agg")

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import metpy
import metpy.calc as mpcalc
from metpy.units import units
from herbie import Herbie
from herbie.toolbox import EasyMap, pc
from herbie import paint
import cartopy
from cartopy import crs as ccrs, feature as cfeature
import scipy
from scipy.ndimage import gaussian_filter
import pandas as pd
from datetime import datetime, timedelta
import json
import matplotlib.colors as mcolors

output_dir = 'assets/maps/q_vector'
os.makedirs(output_dir, exist_ok=True)

#function that builds the images with ds and H as arguments
def extract_build_data (H, fxx):
    ###### 700mb #######

    ds_700 = H.xarray(":(?:HGT|UGRD|VGRD|TMP):700 mb:", remove_grib=True)
    print(ds_700)

    gh_700 = ds_700['gh'] / 10
    u_700 = ds_700['u']
    v_700 = ds_700['v']
    lat = ds_700['latitude']
    lon = ds_700['longitude']
    temp_700 = ds_700['t']
    u_kts_700 = ds_700['u'] * 1.944
    v_kts_700 = ds_700['v'] * 1.944
    level_700 = 700 * units.hPa

    # strip time (only need this once — same valid time across levels)
    time_valid = pd.Timestamp(ds_700['time'].values)
    time_str = time_valid.strftime('%HZ')
    valid_time_moving = pd.Timestamp(ds_700['valid_time'].values)
    valid_time = valid_time_moving.strftime('%HZ %a %b %d %Y')

    dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat)

    uqvect_700, vqvect_700 = mpcalc.q_vector(u_700, v_700, temp_700, level_700, dx=dx, dy=dy)
    q_div_700 = 2 * mpcalc.divergence(uqvect_700, vqvect_700, dx=dx, dy=dy)

    ###### 600mb #######
    ds_600 = H.xarray(":(?:HGT|UGRD|VGRD|TMP):600 mb:", remove_grib=True)
    print(ds_600)

    gh_600 = ds_600['gh'] / 10
    u_600 = ds_600['u']
    v_600 = ds_600['v']
    temp_600 = ds_600['t']
    u_kts_600 = ds_600['u'] * 1.944
    v_kts_600 = ds_600['v'] * 1.944
    level_600 = 600 * units.hPa

    uqvect_600, vqvect_600 = mpcalc.q_vector(u_600, v_600, temp_600, level_600, dx=dx, dy=dy)
    q_div_600 = 2 * mpcalc.divergence(uqvect_600, vqvect_600, dx=dx, dy=dy)

    ###### 500mb #######
    ds_500 = H.xarray(":(?:HGT|UGRD|VGRD|TMP):500 mb:", remove_grib=True)
    print(ds_500)

    gh_500 = ds_500['gh'] / 10
    u_500 = ds_500['u']
    v_500 = ds_500['v']
    temp_500 = ds_500['t']
    u_kts_500 = ds_500['u'] * 1.944
    v_kts_500 = ds_500['v'] * 1.944
    level_500 = 500 * units.hPa

    uqvect_500, vqvect_500 = mpcalc.q_vector(u_500, v_500, temp_500, level_500, dx=dx, dy=dy)
    q_div_500 = 2 * mpcalc.divergence(uqvect_500, vqvect_500, dx=dx, dy=dy)
    gh_smooth_500 = gaussian_filter(gh_500, 2)

    ###### 400mb #######
    ds_400 = H.xarray(":(?:HGT|UGRD|VGRD|TMP):400 mb:", remove_grib=True)
    print(ds_400)

    gh_400 = ds_400['gh'] / 10
    u_400 = ds_400['u']
    v_400 = ds_400['v']
    temp_400 = ds_400['t']
    u_kts_400 = ds_400['u'] * 1.944
    v_kts_400 = ds_400['v'] * 1.944
    level_400 = 400 * units.hPa

    uqvect_400, vqvect_400 = mpcalc.q_vector(u_400, v_400, temp_400, level_400, dx=dx, dy=dy)
    q_div_400 = 2 * mpcalc.divergence(uqvect_400, vqvect_400, dx=dx, dy=dy)

    #### Averaging layer ######
    uqvect_mean = (uqvect_700 + uqvect_600 + uqvect_500 + uqvect_400) / 4
    vqvect_mean = (vqvect_700 + vqvect_600 + vqvect_500 + vqvect_400) / 4
    q_div_mean = 2 * mpcalc.divergence(uqvect_mean, vqvect_mean, dx=dx, dy=dy)
    q_div_smooth = gaussian_filter(q_div_mean * 1e18, 2)

    # plot settings s
    proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35, standard_parallels=(30, 60))
    fig = plt.figure(figsize=(20, 12))
    ax = plt.axes(projection=proj)
    # set bounds
    ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())

    # qvector stuff
    clevs_qdiv = list(range(-35, -4, 5)) + list(range(5, 36, 5))
    norm=mcolors.BoundaryNorm(clevs_qdiv,ncolors=256,extend='both')
    cf = ax.contourf(lon, lat, q_div_smooth, clevs_qdiv, cmap='bwr_r',norm=norm, extend='both', transform=ccrs.PlateCarree())
    line1 = ax.contour(lon, lat, q_div_smooth, colors='black', linewidths=1.5,
                       transform=ccrs.PlateCarree())
    ax.clabel(line1, inline=True, colors='black', fontsize=14, fmt='%d')

    cb = plt.colorbar(cf, ax=ax, orientation='horizontal', shrink=.4, pad=.02, aspect=25)
    cb.set_label('m^2kg^-1s^-1', fontsize=16)
    cb.ax.tick_params(labelsize=16)

    # GH contours
    line = ax.contour(lon, lat, gh_smooth_500, levels=list(range(500, 600, 6)), colors='black', linewidths=3,
                      transform=ccrs.PlateCarree())
    ax.clabel(line, inline=True, colors='black', fontsize=14, fmt='%d')

    # vectors
    wind_slice = (slice(None, None, 5), slice(None, None, 5))
    ax.quiver(lon[wind_slice].values, lat[wind_slice].values,
              uqvect_mean[wind_slice].values,
              vqvect_mean[wind_slice].values,
              pivot='mid', color='black',
              scale=5e-12, scale_units='inches',
              transform=ccrs.PlateCarree())

    # additional plot settings
    ax.coastlines('50m')
    ax.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax.add_feature(cfeature.STATES.with_scale('50m'))
    ax.set_title(f'Model: {time_str} {H.model.upper()} | F{H.fxx:03d}', fontsize=15, loc='left')
    ax.set_title(f'400-700mb Q-Vectors\n (m^2kg^-1s^-1)', fontsize=20, loc='center')
    ax.set_title(f'\nValid: {valid_time}', fontsize=15, loc='right')

    out_path = os.path.join(output_dir, f'{fxx:02d}.png')
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)

    return out_path, valid_time


#pull from model
target_time = (pd.Timestamp.now(tz='UTC') - pd.Timedelta(minutes=60)).tz_localize(None).floor('h')
max_hours=51

#determine if 21 frames or 51
run_number=target_time.hour

if run_number % 6==3:
    max_hours=52
else:
    max_hours=22

#create empty list for frames
frames=[]

for fxx in range(0,max_hours,1):
    try:
        H = Herbie(target_time, model='rap', product='awp236pgrb', fxx=fxx)
        out_path, valid_time = extract_build_data(H, fxx)
        frames.append({"fxx": fxx, "file":os.path.basename(out_path), "valid_time":valid_time})
        print(f"Saved frame f{fxx:02d}")
    except Exception as e:
        print(f"Skipping fxx=P{fxx}: {e}")
        continue

#write document dictionary so the website JS knows what frames exist
document={"run_time":target_time.strftime("%Y-%m-%d %H:%M"), "generated_at": datetime.utcnow().isoformat() + 'Z',
          'frames':frames}
with open(os.path.join(output_dir, "document.json"), 'w') as f:
    json.dump(document, f,indent=2)
print(f"Done {len(frames)} frames saved")

