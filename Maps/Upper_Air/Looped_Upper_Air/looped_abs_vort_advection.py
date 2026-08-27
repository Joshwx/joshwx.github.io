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

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..', '..'))
output_dir = os.path.join(REPO_ROOT, 'assets', 'maps', 'abs_vort_advection')
os.makedirs(output_dir, exist_ok=True)

def extract_build_data (ds,H, fxx):
    gh = ds['gh'] / 10
    u = ds['u']
    v = ds['v']
    lat = ds['latitude']
    lon = ds['longitude']
    time = ds['time']
    u_kts = ds['u'] * 1.944
    v_kts = ds['v'] * 1.944
    speed_500 = np.sqrt(u_kts ** 2 + v_kts ** 2)
    # (print(speed_500))
    print(time)
    # strip time
    time_valid = pd.Timestamp(ds['time'].values)
    time_str = time_valid.strftime('%HZ')
    valid_time_moving = pd.Timestamp(ds['valid_time'].values)
    valid_time = valid_time_moving.strftime('%HZ %a %b %d %Y')

    gh_smooth = gaussian_filter(gh, 1)
    print(lat, lon)

    # vorticity calc
    dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat)
    abs_vort = mpcalc.absolute_vorticity(u, v, dx=dx, dy=dy)
    abs_advection = mpcalc.advection(abs_vort, u, v, dx=dx, dy=dy)
    synoptic_scale = abs_advection * 1e9

    vort_smooth = gaussian_filter(synoptic_scale, 1)

    # plot settings s
    proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35,
                                 standard_parallels=(30, 60))
    fig = plt.figure(figsize=(20, 12))
    ax = plt.axes(projection=proj)
    ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())

    # GH contours
    line = ax.contour(lon, lat, gh_smooth[:, :], levels=list(range(500, 600, 6)), colors='black', linewidths=3,
                      transform=ccrs.PlateCarree())
    ax.clabel(line, inline=True, colors='black', fontsize=12, fmt='%d')

    # vort contour
    cf = ax.contourf(lon, lat, vort_smooth, levels=list(range(-52, 52, 2)), extend='both',
                     cmap='bwr', transform=ccrs.PlateCarree())
    cb = plt.colorbar(cf, ax=ax, orientation='horizontal', shrink=.4, pad=.02, aspect=25)
    cb.set_label('10^-9 s^-2', fontsize=16)
    cb.ax.tick_params(labelsize=16)

    cf_red = ax.contour(lon, lat, vort_smooth, levels=(range(-50, -5, 5)), colors='blue', linewidths=1.5,
                        linestyles='dashed', transform=ccrs.PlateCarree())
    ax.clabel(cf_red, inline=True, colors='blue', fontsize=15, fmt='%d')
    ax.clabel(cf_red, inline=True, colors='red', fontsize=15, fmt='%d')

    cf_blue = ax.contour(lon, lat, vort_smooth, levels=(range(5, 50, 5)), colors='red', linewidths=1.5,
                         linestyles='dashed', transform=ccrs.PlateCarree())
    ax.clabel(cf_blue, inline=True, colors='red', fontsize=15, fmt='%d')
    # wind stuff
    wnd_speed = np.arange(40, 160, 20)
    ax.barbs(lon, lat, u_kts, v_kts,
             length=8, regrid_shape=12, pivot='middle', transform=ccrs.PlateCarree())

    # additional plot settings
    ax.coastlines('50m')
    ax.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax.add_feature(cfeature.STATES.with_scale('50m'))
    ax.set_title(f'Model: {time_str} {H.model.upper()} | F{H.fxx:03d}', fontsize=15, loc='left')
    ax.set_title(f'500mb Height (dam),\n Absolute Vorticity Advection (10^-9 s^-2) ', fontsize=20, loc='center')
    ax.set_title(f'\nValid: {valid_time}', fontsize=15, loc='right')

    out_path = os.path.join(output_dir, f'{fxx:02d}.png')
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)

    return  out_path, valid_time

#pull from model
target_time = (pd.Timestamp.now(tz='UTC') - pd.Timedelta(minutes=60)).tz_localize(None).floor('h')
print('looking for run', target_time)
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
        ds = H.xarray(":(?:HGT|UGRD|VGRD):500 mb:", remove_grib=True)
        out_path, valid_time = extract_build_data(ds, H, fxx)
        frames.append({"fxx": fxx, "file": os.path.basename(out_path), "valid_time": valid_time})
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
