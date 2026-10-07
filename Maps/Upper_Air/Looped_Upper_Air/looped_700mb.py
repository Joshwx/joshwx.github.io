from matplotlib.colors import LinearSegmentedColormap, ListedColormap
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

script_dir=os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.abspath(os.path.join(script_dir, '..', '..', '..'))
output_dir = os.path.join(repo_root, 'assets', 'maps', 'upper_level_700mb')
os.makedirs(output_dir, exist_ok=True)

cmap_path=os.path.join(repo_root,'assets','cmap','MPL_Greens.rgb')


def load_ncl_rgb(filepath, smooth=True):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    rgb_vals = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith(';') or line.startswith('#') or line.startswith('!') or '=' in line:
            continue
        parts = line.split()
        if len(parts) >= 3:
            rgb_vals.append([float(x) for x in parts[:3]])

    rgb = np.array(rgb_vals)

    # normalize only if values are in 0-255 range
    if rgb.max() > 1.0:
        rgb = rgb / 255.0

    if smooth:
        return LinearSegmentedColormap.from_list('custom_cmap', rgb, N=256)
    else:
        return ListedColormap(rgb)

cmap_greens=load_ncl_rgb(cmap_path)

#function that builds the images with ds and H as arguments
def extract_build_data (ds,ds2,H, fxx):
    gh = ds['gh']
    u = ds['u']
    v = ds['v']
    lat = ds['latitude']
    lon = ds['longitude']
    temp = ds['t']
    time = ds['time']
    temp_c = ds['t'] - 273.15
    u_kts = ds['u'] * 1.944
    v_kts = ds['v'] * 1.944
    # extract rhs
    rh_700 = ds2['r'].sel(isobaricInhPa=700)
    rh_650 = ds2['r'].sel(isobaricInhPa=650)
    rh_600 = ds2['r'].sel(isobaricInhPa=600)
    rh_550 = ds2['r'].sel(isobaricInhPa=550)
    rh_500 = ds2['r'].sel(isobaricInhPa=500)

    # mean
    # Stack along a new axis, then take the mean
    mean_rh = np.mean(np.stack([rh_700, rh_650, rh_600, rh_550, rh_500], axis=0), axis=0)

    # strip time
    time_valid = pd.Timestamp(ds['time'].values)
    time_str = time_valid.strftime('%HZ')
    valid_time = time_valid.strftime('%HZ %a %b %d %Y')

    # get the whole wind vector
    speed_700 = np.sqrt(u_kts ** 2 + v_kts ** 2)
    # smooth gh
    gh_smooth = gaussian_filter(gh, 2)
    temp_smooth = gaussian_filter(temp_c, 2)
    smooth_rh = gaussian_filter(mean_rh, 2)



    # plot settings s
    proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35,
                                 standard_parallels=(30, 60))
    fig = plt.figure(figsize=(20, 12))
    ax = plt.axes(projection=proj)

    # set bounds
    ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())
    # GH contours
    line = ax.contour(lon, lat, gh_smooth, levels=list(range(2340, 3780, 30)), colors='black', linewidths=3,
                      transform=ccrs.PlateCarree())
    # rh fill
    cf_rh = ax.contourf(lon, lat, smooth_rh, levels=np.arange(70, 105, 5), cmap=cmap_greens,
                        transform=ccrs.PlateCarree())

    # wind stuff
    ax.barbs(lon, lat, u_kts, v_kts,
             length=8, regrid_shape=12, pivot='middle', transform=ccrs.PlateCarree())
    ax.clabel(line, inline=True, colors='black', fontsize=12, fmt='%d')

    # temp contours
    cf_red = ax.contour(lon, lat, temp_smooth, levels=(range(5, 40, 5)), colors='red', linewidths=1.5,
                        linestyles='dashed', transform=ccrs.PlateCarree())
    ax.clabel(cf_red, inline=True, colors='red', fontsize=15, fmt='%d')
    cf_red_bold = ax.contour(lon, lat, temp_smooth, levels=(range(10, 30, 5)), colors='red', linewidths=2,
                             linestyles='dashed', transform=ccrs.PlateCarree())
    ax.clabel(cf_red_bold, inline=True, colors='red', fontsize=15, fmt='%d')
    cf_blue = ax.contour(lon, lat, temp_smooth, levels=(range(-40, 0, 5)), colors='blue', linewidths=1.5,
                         linestyles='dashed', transform=ccrs.PlateCarree())
    ax.clabel(cf_blue, inline=True, colors='blue', fontsize=15, fmt='%d')

    cf_freez = ax.contour(lon, lat, temp_smooth, levels=[0], colors='blue', linewidths=2,
                          linestyles='solid', transform=ccrs.PlateCarree())
    ax.clabel(cf_freez, inline=True, colors='blue', fontsize=15, fmt='%d')

    cb = plt.colorbar(cf_rh, ax=ax, orientation='horizontal', shrink=.4, pad=.02, aspect=25)
    cb.set_label('Relative Humidity (%)', fontsize=16)
    cb.ax.tick_params(labelsize=16)

    # additional plot settings
    ax.coastlines('50m')
    ax.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax.add_feature(cfeature.STATES.with_scale('50m'))
    ax.set_title(f'Model: {time_str} {H.model.upper()} | F{H.fxx:03d}', fontsize=15, loc='left')
    ax.set_title(f'700mb Height (dam), Wind (kts), Temp (C),\n 700-500mb Mean RH >= 70%', fontsize=20,
                 loc='center')
    ax.set_title(f'\nValid: {valid_time}', fontsize=15, loc='right')

    out_path = os.path.join(output_dir, f'{fxx:02d}.png')
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)

    return out_path, valid_time


#function that subsets the GFS data so you dont plot all the points in the global run
def subset(ds, lon_min=-130, lon_max=-60,lat_min=15, lat_max=58):
    lon_min_convert= lon_min%360
    lon_max_convert= lon_max%360
    return ds.sel(latitude=slice(lat_max, lat_min),longitude=slice(lon_min_convert, lon_max_convert),
    )

#pull from model
target_time = (pd.Timestamp.now(tz='UTC') - pd.Timedelta(hours=4)).tz_localize(None).floor('6h')
print('looking for run', target_time)

fxx_list = list(range(0, 25, 1)) + list(range(27, 121, 3))

#create empty list for frames
frames=[]
for fxx in fxx_list:
    try:
        H = Herbie(target_time, model="gfs", product="pgrb2.0p25", fxx=fxx)
        ds = H.xarray(":(?:HGT|UGRD|VGRD|TMP|RH):700 mb:", remove_grib=True)
        ds2 = H.xarray(":RH:(700|650|600|550|500) mb")
        ds = subset(ds).load()
        ds2 = subset(ds2).load()
        out_path, valid_time = extract_build_data(ds,ds2, H, fxx)
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
