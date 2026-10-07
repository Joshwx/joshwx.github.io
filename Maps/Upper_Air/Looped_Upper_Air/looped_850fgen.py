import os

import matplotlib
matplotlib.use("Agg")

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

output_dir = 'assets/maps/upper_level_850_fgen'
os.makedirs(output_dir, exist_ok=True)

#function that builds the images with ds and H as arguments
def extract_build_data (ds,H):
    gh = ds['gh'] / 10
    u = ds['u']
    v = ds['v']
    lat = ds['latitude']
    lon = ds['longitude']
    temp = ds['t']
    time = ds['time']

    u_kts = ds['u'] * 1.944
    v_kts = ds['v'] * 1.944
    temp_850 = temp.metpy.convert_units('degC')
    level = 850 * units.hPa

    # fgen calculation
    theta_850 = mpcalc.potential_temperature(level, temp_850)
    print(theta_850)
    dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat)
    fronto_850 = mpcalc.frontogenesis(theta_850, u, v, dx, dy)
    convert_to_per_100km_3hr = 1000 * 100 * 3600 * 3

    converted_fgen = fronto_850 * convert_to_per_100km_3hr
    print(converted_fgen)
    # strip time
    time_valid = pd.Timestamp(ds['time'].values)
    time_str = time_valid.strftime('%HZ')
    valid_time_moving = pd.Timestamp(ds['valid_time'].values)
    valid_time = valid_time_moving.strftime('%HZ %a %b %d %Y')

    speed_850 = np.sqrt(u_kts ** 2 + v_kts ** 2)

    # mask settings
    mask = speed_850 < 5
    u_masked = np.ma.masked_array(u_kts, mask)
    v_masked = np.ma.masked_array(v_kts, mask)

    # smoothed stuff
    gh_smooth = gaussian_filter(gh, 2)
    fgen_smoothed = gaussian_filter(converted_fgen, 4)

    fgen_masked = np.ma.masked_where(np.abs(fgen_smoothed) <= .99, fgen_smoothed)

    print(np.nanmin(fgen_smoothed), np.nanmax(fgen_smoothed))

    # print(converted_fgen)
    # plot settings s
    proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35, standard_parallels=(30, 60))
    fig = plt.figure(figsize=(20, 12))
    ax = plt.axes(projection=proj)
    # set bounds
    ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())

    # fgen stuff
    cf = ax.contourf(lon, lat, fgen_masked , np.arange(-8, 8.5, 1), cmap='bwr',
                     extend='both', transform=ccrs.PlateCarree())
    line1 = ax.contour(lon, lat, fgen_masked , colors='black', linewidths=.5,
                       transform=ccrs.PlateCarree())
    ax.clabel(line1, inline=True, colors='black', fontsize=12, fmt='%d')

    # GH contours
    line = ax.contour(lon, lat, gh_smooth, levels=list(range(90, 300, 3)), colors='black', linewidths=3,
                      transform=ccrs.PlateCarree())

    # wind stuff
    wnd_speed = np.arange(25, 85, 5)
    speed_850 = np.sqrt(u_kts ** 2 + v_kts ** 2)
    # wind=ax.contourf(lon,lat,speed_850, wnd_speed, cmap='turbo', transform=ccrs.PlateCarree())
    ax.barbs(lon, lat, u_masked, v_masked,
             length=8, regrid_shape=15, fill_empty=False, pivot='middle', transform=ccrs.PlateCarree())
    ax.clabel(line, inline=True, colors='black', fontsize=12, fmt='%d')

    # additional plot settings
    ax.coastlines('50m')
    ax.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax.add_feature(cfeature.STATES.with_scale('50m'))
    ax.set_title(f'Model: {time_str} {H.model.upper()} | F{H.fxx:03d}', fontsize=22, loc='left')
    ax.set_title(f'850mb Frontogenisis (K/100km/3hr)', fontsize=26, loc='center')
    ax.set_title(f'\nValid: {valid_time}', fontsize=22, loc='right')

    plt.tight_layout()
    plt.show()
    plt.close()
    return fig,ax

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
        ds = H.xarray(":(?:HGT|UGRD|VGRD|TMP|RH):850 mb:", remove_grib=False)
        ds = subset(ds).load()
        out_path, valid_time = extract_build_data(ds, H, fxx)
        frames.append({"fxx": fxx, "file": os.path.basename(out_path), "valid_time": valid_time})
    except Exception as e:
        print(f"Skipping fxx=P{fxx}: {e}")
        continue

#write document dictionary so the website JS knows what frames exist
document={"run_time":target_time.strftime("%Y-%m-%d %H:%M"), "generated_at": datetime.utcnow().isoformat() + 'Z',
          'frames':frames}
with open(os.path.join(output_dir, "document.json"), 'w') as f:
    json.dump(document, f,indent=2)
print(f"Done {len(frames)} frames saved")