import sounderpy as spy
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

#function that builds the images with ds and H as arguments
def extract_build_data (ds,H):
    gh = ds['gh'] / 10
    u = ds['u']
    v = ds['v']
    lat = ds['latitude']
    lon = ds['longitude']
    time = ds['time']
    # strip time
    time_valid = pd.Timestamp(ds['time'].values)
    time_str = time_valid.strftime('%HZ')
    valid_time_moving = pd.Timestamp(ds['valid_time'].values)
    valid_time = valid_time_moving.strftime('%HZ %a %b %d %Y')

    # divergence calc
    dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat)

    divergence = mpcalc.divergence(u, v, dx=dx, dy=dy)
    synoptic_scale = divergence * 1e5
    smoothed_divergence = gaussian_filter(synoptic_scale, 2)

    u_kts = ds['u'] * 1.944
    v_kts = ds['v'] * 1.944
    speed_300 = np.sqrt(u_kts ** 2 + v_kts ** 2)
    (print(speed_300))

    # plot settings
    proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35,
                                 standard_parallels=(30, 60))
    fig = plt.figure(figsize=(20, 12))
    ax = plt.axes(projection=proj)

    # set bounds
    ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())
    # GH contours
    line = ax.contour(lon, lat, gh[:, :], levels=list(range(1020, 1170, 60)), colors='black', linewidths=3,
                      transform=ccrs.PlateCarree())
    # divergence
    # div=ax.contour(lon, lat, synoptic_scale, levels=list(range(2,20,2)), colors='red',fontsize=12, linewidths=2, fmt='%d',
    # transform=ccrs.PlateCarree())
    # wind stuff
    wnd_speed = np.arange(60, 180, 20)
    cf = ax.contourf(lon, lat, speed_300, wnd_speed, cmap='BuPu', transform=ccrs.PlateCarree())
    ax.barbs(lon, lat, u_kts, v_kts,
             length=8, regrid_shape=15, pivot='middle', transform=ccrs.PlateCarree())
    ax.clabel(line, inline=True, colors='black', fontsize=12, fmt='%d')

    # additional plot settings
    ax.coastlines('50m')
    ax.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax.add_feature(cfeature.STATES.with_scale('50m'))
    ax.set_title(f'Model: {time_str} {H.model.upper()} | F{H.fxx:03d}', fontsize=22, loc='left')
    ax.set_title(f'250mb Height (dam), Wind (kts)', fontsize=26, loc='center')
    ax.set_title(f'\nValid: {valid_time}', fontsize=22, loc='right')

    plt.tight_layout()
    plt.show()
    plt.close()
    return fig, ax

#for loop that is bounded to the RAPs hours to plot every hour
for fxx in range(0,52,1):
    #H and ds are the two products that are going into every figure so loop through them between the forecast hours
    #they also are function calls
    H = Herbie("2026-01-24 00:00", model='rap', product='awp236pgrb', fxx=fxx)
    ds = H.xarray(":(?:HGT|UGRD|VGRD):250 mb:", remove_grib=False)
    #you are calling the return products and using the func to use ds and H to plot data
    fig, ax = extract_build_data(ds,H)

