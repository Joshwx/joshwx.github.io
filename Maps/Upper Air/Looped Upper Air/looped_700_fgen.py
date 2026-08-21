import sounderpy as spy
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import metpy
import metpy.calc as mpcalc
from matplotlib.lines import lineStyles
from metpy.units import units
from herbie import Herbie
from herbie.toolbox import EasyMap, pc
from herbie import paint
import cartopy
from cartopy import crs as ccrs, feature as cfeature
import scipy
from scipy.ndimage import gaussian_filter
import pandas as pd
import numpy as np


#function that builds the images with ds and H as arguments
def extract_build_data (ds,H):
    gh = ds['gh']
    u = ds['u']
    v = ds['v']
    lat = ds['latitude']
    lon = ds['longitude']
    temp = ds['t']
    time = ds['time']

    u_kts = ds['u'] * 1.944
    v_kts = ds['v'] * 1.944
    temp_850 = temp.metpy.convert_units('degC')
    level = 700 * units.hPa

    # fgen calculation
    theta_700 = mpcalc.potential_temperature(level, temp_850)
    print(theta_700)
    dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat)
    fronto_850 = mpcalc.frontogenesis(theta_700, u, v, dx, dy)
    convert_to_per_100km_3hr = 1000 * 100 * 3600 * 3

    # strip time
    time_valid = pd.Timestamp(ds['time'].values)
    time_str = time_valid.strftime('%HZ')
    valid_time = time_valid.strftime('%HZ %a %b %d %Y')

    # get the whole wind vector
    speed_700 = np.sqrt(u_kts ** 2 + v_kts ** 2)
    # mask settings
    mask = speed_700 < 5
    u_masked = np.ma.masked_array(u_kts, mask)
    v_masked = np.ma.masked_array(v_kts, mask)

    # smoothed stuff
    gh_smooth = gaussian_filter(gh, 2)
    fgen_smoothed = gaussian_filter(fronto_850, 2)

    fgen_masked = np.ma.masked_where(np.abs(fgen_smoothed) < 0.5, fgen_smoothed)

    print(np.nanmin(fgen_smoothed), np.nanmax(fgen_smoothed))

    # plot settings s
    proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35, standard_parallels=(30, 60))
    fig = plt.figure(figsize=(20, 12))
    ax = plt.axes(projection=proj)
    # set bounds
    ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())

    # fgen stuff
    cf = ax.contourf(lon, lat, fgen_smoothed * convert_to_per_100km_3hr, np.arange(-8, 8.5, .5), cmap='bwr',
                     extend='both', transform=ccrs.PlateCarree())
    line1 = ax.contour(lon, lat, fgen_smoothed * convert_to_per_100km_3hr, colors='black', linewidths=.5,
                       transform=ccrs.PlateCarree())
    ax.clabel(line1, inline=True, colors='black', fontsize=12, fmt='%d')

    # GH contours
    line = ax.contour(lon, lat, gh_smooth, levels=list(range(2340, 3780, 30)), colors='black', linewidths=3,
                      transform=ccrs.PlateCarree())
    ax.clabel(line, inline=True, colors='black', fontsize=12, fmt='%d')

    # wind stuff
    wnd_speed = np.arange(25, 85, 5)
    # wind=ax.contourf(lon,lat,speed_850, wnd_speed, cmap='turbo', transform=ccrs.PlateCarree())
    ax.barbs(lon, lat, u_masked, v_masked,
             length=8, regrid_shape=15, fill_empty=False, pivot='middle', transform=ccrs.PlateCarree())
    ax.clabel(line, inline=True, colors='black', fontsize=12, fmt='%d')

    # additional plot settings
    ax.coastlines('50m')
    ax.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax.add_feature(cfeature.STATES.with_scale('50m'))
    ax.set_title(f'Model: {time_str} {H.model.upper()} | F{H.fxx:03d}', fontsize=22, loc='left')
    ax.set_title(f'700mb Frontogenisis (K/100km/3hr)', fontsize=26, loc='center')
    ax.set_title(f'\nValid: {valid_time}', fontsize=22, loc='right')

    plt.tight_layout()
    plt.show()
    plt.close()

    return fig, ax


for fxx in range(0,51,1):
    H = Herbie("2026-01-24 03:00", model='rap', product='awp236pgrb', fxx=fxx)
    ds = H.xarray(":(?:HGT|UGRD|VGRD|TMP|RH):700 mb:", remove_grib=False)
    fig, ax = extract_build_data(ds,H)