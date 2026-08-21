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

#pull from model
H= Herbie("2026-08-12 13:00", model='rap', product='awp236pgrb', fxx=0)

#print(H.inventory(":850 mb"))

ds=H.xarray(":(?:HGT|UGRD|VGRD|TMP|RH):850 mb:", remove_grib=False)
#print(ds)

#create a function with H and ds as parameters since thats wehre you are pulling from
#a loop can just call this and a loop can loop over hours
def get_and_plot(H, ds):
    gh = ds['gh']
    u = ds['u']
    v = ds['v']
    lat = ds['latitude']
    lon = ds['longitude']
    temp = ds['t']
    time = ds['time']
    rh = ds['r']
    temp_c = ds['t'] - 273.15
    u_kts = ds['u'] * 1.944
    v_kts = ds['v'] * 1.944
    press = 850 * units.hPa

    # product of the wind speed m/s
    speed_850 = np.sqrt(u ** 2 + v ** 2)

    # strip time
    time_valid = pd.Timestamp(ds['time'].values)
    time_str = time_valid.strftime('%HZ')
    valid_time = time_valid.strftime('%HZ %a %b %d %Y')

    # calculate mixing ratio
    mr = metpy.calc.mixing_ratio_from_relative_humidity(press, temp.values * units.kelvin,
                                                        (rh.values / 100) * units.dimensionless,
                                                        phase='liquid')
    # calculate theta-e
    dp = metpy.calc.dewpoint_from_relative_humidity(temp.values * units.kelvin, (rh.values / 100) * units.dimensionless)
    te = metpy.calc.equivalent_potential_temperature(press, temp.values * units.kelvin, dp)
    # moisture transport (multiply wind speed (m/s) and mixing ratio (g/g) and scale by 100)
    mt = (speed_850 * mr.magnitude) * 100

    # smooth
    mt_smooth = gaussian_filter(mt, 2)
    gh_smooth = gaussian_filter(gh, 2)
    te_smooth = gaussian_filter(te.magnitude, 5)
    # plot settings s
    proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35, standard_parallels=(30, 60))
    fig = plt.figure(figsize=(20, 12))
    ax = plt.axes(projection=proj)

    # set bounds
    ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())
    # moisture transport stuff
    cf = ax.contourf(lon, lat, mt_smooth, levels=list(np.arange(10, 40, 5)), cmap='RdYlGn_r',
                     transform=ccrs.PlateCarree())
    line1 = ax.contour(lon, lat, mt_smooth, levels=list(np.arange(10, 40, 5)), colors='black', linewidths=1,
                       transform=ccrs.PlateCarree())
    ax.clabel(line1, inline=True, colors='black', fontsize=12, fmt='%d')
    # transport vectors
    ax.quiver(lon.values[::10, ::10], lat.values[::10, ::10], u.values[::10, ::10], v.values[::10, ::10],
              transform=ccrs.PlateCarree(), scale=500, color='maroon')
    # GH contours
    line2 = ax.contour(lon, lat, gh_smooth, levels=list(range(1008, 1800, 24)), colors='black', linewidths=3,
                       transform=ccrs.PlateCarree())
    ax.clabel(line2, inline=True, colors='black', fontsize=12, fmt='%d')

    # theta-e contours
    line2 = ax.contour(lon, lat, te_smooth, levels=list(range(250, 350, 5)), colors='green', linewidths=1,
                       linestyles='dashed',
                       transform=ccrs.PlateCarree())
    ax.clabel(line2, inline=True, colors='black', fontsize=12, fmt='%d')

    # cbar=plt.colorbar(cf, ax=ax, orientation='horizontal', shrink=0.5,)

    # additional plot settings
    ax.coastlines('50m')
    ax.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax.add_feature(cfeature.STATES.with_scale('50m'))
    ax.set_title(f'Model: {time_str} {H.model.upper()} | F{H.fxx:03d}', fontsize=26, loc='left')
    ax.set_title(f'850mb Moisture Transport (Shaded and Vector), 850mb Height (m),\n Theta-e (dashed)',
                 fontsize=18, loc='center')
    ax.set_title(f'\nValid: {valid_time}', fontsize=22, loc='right')

    plt.tight_layout()
    plt.show()
    return fig,ax
#loop thru forecast hours from Herbie, slap variables to ds, call function to extract
#variables and plot map, repeat until done
for fxx in range(0,24,1):
    H = Herbie("2026-08-12 13:00", model='rap', product='awp236pgrb', fxx=5)
    ds = H.xarray(":(?:HGT|UGRD|VGRD|TMP|RH):850 mb:", remove_grib=False)
    fig, ax=get_and_plot(H,ds)
    # plt.savefig(f"850_moisture_transport_f{fxx:03d}.png", dpi=300)
    plt.close(fig)

