from enum import verify

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



H = Herbie("2026-10-17 06:00", model="gfs", product="pgrb2.0p25", overwrite=True)

ds = H.xarray(r":(?:TMP|UGRD|VGRD|RH):(?:[1-9]\d{2}|1000) mb:", remove_grib=False)
ds = ds.sel(latitude=slice(58, 15), longitude=slice(230, 300)).load()

# print(ds.dims, ds['t'].dims)


#extract varaibles
lat, lon=ds.latitude, ds.longitude
pressure=ds.isobaricInhPa.values *units.hPa
time = ds['time']
relh=ds['r']
temps=ds['t']
u=ds['u']
v=ds['v']
u_kts = ds['u'] * 1.944
v_kts = ds['v'] * 1.944

print('yeet')

isentprs, isentrh, isentu, isentv = mpcalc.isentropic_interpolation(
    np.array([300]) * units.kelvin, pressure, temps, relh, u, v, vertical_dim=0)

lon2d, lat2d = np.meshgrid(lon.values, lat.values)
dx, dy = mpcalc.lat_lon_grid_deltas(lon2d, lat2d)
adv = mpcalc.advection(isentprs[0], u=isentu[0], v=isentv[0], dx=dx, dy=dy)

# strip time
time_valid = pd.Timestamp(ds['time'].values)
time_str = time_valid.strftime('%HZ')
valid_time_moving = pd.Timestamp(ds['valid_time'].values)
valid_time = valid_time_moving.strftime('%HZ %a %b %d %Y')

#calculate all the variables on the surface
isentlevs = np.array([290]) * units.kelvin
# print(isentlevs)
isentropic_analysis=mpcalc.isentropic_interpolation(isentlevs,pressure,temps,relh,
                                                    u, v, vertical_dim=0)
isentprs, isentrh, isentu, isentv = isentropic_analysis




# print(isentprs.shape)

level=0

#calc isentropic wind level
isent_u_kts = isentu[level] * 1.944
isent_v_kts = isentv[level] * 1.944

#lat lon 2d
lon2d,lat2d=np.meshgrid(lon.values,lat.values)

#calculate pressure advection
dx,dy=mpcalc.lat_lon_grid_deltas(lon2d,lat2d)
# print(lon2d.shape, dx.shape, dy.shape)
# pressure_advection = mpcalc.advection(isentprs[level], u=isentu[level], v=isentv[level], dx=dx, dy=dy)
# adv_hr=pressure_advection.to('hPa/s').magnitude*3600

# print(adv_hr)


#smoother
# isentprs_advection_smooth = gaussian_filter(adv_hr, sigma=[2, 2])
isentprs_smooth=gaussian_filter(isentprs.m, sigma=[0,2,2])
isentrh_smooth = gaussian_filter(isentrh.m, sigma=[0, 2, 2])

# plot settings s
proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35,
                             standard_parallels=(30, 60))
fig = plt.figure(figsize=(20, 12))
ax = plt.axes(projection=proj)

#isentropic stuff
cs=ax.contour(lon.values,lat.values,isentprs_smooth[level],np.arange(0,1000,25),transform=ccrs.PlateCarree(),colors='k')
plt.clabel(cs, fmt='%i')

cf= ax.contourf(lon.values,lat.values,isentrh_smooth[level],np.arange(70,101,5), transform=ccrs.PlateCarree(),
                cmap='summer_r')
cb=plt.colorbar(cf, ax=ax, orientation='horizontal',shrink=.4, pad=.02, aspect=25)
cb.set_label('Relative Humidity (%)', fontsize=16)
cb.ax.tick_params(labelsize=16)

#wind barbs
ax.barbs(lon.values, lat.values, isent_u_kts.m, isent_v_kts.m,
         length=8, regrid_shape=17, fill_empty=False, pivot='middle',
         transform=ccrs.PlateCarree())

# set bounds
ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())

# additional plot settings
ax.coastlines('50m')
ax.add_feature(cfeature.BORDERS.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.set_title(f'Model: {time_str} {H.model.upper()} | F{H.fxx:03d}', fontsize=15, loc='left')
ax.set_title(f'290K Pressure (mb), Winds (kts),\n RH (%) ', fontsize=20, loc='center')
ax.set_title(f'\nValid: {valid_time}', fontsize=15, loc='right')

plt.tight_layout()
plt.show()