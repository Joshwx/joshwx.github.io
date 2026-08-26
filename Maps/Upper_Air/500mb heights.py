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

#pull from model
H= Herbie("2026-01-24 00:00", model='rap', product='awp236pgrb', fxx=5)

# print(H.inventory(":500 mb"))

ds=H.xarray(":(?:HGT|UGRD|VGRD):500 mb:", remove_grib=False)
print(ds)

gh=ds['gh']/10
u=ds['u']
v=ds['v']
lat=ds['latitude']
lon=ds['longitude']
time=ds['time']
u_kts=ds['u']*1.944
v_kts=ds['v']*1.944
speed_500=np.sqrt(u_kts**2+v_kts**2)
# (print(speed_500))
print(time)
# strip time
time_valid = pd.Timestamp(ds['time'].values)
time_str = time_valid.strftime('%HZ')
valid_time_moving = pd.Timestamp(ds['valid_time'].values)
valid_time = valid_time_moving.strftime('%HZ %a %b %d %Y')

gh_smooth=gaussian_filter(gh,1)
print(lat,lon)


#plot settings s
proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35,
                              standard_parallels=(30, 60))
fig=plt.figure(figsize=(20,12))
ax=plt.axes(projection=proj)

#set bounds
ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())
#GH contours
line=ax.contour(lon, lat,gh_smooth[:,:],levels=list(range(500, 600, 6)),colors='black',linewidths=3,
                transform=ccrs.PlateCarree())
#wind stuff
wnd_speed= np.arange(40, 160, 20)
cf=ax.contourf(lon,lat, speed_500, wnd_speed, cmap='BuPu', transform=ccrs.PlateCarree())
ax.barbs(lon, lat, u_kts,v_kts,
         length=8, regrid_shape=15, pivot='middle', transform=ccrs.PlateCarree())
ax.clabel(line,inline=True,colors='black',fontsize=12, fmt='%d')
cb=plt.colorbar(cf,ax=ax, orientation='horizontal', shrink=.4,pad=.02, aspect=25)
cb.set_label('Wind Speed (kts)', fontsize=16)
cb.ax.tick_params(labelsize=16)


#additional plot settings
ax.coastlines('50m')
ax.add_feature(cfeature.BORDERS.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.set_title(f'Model: {time_str} {H.model.upper()} | F{H.fxx:03d}', fontsize=22, loc='left')
ax.set_title(f'500mb Height (dam), Winds (kts) ', fontsize=26, loc='center')
ax.set_title(f'\nValid: {valid_time}', fontsize=22, loc='right')


plt.tight_layout()
plt.show()


#vorticity calc
dx,dy=mpcalc.lat_lon_grid_deltas(lon, lat)
vort=mpcalc.absolute_vorticity(u,v,dx=dx, dy=dy)
synoptic_scale=vort*1e5


#plot settings s
proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35,
                              standard_parallels=(30, 60))
fig=plt.figure(figsize=(20,12))
ax=plt.axes(projection=proj)
ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())
#vort contour
cf=ax.contourf(lon,lat,synoptic_scale, levels=list(range(16,44,4)),
               cmap='plasma', transform=ccrs.PlateCarree())
cb=plt.colorbar(cf,ax=ax, orientation='horizontal', shrink=.4,pad=.02, aspect=25)
cb.set_label('10^-5 s^-1', fontsize=16)
cb.ax.tick_params(labelsize=16)
#gh heights
#GH contours
line=ax.contour(lon, lat,gh_smooth[:,:],levels=list(range(500, 600, 6)),colors='black',linewidths=3,
                transform=ccrs.PlateCarree())
ax.clabel(line,inline=True,colors='black',fontsize=12, fmt='%d')
#additional plot settings
ax.coastlines('50m')
ax.add_feature(cfeature.BORDERS.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.set_title(f'Model: {time_str} {H.model.upper()} | F{H.fxx:03d}', fontsize=22, loc='left')
ax.set_title(f'500mb Height (dam),\n Absolute Vorticity (10^-5 s^-1) ', fontsize=26, loc='center')
ax.set_title(f'\nValid: {valid_time}', fontsize=22, loc='right')

plt.tight_layout()
plt.show()
