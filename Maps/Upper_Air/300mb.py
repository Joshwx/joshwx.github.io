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
H= Herbie("2026-08-22 01:00", model='rap', product='awp236pgrb', fxx=0)

print(H.inventory(":300 mb"))

ds=H.xarray(":(?:HGT|UGRD|VGRD):300 mb:", remove_grib=False)
print(ds)

gh=ds['gh']/10
u=ds['u']
v=ds['v']
lat=ds['latitude']
lon=ds['longitude']
time=ds['time']
u_kts=ds['u']*1.944
v_kts=ds['v']*1.944
speed_300=np.sqrt(u_kts**2+v_kts**2)
(print(speed_300))

# strip time
time_valid = pd.Timestamp(ds['time'].values)
time_str = time_valid.strftime('%HZ')
valid_time_moving = pd.Timestamp(ds['valid_time'].values)
valid_time = valid_time_moving.strftime('%HZ %a %b %d %Y')



barb_skip = (slice(None, None, 8), slice(None, None, 8))

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
line=ax.contour(lon, lat,gh_smooth[:,:],levels=list(range(798, 996, 6)),colors='black',linewidths=3,
                transform=ccrs.PlateCarree(), cbar=True, cbar_kwargs={'orientation':'vertical'})
#wind stuff
wnd_speed= np.arange(60, 160, 20)
cf=ax.contourf(lon,lat, speed_300, wnd_speed, cmap='BuPu', transform=ccrs.PlateCarree())
ax.barbs(lon, lat, u_kts,v_kts,
         length=8, regrid_shape=15, pivot='middle', transform=ccrs.PlateCarree())
ax.clabel(line,inline=True,colors='black',fontsize=12, fmt='%d')

#additional plot settings
ax.coastlines('50m')
ax.add_feature(cfeature.BORDERS.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.set_title(f'Model: {time_str} {H.model.upper()} | F{H.fxx:03d}', fontsize=22, loc='left')
ax.set_title(f'300mb Height (dam), Winds (kts) ', fontsize=22, loc='center')
ax.set_title(f'\nValid: {valid_time}', fontsize=22, loc='right')

plt.tight_layout()
plt.show()