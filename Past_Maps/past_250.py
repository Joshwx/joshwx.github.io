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
from datetime import datetime, timedelta, timezone
#initalize past timeframe
target_time = (pd.Timestamp.now(tz='UTC') - pd.Timedelta(hours=24)).tz_localize(None).floor('h')
#pull from model
H= Herbie(target_time, model='rap', product='awp236pgrb')

print(H.inventory())

ds=H.xarray(":(?:HGT|UGRD|VGRD):250 mb:", remove_grib=False)
print(ds)

gh=ds['gh']/10
u=ds['u']
v=ds['v']
lat=ds['latitude']
lon=ds['longitude']
time=ds['time']
#strip time
time_valid=pd.Timestamp(ds['time'].values)
time_str=time_valid.strftime('%Y-%m-%d %H:%M UTC')
fhr=int(ds.step.values/np.timedelta64(1,'h'))

#divergence calc
dx,dy=mpcalc.lat_lon_grid_deltas(lon, lat)

divergence=mpcalc.divergence(u,v, dx=dx, dy=dy)
synoptic_scale=divergence*1e5
smoothed_divergence=gaussian_filter(synoptic_scale,2)


u_kts=ds['u']*1.944
v_kts=ds['v']*1.944
speed_300=np.sqrt(u_kts**2+v_kts**2)
(print(speed_300))

#plot settings
proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35,
                              standard_parallels=(30, 60))
fig=plt.figure(figsize=(20,12))
ax=plt.axes(projection=proj)

#set bounds
ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())
#GH contours
line=ax.contour(lon, lat,gh[:,:],levels=list(range(796, 1148, 6)),colors='black',linewidths=3,
                transform=ccrs.PlateCarree())
#divergence
#div=ax.contour(lon, lat, synoptic_scale, levels=list(range(2,20,2)), colors='red',fontsize=12, linewidths=2, fmt='%d',
               #transform=ccrs.PlateCarree())
#wind stuff
wnd_speed= np.arange(60, 180, 20)
cf=ax.contourf(lon,lat, speed_300, wnd_speed, cmap='turbo', transform=ccrs.PlateCarree())
ax.barbs(lon, lat, u_kts,v_kts,
         length=8, regrid_shape=15, pivot='middle', transform=ccrs.PlateCarree())
ax.clabel(line,inline=True,colors='black',fontsize=12, fmt='%d')

#additional plot settings
ax.coastlines('50m')
ax.add_feature(cfeature.BORDERS.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.set_title(f'Model: RAP', fontsize=16, loc='left')
ax.set_title(f'250mb Height (dam), Wind (kts)', fontsize=22, loc='center')
ax.set_title(f'\nValid {time_str}, FHR: {fhr}', fontsize=16, loc='right')

plt.tight_layout()
plt.show()