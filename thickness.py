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
H= Herbie("2026-08-26 12:00", model='rap', product='awp236pgrb',fxx=1)

print(H.inventory(":APCP:"))

ds_1 = H.xarray(":HGT:(1000) mb:", remove_grib=False)
ds_2 = H.xarray(":HGT:(500) mb:", remove_grib=False)
print(ds_1)
ds3=H.xarray(":MSLMA:mean sea level:", remove_grib=False)
print(ds3)
ds_4=H.xarray(":APCP:surface:0-1 hour acc fcst:", remove_grib=False)
print(ds_4)

gh_1000=ds_1['gh']
gh_500=ds_2['gh']
lat=ds_1['latitude']
lon=ds_1['longitude']
time=ds_1['time']
mslp=ds3['mslma']/100
apcp=ds_4['tp']/25.4
print(apcp)



#calc thick
thick=(gh_500-gh_1000)/10
smooth=gaussian_filter(thick, sigma=2)
#smooth mslp
mslp_smooth=gaussian_filter(mslp,2)

# strip time
time_valid = pd.Timestamp(ds_1['time'].values)
time_str = time_valid.strftime('%HZ')
valid_time_moving = pd.Timestamp(ds_1['valid_time'].values)
valid_time = valid_time_moving.strftime('%HZ %a %b %d %Y')



#plot settings
proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35,
                              standard_parallels=(30, 60))
fig=plt.figure(figsize=(20,12))
ax=plt.axes(projection=proj)

#set bounds
ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())
#GH contours

line=ax.contour(lon, lat,smooth[:,:],levels=list(range(546, 600, 6)),colors='red',linewidths=2, linestyles='dashed',
                transform=ccrs.PlateCarree())
ax.clabel(line,inline=True,colors='red',fontsize=12, fmt='%d')

line2=ax.contour(lon, lat,thick[:,:],levels=list(range(498, 540, 6)),colors='blue',linewidths=2, linestyles='dashed',
                transform=ccrs.PlateCarree())
ax.clabel(line2,inline=True,colors='blue',fontsize=12, fmt='%d')

line3=ax.contour(lon, lat,mslp_smooth,levels=list(np.arange(940, 1072, 2)),colors='black',linewidths=2,
                transform=ccrs.PlateCarree())
ax.clabel(line3,inline=True,colors='black',fontsize=12, fmt='%d')

extent_precip=[.01,.05,.1,.15,.25,.35,.5,.75,1,1.5,2,3,4]
precip=ax.contourf(lon,lat,apcp[:,:],extent_precip ,cmap='NWSRef', transform=ccrs.PlateCarree())
cbar = plt.colorbar(precip, ax=ax, orientation='horizontal', pad=0.03, shrink=0.7)
cbar.set_label('1-hr Precip (in)')
cbar.ax.tick_params(labelsize=16)

#additional plot settings
ax.coastlines('50m')
ax.add_feature(cfeature.BORDERS.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.set_title(f'Model: {time_str} {H.model.upper()} | F{H.fxx:03d}', fontsize=15, loc='left')
ax.set_title(f'1000-500 mb Thickness (dam), MSLP (mb),\n 1hr Precip (in)', fontsize=20, loc='center')
ax.set_title(f'\nValid: {valid_time}', fontsize=15, loc='right')
plt.tight_layout()
plt.show()