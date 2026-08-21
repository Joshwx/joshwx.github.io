import numpy
import sounderpy as spy
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import metpy
import metpy.calc as mpcalc
from cmweather.cm_colorblind import cmap
from metpy.units import units
from herbie import Herbie
from herbie.toolbox import EasyMap, pc
from herbie import paint
import cartopy
from cartopy import crs as ccrs, feature as cfeature
import scipy
from scipy.ndimage import gaussian_filter
import pandas as pd
from matplotlib.colors import ListedColormap
from metpy.plots.ctables import registry
from metpy.plots import ctables

#pull from model
H= Herbie("2026-08-08 09:00", model='rap', product='awp236pgrb')

print(H.inventory(":MSLMA:mean sea level:"))

ds=H.xarray(":PWAT:entire atmosphere", remove_grib=False)
print(ds)

ds2=H.xarray(":MSLMA:mean sea level:", remove_grib=False)
print(ds2)


pwat=ds["pwat"]/25.4
mslp=ds2['mslma']/100
lat=ds['latitude']
lon=ds['longitude']
time=ds['time']

#strip time
time_valid=pd.Timestamp(ds['time'].values)
time_str=time_valid.strftime('%Y-%m-%d %H:%M UTC')
fhr=int(ds.step.values/np.timedelta64(1,'h'))

pwat_smooth=gaussian_filter(pwat,2)
mslp_smooth=gaussian_filter(mslp,2)


#plot settings s
proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35,
                              standard_parallels=(30, 60))
fig=plt.figure(figsize=(20,12))
ax=plt.axes(projection=proj)


#set bounds
ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())

cf=ax.contourf(lon,lat,pwat_smooth,cmap='BrBG',transform=ccrs.PlateCarree())

#mslp contours
line1=ax.contour(lon, lat,mslp_smooth,levels=list(np.arange(940, 1072, 2)),colors='black',linewidths=2,
                transform=ccrs.PlateCarree())
ax.clabel(line1,inline=True,colors='black',fontsize=12, fmt='%d')

#cbar=plt.colorbar(cf, ax=ax, orientation='horizontal', shrink=0.5,)


#additional plot settings
ax.coastlines('50m')
ax.add_feature(cfeature.BORDERS.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.set_title(f'Model: RAP', fontsize=16, loc='left')
ax.set_title(f'MSL Pressure (mb), Precipitable Water (in)', fontsize=22, loc='center')
ax.set_title(f'\nValid {time_str}, FHR: {fhr}', fontsize=16, loc='right')


plt.tight_layout()
plt.show()
