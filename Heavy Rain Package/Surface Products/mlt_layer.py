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
H= Herbie("2026-06-06 18:00", model='nam', product='awphys', fxx=54)

print(H.inventory(":HGT:0C"))

ds=H.xarray(":HGT:0C isotherm:", remove_grib=False)
ds2=H.xarray(":MSLMA:mean sea level:", remove_grib=False)
ds3 = H.xarray(":[UV]GRD:10 m above ground:", remove_grib=False)



mlt=ds['gh']
mslp=ds2['mslma']/100
lat=ds['latitude']
lon=ds['longitude']
time=ds['time']
u_kts=ds3['u10']*1.944
v_kts=ds3['v10']*1.944


#strip time
time_valid=pd.Timestamp(ds['time'].values)
time_str=time_valid.strftime('%Y%m%d %H%MZ')

#smoothing
mslp_smooth=gaussian_filter(mslp,2)
mlt_smooth=gaussian_filter(mlt,2)
#plot settings s
proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35,standard_parallels=(30, 60))
fig=plt.figure(figsize=(20,12))
ax=plt.axes(projection=proj)

#set bounds
ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())
cf=ax.contourf(lon,lat,mlt_smooth,cmap='viridis',levels=list(np.arange(0, 6000, 400)),transform=ccrs.PlateCarree())
#mslp contours
line1=ax.contour(lon, lat,mslp_smooth,levels=list(np.arange(940, 1072, 2)),colors='black',linewidths=2,
                transform=ccrs.PlateCarree())
ax.clabel(line1,inline=True,colors='black',fontsize=12, fmt='%d')
#wind stuff
ax.barbs(lon, lat, u_kts,v_kts,
         length=6, regrid_shape=20, fill_empty=False, pivot='middle', transform=ccrs.PlateCarree())


#cbar=plt.colorbar(cf, ax=ax, orientation='horizontal', shrink=0.5,)

#additional plot settings
ax.coastlines('50m')
ax.add_feature(cfeature.BORDERS.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.set_title(f'Model:{H.model.upper()} | F{H.fxx:03d}', fontsize=18, loc='left')
ax.set_title(f'Freezing Level (m; shaded), MSLP (hPa), 10-m Winds (kts)', fontsize=22, loc='center')
ax.set_title(f'\nValid {time_str}', fontsize=22, loc='right')


plt.tight_layout()
plt.show()
