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
H= Herbie("2026-05-23 19:00", model='rap', product='awp130pgrb')

#print(H.inventory(":850 mb"))

ds=H.xarray("DPT:2 m above", remove_grib=False)
ds2 = H.xarray(":[UV]GRD:10 m above ground:", remove_grib=False)

dpt=ds['d2m']
lat=ds['latitude']
lon=ds['longitude']
time=ds['time']
u_kts=ds2['u10']*1.944
v_kts=ds2['v10']*1.944
#strip time
time_valid=pd.Timestamp(ds['time'].values)
time_str=time_valid.strftime('%Y%m%d %H%MZ')

dpt_c=dpt-273.15
deg_f=(dpt_c*(9/5))+32

dpt_smooth=gaussian_filter(deg_f, 2)

#plot settings s
proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35,
                              standard_parallels=(30, 60))
fig=plt.figure(figsize=(20,12))
ax=plt.axes(projection=proj)

#set bounds
ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())

#contours
cf=ax.contourf(lon, lat, dpt_smooth, levels=(range(10,80,2)), cmap='BrBG',linewidths=2,
                  linestyles='dashed',transform=ccrs.PlateCarree())
line1=ax.contour(lon, lat,dpt_smooth,levels=list(np.arange(55, 80, 5)),colors='black',linewidths=2,
                transform=ccrs.PlateCarree())
ax.clabel(line1,inline=True,colors='black',fontsize=12, fmt='%d')
line2=ax.contour(lon, lat,dpt_smooth,levels=list(np.arange(45, 55, 5)),colors='black',linewidths=1,linestyles='dashed',
                transform=ccrs.PlateCarree())
ax.clabel(line1,inline=True,colors='black',fontsize=12, fmt='%d')
#wind stuff
ax.barbs(lon, lat, u_kts,v_kts,
         length=6, regrid_shape=20, fill_empty=False, pivot='middle', transform=ccrs.PlateCarree())


#additional plot settings
ax.coastlines('50m')
ax.add_feature(cfeature.BORDERS.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.set_title(f'Model: RAP', fontsize=22, loc='left')
ax.set_title(f'Surface Dewpoint (F), Wind (kts)', fontsize=22, loc='center')
ax.set_title(f'\nValid {time_str}', fontsize=22, loc='right')

plt.tight_layout()
plt.show()
