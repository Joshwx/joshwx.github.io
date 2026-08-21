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
H= Herbie("2026-05-23 12:00", model='rap', product='awp236pgrb')

print(H.inventory(":CAPE"))

ds=H.xarray(":CAPE:surface:anl:", remove_grib=False)
print(ds)

ds2=H.xarray(":(?:UGRD|VGRD|):850 mb:", remove_grib=False)
print(ds2)

ds3=H.xarray(":(?:UGRD|VGRD|):500 mb:", remove_grib=False)
print(ds3)
#pull u and v
u_850=ds2['u']
v_850=ds2['v']
u_850_kts=u_850*1.944
v_850_kts=v_850*1.944
speed_850=np.sqrt(u_850_kts**2+v_850_kts**2)

u_500=ds3['u']
v_500=ds3['v']
u_500_kts=u_500*1.944
v_500_kts=v_500*1.944
speed_500=np.sqrt(u_500_kts**2+v_500_kts**2)


CAPE=ds['cape']
lat=ds['latitude']
lon=ds['longitude']

#strip time
time_valid=pd.Timestamp(ds['time'].values)
time_str=time_valid.strftime('%Y%m%d %H%MZ')

#plot settings
proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35,
                              standard_parallels=(30, 60))
fig=plt.figure(figsize=(20,12))
ax=plt.axes(projection=proj)

#set bounds
ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())
#cape shading
cf=ax.contourf(lon, lat,CAPE[:,:],levels=list(range(100, 10000, 500)),cmap='hsv',linewidths=2, linestyles='dashed',
                transform=ccrs.PlateCarree())
#wind stuff 850
ax.barbs(lon, lat, u_850_kts,v_850_kts,
         length=6, regrid_shape=20, fill_empty=False, pivot='middle', transform=ccrs.PlateCarree())
#wind stuff 500
ax.barbs(lon, lat, u_500_kts,v_500_kts,
         length=6, regrid_shape=20, fill_empty=False, pivot='middle', transform=ccrs.PlateCarree())

#additional plot settings
ax.coastlines('50m')
ax.add_feature(cfeature.BORDERS.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.set_title(f'Model: RAP', fontsize=22, loc='left')
ax.set_title(f'SBCAPE (contour), ', fontsize=22, loc='center')
ax.set_title(f'\nValid {time_str}', fontsize=22, loc='right')

plt.tight_layout()
plt.show()