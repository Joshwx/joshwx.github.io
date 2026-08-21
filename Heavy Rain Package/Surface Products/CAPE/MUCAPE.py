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
H= Herbie("2026-05-29 14:00", model='rap', product='awp236pgrb')

#print(H.inventory(":850 mb"))

ds=H.xarray(":CAPE:surface:", remove_grib=False)
ds2 = H.xarray(":[UV]GRD:10 m above ground:", remove_grib=False)
ds3=H.xarray(":CIN:surface:", remove_grib=False)
print(ds)


cape=ds['cape']
cin=ds3['cin']
lat=ds['latitude']
lon=ds['longitude']
time=ds['time']
u_kts=ds2['u10']*1.944
v_kts=ds2['v10']*1.944

#strip time
time_valid=pd.Timestamp(ds['time'].values)
time_str=time_valid.strftime('%Y%m%d %H%MZ')

cape_mask=np.ma.masked_where(cape<=100,cape)

#smooth
cape_smooth=gaussian_filter(cape,1)
cin_smooth=gaussian_filter(cin,2)

#plot settings s
proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35,
                              standard_parallels=(30, 60))
fig=plt.figure(figsize=(20,12))
ax=plt.axes(projection=proj)

#set bounds
ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())

#contours
line=ax.contour(lon,lat,cape_smooth,levels=list(range(100,250,150)),colors='red',transform=ccrs.PlateCarree())
ax.clabel(line,inline=True,colors='black',fontsize=12, fmt='%d')

#wind stuff
ax.barbs(lon, lat, u_kts,v_kts,
         length=6, regrid_shape=20, fill_empty=False, pivot='middle', transform=ccrs.PlateCarree())


#additional plot settings
ax.coastlines('50m')
ax.add_feature(cfeature.BORDERS.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.set_title(f'Model: RAP', fontsize=22, loc='left')
ax.set_title(f'SBCAPE (J/kg), Wind (kts)', fontsize=22, loc='center')
ax.set_title(f'\nValid {time_str}', fontsize=22, loc='right')

plt.tight_layout()
plt.show()
