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
H= Herbie("2026-08-12 13:00", model='rap', product='awp236pgrb', fxx=0)



ds = H.xarray(":(?:TMP|RH):2 m above ground:", remove_grib=False)
ds2 = H.xarray(":(?:UGRD|VGRD):10 m above ground:", remove_grib=False)
ds3= H.xarray(":PRES:surface:", remove_grib=False)
ds4=H.xarray(":(?:TMP|RH|HGT):925 mb", remove_grib=False)
print(ds, ds2,ds3)

u=ds2['u10']
v=ds2['v10']
lat=ds['latitude']
lon=ds['longitude']
temp=ds['t2m']
time=ds['time']
rh=ds['r2']
pressure=ds3['sp']

#925
temp2=ds4['t']
rh2=ds4['r']
pressure_925=92500*units.Pa

#conversions
temp_c=ds['t2m']-273.15
u_kts=ds2['u10']*1.944
v_kts=ds2['v10']*1.944


#strip time
time_valid=pd.Timestamp(ds['time'].values)
time_str=time_valid.strftime('%HZ')
valid_time = time_valid.strftime('%HZ %a %b %d %Y')


#calculate mixing ratio sfc
mr=metpy.calc.mixing_ratio_from_relative_humidity(pressure.values*units.Pa,temp.values*units.kelvin,
                                                  (rh.values/100)*units.dimensionless,
                                                  phase='liquid')

#calculate mixing ratio 925mb
mr925=metpy.calc.mixing_ratio_from_relative_humidity(pressure_925,temp2.values*units.kelvin,
                                                  (rh2.values/100)*units.dimensionless,
                                                  phase='liquid')

mrgkg=mr.magnitude*1000
mrgkg_925=mr925.magnitude*1000

avg_mr=(mrgkg+mrgkg_925)/2

mrgkg_smooth=gaussian_filter(avg_mr,2)

#plot settings s
proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35,standard_parallels=(30, 60))
fig=plt.figure(figsize=(20,12))
ax=plt.axes(projection=proj)

#set bounds
ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())
cf=ax.contourf(lon,lat,mrgkg_smooth,cmap='BrBG',transform=ccrs.PlateCarree())
#mslp contours
line1=ax.contour(lon, lat,mrgkg_smooth,levels=list(np.arange(0, 20, 2)),colors='black',linewidths=2,
                transform=ccrs.PlateCarree())
ax.clabel(line1,inline=True,colors='black',fontsize=12, fmt='%d')

#streamlines
ax.streamplot(lon,lat, u, v, transform=ccrs.PlateCarree(),color='black', linewidth=1, density=5)

#cbar=plt.colorbar(cf, ax=ax, orientation='horizontal', shrink=0.5,)

#additional plot settings
ax.coastlines('50m')
ax.add_feature(cfeature.BORDERS.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.set_title(f'Model: {time_str} {H.model.upper()} | F{H.fxx:03d}', fontsize=22, loc='left')
ax.set_title(f'100mb Mean Mixing Ratio and Surface Streamlines', fontsize=22, loc='center')
ax.set_title(f'\nValid: {valid_time}', fontsize=22, loc='right')


plt.tight_layout()
plt.show()