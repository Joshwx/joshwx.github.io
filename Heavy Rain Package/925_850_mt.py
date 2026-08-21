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
H= Herbie("2026-05-25 16:00", model='rap', product='awp130pgrb')

#850mb
ds=H.xarray(":(?:HGT|UGRD|VGRD|TMP|RH):850 mb:", remove_grib=False)

gh=ds['gh']
u=ds['u']
v=ds['v']
lat=ds['latitude']
lon=ds['longitude']
temp=ds['t']
time=ds['time']
rh=ds['r']
temp_c=ds['t']-273.15
u_kts=ds['u']*1.944
v_kts=ds['v']*1.944
press=850*units.hPa

#product of the 850 wind speed m/s
speed_850=np.sqrt(u**2+v**2)

#925
ds2=H.xarray(":(?:HGT|UGRD|VGRD|TMP|RH):925 mb:", remove_grib=False)
#print(ds)

gh2=ds2['gh']
u2=ds2['u']
v2=ds2['v']
temp2=ds2['t']
rh2=ds2['r']
temp_c2=ds2['t']-273.15
u_kts2=ds2['u']*1.944
v_kts2=ds2['v']*1.944
press2=925*units.hPa

#product of the wind speed m/s
speed_925=np.sqrt(u2**2+v2**2)

#strip time
time_valid=pd.Timestamp(ds['time'].values)
time_str=time_valid.strftime('%Y%m%d %H%MZ')

#calculate mixing ratio 850
mr=metpy.calc.mixing_ratio_from_relative_humidity(press,temp.values*units.kelvin,
                                                  (rh.values/100)*units.dimensionless,
                                                  phase='liquid')
#calculate mixing ratio 925
mr2=metpy.calc.mixing_ratio_from_relative_humidity(press2,temp2.values*units.kelvin,
                                                  (rh2.values/100)*units.dimensionless,
                                                  phase='liquid')

#calculate theta-e 850
dp=metpy.calc.dewpoint_from_relative_humidity(temp.values*units.kelvin,(rh.values/100)*units.dimensionless)
te=metpy.calc.equivalent_potential_temperature(press, temp.values*units.kelvin, dp)

#calculate theta-e 925
dp2=metpy.calc.dewpoint_from_relative_humidity(temp2.values*units.kelvin,(rh2.values/100)*units.dimensionless)
te2=metpy.calc.equivalent_potential_temperature(press2, temp2.values*units.kelvin, dp2)

#moisture transport (multiply wind speed (m/s) and mixing ratio (g/g) and scale by 100)
mt=(speed_850*mr.magnitude)*100

#moisture transport (multiply wind speed (m/s) and mixing ratio (g/g) and scale by 100)
mt2=(speed_925*mr2.magnitude)*100

#average
total_mt=(mt+mt2)/2
total_te=(te+te2)/2
u_avg=(u+u2)/2
v_avg=(v+v2)/2

#smooth
mt_smooth=gaussian_filter(total_mt,2)
te_smooth = gaussian_filter(total_te.magnitude, 5)

#plot settings s
proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35,standard_parallels=(30, 60))
fig=plt.figure(figsize=(20,12))
ax=plt.axes(projection=proj)

#set bounds
ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())
#moisture transport stuff
cf=ax.contourf(lon,lat,mt_smooth,levels=list(np.arange(10, 40, 5)),cmap='RdYlGn_r',transform=ccrs.PlateCarree())
line1=ax.contour(lon, lat,mt_smooth,levels=list(np.arange(10, 40, 5)),colors='black',linewidths=1,
                transform=ccrs.PlateCarree())
ax.clabel(line1,inline=True,colors='black',fontsize=12, fmt='%d')
#transport vectors
ax.quiver(lon.values[::10, ::10], lat.values[::10, ::10],u_avg.values[::10, ::10], v_avg.values[::10, ::10],
          transform=ccrs.PlateCarree(), scale=500, color='maroon')

#theta-e contours
line2=ax.contour(lon, lat,te_smooth,levels=list(range(250, 350, 5)),colors='green',linewidths=1,linestyles='dashed',
                transform=ccrs.PlateCarree())
ax.clabel(line2,inline=True,colors='black',fontsize=12, fmt='%d')




#cbar=plt.colorbar(cf, ax=ax, orientation='horizontal', shrink=0.5,)

#additional plot settings
ax.coastlines('50m')
ax.add_feature(cfeature.BORDERS.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.set_title(f'Model: RAP', fontsize=18, loc='left')
ax.set_title(f'925-850 mb Moisture Transport (Shaded and Vector), Layer Averaged Theta-e (dashed)',
             fontsize=18, loc='center')
ax.set_title(f'\nValid {time_str}', fontsize=18, loc='right')


plt.tight_layout()
plt.show()