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
H= Herbie("2026-05-29 14:00", model='rap', product='awp236pgrb')

print(H.inventory(":850 mb"))

ds=H.xarray(":(?:HGT|UGRD|VGRD|TMP|RH):850 mb:", remove_grib=False)
print(ds)

gh=ds['gh']/10
u=ds['u']
v=ds['v']
lat=ds['latitude']
lon=ds['longitude']
temp=ds['t']
time=ds['time']
temp_c=ds['t']-273.15
u_kts=ds['u']*1.944
v_kts=ds['v']*1.944
rh=ds['r']
#strip time
time_valid=pd.Timestamp(ds['time'].values)
time_str=time_valid.strftime('%Y%m%d %H%MZ')
#get the whole wind vector
speed_850=np.sqrt(u_kts**2+v_kts**2)

dpt=metpy.calc.dewpoint_from_relative_humidity(temp_c.values*units.celsius,
                                               (rh.values/100)*units.dimensionless)
theta_e_calc=metpy.calc.equivalent_potential_temperature(850*units.hPa, temp.values*units.kelvin, dpt)
print(theta_e_calc)

#smooth gh
gh_smooth=gaussian_filter(gh,2)
temp_c_smooth=gaussian_filter(temp_c,1)
dpt_smooth=gaussian_filter(dpt, 1)
theta_e_smooth=gaussian_filter(theta_e_calc, 1)

#plot settings s
proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35,
                              standard_parallels=(30, 60))
fig=plt.figure(figsize=(20,12))
ax=plt.axes(projection=proj)

#set bounds
ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())
#GH contours
line=ax.contour(lon, lat,gh_smooth,levels=list(range(90, 300, 3)),colors='black',linewidths=3,
                transform=ccrs.PlateCarree())

#theta-e contours
cf=ax.contourf(lon, lat, theta_e_smooth, cmap='Spectral_r',linewidths=2,
                  linestyles='dashed',transform=ccrs.PlateCarree())
#wind stuff
ax.barbs(lon, lat, u_kts,v_kts,
         length=6, regrid_shape=20, fill_empty=False, pivot='middle', transform=ccrs.PlateCarree())
ax.clabel(line,inline=True,colors='black',fontsize=12, fmt='%d')


#additional plot settings
ax.coastlines('50m')
ax.add_feature(cfeature.BORDERS.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.set_title(f'Model: RAP', fontsize=22, loc='left')
ax.set_title(f'850mb Height (dam; black), Theta-e (K;shaded), Wind (kts)', fontsize=18, loc='center')
ax.set_title(f'\nValid {time_str}', fontsize=22, loc='right')

plt.tight_layout()
plt.show()
