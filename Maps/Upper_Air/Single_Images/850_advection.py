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
import numpy as np

#pull from model
H= Herbie("2026-01-24 00:00", model='rap', product='awp236pgrb', fxx=0)

# print(H.inventory(":850 mb"))

ds=H.xarray(":(?:HGT|UGRD|VGRD|TMP|RH):850 mb:", remove_grib=False)
# print(ds)

gh=ds['gh']/10
u=ds['u']
v=ds['v']
lat=ds['latitude']
lon=ds['longitude']
temp=ds['t']
time=ds['time']

u_kts=ds['u']*1.944
v_kts=ds['v']*1.944
temp_850=temp.metpy.convert_units('degC')
# print(temp_850)

# strip time
time_valid = pd.Timestamp(ds['time'].values)
time_str = time_valid.strftime('%HZ')
valid_time_moving = pd.Timestamp(ds['valid_time'].values)
valid_time = valid_time_moving.strftime('%HZ %a %b %d %Y')
#get the whole wind vector
speed_850=np.sqrt(u_kts**2+v_kts**2)
#mask settings
mask = speed_850<5
u_masked=np.ma.masked_array(u_kts,mask)
v_masked=np.ma.masked_array(v_kts,mask)

#advection calculation and conversion to C/3hr
dx,dy=metpy.calc.lat_lon_grid_deltas(lon, lat)
advection_calc=metpy.calc.advection(temp_850,u,v,dx=dx,dy=dy)
adv=advection_calc.metpy.convert_units('degC/s')
three_hr=(adv.metpy.convert_units('degC/hour')*3)
advection=three_hr

#smooth gh
gh_smooth=gaussian_filter(gh,2)
advection_smooth=gaussian_filter(advection,2)
temp_c_smooth=gaussian_filter(temp_850,2)

advection_masked = np.ma.masked_where(np.abs(advection_smooth) < 0.5, advection_smooth)

#plot settings s
proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35,standard_parallels=(30, 60))
fig=plt.figure(figsize=(20,12))
ax=plt.axes(projection=proj)
#set bounds
ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())

#advection stuff
cint=np.arange(-8,9)
cf=ax.contourf(lon,lat,advection_masked,cint[cint!=0],extend='both',cmap='bwr',transform=ccrs.PlateCarree())
line1=ax.contour(lon, lat,advection_masked,colors='black',linewidths=.5,
                transform=ccrs.PlateCarree())
ax.clabel(line1,inline=True,colors='black',fontsize=12, fmt='%d')

#GH contours
line=ax.contour(lon, lat,gh_smooth,levels=list(range(90, 300, 3)),colors='black',linewidths=3,
                transform=ccrs.PlateCarree())
#wind stuff
wnd_speed= np.arange(25, 85, 5)
speed_850=np.sqrt(u_kts**2+v_kts**2)
# wind=ax.contourf(lon,lat,speed_850, wnd_speed, cmap='turbo', transform=ccrs.PlateCarree())
ax.barbs(lon, lat, u_masked,v_masked,
         length=8, regrid_shape=15, fill_empty=False, pivot='middle', transform=ccrs.PlateCarree())
ax.clabel(line,inline=True,colors='black',fontsize=12, fmt='%d')




#additional plot settings
ax.coastlines('50m')
ax.add_feature(cfeature.BORDERS.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.set_title(f'Model: {time_str} {H.model.upper()} | F{H.fxx:03d}', fontsize=22, loc='left')
ax.set_title(f'850mb Temperature Advection (10^-5 C/3hr)', fontsize=26, loc='center')
ax.set_title(f'\nValid: {valid_time}', fontsize=22, loc='right')

plt.tight_layout()
plt.show()