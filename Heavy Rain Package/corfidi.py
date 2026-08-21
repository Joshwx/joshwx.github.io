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
H= Herbie("2026-08-19 01:00", model='rap', product='awp236pgrb', fxx=0)
ds = H.xarray(":(UGRD|VGRD):(850|700|500|300) mb")

#strip time
time_valid=pd.Timestamp(ds['time'].values)
time_str=time_valid.strftime('%HZ')
valid_time = time_valid.strftime('%HZ %a %b %d %Y')

MS_TO_KTS = 1.94384

u_850 = ds['u'].sel(isobaricInhPa=850) * MS_TO_KTS
u_700 = ds['u'].sel(isobaricInhPa=700) * MS_TO_KTS
u_500 = ds['u'].sel(isobaricInhPa=500) * MS_TO_KTS
u_300 = ds['u'].sel(isobaricInhPa=300) * MS_TO_KTS

v_850 = ds['v'].sel(isobaricInhPa=850) * MS_TO_KTS
v_700 = ds['v'].sel(isobaricInhPa=700) * MS_TO_KTS
v_500 = ds['v'].sel(isobaricInhPa=500) * MS_TO_KTS
v_300 = ds['v'].sel(isobaricInhPa=300) * MS_TO_KTS

lon=ds['longitude']
lat=ds['latitude']
#mean wind
u_mean = (u_850 + u_700 + u_500 + u_300) / 4
v_mean = (v_850 + v_700 + v_500 + v_300) / 4
mean_wind_speed = np.sqrt(u_mean ** 2 + v_mean ** 2)


#llj

llj=np.sqrt(u_850**2 + v_850**2)


u_vmbe=u_mean-u_850
v_vmbe=v_mean-v_850

vmbe_mag=np.sqrt(u_vmbe ** 2 + v_vmbe ** 2)

vmbe=mean_wind_speed - llj



print(vmbe)

# #plot settings s
proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35,
                              standard_parallels=(30, 60))
fig=plt.figure(figsize=(20,12))
ax=plt.axes(projection=proj)

#set bounds
ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())

#wind stuff
wnd_speed= np.arange(25, 85, 5)
speed_850=np.sqrt(u_850**2+v_850**2)
# wind=ax.contourf(lon,lat,speed_850, wnd_speed, cmap='turbo', transform=ccrs.PlateCarree())
# ax.barbs(lon, lat, u_vmbe,v_vmbe,length=8, regrid_shape=15, fill_empty=False, pivot='middle', transform=ccrs.PlateCarree())

ax.barbs(lon, lat, u_vmbe,v_vmbe,
         length=8, regrid_shape=15, fill_empty=False, pivot='middle', transform=ccrs.PlateCarree(), color='red')



#additional plot settings
ax.coastlines('50m')
ax.add_feature(cfeature.BORDERS.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.set_title(f'Model: Rapid Refresh (RAP)', fontsize=15, loc='left')
ax.set_title(f'Upwind Propagation Vectors\n (V_MCL: Blue, V_LLJ: Green, V_MBE: Yellow) ', fontsize=20, loc='center')
ax.set_title(f'\nValid {time_str}, FHR: {valid_time}', fontsize=15, loc='right')

plt.tight_layout()
plt.show()