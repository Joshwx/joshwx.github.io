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
from datetime import datetime, timedelta, timezone
#initalize past timeframe
target_time = (pd.Timestamp.now(tz='UTC') - pd.Timedelta(hours=24)).tz_localize(None).floor('h')
#pull from model
H= Herbie(target_time, model='rap', product='awp236pgrb')

print(H.inventory(":700 mb"))

ds=H.xarray(":(?:HGT|UGRD|VGRD|TMP|RH):700 mb:", remove_grib=False)
print(ds)
ds2=H.xarray(":RH:(700|650|600|550|500) mb")
print(ds2)

gh=ds['gh']
u=ds['u']
v=ds['v']
lat=ds['latitude']
lon=ds['longitude']
temp=ds['t']
time=ds['time']
temp_c=ds['t']-273.15
u_kts=ds['u']*1.944
v_kts=ds['v']*1.944
#extract rhs
rh_700=ds2['r'].sel(isobaricInhPa=700)
rh_650=ds2['r'].sel(isobaricInhPa=650)
rh_600=ds2['r'].sel(isobaricInhPa=600)
rh_550=ds2['r'].sel(isobaricInhPa=550)
rh_500=ds2['r'].sel(isobaricInhPa=500)

#mean
# Stack along a new axis, then take the mean
mean_rh = np.mean(np.stack([rh_700, rh_650, rh_600, rh_550, rh_500], axis=0), axis=0)



#strip time
time_valid=pd.Timestamp(ds['time'].values)
time_str=time_valid.strftime('%HZ')
valid_time = time_valid.strftime('%HZ %a %b %d %Y')

#get the whole wind vector
speed_700=np.sqrt(u_kts**2+v_kts**2)
#smooth gh
gh_smooth=gaussian_filter(gh,2)
temp_smooth=gaussian_filter(temp_c,2)
smooth_rh=gaussian_filter(mean_rh,2)

from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import numpy as np

def load_ncl_rgb(filepath, smooth=True):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    rgb_vals = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith(';') or line.startswith('#') or line.startswith('!') or '=' in line:
            continue
        parts = line.split()
        if len(parts) >= 3:
            rgb_vals.append([float(x) for x in parts[:3]])

    rgb = np.array(rgb_vals)

    # normalize only if values are in 0-255 range
    if rgb.max() > 1.0:
        rgb = rgb / 255.0

    if smooth:
        return LinearSegmentedColormap.from_list('custom_cmap', rgb, N=256)
    else:
        return ListedColormap(rgb)


cmap_greens = load_ncl_rgb('C:/Users/sherm/Desktop/CMAP/MPL_Greens.rgb')


#plot settings s
proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35,
                              standard_parallels=(30, 60))
fig=plt.figure(figsize=(20,12))
ax=plt.axes(projection=proj)

#set bounds
ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())
#GH contours
line=ax.contour(lon, lat,gh_smooth,levels=list(range(2340, 3780, 30)),colors='black',linewidths=3,
                transform=ccrs.PlateCarree())
#rh fill
cf_rh=ax.contourf(lon, lat, smooth_rh,levels=np.arange(70,100,5), cmap=cmap_greens, transform=ccrs.PlateCarree())

#wind stuff
ax.barbs(lon, lat, u_kts,v_kts,
         length=8, regrid_shape=15, pivot='middle', transform=ccrs.PlateCarree())
ax.clabel(line,inline=True,colors='black',fontsize=12, fmt='%d')

#temp contours
cf_red=ax.contour(lon, lat, temp_smooth, levels=(range(5,40,5)), colors='red',linewidths=1.5,
                  linestyles='dashed',transform=ccrs.PlateCarree())
ax.clabel(cf_red,inline=True,colors='red',fontsize=15, fmt='%d')
cf_red_bold=ax.contour(lon, lat, temp_smooth, levels=(range(10,30,5)), colors='red',linewidths=2,
                  linestyles='dashed',transform=ccrs.PlateCarree())
ax.clabel(cf_red_bold,inline=True,colors='red',fontsize=15, fmt='%d')
cf_blue=ax.contour(lon, lat, temp_smooth, levels=(range(-40,0,5)), colors='blue',linewidths=1.5,
                   linestyles='dashed',transform=ccrs.PlateCarree())
ax.clabel(cf_blue,inline=True,colors='blue',fontsize=15, fmt='%d')

cf_freez=ax.contour(lon, lat, temp_smooth, levels=[0], colors='blue',linewidths=2,
                   linestyles='solid',transform=ccrs.PlateCarree())
ax.clabel(cf_freez,inline=True,colors='blue',fontsize=15, fmt='%d')



#additional plot settings
ax.coastlines('50m')
ax.add_feature(cfeature.BORDERS.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.set_title(f'Model: {time_str} {H.model.upper()} | F{H.fxx:03d}', fontsize=22, loc='left')
ax.set_title(f'700mb Height (dam), Wind (kts), Temp (C, red/blue),\n 700-500mb Mean RH >= 70%', fontsize=26, loc='center')
ax.set_title(f'\nValid: {valid_time}', fontsize=22, loc='right')

plt.tight_layout()
plt.show()



#vorticity calc
dx,dy=mpcalc.lat_lon_grid_deltas(lon, lat)
vort=mpcalc.absolute_vorticity(u,v,dx=dx, dy=dy)
synoptic_scale=vort*1e5

#plot settings s
proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35,
                              standard_parallels=(30, 60))
fig=plt.figure(figsize=(20,12))
ax=plt.axes(projection=proj)
ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())
#vort contour
cf=ax.contourf(lon,lat,synoptic_scale, levels=list(range(16,44,4)),
               cmap='plasma', transform=ccrs.PlateCarree())
#gh heights
#GH contours
line=ax.contour(lon, lat,gh_smooth[:,:],levels=list(range(2340, 3780, 30)),colors='black',linewidths=3,
                transform=ccrs.PlateCarree())
ax.clabel(line,inline=True,colors='black',fontsize=12, fmt='%d')
#additional plot settings
ax.coastlines('50m')
ax.add_feature(cfeature.BORDERS.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.set_title(f'Model: {time_str} {H.model.upper()} | F{H.fxx:03d}', fontsize=22, loc='left')
ax.set_title(f'700mb Height (dam), Absolute Vorticity (10^-5 s^-1) ', fontsize=22, loc='center')
ax.set_title(f'\nValid: {valid_time}', fontsize=22, loc='right')

plt.tight_layout()
plt.show()