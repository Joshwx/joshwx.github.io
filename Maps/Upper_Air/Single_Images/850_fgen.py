import matplotlib.pyplot as plt
import numpy as np
import metpy.calc as mpcalc
from metpy.units import units
from herbie import Herbie
from cartopy import crs as ccrs, feature as cfeature
from scipy.ndimage import gaussian_filter
import pandas as pd

# pull from model
H = Herbie("2026-10-10 06:00", model="gfs", product="pgrb2.0p25")
ds = H.xarray(":(?:HGT|UGRD|VGRD|TMP|RH):850 mb:", remove_grib=False)

# subset BEFORE computing anything
ds = ds.sel(latitude=slice(58, 15), longitude=slice(230, 300)).load()

gh = ds['gh'] / 10
u = ds['u']
v = ds['v']
lat = ds['latitude']
lon = ds['longitude']
temp_k = ds['t']
level = 850 * units.hPa

u_kts = u * 1.944
v_kts = v * 1.944

dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat)

# frontogenesis
theta_850 = mpcalc.potential_temperature(level, temp_k)
fronto_850 = mpcalc.frontogenesis(theta_850, u, v, dx, dy)
converted_fgen = fronto_850 * (1000 * 100 * 3600 * 3)

#temp adv
adv = mpcalc.advection(temp_k, u, v, dx=dx, dy=dy)
adv_3hr = adv.metpy.convert_units('kelvin/hour').metpy.magnitude * 3

#timestuff
time_str = pd.Timestamp(ds['time'].values).strftime('%HZ')
valid_time = pd.Timestamp(ds['valid_time'].values).strftime('%HZ %a %b %d %Y')

#wind and masking
speed_850 = np.sqrt(u_kts**2 + v_kts**2)
mask = np.asarray(speed_850) < 5
u_masked = np.ma.masked_array(np.asarray(u_kts), mask)
v_masked = np.ma.masked_array(np.asarray(v_kts), mask)

#smoothing
gh_smooth = gaussian_filter(np.asarray(gh), 2)
adv_smooth = gaussian_filter(np.asarray(adv_3hr), 2)
fgen_smooth = gaussian_filter(np.asarray(converted_fgen), 2)

adv_masked = np.ma.masked_where(~np.isfinite(adv_smooth) | (np.abs(adv_smooth) < 0.5), adv_smooth)
fgen_masked = np.ma.masked_where(~np.isfinite(fgen_smooth) | (np.abs(fgen_smooth) <= 0.99), fgen_smooth)

#plot settings
proj = ccrs.LambertConformal(central_longitude=-96, central_latitude=35, standard_parallels=(30, 60))
fig = plt.figure(figsize=(20, 12))
ax = plt.axes(projection=proj)
ax.set_extent([-120., -72., 22., 50.], crs=ccrs.PlateCarree())

#temp adv
adv_levels = np.arange(-8, 8.5, 1)
cf = ax.contourf(lon, lat, adv_masked, adv_levels, cmap='bwr', extend='both',
                 transform=ccrs.PlateCarree())
cbar = plt.colorbar(cf, ax=ax, orientation='horizontal', pad=0.03, shrink=0.6, aspect=40)
cbar.set_label('850mb Temperature Advection (°C/3hr)', fontsize=14)

# frontogenesis
fgen_levels = list(range(1, 25, 1))
line1 = ax.contour(lon, lat, fgen_masked, levels=fgen_levels, colors='purple',
                   linewidths=2, transform=ccrs.PlateCarree())
ax.clabel(line1, inline=True, colors='purple', fontsize=14, fmt='%d')

#gh
line = ax.contour(lon, lat, gh_smooth, levels=list(range(90, 300, 3)), colors='black',
                  linewidths=3, transform=ccrs.PlateCarree())
ax.clabel(line, inline=True, colors='black', fontsize=14, fmt='%d')

#barbs
ax.barbs(lon, lat, u_masked, v_masked, length=8, regrid_shape=15, fill_empty=False,
         pivot='middle', transform=ccrs.PlateCarree())

#map features
ax.coastlines('50m')
ax.add_feature(cfeature.BORDERS.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.set_title(f'Model: {time_str} {H.model.upper()} | F{H.fxx:03d}', fontsize=15, loc='left')
ax.set_title('850mb Temp Advection (shaded, °C/3hr) & \n Frontogenesis (purple, K/100km/3hr)',
             fontsize=20, loc='center')
ax.set_title(f'\nValid: {valid_time}', fontsize=15, loc='right')

plt.tight_layout()
plt.show()