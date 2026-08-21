import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
#connect to the dataset that I downloaded
ds=xr.open_dataset('/Users/joshsherman/Desktop/cape.1993.nc')
print(ds)

cape=ds['cape']
lat=ds['lat'].values
lon=ds['lon'].values
print(cape)

cape_single=cape.isel(time=0)
#mask values lower than 100 to show more significant cape values, dont really care about 1 J/kg
cape_masked=cape_single.where(cape_single>100)
print(cape_single)

#plot settings along with contour settings
fig=plt.figure(figsize=(20,20))
ax=plt.axes(projection=ccrs.PlateCarree())
s=cape_masked.plot(ax=ax, transform=ccrs.PlateCarree(), cmap='YlOrRd', add_colorbar=False)
cbar=plt.colorbar(s, ax=ax,orientation='horizontal')
cbar.set_label('SBCAPE (J/kg)', fontsize=25)
cbar.ax.tick_params(labelsize=25)

#set bounds, contours, and additional settings
bounds = [(-100., -60., 23., 50.)]
line=ax.contour(lon, lat,cape_single[:,:],levels=list(range(100, 10000, 100)),colors='black',linewidths=3, transform=ccrs.PlateCarree())
ax.clabel(line,inline=True,colors='black',fontsize=20)
ax.set_extent(*bounds,crs=ccrs.PlateCarree())
ax.coastlines('50m')
ax.add_feature(cfeature.BORDERS.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))
ax.set_title('3-13-93 12Z Surface Based CAPE (J/kg) (Shaded, Contour)', fontsize=30)


plt.tight_layout()
plt.show()