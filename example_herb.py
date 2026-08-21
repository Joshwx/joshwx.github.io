import sounderpy as spy
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import metpy
from herbie import Herbie
from herbie.toolbox import EasyMap, pc
from herbie import paint
#pull from model
H= Herbie("2025-05-10 12:00", model='rap', product='awp130pgrb')
#search for feature you want and put into ds
ds=H.xarray("TMP:2 m above")

print(ds)

lat=ds['latitude'].values
lon=ds['longitude'].values
#plot settings
ax = EasyMap(crs=ds.herbie.crs, figsize=[8, 8]).STATES().ax
p = ax.pcolormesh(
    ds.longitude,
    ds.latitude,
    ds.t2m - 273.15,
    transform=pc,
    **paint.NWSTemperature.kwargs2,
)
plt.colorbar(
    p, ax=ax, orientation="horizontal", pad=0.05, **paint.NWSTemperature.cbar_kwargs2
)

ax.set_title(ds.t2m.GRIB_name, loc="right")
ax.set_title(f"{ds.model.upper()}: {H.product_description}", loc="left")

plt.show()