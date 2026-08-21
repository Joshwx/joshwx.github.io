# Core packages
import gzip
import tempfile

# File handling (if you're downloading MRMS .grib2.gz files manually)
import urllib.request
from datetime import datetime, timedelta
from io import StringIO

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmweather  # noqa: F401
import matplotlib.colors as mcolors

# Plotting
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd
import requests
import s3fs
import xarray as xr
from IPython.display import HTML  # To display the animation

# Animation
from matplotlib.animation import ArtistAnimation, PillowWriter
from metpy.plots import ctables, StationPlot  # For NWS reflectivity colormap
from scipy.interpolate import RegularGridInterpolator
import os
import tempfile
import pyart
from siphon.catalog import TDSCatalog
from datetime import timedelta, datetime
import metpy.plots as mpplots
import matplotlib.pyplot as plt
import cartopy as ccrs
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from metpy.calc import reduce_point_density
from metpy.cbook import get_test_data
from metpy.io import metar
from metpy.plots import add_metpy_logo, current_weather, sky_cover, StationPlot
from metpy.units import units
from metpy.io import parse_metar_file

now=datetime.now()
year=now.year
month=now.month
day=now.day
hour=now.hour
minute=now.minute


### Metar Block ###
cat=TDSCatalog('https://thredds.ucar.edu/thredds/catalog/noaaport/text/metar/catalog.xml')
metar_lst=list(cat.datasets.values())[-1]
print(metar_lst)



metar_ds=cat.datasets[0]
metar_ds.download()
data=parse_metar_file(metar_ds.name)
print(data)

# Drop rows with missing winds
data = data.dropna(how='any', subset=['wind_direction', 'wind_speed'])

proj = ccrs.LambertConformal(central_longitude=-95, central_latitude=35,
                             standard_parallels=[35])
# Use the Cartopy map projection to transform station locations to the map and
# then refine the number of stations plotted by setting a 300km radius
point_locs = proj.transform_points(ccrs.PlateCarree(), data['longitude'].values,
                                   data['latitude'].values)
data = data[reduce_point_density(point_locs, 75000.)]



tempf=(data['air_temperature'].values* units.degC).to('degF')
dpt=(data['dew_point_temperature'].values* units.degC).to('degF')


# fig = plt.figure(figsize=(20, 10))
# ax = fig.add_subplot(1, 1, 1, projection=proj)
#
# ax.add_feature(cfeature.COASTLINE)
# ax.add_feature(cfeature.STATES)
# ax.add_feature(cfeature.BORDERS)
#
# ax.set_extent([-91, -80, 34.5, 41])
# stationplot = StationPlot(ax, data['longitude'].values, data['latitude'].values,
#                           clip_on=True, transform=ccrs.PlateCarree(), fontsize=12)
#
# stationplot.plot_parameter('NW', tempf , color='red', formatter=lambda v: format(v, '.0f'))
# stationplot.plot_parameter('SW', dpt,formatter=lambda v: format(v, '.0f'),
#                            color='darkgreen')
#
# stationplot.plot_parameter('NE', data['air_pressure_at_sea_level'].values,
#                            formatter=lambda v: format(10 * v, '.0f')[-3:])
# stationplot.plot_symbol('C', data['cloud_coverage'].fillna(0), sky_cover)
#
# stationplot.plot_symbol('W', data['current_wx1_symbol'].fillna(0), current_weather)
#
# stationplot.plot_barb(data['eastward_wind'].values, data['northward_wind'].values)
#
# # stationplot.plot_text((2, 0), data['station_id'].values)
# plt.tight_layout()
# plt.show()


################################ MRMS Block ################################

# Define the URL to the compressed MRMS GRIB2 file for a specific timestamp
url = "https://noaa-mrms-pds.s3.amazonaws.com/CONUS/LayerCompositeReflectivity_Low_00.50/20260814/MRMS_LayerCompositeReflectivity_Low_00.50_20260814-195641.grib2.gz	"

# Download the file as bytes
response = urllib.request.urlopen(url)
compressed_file = response.read()

# Decompress and load into xarray using a temporary file


# Decompress and load into xarray using a temporary file
tmp_path = None
try:
    with tempfile.NamedTemporaryFile(suffix=".grib2", delete=False) as f:
        f.write(gzip.decompress(compressed_file))
        tmp_path = f.name  # file handle closes when the `with` block exits

    # Now the file is closed, so cfgrib can open it on Windows
    data_in = xr.load_dataarray(tmp_path, engine="cfgrib", decode_timedelta=True)
finally:
    if tmp_path and os.path.exists(tmp_path):
        os.remove(tmp_path)

refl_cmap = ctables.registry.get_colortable('NWSReflectivity')

# 2. Extract coords & data
lons = data_in.longitude.values
lats = data_in.latitude.values
refl = data_in.values

# If coords are 1D, make them 2D
if lons.ndim == 1 and lats.ndim == 1:
    lons, lats = np.meshgrid(lons, lats)


# 3. Plot
fig = plt.figure(figsize=(20, 16))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-91, -80, 34.5, 41], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE, linewidth=1)
ax.add_feature(cfeature.BORDERS, linewidth=1)
ax.add_feature(cfeature.STATES, linewidth=0.5)

### MRMS Settings ###
mesh = ax.pcolormesh(
    lons,
    lats,
    ma.masked_where(refl < 5, refl),
    cmap=refl_cmap,
    vmin=5,
    vmax=60,
    transform=ccrs.PlateCarree(), alpha=.5,
)

cb = plt.colorbar(mesh, ax=ax, orientation="horizontal", pad=0.05, aspect=50)
cb.set_label("Reflectivity (dBZ)")

### Metar Settings ###
stationplot = StationPlot(ax, data['longitude'].values, data['latitude'].values,
                          clip_on=True, transform=ccrs.PlateCarree(), fontsize=12)

stationplot.plot_parameter('NW', tempf , color='red', formatter=lambda v: format(v, '.0f'))
stationplot.plot_parameter('SW', dpt,formatter=lambda v: format(v, '.0f'),
                           color='darkgreen')

stationplot.plot_parameter('NE', data['air_pressure_at_sea_level'].values,
                           formatter=lambda v: format(10 * v, '.0f')[-3:])
stationplot.plot_symbol('C', data['cloud_coverage'].fillna(0), sky_cover)

stationplot.plot_symbol('W', data['current_wx1_symbol'].fillna(0), current_weather)

stationplot.plot_barb(data['eastward_wind'].values, data['northward_wind'].values)


plt.title("MRMS Layer Composite Reflectivity – Kentucky", fontsize=14)
plt.tight_layout()

plt.show()