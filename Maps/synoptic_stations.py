import matplotlib
matplotlib.use("Agg")

import os
import gzip
import tempfile
import warnings
import urllib.request
from datetime import datetime, timedelta, timezone

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import fsspec
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd
import xarray as xr
import cfgrib

from metpy.calc import reduce_point_density
from metpy.io import parse_metar_to_dataframe
from metpy.plots import StationPlot, sky_cover, current_weather, ctables
from metpy.units import units

from synoptic import TimeSeries

# Animation
from matplotlib.animation import ArtistAnimation, PillowWriter
from metpy.plots import ctables  # For NWS reflectivity colormap
from scipy.interpolate import RegularGridInterpolator
import os
from datetime import datetime, timezone
import urllib.request, gzip, tempfile, os
from datetime import datetime, timezone
import pyart
import json


from synoptic import Latest
output_dir= '/assets/maps/metar_mrms'
os.makedirs(output_dir, exist_ok=True)
######################### METAR BLOCK #########################

synoptic_token = os.environ["SYNOPTIC_TOKEN"]

#pull metar data from listed bounds (OH Valley and surrounding areas), put into df
bounds=[-91, -80, 34.5, 41]
df = TimeSeries(bbox=bounds,vars='metar',recent=timedelta(minutes=90),token=synoptic_token,).df()
metars=df

#extract raw code from metar df with nulls dropped
df_metar=metars['value_sting'].drop_nulls()
print(df_metar.head())
#parsing metar data to be readable and shoving it into a data frame also dropping any goofy metar strings
parse_data = []
print('parsing...')
for m in df_metar:
    text = m if m.strip().startswith(("METAR", "SPECI")) else "METAR " + m.strip()
    try:
        parse_data.append(parse_metar_to_dataframe(text))
    except Exception as e:
        print("skip:", m[:40], e)

data= pd.concat(parse_data, ignore_index=True)


print('parsed, now getting mrms data')

#drop any NA values in wind speed/dir and station IDS
data=data.dropna(how='any', subset=['wind_direction', 'wind_speed','station_id'])

#set projection for station plots
proj = ccrs.LambertConformal(central_longitude=-95, central_latitude=35,standard_parallels=[35])
point_locs = proj.transform_points(ccrs.PlateCarree(), data['longitude'].values,data['latitude'].values)

#change station plot density (meters)
data = data[reduce_point_density(point_locs, 35000.)]

#change deg C to F
tempf=(data['air_temperature'].values* units.degC).to('degF')
dpt=(data['dew_point_temperature'].values* units.degC).to('degF')

######################### MRMS BLOCK #########################

#ignore all RuntimeWarnings
warnings.filterwarnings("ignore", category=RuntimeWarning)


#Setup the AWS S3 filesystem
fs = fsspec.filesystem("s3", anon=True)

#extract current time
DATE = datetime.now(timezone.utc).strftime('%Y%m%d')

#grab url without having problems with aws
url = "https://mrms.ncep.noaa.gov/2D/LayerCompositeReflectivity_Low/MRMS_LayerCompositeReflectivity_Low.latest.grib2.gz"
req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
with urllib.request.urlopen(req) as response:
    compressed_file = response.read()
tmp_path = None

try:
    with tempfile.NamedTemporaryFile(suffix=".grib2", delete=False) as f:
        f.write(gzip.decompress(compressed_file))
        tmp_path = f.name
    data_in = xr.load_dataarray(tmp_path, engine="cfgrib", decode_timedelta=True)
finally:
    if tmp_path and os.path.exists(tmp_path):
        os.remove(tmp_path)

# refl_cmap = ctables.registry.get_colortable('NWSRef')

#extract lat lons and reflectvity values
lons = data_in.longitude.values
lats = data_in.latitude.values
refl = data_in.values

#strip time
time_valid=pd.Timestamp(data_in['valid_time'].values)
time_str=time_valid.strftime('%HZ')
valid_time = time_valid.strftime('%H:%MZ %a %b %d %Y')


# If coords are 1D, make them 2D
if lons.ndim == 1 and lats.ndim == 1:
    lons, lats = np.meshgrid(lons, lats)

#fig settings
fig = plt.figure(figsize=(20, 16))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-91, -80, 34.5, 41], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE, linewidth=1)
ax.add_feature(cfeature.BORDERS, linewidth=1)
ax.add_feature(cfeature.STATES, linewidth=1.25)

### Station Plot Settings ###
stationplot = StationPlot(ax, data['longitude'].values, data['latitude'].values,
                          clip_on=True, transform=ccrs.PlateCarree(), fontsize=15)
stationplot.plot_parameter('NW', tempf , color='red', formatter=lambda v: format(v, '.0f'))
stationplot.plot_parameter('SW', dpt,formatter=lambda v: format(v, '.0f'),color='darkgreen')
stationplot.plot_parameter('NE', data['air_pressure_at_sea_level'].values,
                           formatter=lambda v: format(10 * v, '.0f')[-3:])
stationplot.plot_symbol('C', data['cloud_coverage'].fillna(0), sky_cover)
stationplot.plot_symbol('W', data['current_wx1_symbol'].fillna(0), current_weather)
stationplot.plot_barb(data['eastward_wind'].values, data['northward_wind'].values)

### MRMS Settings ###
#plot mrms data with lowered opacity to see station plot obs if reflec and obs overlap
mesh = ax.pcolormesh(lons,lats,ma.masked_where(refl < 5, refl),cmap='NWSRef',vmin=-1,vmax=80,
                     transform=ccrs.PlateCarree(), alpha=.5)

#cbar settings, set alpha to 1 so the opacity doesnt change for cbar
cb = plt.colorbar(mesh, ax=ax, orientation="horizontal", pad=0.05, aspect=50)
cb.solids.set_alpha(1)
cb.set_label("Reflectivity (dBZ)",fontsize=22)
cb.ax.tick_params(labelsize=16)
#title settings
ax.set_title(f'MRMS & METAR',  fontsize=22, loc='left')
plt.title("Composite Reflectivity\n Surface Observations", fontsize=26)
ax.set_title(f'\nValid: {valid_time}', fontsize=22, loc='right')

out_path=os.path.join(output_dir, f'latest.png')
plt.savefig(out_path,dpi=150,bbox_inches='tight')
plt.close()

with open(os.path.join(output_dir, f'meta.json'), 'w') as f:
    json.dump({"generated_at":datetime.now(timezone.utc).isoformat() + "Z"}, f)
