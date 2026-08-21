import matplotlib.pyplot as plt
import urllib.request

from metpy.io import parse_wpc_surface_bulletin
from io import BytesIO
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from metpy.plots import ColdFront, WarmFront, OccludedFront, StationaryFront
import numpy as np
import geopandas as gpd
url='https://www.wpc.ncep.noaa.gov/discussions/codsus'

with urllib.request.urlopen(url) as response:
    content=response.read()

print(content)

df= parse_wpc_surface_bulletin(BytesIO(content))
print(df)

df['feature'].unique()
print(df['feature'].unique())

map_crs = ccrs.LambertConformal(central_latitude=40, central_longitude=-85)

fig = plt.figure(figsize=(15, 15), dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=map_crs)
ax.set_extent((-95, -75, 34, 48))
ax.add_feature(cfeature.STATES)
ax.add_feature(cfeature.BORDERS)

# set lat/lon bounds to match the map extent for filtering
lon_min, lon_max, lat_min, lat_max = -95, -75, 34, 48
#include if statement that masks out any features outside bounds
features = df[df['feature'] == 'LOW']
for f in features['geometry']:
    if lon_min <= f.x <= lon_max and lat_min <= f.y <= lat_max:
        ax.text(f.x, f.y, 'L', transform=ccrs.PlateCarree(), color='red',
                 fontsize=25, fontweight='bold')

features = df[df['feature'] == 'HIGH']
for f in features['geometry']:
    if lon_min <= f.x <= lon_max and lat_min <= f.y <= lat_max:
        ax.text(f.x, f.y, 'H', transform=ccrs.PlateCarree(), color='blue',
                 fontsize=25, fontweight='bold')
s=10
feature_names=['WARM','COLD','STNRY','OCFNT','TROF']
feature_styles=[{'linewidth':1,'path_effects': [WarmFront(size=s)]},
                {'linewidth':1,'path_effects': [ColdFront(size=s)]},
                {'linewidth': 1, 'path_effects': [StationaryFront(size=s)]},
                {'linewidth': 1, 'path_effects': [OccludedFront(size=s)]},
                {'linewidth': 2, 'linestyle': 'dashed', 'edgecolor': 'darkorange'},
                ]


for name, style in zip(feature_names, feature_styles):
    f=df[df['feature']==name]
    ax.add_geometries(f.geometry, crs=ccrs.PlateCarree(), **style, facecolor='none')
plt.tight_layout()
plt.show()
