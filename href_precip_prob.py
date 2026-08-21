import numpy
import sounderpy as spy
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import metpy
import metpy.calc as mpcalc
from cmweather.cm_colorblind import cmap
from metpy.units import units
from herbie import Herbie
from herbie.toolbox import EasyMap, pc
from herbie import paint
import cartopy
from cartopy import crs as ccrs, feature as cfeature
import scipy
from scipy.ndimage import gaussian_filter
import pandas as pd
from matplotlib.colors import ListedColormap
from metpy.plots.ctables import registry
from metpy.plots import ctables


dt = pd.Timestamp("now", tz="utc").floor('12h').replace(tzinfo=None)
H = Herbie(dt, model="href", product="mean", domain="conus", fxx=1)

print(H.inventory(r"APCP:.+:prob\s>12.7:"))

ds = H.xarray(r"APCP:surface:0-1.+:prob\s>12.7:")
print(ds)
