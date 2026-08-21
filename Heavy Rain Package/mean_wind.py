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

from basic_severe_wx import speed_850

#pull from model
H = Herbie("2026-08-09 17:00", model='rap', product='awp236pgrb', fxx=2)

ds = H.xarray(":(?:UGRD|VGRD):\\d+ mb:", remove_grib=False)
print(ds)
