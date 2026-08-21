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

#pull from model
H= Herbie("2026-05-10 01:00", model='rap', product='awp236pgrb')

print(H.inventory(":700 mb"))

ds=H.xarray(":(?:HGT|UGRD|VGRD|TMP|RH):700 mb:", remove_grib=False)
ds_rh = H.xarray(":RH:(?:700|650|600|550|500) mb:", remove_grib=False)
print(ds_rh)