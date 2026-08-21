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

#pull from model
H= Herbie("2026-05-24 14:00", model='rap', product='awp130pgrb')

#print(H.inventory(":850 mb"))

ds=H.xarray(":(?:HGT|UGRD|VGRD|TMP|RH):850 mb:", remove_grib=False)
ds=ds.metpy.parse_cf()

gh=ds['gh']
u=ds['u']
v=ds['v']
lat=ds['latitude']
lon=ds['longitude']
temp=ds['t']
time=ds['time']
rh=ds['r']
temp_c=ds['t']-273.15
u_kts=ds['u']*1.944
v_kts=ds['v']*1.944

#strip time
time_valid=pd.Timestamp(ds['time'].values)
time_str=time_valid.strftime('%Y%m%d %H%MZ')







