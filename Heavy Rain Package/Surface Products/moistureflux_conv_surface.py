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
import geocat.comp as gc
#pull from model
H= Herbie("2026-05-24 14:00", model='rap', product='awp130pgrb')
ds=H.xarray(":SPFH:2 m above ground", remove_grib=False)
ds2= H.xarray(":[UV]GRD:10 m above ground:", remove_grib=False)

print(ds)
sh=ds['sh2']
u=ds2['u10']
v=ds2['v10']
u_kts=ds2['u10']*1.944
v_kts=ds2['v10']*1.944
lat = ds['latitude']
lon = ds['longitude']

#mfc=-u dq/dx-v dq/dy-q(du/dx+dv/dy)
dqdx = mpcalc.first_derivative(sh, axis=1, x=lon)
dqdy = mpcalc.first_derivative(sh, axis=0, x=lat)
dudx = mpcalc.first_derivative(u, axis=1, x=lon)
dvdy = mpcalc.first_derivative(v, axis=0, x=lat)

advection=(-u*dqdx)-(v*dqdy)
conv=sh*np.add(dudx,dvdy)
mfc=advection-conv

print(mfc)