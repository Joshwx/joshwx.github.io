import sounderpy as spy
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import metpy
import pandas as pd
import datetime
from datetime import timedelta

# pull current time
time_now=datetime.datetime.utcnow()
sounding_time_morning=12
sounding_time_evening=00



site_id='BNA'
year=time_now.strftime('%Y')
month=time_now.strftime('%m')
day=time_now.strftime('%d')

if time_now.hour>=12 and time_now.hour<=23:
    hour=sounding_time_morning
else:
    hour=sounding_time_evening


obs_data=spy.get_obs_data(site_id,year,month,day,hour)
print('plotting...')
sounding=spy.build_sounding(obs_data, color_blind=False,save=True, filename='sounding')

plt.tight_layout()
plt.show()

def get_sounding ():
    # pull current time
    time_now = datetime.datetime.utcnow()
    sounding_time_morning = 12
    sounding_time_evening = 00

