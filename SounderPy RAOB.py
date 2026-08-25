import sounderpy as spy
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import metpy
import pandas as pd
import datetime
from datetime import datetime, timedelta

#get utc time for the most recent one
target_time = (pd.Timestamp.now(tz='UTC') - pd.Timedelta(minutes=60)).tz_localize(None).floor('h')

site_id='BNA'
year=target_time.strftime('%Y')
month=target_time.strftime('%m')
day=target_time.strftime('%d')
hour=target_time.strftime('%H')

obs_data=spy.get_obs_data(site_id,year,month,day,hour)
print('plotting...')
sounding=spy.build_sounding(obs_data, color_blind=False,save=True, filename='sounding')

plt.tight_layout()
plt.show()

site_id='ILN'
year=target_time.strftime('%Y')
month=target_time.strftime('%m')
day=target_time.strftime('%d')
hour=target_time.strftime('%H')

obs_data=spy.get_obs_data(site_id,year,month,day,hour)
print('plotting...')
sounding=spy.build_sounding(obs_data, color_blind=False,save=True, filename='sounding')

plt.tight_layout()
plt.show()


site_id='ILX'
year=target_time.strftime('%Y')
month=target_time.strftime('%m')
day=target_time.strftime('%d')
hour=target_time.strftime('%H')

obs_data=spy.get_obs_data(site_id,year,month,day,hour)
print('plotting...')
sounding=spy.build_sounding(obs_data, color_blind=False,save=True, filename='sounding')

plt.tight_layout()
plt.show()

site_id='SGF'
year=target_time.strftime('%Y')
month=target_time.strftime('%m')
day=target_time.strftime('%d')
hour=target_time.strftime('%H')

obs_data=spy.get_obs_data(site_id,year,month,day,hour)
print('plotting...')
sounding=spy.build_sounding(obs_data, color_blind=False,save=True, filename='sounding')

plt.tight_layout()
plt.show()