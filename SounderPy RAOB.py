import sounderpy as spy
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import metpy


site_id='BNA'
year='2026'
month='08'
day='16'
hour=12

obs_data=spy.get_obs_data(site_id,year,month,day,hour)
print('plotting...')
sounding=spy.build_sounding(obs_data, color_blind=False,save=True, filename='sounding')

plt.tight_layout()
plt.show()