import sounderpy as spy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sounderpy as spy
import os
import json
import datetime

script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.abspath(os.path.join(script_dir, '..', '..'))
output_dir = os.path.join(repo_root, 'assets', 'soundings')
os.makedirs(output_dir, exist_ok=True)

def get_sounding_time():
    time_now = datetime.datetime.utcnow()
    sounding_time_morning = 12
    sounding_time_evening = 0

    year = time_now.strftime('%Y')
    month = time_now.strftime('%m')
    day = time_now.strftime('%d')

    hour = sounding_time_morning if 12 <= time_now.hour <= 23 else sounding_time_evening

    return hour, year, month, day


hour, year, month, day = get_sounding_time()

sites = ['BNA','ILN','ILX','SGF']
#'ILN','ILX','SGF'
frames = []

for site in sites:
    try:
        print(f'pulling data for {site}')
        obs_data = spy.get_obs_data(site, year, month, day, hour)

        print(f'plotting {site}')
        base = os.path.join(output_dir, f'sounding_{site}')
        spy.build_sounding(obs_data, color_blind=False, save=True,
                           filename=base, radar=None, map_zoom=0, dpi=150)
        plt.close('all')

        saved = f'sounding_{site}.png'
        frames.append({"site": site, "file": saved})
        print(f"Saved sounding for {site}")

    except Exception as e:
        print(f"Skipping {site}: {e}")
        continue

document = {
    "run_time": f'{year}-{month}-{day} {hour:02d}:00',
    "generated_at": datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
    "frames": frames,
}

with open(os.path.join(output_dir, "document.json"), 'w') as f:
    json.dump(document, f, indent=2)