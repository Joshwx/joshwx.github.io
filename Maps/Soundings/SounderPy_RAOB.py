import sounderpy as spy
import matplotlib.pyplot as plt
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

sites = ['BNA', 'ILN']
frames = []

for site in sites:
    try:
        print(f'pulling data for {site}...')
        obs_data = spy.get_obs_data(site, year, month, day, hour)

        print(f'plotting {site}...')
        filename = f'sounding_{site}'
        sounding = spy.build_sounding(obs_data, color_blind=False, save=True, filename=filename)
        plt.tight_layout()
        plt.close()

        frames.append({
            "site": site,
            "file": f'{filename}.png',
        })
        print(f"Saved sounding for {site}")

    except Exception as e:
        print(f"Skipping {site}: {e}")
        continue

document = {"run_time": f'{year}-{month}-{day} {hour:02d}:00',"generated_at": datetime.datetime.utcnow().isoformat() + 'Z',
            "frames": frames}

with open(os.path.join(output_dir, "document.json"), 'w') as f:
    json.dump(document, f, indent=2)

