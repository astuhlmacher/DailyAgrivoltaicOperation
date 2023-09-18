# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 10:47:17 2023

"""

import numpy as np
from pvlib import solarposition
import pandas as pd
import matplotlib.pyplot as plt

tz = 'US/Eastern' # Time Zone
lat, lon = 42.3, -83.743  # Latitude, Longitude

times = pd.date_range('2021-07-01 00:00:00', '2021-07-31', freq="30min", tz=tz)
solpos = solarposition.get_solarposition(times, lat, lon)
# remove nighttime
#solpos = solpos.loc[solpos['apparent_elevation'] > 0, :]

spa_data_july = solpos.filter(['azimuth','apparent_zenith'],axis=1)
#solpos.apparent_zenith
#solar_altitude = 90 - solpos.apparent_zenith

july15 = spa_data_july[(spa_data_july.index < '2021-7-16 00:00:00') & (spa_data_july.index >= '2021-7-15 00:00:00')]
july14 = spa_data_july[(spa_data_july.index < '2021-7-15 00:00:00') & (spa_data_july.index >= '2021-7-14 00:00:00')]


filepath1 = 'july2021data_halfhour.csv'
filepath2  = 'july-15-2021data_halfhour.csv'
filepath3  = 'july-14-2021data_halfhour.csv'


spa_data_july.to_csv(filepath1)
#july15.to_csv(filepath2)
#july14.to_csv(filepath3)


# remove nighttime
solpos_day = july14.loc[solpos['apparent_zenith'] < 90, :]

fig, ax = plt.subplots()

for hour in np.unique(solpos_day.index.hour):
    # choose label position by the largest elevation for each hour
    subset = solpos_day.loc[solpos_day.index.hour == hour, :]
    height = 90-subset.apparent_zenith
    pos = solpos_day.loc[height.idxmax(), :]
    ax.text(180-pos['azimuth'], 90-pos['apparent_zenith'], str(hour))
    
date = '2021-7-14'
ax.plot(180-solpos_day.azimuth, 90-solpos_day.apparent_zenith, label=date)

ax.figure.legend(loc='upper left')
ax.set_xlabel('Solar Azimuth (degrees East of South)')
ax.set_ylabel('Solar Altitude (degrees)')

plt.show()