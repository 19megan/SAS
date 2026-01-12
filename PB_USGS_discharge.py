#This script will read the USGS discharge data from a csv file obtained on HydroShare
import pandas as pd
import numpy as np
import statistics as stats
import matplotlib.pyplot as plt
from mesas.sas.model import Model
from mesas.utils import vis

# Load the USGS discharge data
usgs_data = pd.read_csv('nwisdv-pond_branch_at_oregon_ridge__md-discharge__cubic_feet_per_second.csv', skiprows=2)
usgs_data = usgs_data[['UTCTimeStamp', 'Value']]
print(usgs_data.head())
print(usgs_data[1367:1371])# The data skips from 1986 to 1998

# Print some statistics about the data
print('\n4.\nmin:', usgs_data['Value'].min())
print('max:', usgs_data['Value'].max())
print('standard deviation:', stats.stdev(usgs_data['Value']))
print('length:', len(usgs_data))

# Plot the USGS discharge data
usgs_data['UTCTimeStamp'] = pd.to_datetime(usgs_data['UTCTimeStamp'])
usgs_data.set_index('UTCTimeStamp', inplace=True)
usgs_data.plot()
plt.title('USGS Discharge Data')
plt.xlabel('Time')
plt.ylabel('Discharge (cubic feet per second)')
plt.show()