# Simple example to run a mesas.py model for a 
# tracer transport through a well-mixed system under
# steady flow.

# Date: 04/15/2025
#%%
# imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mesas.sas.model import Model
from mesas.utils import vis
# from scipy.stats import beta, gamma #for beta and gamma distributions

# First need timeseries data
timeseries_duration = 3.
timeseries_length = 300
dt = timeseries_duration/timeseries_length

pulse_start = 0.05
pulse_end = 0.15

C_tracer_input = 0.
Q_steady = 1.
Storage_vol = 0.1

# Create datafram in pandas
data_df = pd.DataFrame(index=np.arange(timeseries_length)*dt)
data_df['Q out [vol/time]'] = Q_steady
data_df['J in [vol/time]'] = Q_steady
data_df['C [conc]'] = 0
data_df['Storage_vol'] = Storage_vol #I needed to add this so Model could access this variable when config.json uses it
data_df.loc[pulse_start:pulse_end, 'C [conc]'] = float(C_tracer_input)
# save data to .csv file
data_df.to_csv('data.csv')

#%%
# data_df.plot()
# plt.show()

#%%
#Now we can run messas.py
model=Model(data_df='data.csv',config='config.json')
model.run()
model.data_df.to_csv('results.csv')
#%%
#Plot the results
data_df = model.data_df
plt.plot(data_df.index, data_df['C [conc]'])
plt.plot(data_df.index,data_df['C [conc] --> Q out [vol/time]'])
plt.show()
#%%
# Plot the breakthrough curve
plt.plot(data_df.index,data_df['C [conc]'])
plt.plot(data_df.index,data_df['C [conc] --> Q out [vol/time]'])
plt.xlim([0,100])
plt.ylim([0,0.51])
#%%
# Plot log vertical axis
plt.plot(data_df.index,data_df['C [conc]'])
plt.plot(data_df.index,data_df['C [conc] --> Q out [vol/time]'])	
plt.yscale('log')
plt.xlim([0,300])
plt.ylim([1E-3,1])
#%%
#Plot using vis
i=100 #i chooses the time step to plot the data
vis.plot_transport_column(model, flux='Q out [vol/time]',sol='C [conc]',
                          i=i, ST_max=0.2)
#can toggle 4*Storage_vol to see more or less
#%%
# Plot the transport column with timeseries
fig=plt.figure(figsize=[18, 6])
i=6
vis.plot_transport_column_with_timeseries(model, flux='Q out [vol/time]',
										   sol='C [conc]', fig=fig, i=i,
										   ST_max=4*Storage_vol)








#%%
#config.json code:
# "sas_specs": {
#  	"Q out [vol/time]":{
# 		"a beta distribution": {
# 		 	"use": "scipy.stats",
# 			"func": "beta",
# 		 	"args" : {
# 		 		"a": 1,
# 		 		"b": 1,
# 		 		"scale": 0.1,
# 		 		"loc": 0
# 		 	}
# 		 }
# 		}
# },

# "sas_specs": {
# 	"Q out [vol/time]":
#                      {"a gamma distribution":
#                           {"scipy.stats": "gamma",
#                            "args": {"a": 1,
#                                     "scale": "Storage_vol",
#                                     "loc": 0},
#                            "nsegment": 200}
# 						   }
# },

# "sas_specs": {
#  	"Q out [vol/time]":{
# 		 "my first SAS func!":{
#  			"ST": [0, "Storage_vol"]
# 		 }
# 		}
# },

#  "options": {
#   "record_state": true
#  }