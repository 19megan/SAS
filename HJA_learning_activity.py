# This script will upload HJA_data.csv and visualize the tracer timeseries

# Date: 06/11/2025

#%%
import pandas as pd
import numpy as np  
import matplotlib.pyplot as plt

data_df = pd.read_csv('HJA_data.csv', index_col=0, parse_dates=[0])
data_df
# %%
data_df.plot(marker=',',xlabel='[day hr:sec]', ylabel='[mg/L]', title='HJA Tracer Timeseries 08/01/2016', figsize=(12, 4))
# %%
# -------------Clean up the data----------------
# fill missing data for upstream concentration
C_background = 10.96
data_df['Upstream conc. [mg/L]'] = data_df['Upstream conc. [mg/L]'].fillna(C_background)
# downsample from 1s to 1min
data_df_1min = data_df.resample('1min').mean()
# add steady flow upstream and downstream
data_df_1min['Q US [L/s]'] = 11.
data_df_1min['Q DS [L/s]'] = 11.

# account for unrecovered tracer by a tracer recovery factor
f_lost = 0.1 # fraction of tracer input pulse that was not recovered
data_df_1min['Recovered tracer conc. [mg/L]'] = (1-f_lost) * ((data_df_1min['Upstream conc. [mg/L]']) - C_background) + C_background

data_df_1min['C_background'] = float(C_background)
#%%
# create model                          
from mesas.sas.model import Model
model = Model(data_df=data_df_1min, config="HJAconfig.json")
model.run()
#%%
# visualize the output
data_df_1min = model.data_df #don't forget to do this
from mesas.utils import vis
fig = plt.figure(figsize=[11,4])
ax1 = plt.subplot2grid((1,2), (0,0))
vis.plot_SAS_cumulative(model, 'Q DS [L/s]', ax=ax1)

ax2 = plt.subplot2grid((1,2), (0,1))
ax2.plot(data_df_1min.index, data_df_1min['Downstream conc. [mg/L]'], label='Observed tracer outflow')
ax2.plot(data_df_1min.index, data_df_1min['Upstream conc. [mg/L] --> Q DS [L/s]'], label='Predicted tracer outflow')
ax2.legend()

# %%
fig=plt.figure(figsize=[18,6])
i = 20
print(f'Plotting transport column at timestep {i}')
vis.plot_transport_column_with_timeseries(model, flux='Q DS [L/s]', sol='Upstream conc. [mg/L]', i=i, fig=fig, ST_max=40000)
# %%

# Evalaute goodness of fit with RMSE
obs = data_df_1min['Downstream conc. [mg/L]']
pred = (1-f_lost) * data_df_1min['Upstream conc. [mg/L] --> Q DS [L/s]']
RMSE = np.sqrt(np.mean((pred-obs)**2))
print(f'RMSE = {RMSE}')
# %%



# Fitting a distribution
from scipy.stats import beta, gamma

def make_beta_model_from(params): # for beta distribution
    a, b, S_min, S_max, f_lost = params
    sas_specs = {'Q DS [L/s]':
                     {'reach':
                          {'scipy.stats': beta,
                           'args': { 'a': a,
                                     'b': b,
                                 'scale': S_max-S_min,
                                 'loc': S_min },
                           'nsegment': 100}}}
    data_df_1min['Recovered tracer conc. [mg/L]'] = (1-f_lost) * ((data_df_1min['Upstream conc. [mg/L]']) - C_background) + C_background
    solute_parameters = {'Recovered tracer conc. [mg/L]': {'C_old': C_background}}
    return Model(data_df_1min, sas_specs=sas_specs, solute_parameters=solute_parameters, dt=60, influx='Q US [L/s]')
# Note: dt must match the data timestep (60 seconds here)
def make_gamma_model_from(params): # for gamma distribution
    a, S_min, S_max, f_lost = params
    sas_specs = {'Q DS [L/s]':
                     {'reach':
                          {'scipy.stats': gamma,
                           'args': { 'a': a,
                                 'scale': S_max-S_min,
                                 'loc': S_min },
                           'nsegment': 100}}}
    data_df_1min['Recovered tracer conc. [mg/L]'] = (1-f_lost) * ((data_df_1min['Upstream conc. [mg/L]']) - C_background) + C_background
    solute_parameters = {'Recovered tracer conc. [mg/L]': {'C_old': C_background}}
    return Model(data_df_1min, sas_specs=sas_specs, solute_parameters=solute_parameters, dt=60, influx='Q US [L/s]')

# Function that builds SAS model  and returns RMSE
def minimize_me(params):
   model = make_beta_model_from(params) #***edit which distribution to minimize***
   model.run()
   obs = model.data_df['Downstream conc. [mg/L]']
   pred = model.data_df['Recovered tracer conc. [mg/L] --> Q DS [L/s]']
   RMSE = np.sqrt(np.mean((pred-obs)**2))
   print(f'RMSE = {RMSE} for params = {params}')
   return RMSE

#%%
# Now let's supply initial estimates of parameters
a = 2.0
b = 2.0
reach_length = 52. # meters
S_max = 0.22 * reach_length * 1000
S_min = 0.18 * reach_length * 1000
f_lost = 0.1
params_init = a, b, S_min, S_max, f_lost #***edit b for distribution type***
#%%
# Then we can optimize by fmin
from scipy.optimize import fmin
params = fmin(minimize_me, params_init)
#%%
# Optimized parameters: (fmin takes a long time)
params = [ 1.21228858e+03,  8.93484406e+02, -2.23029148e+05,  1.95329691e+05,
        1.68929198e-01] #for beta dist
# params = [2.28521001e+00, 7.29019741e+03, 1.37865627e+04, 8.54532897e-02] #for gamma dist
#%%
# Now build a model with optimized parameters
model = make_beta_model_from(params) #***edit which distribution***
model.run()

# Visualize result
fig=plt.figure(figsize=[11,4])
ax1 = plt.subplot2grid((1,2), (0,0))
ax2 = plt.subplot2grid((1,2), (0,1))
vis.plot_SAS_cumulative(model, 'Q DS [L/s]', ax=ax1)
ax2.plot(model.data_df.index, model.data_df['Downstream conc. [mg/L]'], label='Observed tracer outflow')
ax2.plot(model.data_df.index, model.data_df['Recovered tracer conc. [mg/L] --> Q DS [L/s]'], label='Predicted tracer outflow')
ax2.legend()
# %%
