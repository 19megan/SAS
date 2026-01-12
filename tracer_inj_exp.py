# This script will run a tracer injection of 
# salt to understand how much of the alluvial
# aquifer is exchanging with the section of a stream.

# Date: 06/05/2025

#%%
# imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mesas.sas.model import Model
from mesas.utils import vis
from scipy.stats import beta, gengamma

#%%
timeseries_duration = 2 # duration in hours
timeseries_length = 300 # number of timesteps in the simulation
dt = timeseries_duration/timeseries_length

C_background = 10.96 # background concentration in mg/L
pulse_start = 0.05 # hr
#pulse_end = pulse_start+timeseries_length/2*dt
pulse_end = pulse_start + dt*2  # tracer released over two timesteps

Q_steady = 11./1000.*(3600) #m3/hr

M_tracer_input = 6000 # mg  - must be >11880mg background for 3hrs injection and <339000mg solubility limit
C_tracer_input = M_tracer_input / (Q_steady * (pulse_end - pulse_start))

print(f'      Input pulse duration: {(pulse_end - pulse_start)*60} min')
print(f'Tracer input concentration: {C_tracer_input} mg/L')
print(f'          Solubility limit: 350000 mg/L')

data_df = pd.DataFrame(index=np.arange(timeseries_length) * dt)

data_df['Q out [vol/time]'] = Q_steady
data_df['J in [vol/time]'] = Q_steady

data_df['C [conc]'] = C_background
data_df.loc[pulse_start:pulse_end, 'C [conc]'] = C_tracer_input

data_df_b = data_df.copy() # a copy for the beta model
data_df_g = data_df.copy() # a copy for the gamma model
#%%
# To specify the beta distribution
vars_b = {"S_min":12.3/4,"a_b":1.3,"b_b":1.9, "S_b":29.2/4}
data_df_b = data_df_b.assign(**vars_b)

# To specify the generalized gamma distribution
vars_g = {"S_min":12.3/4, "a_g":0.562, "c_g":2.2, "S_g":20.0/4}
data_df_g = data_df_g.assign(**vars_g)

# Run mesas.py
model_b = Model(data_df=data_df_b, config='config_beta.json')
model_b.run()
model_g = Model(data_df=data_df_g, config='config_gamma.json')
model_g.run()
#%%
# Plot the breakthrough curve
# plt.figure(1)
# plt.plot(model_b.data_df.index, model_b.data_df['C [conc]'])
# plt.plot(model_b.data_df.index, model_b.data_df['C [conc] --> Q out [vol/time]'])
# plt.title('Beta Model')
# plt.xlim([0,.5])
# plt.ylim([0,70000])

# plt.figure(2)
# plt.plot(model_g.data_df.index, model_g.data_df['C [conc]'])
# plt.plot(model_g.data_df.index, model_g.data_df['C [conc] --> Q out [vol/time]'])
# plt.title('Gamma Model')
# plt.xlim([0,.5])
# plt.ylim([0,70000])
# plt.show()

plt.figure(3)
plt.plot(model_b.data_df.index, model_b.data_df['C [conc] --> Q out [vol/time]'], label='Beta Model')
plt.plot(model_g.data_df.index, model_g.data_df['C [conc] --> Q out [vol/time]'], label='Gamma Model')
plt.title('Tracer Input Comparison')
plt.xlim([0,.5])
plt.legend()
plt.show()
# %%
