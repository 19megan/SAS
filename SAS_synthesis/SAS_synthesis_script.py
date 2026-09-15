# This script will compute the SAS-derived metrics for the SAS synthesis
# Each synthesis study will be run here to compute the SAS-derived metrics

# Date: 06/24/2026

#%%
# ------------------Import dataset-------------------
from random import random

import pandas as pd
import numpy as np  
import matplotlib.pyplot as plt
from mesas.sas.model import Model
from permetrics.regression import RegressionMetric
from SAS_synthesis.models.functions import make_model, load_data
np.random.seed(42)

# SET FOLDERS-----------------------------------------
data_file_path = "SAS_synthesis/data"
model_file_path = "SAS_synthesis/SAS models"
#-----------------------------------------------------

#CHOOSE CATCHMENT TO RUN-------------------------------
dataset = ['Lower Hafren',
           'Bruntland Burn', #no data yet
           'Providence Creek', #no data yet
           'Weierbach (2019)', #no et data yet
           'Corin', #no data yet
           'Weierbach (2021)', #no et data yet
           'Chenqi', 
           'Can Villa', 
           'Dry Creek', 
           'Selke', 
           'Hemuqiao']

catchment_to_run = dataset[6]

#%% #-----------------------RUN MODEL-----------------------
if catchment_to_run == 'Lower Hafren':
    spin_up = False
    MC = False
    model = make_model("Lower Hafren", data_file_path, spin_up=spin_up, MC=MC)
    model.run()
    df, spinup, df_val, issample, influx, et, discharge, iso_out, iso_in, age_unit = load_data(catchment_to_run, data_file_path)
elif catchment_to_run == 'Bruntland Burn':
    #need data
    pass
elif catchment_to_run == 'Providence Creek':
    #need data
    # spinup 2004-2013 repeated twice
    # simulation period 2004-2016
    pass
elif catchment_to_run == 'Weierbach (2019)':
    spin_up = True
    MC = False
    model = make_model("Weierbach (2019)", data_file_path, spin_up=spin_up, MC=MC)
    model.run()
    df, spinup, df_val, issample, influx, et, discharge, iso_out, iso_in, age_unit = load_data(catchment_to_run, data_file_path)
elif catchment_to_run == 'Corin':
    #need data
    pass
elif catchment_to_run == 'Weierbach (2021)':
    spin_up = True
    MC = False
    model = make_model("Weierbach (2021)", data_file_path, spin_up=spin_up, MC=MC)
    model.run()
    df, spinup, df_val, issample, influx, et, discharge, iso_out, iso_in, age_unit = load_data(catchment_to_run, data_file_path)
elif catchment_to_run == 'Chenqi':
    spin_up = True
    MC = True
    model = make_model("Chenqi", data_file_path, spin_up=spin_up, MC=MC, validate=True)
    model.run()
    df, spinup, df_val, issample, influx, et, discharge, iso_out, iso_in, age_unit = load_data(catchment_to_run, data_file_path)
    issample = (model.data_df['outlet-D (‰)'].notna()) & (model.data_df['Q (mm/day)']>0)
elif catchment_to_run == 'Can Villa':
    spin_up = True
    MC = False
    model = make_model("Can Villa", data_file_path, spin_up=spin_up, MC=MC, validate=True)
    model.run()
    df, spinup, df_val, issample, influx, et, discharge, iso_out, iso_in, age_unit = load_data(catchment_to_run, data_file_path)
    issample = (model.data_df['measC_Q [-]'].notna()) & (model.data_df['Q [mm/h]']>0)
elif catchment_to_run == 'Dry Creek':
    spin_up = True
    MC = True
    model = make_model("Dry Creek", data_file_path, spin_up=spin_up, MC=MC, validate=True)
    model.run()
    df, spinup, df_val, issample, influx, et, discharge, iso_out, iso_in, age_unit = load_data(catchment_to_run, data_file_path)
    issample = (model.data_df['measC_Q [-]'].notna()) & (model.data_df['Q [mm/4h]']>0)
elif catchment_to_run == 'Selke':
    spin_up = True
    MC = True
    model = make_model("Selke", data_file_path, spin_up=spin_up, MC=MC)
    model.run()
    df, spinup, df_val, issample, influx, et, discharge, iso_out, iso_in, age_unit = load_data(catchment_to_run, data_file_path)
elif catchment_to_run == 'Hemuqiao':
    #need data
    pass
else:
    raise ValueError(f"Unknown catchment: {catchment_to_run}. Please provide a valid catchment.")






# %% ----------------------Visualize results---------------------------------
from mesas.utils import vis
fig = plt.figure(figsize=[14,4])
plt.plot(model.data_df.index[issample], model.data_df[iso_out][issample],'.', color='grey', label=f'Observed {iso_out} outflow')
plt.plot(model.data_df.index[issample], model.data_df[f'{iso_in} --> {discharge}'][issample], '.', color='orange', alpha=0.9,label=f'Predicted {iso_out} outflow')
# plt.axvspan(pd.Timestamp('1993-01-01'), pd.Timestamp('2003-01-01'), color='lightgrey', alpha=0.2)
plt.legend()
plt.title(f'Isotope outflow at {catchment_to_run}')
# %%
# check accuracy
from permetrics.regression import RegressionMetric
import hydroeval as he
obs = model.data_df[iso_out][issample].to_numpy()
pred = model.data_df[f'{iso_in} --> {discharge}'][issample].to_numpy()
nse = he.evaluator(he.nse, pred, obs)
print(f'NSE = {nse[0]}')
RMSE = np.sqrt(np.mean((pred-obs)**2))
print(f'RMSE = {RMSE}')

#%%
# Plot TTD
import matplotlib.cm as cm
cmap = plt.get_cmap('viridis')
colors = [cmap(i) for i in np.linspace(0,1,len(model.data_df))]
pq = model.get_pQ(discharge)
T = model.options['dt'] *np.arange(model.options['max_age'])
PQ = np.cumsum(pq, axis=0) * model.options['dt']

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[11,4])
for i in range(0, len(model.data_df)):
    ax[0].plot(T, PQ[:,i], color=colors[i])
    ax[1].plot(T, pq[:,i], color=colors[i])

ax[0].set_ylim([0, 1.1])
ax[0].set_xlim(xmin=0)
ax[0].axhline(1, color='0.1', lw=0.8, ls=':')
ax[0].axhline(0, color='0.1', lw=0.8, ls=':')
ax[0].set_xlabel(f'Age ({age_unit})')
ax[0].set_ylabel('$P_Q$')
ax[0].set_title('Cumulative TTD')
ax[1].set_xlabel(f'Age ({age_unit})')
ax[1].set_ylabel('$p_Q$')
ax[1].set_title('TTD')
sm = cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0,
vmax=len(model.data_df)-1))
sm.set_array([])
fig.colorbar(sm, ax=ax, label='Time Index', fraction=0.046, pad=0.04)

# %%

#-----------------------SAS-derived metrics-----------------------
# 1. YW runoff ratio
pq = model.get_pQ(discharge)
PQ = np.cumsum(pq, axis=0) * model.options['dt']
YWq = PQ[90, :]*model.data_df[discharge]
P = model.data_df[influx]
YW_RR = YWq / P
RR = model.data_df[discharge]/model.data_df[influx]

plt.figure(figsize=[14,4])
plt.plot(YW_RR, label='YW runoff ratio')
plt.plot(RR, alpha=0.5, zorder=0.5, label='Total runoff ratio')
plt.title(f'YW (Age <=90 {age_unit}) Runoff Ratio')
plt.xlabel('Time')
plt.ylabel('Runoff Ratio')
plt.legend()
plt.show()

# %%
# 2. YW discharge vs (inverse) wetness (ISE)
top = PQ[90, :]
bottom = (model.data_df[discharge].max()-model.data_df[discharge])/(model.data_df[discharge].max()-model.data_df[discharge].min())
ISE = top/bottom

 #Plot ISE over time
plt.figure(figsize=[14,4])
plt.plot(ISE, label='Storage Effect (< 1: Direct Storage Effect, > 1: Inverse Storage Effect)')
plt.title('Inverse Storage Effect (>1), Direct Storage Effect (<1)')
plt.xlabel('Time')
plt.ylabel('ISE')
plt.ylim([0, 2])
plt.axhline(1, color='0.1', lw=0.8, ls=':', label='No Storage Effect')
plt.legend()
plt.show()
#%%
# plot YWF vs wetness
plt.figure(figsize=[14,4])
plt.plot(1-bottom, top, '.', label='>1:1 Inverse Storage Effect, <1:1 Direct Storage Effect')
plt.title(f'YWF (Age <=90 {age_unit}) of discharge vs (inverse) Wetness')
plt.xlabel('Wetness (function of Q) (1=High wetness, 0=Low wetness)')
plt.ylabel(f'YWF (Age <=90 {age_unit})')
plt.legend()
plt.show()
#%%
# plot YWF vs dryness
plt.figure(figsize=[14,4])
plt.plot(bottom, top, '.', label='slope<1: Inverse Storage Effect, slope>1: Direct Storage Effect')
plt.title(f'YWF (Age <=90 {age_unit}) of discharge vs (inverse) Wetness')
plt.xlabel('Dryness (function of Q) (0=High wetness, 1=Low wetness)')
plt.ylabel(f'YWF (Age <=90 {age_unit})')
plt.plot(bottom, 1-bottom, '--', color='r', alpha=.3, label='Inverse Storage Effect')
plt.plot(bottom, bottom, '--', color='b', alpha=.3, label='Direct Storage Effect')
plt.legend()
plt.show()


# %%
# 3. Average age of discharge vs storage
mean_pq = np.mean(pq, axis=0)
abs_diff = np.abs(pq-mean_pq)
idx = np.argmin(abs_diff, axis=0)
idx = idx.reshape(-1,1)
storage = model.data_df[influx]-model.data_df[discharge]-model.data_df[et]

plt.figure(figsize=[14,4])
plt.plot(T[idx], storage, '.', label='Average Age of Discharge')
plt.title('Average Age of Discharge vs Storage')
plt.xlabel('Average Age of Discharge')
plt.ylabel('Storage')
plt.legend()
plt.show()
# %%
# 4. Age-ranked storage < average annual precip vs catchment storage
ST = model.get_sT()
aap = model.data_df[influx].groupby(model.data_df.index.year).sum().mean()
i=(ST>=aap).argmax(axis=0) #takes the first index that meets criteria
r = PQ[i[1:], np.arange(0,len(ST))]

# storage = df_harman['S_scale']/(-103)+48 # for harman data
storage = model.data_df[influx]-model.data_df[discharge]-model.data_df[et] # for benettin data

# plt.figure(figsize=[12,3])
plt.plot(storage,r, '.')
plt.ylabel(f'Age-ranked storage < average annual rainfall ({aap:.1f} mm)')
plt.xlabel('storage')
# plt.ylim([0, .4])
plt.title(f'Age-ranked storage < average annual rainfall vs catchment storage at {catchment_to_run}')

print(f'Average fraction of discharge from storage < average annual rainfall volume: {np.mean(r):.3f}')
print('Tells you how much discharge is coming from storage added within the average year')
# %%
