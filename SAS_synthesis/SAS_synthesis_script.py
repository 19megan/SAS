# This script will compute the SAS-derived metrics for the SAS synthesis
# Each synthesis study will be run here to compute the SAS-derived metrics

# Date: 06/24/2026

#%%
# ------------------Import dataset-------------------
import sys
from pathlib import Path
# find the repo root from wherever this is run, so cwd does not matter
ROOT = next(p for p in [Path.cwd(), *Path.cwd().parents] if (p/'SAS_synthesis'/'data').is_dir())
sys.path.insert(0, str(ROOT))

from random import random

import pandas as pd
import numpy as np  
import matplotlib.pyplot as plt
from mesas.sas.model import Model
from permetrics.regression import RegressionMetric
from SAS_synthesis.models.functions import make_model, load_data, save_model_run, load_model_run
np.random.seed(42)

# SET FOLDERS-----------------------------------------
data_file_path = str(ROOT/'SAS_synthesis'/'data')
model_file_path = str(ROOT/'SAS_synthesis'/'SAS models')
model_run_path = ROOT/'SAS_synthesis'/'model_runs' #saved model runs (*.pkl, not tracked by git)
#-----------------------------------------------------

#CHOOSE CATCHMENT TO RUN-------------------------------
dataset = ['Lower Hafren', #need spinup?
           'Bruntland Burn', #done
           'Providence Creek', #no data yet
           'Weierbach (2019)', #Rodriguez and Klaus (2019): composite SAS, d2H, Sref calibrated
           'Corin', #no data yet
           'Weierbach (2021)', #Rodriguez et al. (2021): same model with Sref = 2000 mm; 3H needs precip 3H data
           'Chenqi', 
           'Can Villa', 
           'Dry Creek', 
           'Selke', #done
           'Hemuqiao']

catchment_to_run = dataset[3]

# SAVE/RELOAD MODEL RUN---------------------------------
load_saved_run = False #True: reload the saved run for this catchment instead of running the model
save_run = True #save the model after running it (ignored when reloading)
run_file = model_run_path/f'{catchment_to_run}.pkl'

#%%
# #-----------------------LOAD/RUN MODEL-----------------------
if load_saved_run:
    df, spinup, df_val, issample, influx, et, discharge, iso_out, iso_in, age_unit = load_data(catchment_to_run, data_file_path)
    run = load_model_run(run_file)
    model, issample = run['model'], run['issample']

elif catchment_to_run == 'Lower Hafren':
    spin_up = False
    MC = False
    model = make_model("Lower Hafren", data_file_path, spin_up=spin_up, MC=MC)
    model.run()
    df, spinup, df_val, issample, influx, et, discharge, iso_out, iso_in, age_unit = load_data(catchment_to_run, data_file_path)
elif catchment_to_run == 'Bruntland Burn':
    spin_up = True
    MC = True
    model = make_model("Bruntland Burn", data_file_path, spin_up=spin_up, MC=MC)
    model.run()
    df, spinup, df_val, issample, influx, et, discharge, iso_out, iso_in, age_unit = load_data(catchment_to_run, data_file_path)
    issample = (model.data_df['measC_Q d2H'].notna()) & (model.data_df['Q [mm/h]']>0)
elif catchment_to_run == 'Providence Creek':
    #need data
    # spinup 2004-2013 repeated twice
    # simulation period 2004-2016
    pass
elif catchment_to_run == 'Weierbach (2019)':
    spin_up = True
    MC = True#False
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


# Save model run (reload it later with load_saved_run = True)
if save_run and not load_saved_run:
    save_model_run(run_file, model, issample, catchment=catchment_to_run, spin_up=spin_up, MC=MC)



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
pq = model.get_pQ(discharge)
T = model.options['dt'] *np.arange(model.options['max_age'])
PQ = np.cumsum(pq, axis=0) * model.options['dt']
# one column per RECORDED timestep, which is not every timestep when the state is recorded sparsely
# (Bruntland Burn records daily, see make_SAS_model), so take the count from the array
n_rec = pq.shape[1]
colors = [cmap(i) for i in np.linspace(0,1,n_rec)]

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[11,4])
for i in range(0, n_rec):
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
vmax=n_rec-1))
sm.set_array([])
fig.colorbar(sm, ax=ax, label='Time Index', fraction=0.046, pad=0.04)

# %%

#-----------------------SAS-derived metrics-----------------------
# The state arrays carry one column per RECORDED timestep, which is not necessarily every timestep,
# so line them up with the rows of data_df they belong to. This is a no-op for the catchments that
# record every step; Bruntland Burn runs hourly and records daily (see make_SAS_model), so it has
# 1/24 as many columns as rows.
rec = getattr(model, '_index_ts', np.arange(len(model.data_df)))
rec_df = model.data_df.iloc[rec]
# the young water threshold is 90 DAYS, so convert it to timesteps instead of assuming a daily model
dt_hours = (model.data_df.index[1]-model.data_df.index[0])/pd.Timedelta('1h')
YW = int(np.round(90*24/dt_hours)) #90 for a daily model, 540 at 4h, 2160 at hourly
age_90d = f'Age <=90 day' #label, was the same as 'Age <=90 {age_unit}' only for the daily models

# 1. YW runoff ratio
pq = model.get_pQ(discharge)
PQ = np.cumsum(pq, axis=0) * model.options['dt']
YWq = PQ[YW, :]*rec_df[discharge]
P = rec_df[influx]
YW_RR = YWq / P
RR = rec_df[discharge]/rec_df[influx]

plt.figure(figsize=[14,4])
plt.plot(YW_RR, label='YW runoff ratio')
plt.plot(RR, alpha=0.5, zorder=0.5, label='Total runoff ratio')
plt.title(f'YW ({age_90d}) Runoff Ratio')
plt.xlabel('Time')
plt.ylabel('Runoff Ratio')
plt.legend()
plt.show()

# %%
# 2. YW discharge vs (inverse) wetness (ISE)
top = PQ[YW, :]
bottom = (rec_df[discharge].max()-rec_df[discharge])/(rec_df[discharge].max()-rec_df[discharge].min())
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
plt.title(f'YWF ({age_90d}) of discharge vs (inverse) Wetness')
plt.xlabel('Wetness (function of Q) (1=High wetness, 0=Low wetness)')
plt.ylabel(f'YWF ({age_90d})')
plt.legend()
plt.show()
#%%
# plot YWF vs dryness
plt.figure(figsize=[14,4])
plt.plot(bottom, top, '.', label='slope<1: Inverse Storage Effect, slope>1: Direct Storage Effect')
plt.title(f'YWF ({age_90d}) of discharge vs (inverse) Wetness')
plt.xlabel('Dryness (function of Q) (0=High wetness, 1=Low wetness)')
plt.ylabel(f'YWF ({age_90d})')
plt.plot(bottom, 1-bottom, '--', color='r', alpha=.3, label='Inverse Storage Effect')
plt.plot(bottom, bottom, '--', color='b', alpha=.3, label='Direct Storage Effect')
plt.legend()
plt.show()

#%%
#check exponential fit
from scipy.optimize import curve_fit
def exponential_model(x, a, b, c):
    """
    a: scale factor / initial value
    b: growth/decay rate constant
    c: vertical offset (shift)
    """
    return a * np.exp(b * x) + c

popt, pcov = curve_fit(exponential_model, bottom, top, p0=[1,0,1])
# Extract optimized parameters
a_opt, b_opt, c_opt = popt
print(f"Fitted exponential parameters:\na = {a_opt:.3f}\nb = {b_opt:.3f}\nc = {c_opt:.3f}")
b_smooth = np.linspace(min(bottom), max(bottom), 100)
t_smooth = exponential_model(b_smooth, a_opt, b_opt, c_opt)
plt.figure(figsize=[14,4])
plt.plot(bottom, top, '.', label='Data points')
plt.plot(b_smooth, t_smooth, color='black', lw=2, label=f'a={a_opt:.2f}, b={b_opt:.2f}, c={c_opt:.2f}')
plt.title(f'YWF ({age_90d}) of discharge vs (inverse) wetness exponential fit')
plt.xlabel('Dryness (function of Q) (0=High wetness, 1=Low wetness)')
plt.ylabel(f'YWF ({age_90d})')
plt.ylim([0, 1.1])
plt.legend()
plt.show()
#%%
#find slope of exponential curve
slope, int = np.polyfit(b_smooth, t_smooth, 1)
print(f'average slope of exponential fit = {slope:.3f}')
#find dFY/dQ
slope, int = np.polyfit(rec_df[discharge], top*rec_df[discharge], 1)
print(f'average dFYWQ/dQ = {slope:.3f}')


# %%
# 3. Average age of discharge vs storage
mean_pq = np.mean(pq, axis=0)
abs_diff = np.abs(pq-mean_pq)
idx = np.argmin(abs_diff, axis=0)
idx = idx.reshape(-1,1)
storage = rec_df[influx]-rec_df[discharge]-rec_df[et]

plt.figure(figsize=[14,4])
plt.plot(T[idx], storage, '.', label='Average Age of Discharge')
plt.title('Average Age of Discharge vs Storage')
plt.xlabel('Average Age of Discharge')
plt.ylabel('Storage')
plt.legend()
plt.show()
# %%
# 4. Age-ranked storage < average annual precip vs catchment storage
ST = model.get_sT() # (max_age, n_recorded+1), column 0 is the initial condition
ST = np.cumsum(ST, axis=0) * model.options['dt']
aap = model.data_df[influx].groupby(model.data_df.index.year).sum().mean()
i=(ST>=aap).argmax(axis=0) #takes the first index that meets criteria
# i[1:] drops the initial condition, so it is one entry per recorded timestep and pairs up with the
# columns of PQ. Indexing the columns by len(ST) instead would use max_age, which is the AGE axis.
r = PQ[i[1:], np.arange(0,PQ.shape[1])]


storage = (rec_df[influx]-rec_df[discharge]-rec_df[et])#.cumsum()

plt.figure(figsize=[12,3])
plt.plot(storage,r, '.')
plt.ylabel(f'Fraction of discharge from storage < average \n annual rainfall ({aap:.1f} mm)')
plt.xlabel('storage')
# plt.ylim([0, .4])
plt.title(f'Fraction of discharge from storage < average annual rainfall \n vs catchment storage at {catchment_to_run}')

print(f'Average fraction of discharge from storage < average annual rainfall volume: {np.mean(r):.3f}')
print('Tells you how much discharge is coming from storage added within the average year')
# %%
