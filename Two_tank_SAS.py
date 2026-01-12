# This script will find the analytical SAS function for a two tank system
#  based on analytical approach in analytical_SAS.py

# This script needs a given SAS function form, J, Q1, Q2, and k1, k2

# Date: 10/13/2025

#%% ------------------Create J & Q1, Q2-------------------
import pandas as pd
import numpy as np  
import matplotlib.pyplot as plt
from mesas.sas.model import Model

date = pd.date_range(start='2018-01-01 00:00:00', end='2020-12-31 23:00:00', freq='h')
df = pd.DataFrame(index=date)
# print(len(date)) # 52608 hours

# Define J by random Poisson process and gamma distribution of rates
# Each rain event has a random start time, duration, and intensity
def simulate_rainfall_rate(
    n_hours=24*30,         # number of hours to simulate
    lambda_per_hour=0.007,   # mean number of rain events per hour (~5 per month)
    mean_intensity=3.0,    # mean intensity of each event (mm/hr)
    duration_mean=3.0,     # mean duration of each event (hours)
    start_time="2020-01-01 00:00",
    random_seed=None
):
    """
    Simulate rainfall rate (mm/hr) as a Poisson process of rain events.
    Each event has a random start time, duration, and intensity.
    Returns a pandas DataFrame with hourly rainfall rate.
    """
    if random_seed is not None:
        np.random.seed(random_seed)

    # Step 1: number of events expected
    total_events = np.random.poisson(lam=lambda_per_hour * n_hours)
    
    # Step 2: event start times (uniform over simulation period)
    start_times = np.random.uniform(0, n_hours, total_events)
    
    # Step 3: random duration and intensity per event
    durations = np.random.exponential(scale=duration_mean, size=total_events)
    intensities = np.random.exponential(scale=mean_intensity, size=total_events)
    
    # Step 4: create time grid and initialize rate array
    t = np.arange(n_hours)
    rainfall_rate = np.zeros_like(t, dtype=float)
    
    # Step 5: for each event, add its intensity over its duration
    for s, d, i in zip(start_times, durations, intensities):
        start_idx = int(np.floor(s))
        end_idx = int(np.floor(s + d))
        end_idx = min(end_idx, n_hours - 1)
        rainfall_rate[start_idx:end_idx+1] += i
    
    return rainfall_rate

df['J (mm/hr)'] = simulate_rainfall_rate(n_hours=len(date), lambda_per_hour=0.007, mean_intensity=3, random_seed=42)

# Plot simulated rainfall
plt.figure(figsize=(10,4))
plt.plot(df.index, df['J (mm/hr)'])
plt.ylabel("Rainfall rate (mm/hr)")
plt.xlabel("Time")
plt.title("Simulated rainfall rate from Poisson process")
plt.tight_layout()
plt.show()

# Define Q by a simple linear reservoir model with random K
def simulate_discharge(
    rainfall_rate,  # array or series of rainfall rate [mm/hr]
    dt=1.0,         # time step in hours
    K=20.0,         # reservoir time constant (hours) time to empty
    S0=0.0          # initial storage
):
    """
    Simulate discharge (mm/hr) from rainfall (mm/hr)
    using a simple linear reservoir model.
    """
    n = len(rainfall_rate)
    S = np.zeros(n+1)
    Smean = np.zeros(n)
    Q = np.zeros(n)

    S[0] = S0
    for t in range(0, n):
        Q[t] = S[t] / K
        S[t+1] = S[t] + (rainfall_rate[t] - Q[t]) * dt
        S[t+1] = max(S[t+1], 0)  # prevent negative storage
        Smean[t] = (S[t] + S[t+1]) / 2

    return Q, Smean
# Q1 = 1/k1*S1
# Q2 = 1/k2*S2
# Convert to discharge
k1 = 500 # reservoir time constant (hours in tank)
k2= 1000
S0_1 = (df['J (mm/hr)']).mean()*k1
print(S0_1)
print(0.105*k1)
Q1, Sm1 = simulate_discharge(df['J (mm/hr)'].values, dt=1, K=k1, S0=S0_1)
df['Q1 (mm/hr)'] = Q1
df['Sm1 (mm)'] = Sm1
S0_2 = df['Q1 (mm/hr)'].mean()*k2
print(S0_2)
Q2, Sm2 = simulate_discharge(df['Q1 (mm/hr)'].values, dt=1, K=k2, S0=S0_2) #Q1 input should not be scaled by C as in the simulate_discharge function, so scale it back
df['Q2 (mm/hr)'] = Q2
df['Sm2 (mm)'] = Sm2
df['Sm12 (mm)'] = df['Sm1 (mm)'] + df['Sm2 (mm)']
    
# Plot rainfall vs discharge
plt.figure(figsize=(10,4))
plt.plot(df.index, df["Q1 (mm/hr)"], color="darkgreen", label="Discharge 1 (mm/hr)")
plt.plot(df.index, df["Q2 (mm/hr)"], color="lightgreen", alpha=0.75, label="Discharge 2 (mm/hr)")
plt.ylabel("Discharge (mm/hr)")
plt.title("Simulated discharge response to Poisson rainfall")
plt.legend()
plt.tight_layout()
plt.show()

# Resample to daily
df = df.resample('D').mean()
dt = 24 # hour in a day
# df['J (mm/hr)'] = df['J (mm/hr)']*.6 # since runoff coefficient is 0.6

#%%
# -------------------Create staggered tracer input-------------------
from scipy import sparse
tracers = [f'C{i}' for i in range(1,len(df)+1)] # list of [C1, C2, C3, ...]
df = pd.concat([df, pd.DataFrame.sparse.from_spmatrix(sparse.eye(len(df), dtype=np.float32, format='csr'), columns=tracers, index=df.index)], axis=1)

#%%
# ------------------Fit a distribution-------------------
from scipy.stats import beta, gamma

def make_uniform_model1_from(params): # for uniform distribution                          
    S1 = params
    tracers = [f'C{i}' for i in range(1,len(df)+1)]
    sas_specs = {'Q1 (mm/hr)':
                        {'tank1':
                          {"ST": [0, 'Sm1 (mm)'],
                           "P": [0.0, 1.0]}
                        }
                        }
    solute_parameters = {tracers[i]:{} for i in range(len(tracers))}
    return Model(df, sas_specs=sas_specs, solute_parameters=solute_parameters, dt=24, influx='J (mm/hr)')#,n_substeps=1, record_state=True)

def make_uniform_model2_from(params): # for uniform distribution                          
    S2 = params
    tracers = [f'C{i}' for i in range(1,len(df)+1)]
    sas_specs = {'Q2 (mm/hr)':
                        {'tank2':
                          {"ST": [0, 'Sm2 (mm)'],
                           "P": [0.0, 1.0]}
                        }
                        }
    solute_parameters = {tracers[i]:{} for i in range(len(tracers))}
    return Model(df2, sas_specs=sas_specs, solute_parameters=solute_parameters, dt=24, influx='Q1 (mm/hr)')#,n_substeps=1, record_state=True)

def make_uniform_model3_from(params): # for uniform distribution                          
    S1 = params
    tracers = [f'C{i}' for i in range(1,len(df)+1)]
    sas_specs = {'Q2 (mm/hr)':
                        {'tank12':
                          {"ST": [0, 'Sm12 (mm)'],
                           "P": [0.0, 1.0]}
                        }
                        }
    solute_parameters = {tracers[i]:{'observations':f'obs Q2 for C{i+1}'} for i in range(len(tracers))}
    return Model(df3, sas_specs=sas_specs, solute_parameters=solute_parameters, dt=24, influx='J (mm/hr)')#,n_substeps=1, record_state=True)

def make_beta_model_from(params): # for beta distribution
    a, b = params
    tracers = [f'C{i}' for i in range(1,len(df)-1)]
    sas_specs = {'Q (mm/hr)':
                        {'ORPB':
                          {'func': "beta",
                           'args': { 'a': a,
                                     'b': b,
                                 'scale': 'abs_storage (mm)',
                                 'loc': 0 },
                           'nsegment': 100}}}
    solute_parameters = {tracers[i]:{} for i in range(len(tracers))}
    return Model(df, sas_specs=sas_specs, solute_parameters=solute_parameters, dt=1, influx='J (mm/hr)')
# Note: dt must match the data timestep (1 hour here)
def make_gamma_model1_from(params): # for gamma distribution
    a = params
    tracers = [f'C{i}' for i in range(1,len(df)+1)]
    sas_specs = {'Q1 (mm/hr)':
                     {'tank1':
                          {'func': 'gamma',
                           'args': { 'a': 1,
                                 'scale': 'Sm1 (mm)',
                                 'loc': 0 },
                           'nsegment': 100}}}
    solute_parameters = {tracers[i]:{} for i in range(len(tracers))} #add other solutes here and c_old can be calibration or mean
    return Model(df, sas_specs=sas_specs, solute_parameters=solute_parameters, dt=24, influx='J (mm/hr)')

def make_gamma_model2_from(params): # for gamma distribution
    a = params
    tracers = [f'C{i}' for i in range(1,len(df)+1)]
    sas_specs = {'Q2 (mm/hr)':
                     {'tank2':
                          {'func': 'gamma',
                           'args': { 'a': 1,
                                 'scale': 'Sm2 (mm)',
                                 'loc': 0 },
                           'nsegment': 100}}}
    solute_parameters = {tracers[i]:{} for i in range(len(tracers))} #add other solutes here and c_old can be calibration or mean
    return Model(df2, sas_specs=sas_specs, solute_parameters=solute_parameters, dt=24, influx='Q1 (mm/hr)')

def make_gamma_model3_from(params): # for gamma distribution
    a = params
    tracers = [f'C{i}' for i in range(1,len(df)+1)]
    sas_specs = {'Q2 (mm/hr)':
                     {'tank12':
                          {'func': 'gamma',
                           'args': { 'a': 1,
                                 'scale': 'Sm12 (mm)',
                                 'loc': 0 },
                           'nsegment': 100}}}
    solute_parameters = {tracers[i]:{'observations':f'obs Q2 for C{i+1}'} for i in range(len(tracers))} #add other solutes here and c_old can be calibration or mean
    return Model(df3, sas_specs=sas_specs, solute_parameters=solute_parameters, dt=24, influx='J (mm/hr)')

# gamma distribution for optimizing S_0 has different assumptions of the affect of storage with the shape of the SAS function
# if using 'scipy.stats', then replace 'func' with 'scipy.stats' and function does not need ""

#%%
# -------------------Set parameters ----------------------------
#***params aren't use currently***
#--------storage----------
S1 = 1500 # (mm)
S2 = 3000 # (mm)
#--------distribution parameters--------
a = 1
b = 1 

#--------set params_init---------
params_init = S1 #***edit for distribution type***
params=params_init

#consider normalizing parameters to obtain better convergence of optimization
# Beta notes: a<1, b=1 young water prefernce, a=1,b<1 old water preference, a=b=1 uniform selection


#%% 
# ------------------Build the 1st model --------------------
from mesas.sas.model import Model
import warnings
warnings.simplefilter('ignore', pd.errors.PerformanceWarning) # ignore sparse dataframe performance warning from mesas

model = make_uniform_model1_from(params) #***edit which distribution***
# model = make_gamma_model1_from(params)
model.run()

#%% Calculate analytical SAS for 1st tank ------------------------------
# calculate bTTD using Neimi's theorem
# Q1
bpq1 = pd.DataFrame(index=range(0,len(model.data_df),1), columns=[f'obs time {i}' for i in range(1,len(tracers)+1)])
cq1 = model.data_df[[col for col in model.data_df.columns if '--> Q1' in col]].copy()
cj1 = df[[col for col in df.columns if 'C' in col]].copy() # input tracer concentrations
for i in range(0, len(bpq1)): # loop over obs time
    bpq1.iloc[:i+1,i] = 1/cj1.iloc[i,i]*cq1.iloc[i,:i+1][::-1]/dt # reverse order of columns C* --> Q
bpq1=bpq1.fillna(0)
# calcuate cumulative TTD
# Q1
Pq1 = (bpq1*dt).cumsum(axis=0) #cumsum over ages
Pq1 = pd.concat([pd.DataFrame([[0]*len(Pq1.columns)], columns=Pq1.columns), Pq1]).reset_index(drop=True) # pad with zeros at age 0
# # We can't observe any tracer at time 0 or age 0 coming out of storage
# calculate sT
sT1 = pd.DataFrame(index=range(0,len(model.data_df),1), columns=[f'obs time {i}' for i in range(len(tracers)+1)])
sT1.iloc[:,0]=0
for j in range(len(sT1.columns)-1):
    sT1prev = sT1[f'obs time {j}'].iloc[:-1].copy().to_list() # shift down by 1
    sT1.iloc[1:,j+1] = sT1prev - dt*model.data_df.iloc[j]['Q1 (mm/hr)']*bpq1.iloc[1:,j] #j on bpq since I didn't pad with zeros at time 0
    sT1.iloc[0,j+1] = model.data_df.iloc[j]['J (mm/hr)'] - dt*model.data_df.iloc[j]['Q1 (mm/hr)']*bpq1.iloc[0,j]
# sT columns are observed times, sT rows are tracer ages

# then average between time interval for ST
sT1mean = sT1.copy()
for j in range(len(sT1.columns)-1):
    sT1mean.iloc[:,j+1] = np.nanmean(np.vstack([pd.concat([pd.Series([0]), sT1.iloc[:-1,j].copy()]).to_numpy(), sT1.iloc[:,j+1].copy().to_numpy()]), axis=0) # average over time interval
sT1mean = sT1mean.iloc[:,1:] # remove time 0 column that was copied over from sT
sT1mean = sT1mean.fillna(0)
sT1 = sT1.fillna(0)
ST1 = (sT1mean*dt).cumsum(axis=0) # cumsum over ages
ST1 = pd.concat([pd.DataFrame([[0]*len(ST1.columns)], columns=ST1.columns), ST1]).reset_index(drop=True) # pad with zeros at age 0





#%%
# ------------------Build the 2nd model --------------------
# Set up 2nd tank with output tracers from 1st tank
df2 = df[['J (mm/hr)', 'Q1 (mm/hr)', 'Sm1 (mm)', 'Q2 (mm/hr)', 'Sm2 (mm)']].copy()
df22 = pd.DataFrame(model.data_df[[col for col in model.data_df.columns if '--> Q1' in col]].copy(), index = df.index)
df22.columns = tracers
df2 = pd.concat([df2, df22], axis=1)

model = make_uniform_model2_from(params) #***edit which distribution***
# model = make_gamma_model2_from(params)
model.run()


#%% Calculate analytical SAS for 2nd tank ------------------------------
# calculate bTTD using Neimi's theorem
# Q2
bpq2 = pd.DataFrame(index=range(0,len(model.data_df),1), columns=[f'obs time {i}' for i in range(1,len(tracers)+1)])
cq2 = model.data_df[[col for col in model.data_df.columns if '--> Q2' in col]].copy()
cj1 = df[[col for col in df2.columns if 'C' in col]].copy() # input tracer concentrations
for i in range(0, len(bpq2)): # loop over obs time
    bpq2.iloc[:i+1,i] = 1/cj1.iloc[i,i]*cq2.iloc[i,:i+1][::-1]/dt # reverse order of columns C* --> Q
bpq2=bpq2.fillna(0)
# bpq columns are observed times, bpq rows are tracer ages
# # calcuate cumulative TTD
# # Q2
Pq2 = (bpq2*dt).cumsum(axis=0) #cumsum over ages
Pq2 = pd.concat([pd.DataFrame([[0]*len(Pq2.columns)], columns=Pq2.columns), Pq2]).reset_index(drop=True) # pad with zeros at age 0
# # We can't observe any tracer at time 0 or age 0 coming out of storage
# calculate sT
sT2 = pd.DataFrame(index=range(0,len(model.data_df),1), columns=[f'obs time {i}' for i in range(len(tracers)+1)])
sT2.iloc[:,0]=0
for j in range(len(sT2.columns)-1):
    sT2prev = sT1[f'obs time {j}'].iloc[:-1].copy().to_list() # shift down by 1
    sT2.iloc[1:,j+1] = sT2prev - dt*model.data_df.iloc[j]['Q2 (mm/hr)']*bpq2.iloc[1:,j] #j on bpq since I didn't pad with zeros at time 0
    sT2.iloc[0,j+1] = model.data_df.iloc[j]['Q1 (mm/hr)'] - dt*model.data_df.iloc[j]['Q2 (mm/hr)']*bpq2.iloc[0,j]
# sT columns are observed times, sT rows are tracer ages

# then average between time interval for ST
sT2mean = sT2.copy()
for j in range(len(sT2.columns)-1):
    sT2mean.iloc[:,j+1] = np.nanmean(np.vstack([pd.concat([pd.Series([0]), sT2.iloc[:-1,j].copy()]).to_numpy(), sT2.iloc[:,j+1].copy().to_numpy()]), axis=0) # average over time interval
sT2mean = sT2mean.iloc[:,1:] # remove time 0 column that was copied over from sT
sT2mean = sT2mean.fillna(0)
sT2 = sT2.fillna(0)
ST2 = (sT2mean*dt).cumsum(axis=0) # cumsum over ages
ST2 = pd.concat([pd.DataFrame([[0]*len(ST2.columns)], columns=ST2.columns), ST2]).reset_index(drop=True) # pad with zeros at age 0





#%%
#------------------Build the 3rd model--------------------
# Set up overall tank with output tracers from 2nd tank
df3 = pd.DataFrame(model.data_df[[col for col in model.data_df.columns if '--> Q2' in col]].copy(), index=df.index)
df3.columns = [f'obs Q2 for C{i}' for i in range(1,len(tracers)+1)] #change column name so model.data_df can have --> Q2 in columns
df3 = pd.concat([df, df3], axis=1)

model = make_uniform_model3_from(params) #***edit which distribution***
# model = make_gamma_model3_from(params)
model.run()


#%%------------------Visualize model fit to data-------------------
from mesas.utils import vis
fig = plt.figure(figsize=[16,4])
for i in range(1,6):
    plt.plot(model.data_df.index, model.data_df[f'obs Q2 for C{i}'],'-', color='grey', lw=4, label=f'Observed C{i} outflow')
    plt.plot(model.data_df.index, model.data_df[f'C{i} --> Q2 (mm/hr)'],'-', color='cyan', zorder=10, label=f'Predicted C{i} outflow')
plt.axvspan(pd.Timestamp('2018-10-01'), pd.Timestamp('2019-09-30'), color='lightgrey', alpha=0.5)
plt.legend()
plt.title('Isotope outflow')
plt.show()

# %% 
# -------------------Calculate backwards TTD using Neimi's theorem-------------------

# Use tracer1 inputs and -->Q2 outputs to get total bpq
bpq = pd.DataFrame(index=range(0,len(model.data_df),1), columns=[f'obs time {i}' for i in range(1,len(tracers)+1)])
cq2 = model.data_df[[col for col in model.data_df.columns if '--> Q2' in col]].copy() # output tracer concentrations
cj1 = df[[col for col in df.columns if 'C' in col]].copy() # input tracer concentrations
for i in range(0, len(bpq)):
    bpq.iloc[:i+1,i] =  1/cj1.iloc[i,i]*cq2.iloc[i,:i+1][::-1]/dt # reverse order of columns C* --> Q2
bpq=bpq.fillna(0)
# bpq columns are observed times, bpq rows are tracer ages

# %%
#  -------------------Calculate cumulative TTD and S_T-------------------

# Weight by discharge to get total Pq
Pq = (bpq*dt).cumsum(axis=0) #cumsum over ages
Pq = pd.concat([pd.DataFrame([[0]*len(Pq.columns)], columns=Pq.columns), Pq]).reset_index(drop=True) # pad with zeros at age 0 # add Pq at time 0 = 0
# We can't observe any tracer at time 0 or age 0 coming out of storage

sT = pd.DataFrame(index=range(0,len(model.data_df),1), columns=[f'obs time {i}' for i in range(len(tracers)+1)])
sT.iloc[:,0]=0
for j in range(0, len(sT.columns)-1):
    sTprev = sT[f'obs time {j}'].iloc[:-1].copy().to_list() # shift down by 1
    sT.iloc[1:,j+1] = sTprev - dt*model.data_df.iloc[j]['Q2 (mm/hr)']*bpq.iloc[1:,j] #j on bpq since I didn't pad with zeros at time 0
    sT.iloc[0,j+1] = model.data_df.iloc[j]['J (mm/hr)'] - dt*model.data_df.iloc[j]['Q2 (mm/hr)']*bpq.iloc[0,j]
# sT columns are observed times, sT rows are tracer ages

# then average between time interval for ST
sTmean = sT.copy()
for j in range(len(sT.columns)-1):
    sTmean.iloc[:,j+1] = np.nanmean(np.vstack([pd.concat([pd.Series([0]), sT.iloc[:-1,j].copy()]).to_numpy(), sT.iloc[:,j+1].copy().to_numpy()]), axis=0) # average over time interval
sTmean = sTmean.iloc[:,1:] # remove time 0 column that was copied over from sT
sTmean = sTmean.fillna(0)
sT = sT.fillna(0)
ST = (sTmean*dt).cumsum(axis=0) # cumsum over ages
ST = pd.concat([pd.DataFrame([[0]*len(ST.columns)], columns=ST.columns), ST]).reset_index(drop=True) # pad with zeros at age 0


# %% -----------------Animate cumTTD and ST-------------------
# from matplotlib.animation import FuncAnimation
# fig, ax = plt.subplots(figsize=(10,4))
# line, = ax.plot([],[],'r-')

# def init():
#     line.set_data([],[])
#     return line,

# def animate(i):
#     x=ST[f'obs time {i+1}']
#     y=Pq[f'obs time {i+1}']
#     line.set_data(x,y)
#     ax.set_title(f'Observed time {i}')
#     return line,

# ani = FuncAnimation(fig, animate, init_func=init, frames=len(model.data_df), interval=100, blit=True)
# plt.show()
# ani.save('animation.gif',writer='pillow')

#%%
# -------------------Plot cumulative TTD and S_T-------------------
plt.figure(figsize=(10,4))
# for i in range(1, len(model.data_df)):
i=len(model.data_df) #not -1 since obs times start at 1
if True:
    plt.plot(ST[f'obs time {i}'], Pq[f'obs time {i}'], label='Analystical SAS', marker='.', alpha=.5) #if i cumsum over ages
plt.plot([0,df['Sm12 (mm)'][-1]], [0,1], label='true uniform SAS', zorder=1) #uniform SAS line
# plt.plot([0,df['Sm1 (mm)'][-1]], [0,1], label='true tank1 SAS', zorder=1) #uniform SAS line
# plt.plot([0,df['Sm2 (mm)'][-1]], [0,1], label='true tank2 SAS', zorder=1) #uniform SAS line
plt.xlabel('$S_T$ (mm)')
plt.ylabel('Cumulative TTD $P_Q$')
plt.title('Analytical SAS')
# plt.ylim([0,1])
# plt.xlim([0, 400])
plt.legend()
plt.tight_layout()
plt.show()


# Should max(df['Sm12 (mm)']) be equal to where analytical SAS reaches 1? or is it a difference of dynamic storage vs total storage?

# from mesas.utils import vis
# plt.figure(figsize=(10,4))
# vis.plot_SAS_cumulative(model, 'Q2 (mm/hr)')
# vis.plot_SAS_cumulative(model, 'Q1 (mm/hr)')
# plt.title('SAS model')


# %%
