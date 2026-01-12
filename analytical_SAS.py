# This script will find SAS functions based on an analytical
# approach by feeding a discretized concentration once in a sequential
# timestep over a long period of time.
# Neimi's theorem will be used to find the backward TTD to sum over
# time to get cumulative TTD.
# S_T will be calculated by summing J(t)-Q(t)*P_Q(t) over time.

# This script needs a given SAS function form, J, and Q

# Date: 10/07/2025

#%% ------------------Create J and Q-------------------
import pandas as pd
import numpy as np  
import matplotlib.pyplot as plt
from mesas.sas.model import Model

date = pd.date_range(start='2018-01-01 00:00:00', end='2020-12-31 23:00:00', freq='D')
df = pd.DataFrame(index=date)
# print(len(date)) # 52608 hours

# use freq = H for this function
# Define J by random Poisson process and gamma distribution of rates
# Each rain event has a random start time, duration, and intensity
# def simulate_rainfall_rate(
#     n_hours=24*30,         # number of hours to simulate
#     lambda_per_hour=0.007,   # mean number of rain events per hour (~5 per month)
#     mean_intensity=5.0,    # mean intensity of each event (mm/hr)
#     duration_mean=3.0,     # mean duration of each event (hours)
#     start_time="2020-01-01 00:00",
#     random_seed=None
# ):
#     """
#     Simulate rainfall rate (mm/hr) as a Poisson process of rain events.
#     Each event has a random start time, duration, and intensity.
#     Returns a pandas DataFrame with hourly rainfall rate.
#     """
#     if random_seed is not None:
#         np.random.seed(random_seed)

#     # Step 1: number of events expected
#     total_events = np.random.poisson(lam=lambda_per_hour * n_hours)
    
#     # Step 2: event start times (uniform over simulation period)
#     start_times = np.random.uniform(0, n_hours, total_events)
    
#     # Step 3: random duration and intensity per event
#     durations = np.random.exponential(scale=duration_mean, size=total_events)
#     intensities = np.random.exponential(scale=mean_intensity, size=total_events)
    
#     # Step 4: create time grid and initialize rate array
#     t = np.arange(n_hours)
#     rainfall_rate = np.zeros_like(t, dtype=float)
    
#     # Step 5: for each event, add its intensity over its duration
#     for s, d, i in zip(start_times, durations, intensities):
#         start_idx = int(np.floor(s))
#         end_idx = int(np.floor(s + d))
#         end_idx = min(end_idx, n_hours - 1)
#         rainfall_rate[start_idx:end_idx+1] += i
    
#     return rainfall_rate

# df['J (mm/hr)'] = simulate_rainfall_rate(n_hours=len(date), lambda_per_hour=0.007, mean_intensity=3, random_seed=42)


# Try another daily method
def simulate_daily_rainfall_rate(
    n_days=30*365,         # number of days to simulate
    prob_rain_per_day=0.3,   # probability of rain each day
    shape=2.0,            # shape parameter for gamma distribution (controls skew)
    scale=1.5,            # scale parameter for gamma distribution  (controls mean intensity)
    hours_per_day=24,
    random_seed=None
):
    """
    Simulate daily rainfall rate (mm/hr) using a Bernoulli process for rain occurrence
    and a gamma distribution for intensity on rainy days.
    Returns a pandas DataFrame with hourly rainfall rate.
    """
    if random_seed is not None:
        np.random.seed(random_seed)

    # Rain days (1 if rain, 0 if no rain)
    rain_days = np.random.rand(n_days) < prob_rain_per_day

    # Rain intensity
    rain_mm_day = np.zeros(n_days)
    rain_mm_day[rain_days] = np.random.gamma(shape, scale, size=rain_days.sum())

    # Convert to average hourly rate (mm/hr)
    rain_mm_hr = rain_mm_day / hours_per_day
    
    return rain_mm_hr

df['J (mm/hr)'] = simulate_daily_rainfall_rate(n_days=len(date), prob_rain_per_day=.3, shape=2.0, scale=1.5, random_seed=42)

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
    K=20.0,         # reservoir time constant (hours)
    C=0.8,          # runoff coefficient (fraction of rainfall that becomes runoff)
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
        S[t+1] = S[t] + (C * rainfall_rate[t] - Q[t]) * dt
        S[t+1] = max(S[t+1], 0)  # prevent negative storage
        Smean[t] = (S[t+1] + S[t]) / 2

    return Q, Smean

# For daily rainfall use dt=24
dt=24
k=2400 # hours (100 days)
Q,Sm = simulate_discharge(df['J (mm/hr)'].values, dt=dt, K=k, C=0.6, S0=df['J (mm/hr)'].mean()*k*0.6)
df['Q (mm/hr)'] = Q
df['Storage (dt-mean) (mm)'] = Sm

#Check with constants
# df['J (mm/hr)'] = df['J (mm/hr)'].mean()
# df['Q (mm/hr)'] = df['J (mm/hr)']
# df['Storage (dt-mean) (mm)'] = df['J (mm/hr)'].mean()*2400*.6

# Plot rainfall vs discharge
plt.figure(figsize=(10,4))
plt.plot(df.index, df["Q (mm/hr)"], color="darkgreen", label="Discharge (mm/hr)")
plt.ylabel("Discharge (mm/hr)", color="green")
plt.title("Simulated discharge response to Poisson rainfall")
plt.tight_layout()
plt.show()

# Resample to daily
# df = df.resample('D').sum()#.mean()
# dt = 24 # hour in a day


#%% 
# -------------------Create staggered tracer input-------------------
from scipy import sparse
tracers = [f'C{i}' for i in range(1,len(df)+1)] # list of [C1, C2, C3, ...]
# C = sparse.eye(len(df), dtype=np.float32, format='csr')
df = pd.concat([df, pd.DataFrame.sparse.from_spmatrix(sparse.eye(len(df), dtype=np.float32, format='csr'), columns=tracers, index=df.index)], axis=1)

#%%
# ------------------Fit a distribution-------------------
from scipy.stats import beta, gamma

def make_uniform_model_from(params): # for uniform distribution                          
    S_max = params
    tracers = [f'C{i}' for i in range(1,len(df)+1)]
    sas_specs = {'Q (mm/hr)':
                        {'ORPB':
                          {"ST": [0, 'Storage (dt-mean) (mm)'],
                           "P": [0.0, 1.0]}
                           }
                           }
    solute_parameters = {tracers[i]:{} for i in range(len(tracers))}
    return Model(df, sas_specs=sas_specs, solute_parameters=solute_parameters, dt=24, influx='J (mm/hr)', n_substeps=1, record_state=True)


def make_beta_model_from(params): # for beta distribution
    a, b = params
    tracers = [f'C{i}' for i in range(1,len(df)+1)]
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
def make_gamma_model_from(params): # for gamma distribution
    a = params
    tracers = [f'C{i}' for i in range(1,len(df)+1)]
    sas_specs = {'Q (mm/hr)':
                     {'ORPB':
                          {"scipy.stats": "gamma",#'func': 'gamma',
                           'args': { 'a': a,
                                 'scale': 'Storage (dt-mean) (mm)',
                                 'loc': 0 },
                           'nsegment': 100}}}
    solute_parameters = {tracers[i]:{} for i in range(len(tracers))} #add other solutes here and c_old can be calibration or mean
    return Model(df, sas_specs=sas_specs, solute_parameters=solute_parameters, dt=24, influx='J (mm/hr)', n_substeps=3)
# gamma distribution for optimizing S_0 has different assumptions of the affect of storage with the shape of the SAS function
# if using 'scipy.stats', then replace 'func' with 'scipy.stats' and function does not need ""



#%%
# -------------------Set parameters ----------------------------

#--------storage----------
S_max = 1500 # (mm)

#--------distribution parameters--------
a = df['Q (mm/hr)'].mean()**2/df['Q (mm/hr)'].var()
b = 1 

#--------set params_init---------
params_init = a #***edit for distribution type***
params=params_init

#consider normalizing parameters to obtain better convergence of optimization
# Beta notes: a<1, b=1 young water prefernce, a=1,b<1 old water preference, a=b=1 uniform selection


#%% 
# ------------------Build the model --------------------
from mesas.sas.model import Model
model = make_uniform_model_from(params) #***edit which distribution***
# model = make_gamma_model_from(params)
model.run()


#%%------------------Visualize the model fit to data-------------------
from mesas.utils import vis
fig = plt.figure(figsize=[16,4])
for i in range(1,6):
    plt.plot(model.data_df.index, model.data_df[f'C{i}'],'-', color='grey', label=f'Observed C{i} outflow')
    plt.plot(model.data_df.index, model.data_df[f'C{i} --> Q (mm/hr)'],'-', color='cyan', label=f'Predicted C{i} outflow')
plt.axvspan(pd.Timestamp('2018-10-01'), pd.Timestamp('2019-09-30'), color='lightgrey', alpha=0.5)
plt.legend()
plt.title('Isotope outflow')
plt.xlim([pd.Timestamp('2018-01-01'), pd.Timestamp('2018-01-31')])
plt.show()


# %% 
# -------------------Calculate backwards TTD using Neimi's theorem-------------------
bpq = pd.DataFrame(index=range(0,len(model.data_df),1), columns=[f'obs time {i}' for i in range(1,len(tracers)+1)])
cq =  model.data_df[[col for col in model.data_df.columns if '--> Q' in col]].copy()
for i in range(0, len(bpq.columns)): # loop over obs time
    # bpq columns = 1/c_in * c_out(row backwards from diagonal) / dt
    # bpq.iloc[:i+1,i] = (1/model.data_df.iloc[i][f'C{i+1}'])*model.data_df.iloc[i,len(bpq)+3:i+len(bpq)+3+1][::-1]/dt # reverse order of columns C* --> Q
    bpq.iloc[:i+1,i] = (1/model.data_df.iloc[i][f'C{i+1}'])*cq.iloc[i,:i+1][::-1]/dt # reverse order of columns C* --> Q
bpq=bpq.fillna(0)
# bpq[np.abs(bpq)<1e-6] = 0  # set very small values to zero
# bpq columns are observed times, bpq rows are tracer ages


# %%
#  -------------------Calculate cumulative TTD and S_T-------------------
Pq = (bpq*dt).cumsum(axis=0) #cumsum over ages
Pq = pd.concat([pd.DataFrame([[0]*len(Pq.columns)], columns=Pq.columns), Pq]).reset_index(drop=True) # pad with zeros at age 0 # add Pq at time 0 = 0
# Pq = pd.concat([pd.DataFrame({'obs time 0': [0]*len(Pq)}), Pq], axis=1) # pad with zeros at time=0
# We can't observe any tracer at time 0 or age 0 coming out of storage

sT = pd.DataFrame(index=range(len(model.data_df)), columns=[f'obs time {i}' for i in range(len(tracers)+1)])
# sT = pd.concat([pd.DataFrame({'obs time 0': [0]*len(sT)}), sT], axis=1) # pad with zeros at time=0
sT.iloc[:,0] = 0
for j in range(len(sT.columns)-1):
    sTprev = sT[f'obs time {j}'].iloc[:-1].copy().to_list() # shift down by 1
    sT.iloc[1:,j+1] = sTprev - dt*model.data_df.iloc[j]['Q (mm/hr)']*bpq.iloc[1:,j] #j on bpq since I didn't pad with zeros at time 0
    sT.iloc[0,j+1] = model.data_df.iloc[j]['J (mm/hr)'] - dt*model.data_df.iloc[j]['Q (mm/hr)']*bpq.iloc[0,j]

# sT columns are observed times, sT rows are tracer ages
# Error bugging: try comparing model.get_sT() and model.get_Pq(flux='Q (mm/hr)') to sT and Pq here
# then average between time interval for ST
sTmean = sT.copy()
for j in range(len(sT.columns)-1):
    sTmean.iloc[:,j+1] = np.nanmean(np.vstack([pd.concat([pd.Series([0]), sT.iloc[:-1,j].copy()]).to_numpy(), sT.iloc[:,j+1].copy().to_numpy()]), axis=0) # average over time interval
sTmean = sTmean.iloc[:,1:] # remove time 0 column that was copied over from sT
sTmean = sTmean.fillna(0)
sT = sT.fillna(0)
ST = (sTmean*dt).cumsum(axis=0) # cumsum over ages #sTmean doesn't match ST_mesas
ST = pd.concat([pd.DataFrame([[0]*len(ST.columns)], columns=ST.columns), ST]).reset_index(drop=True) # pad with zeros at age 0 # add ST at time 0 = 0


# %%------------------Animate cumTTD and S_T-------------------
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
i=len(model.data_df)-1
if True:
    # plt.plot(Pq.iloc[:,i], label='Analystical SAS')
    # plt.plot(ST[f'obs time {i}']/df.iloc[i]['Storage (dt-mean) (mm)'], Pq[f'obs time {i}'], label='Analystical SAS', marker='.') #if i cumsum over ages
    plt.plot(ST[f'obs time {i+1}'], Pq[f'obs time {i+1}'], label='Analystical SAS', marker='.') #if i cumsum over ages
plt.xlabel('$S_T$ (mm)')
plt.ylabel('Cumulative TTD $P_Q$')
plt.title('Analytical SAS')
# plt.ylim([0,1])
# plt.xlim([0, 1500])
# plt.legend()
plt.plot([0,df['Storage (dt-mean) (mm)'][-1]], [0,1], label='true SAS', zorder=1) #uniform SAS line
pq_mesas=model.get_pQ(flux='Q (mm/hr)')
Pq_mesas=(pq_mesas*dt).cumsum(axis=0)
ST_mesas=model.get_ST()
plt.plot(ST_mesas[:,i], Pq_mesas[:,i], label='MESAS SAS', marker='.',alpha=.5, color='purple', zorder=.1) #MESAS SAS
plt.xlim([0, df['Storage (dt-mean) (mm)'].max()*1.1])
plt.legend()
plt.tight_layout()
plt.show()

# from mesas.utils import vis
# plt.figure(figsize=(10,4))
# # for i in range(1, len(model.data_df)):
# vis.plot_SAS_cumulative(model, 'Q (mm/hr)',i=len(model.data_df)-1)
# plt.title('SAS model')


# %%
plt.figure(figsize=(10,4))
plt.plot(Pq, label='Pq')
plt.title('Pq')


plt.figure(figsize=(10,4))
plt.plot(ST, label='ST')
plt.title('ST')
# %%
