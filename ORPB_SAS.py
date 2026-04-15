# This script will upload ORPB_isotope_data.csv and 
# visualize the tracer timeseries along with SAS functions

# Date: 06/17/2025

#%% 
# ------------------Import dataset-------------------
import pandas as pd
import numpy as np  
import matplotlib.pyplot as plt
from mesas.sas.model import Model
from permetrics.regression import RegressionMetric

data_df = pd.read_csv('ORPB_isotope_data.csv', index_col=0, parse_dates=[0])

#%% plot timeseries of full data (Figure 1 in my dissertation proposal)
fig,(ax1,ax2, ax3)=plt.subplots(nrows=3,ncols=1,figsize=[10,9])
ax1.plot(data_df.index, data_df['discharge (mm/hr)'], color='blue', label='discharge (mm/hr)')
ax1.set_xlabel('Date')
ax1.set_ylabel('Discharge (mm/hr)', color='blue')
ax1.tick_params(axis='y', labelcolor='blue')
ax12 = ax1.twinx()
ax12.bar(data_df.index, data_df['rainfall (mm/hr)'], color='gray', alpha=0.6, width=0.01, label='rainfall (mm/hr)')
ax12.bar(data_df.index, data_df['snowfall SWE (mm/hr)'], color='black', alpha=0.6, width=0.01, label='snowfall SWE (mm)')
ax12.set_ylabel('(mm/hr)', color='gray')
ax12.tick_params(axis='y', labelcolor='gray')
ax12.set_ylim([0,0.75])#max(max(data_df['rainfall (mm/hr)']), max(data_df['snowfall SWE (mm/hr)']))*1.1])
ax12.invert_yaxis()
ax1.set_title('ORPB Discharge and Rainfall')
ax1.legend()
ax12.legend()

issample = np.logical_not(np.isnan(data_df['ORPB 18O']))
ax2.plot(data_df.index[issample], 
         data_df[['ORPB 2H',
             'ORPB 18O', 
             'ORPB 17O',]][issample],'.-', label=['2H', '18O', '17O'])
ax2.legend()
ax2.set_title('ORPB Tracer at stream')
ax3.plot(data_df.index,
         data_df[['precip 2H', 
             'precip 18O', 
             'precip 17O']], label=['precip 2H','precip 18O','precip 17O'])
ax3.legend()
ax3.set_title('ORPB Tracer in precipitation')
plt.tight_layout()
plt.show()

#%% Plot hourly, weekly, monthly timeseries of ORPB 18O
ORPB18O_hourly = data_df['ORPB 18O'].resample('h').mean()
ORPB18O_weekly = data_df['ORPB 18O'].resample('W').mean()
ORPB18O_monthly = data_df['ORPB 18O'].resample('M').mean()
fig,(ax1,ax2, ax3)=plt.subplots(nrows=3,ncols=1,figsize=[10,9], sharey=True)
ax1.plot(ORPB18O_hourly.index, ORPB18O_hourly,'.')
ax1.set_title('ORPB $^{18}O$ hourly')
ax2.plot(ORPB18O_weekly.index, ORPB18O_weekly,'.')
ax2.set_title('ORPB $^{18}O$ weekly')
ax3.plot(ORPB18O_monthly.index, ORPB18O_monthly,'.')
ax3.set_title('ORPB $^{18}O$ monthly')
plt.tight_layout()
plt.show()


#%%
data_df = data_df.loc[pd.Timestamp('2014-01-01'): pd.Timestamp('2018-09-21')] # 80% training data (data_df.iloc[int(len(data_df)*0.8)])
# data_df = data_df.loc[pd.Timestamp('2014-08-01'): pd.Timestamp('2014-08-31')] #subset to Putnam's data range 2014-08-01 - 2016-08-31
issample = np.logical_not(np.isnan(data_df['ORPB 18O']))
#--------influx----------
# df['influx (mm/hr)'] = df['rainfall (mm/hr)']
# df['influx (mm/hr)'] = df[['rainfall (mm/hr)','snowfall SWE (mm/hr)','snowmelt (mm/hr)']].sum(axis=1)
data_df['influx (mm/hr)'] = data_df[['rainfall (mm/hr)','snowmelt (mm/hr)']].sum(axis=1)

data_df['quickflow (mm/hr)'] = data_df['discharge (mm/hr)'] - data_df['baseflow 1 (mm/hr)']
data_df['bf1_weight'] = data_df['baseflow 1 (mm/hr)'] / data_df['discharge (mm/hr)']
data_df['qf_weight'] = data_df['quickflow (mm/hr)'] / data_df['discharge (mm/hr)']

# Find data where quickflow is small - from baseflow separation code in GenerateCleanData_v2.ipyb
# data_df["rain+melt (mm/hr)"] = data_df["rainfall (mm/hr)"] + data_df["snowmelt (mm/hr)"]
# data_df["inputs in last 2d?"] = data_df["rain+melt (mm/hr)"].rolling('2d').sum() > 0
# data_df["inputs in next 3 hr?"] = data_df["rain+melt (mm/hr)"].rolling('3h').sum().shift(-3) > 0
# data_df['discharge censored (mm/hr)'] = data_df['discharge (mm/hr)']
# data_df.loc[data_df["inputs in last 2d?"] | data_df["inputs in next 3 hr?"], 'discharge censored (mm/hr)'] = np.NaN
# isbaseflow = data_df.loc[(data_df['discharge censored (mm/hr)'].notna()) & issample].index
# # or try finding based on quickflow < threshold
# plt.hist(data_df['quickflow (mm/hr)'], bins=2500)
# plt.xlim([0,0.04])
# plt.ylabel('frequency')
# plt.xlabel('quickflow (mm/hr)')
# print(len(data_df.loc[data_df['quickflow (mm/hr)']<0.001, 'quickflow (mm/hr)']), 'meet criteria out of ', len(data_df['quickflow (mm/hr)']))
# isbaseflow = data_df.loc[(data_df['quickflow (mm/hr)']<0.001) & issample].index
# isquickflow = data_df.loc[(data_df['quickflow (mm/hr)']>=0.001) & issample].index
# print(data_df.columns)

#%% # ------------------Check for nans-------------------
# ------------------Check for nans-------------------
#Just solving for 18O is sufficient since they all vary similarly
#also fit C_old adn S_0 by minimizing RMSE

# Check for nans
print('Number of nans: ',len(data_df.loc[data_df['precip 18O'].isna()==True])) #1533 nan rows
# diffs=data.loc[data['precip 18O'].isna()==True].index.to_series().diff()
# diffs.value_counts() # check for frequency of gaps
# decide on what to fill nans with
#fill nas with either mean of data or nearest neighbor estimate

print(data_df['precip 18O'].describe())
plt.hist(data_df['precip 18O'])
plt.title('precip 18O histogram')

#%% # ----------------Fill nans -----------------------
# fill nans with mean
mean = data_df['precip 18O'].mean()
df= data_df.copy() #make a copy of the data_df
# df.loc[df['precip 18O'].isna()==True, 'precip 18O']=mean
# Instead of filling with mean, use forward/backward fill or interpolation
df['precip 18O'] = df['precip 18O'].ffill().bfill()

# assert positive ET values
# df.loc[df['ET (mm/hr)']<0, 'ET (mm/hr)']=0

# Putnam's data subset
# df = df.loc['2014-09-01':'2016-08-31']
# issample = np.logical_not(np.isnan(df['ORPB 18O']))
#%% # ------------------Check for filled nans and frequency-------------------
# Check for filled nas and frequency

print('Number of rows turned from nan: ', len(df.loc[df['precip 18O']==mean])) #should be 1533 rows unless some obs happen to be the mean value
print('Number of nans: ', len(df.loc[df['precip 18O'].isna()])) #should be 0 rows
# double check that the data is hourly or has constant frequency for SAS
print('Inferred timeseries frequency: ', pd.infer_freq(df.index)) # should be h for hourly



#%% 
# ------------------Fit a distribution-------------------
from scipy.stats import beta, gamma

def make_uniform_model_from(params): # for uniform distribution                          
    S_0, c18O_old = params
    df['abs_storage (mm)'] = df['storage (mm)'] + S_0
    sas_specs = {'discharge (mm/hr)':
                        {'ORPB':
                          {"ST": [0, 'abs_storage (mm)'],
                           "P": [0.0, 1.0]}
                           }
                           }
    solute_parameters = {'precip 18O': {'C_old': c18O_old}}
    return Model(df, sas_specs=sas_specs, solute_parameters=solute_parameters, dt=1, influx='influx (mm/hr)')


def make_beta_model_from(params): # for beta distribution
    S_0, c18O_old, a, b = params
    df['abs_storage (mm)'] = df['storage (mm)'] + S_0
    sas_specs = {'discharge (mm/hr)':
                        {'ORPB':
                          {'func': "beta",
                           'args': { 'a': a,
                                     'b': b,
                                 'scale': 'abs_storage (mm)',
                                 'loc': 0 },
                           'nsegment': 100}}}
    solute_parameters = {'precip 18O': {'C_old': c18O_old}}
    return Model(df, sas_specs=sas_specs, solute_parameters=solute_parameters, dt=1, influx='influx (mm/hr)')
# Note: dt must match the data timestep (1 hour here)
def make_gamma_model_from(params): # for gamma distribution
    c18O_old, a, lamda, S_c, et_scale = params
    normalized_storage = df['storage (mm)']-df['storage (mm)'].mean()
    df['S_scale'] = lamda*(normalized_storage-S_c) #slope and intercept 
    sas_specs = {'discharge (mm/hr)':
                     {'ORPB':
                          {'func': 'gamma',
                           'args': { 'a': a,
                                 'scale': 'S_scale', #'abs_storage (mm)', #mean travel time
                                 'loc': 0 }}},
                 'ET (mm/hr)':
                     {'ET':
                          {'func': 'kumaraswamy',
                           'args':{
                                'a': 1.0,
                                'b': 1.0,
                                'loc': 0.0,
                                'scale': et_scale}}}
                }
    solute_parameters = {'precip 18O': {'C_old': c18O_old, 'observations': 'ORPB 18O'}}#{'discharge (mm/hr)': 'ORPB 18O'}}} #add other solutes here and c_old can be calibration or mean
    return Model(df, sas_specs=sas_specs, solute_parameters=solute_parameters, dt=1, influx='influx (mm/hr)', record_state=True, verbose=True, n_substeps=1)

def make_gamma_split_model_from(params): # split quickflow and baseflow with different SAS functions
    c18O_old, a_qf, t_qf, a_bf, lamda, S_c, et_scale = params
    df['S_scale'] = lamda*(df['storage (mm)']-S_c) #slope and intercept 
    sas_specs = {'discharge (mm/hr)':
                     {'qf_weight':
                          {'func': 'gamma',
                           'args': { 'a': a_qf,
                                 'scale': t_qf,
                                 'loc': 0 }},
                      'bf1_weight':
                          {'func': 'gamma',
                           'args': { 'a': a_bf,
                                 'scale': 'S_scale',
                                 'loc': 0 }}},
                 'ET (mm/hr)':
                     {'ET':
                          {'func': 'kumaraswamy',
                           'args':{
                                'a': 1.0,
                                'b': 1.0,
                                'loc': 0.0,
                                'scale': et_scale}}}
                }
    solute_parameters = {'precip 18O': {'C_old': c18O_old, 'observations': 'ORPB 18O'}}
    return Model(df, sas_specs=sas_specs, solute_parameters=solute_parameters, dt=1, influx='influx (mm/hr)', record_state=False, verbose=True, n_substeps=1)
# gamma distribution for optimizing S_0 has different assumptions of the affect of storage with the shape of the SAS function
# if using 'scipy.stats', then replace 'func' with 'scipy.stats' and function does not need ""
def make_kumar_model_from(params): # for kumaraswamy distribution
    c18O_old, a, b, et_storage = params
    sas_specs = {'discharge (mm/hr)':
                     {'ORPB':
                          {'func': 'kumaraswamy',
                           'args': { 'a': a,
                                     'b': b,
                                 'scale': 'abs_storage (mm)',
                                 'loc': 0 }
                        #    'nsegment': 100
                           }},
                 'ET (mm/hr)':
                     {'ET':
                          {'ST': [0, et_storage],
                           'P': [0.0, 1.0]}}}
    solute_parameters = {'precip 18O': {'C_old': c18O_old}} #add other solutes here and c_old can be calibration or mean
    return Model(df, sas_specs=sas_specs, solute_parameters=solute_parameters, dt=1, influx='influx (mm/hr)')

def make_Putnam_model_from(params): # from Putnam Chapter 3
    a_bf, t_bf, LnS_Tet = params
    #Normalize parameters
    S_Tet = np.exp(LnS_Tet)*1530.484 #ensures S_Tet is always positive
    # c18O_old = c18O_old*(-7.6) #normalize to -7.6
    t_bf = t_bf*1423.08 # normalize to corrected with >0 ET, storage calc of storage.max=storage.min * bf1_weight.mean
    # make new column that is a parameter * wetness and put into ST max for quickflow
    # df['wwetness'] = df['wetness']*q_max #normalize to wetness
    sas_specs = {'quickflow (mm/hr)':
                        {'ORPB qf':
                            {"ST": [0, 0.254],
                             "P": [0.0, 1.0]}
                        }, #if >1 dict for 'quickflow (mm/hr)', then mesas looks for column named 'ORPB qf' and other dict key(s) -- >1 dict will allow for weighted SAS functions, these columns in df will be the weights [0,1]
                 'baseflow 1 (mm/hr)':
                        {'ORPB bf':
                                {'func': 'gamma',
                                 'args': { 'a': a_bf,
                                       'scale': t_bf,
                                         'loc': 0 }}
                        },
                 'ET (mm/hr)': 
                        {'ET':
                            {'ST': [0, S_Tet],
                             'P': [0.0, 1.0]}
                        }
                }
    solute_parameters = {'precip 18O': {'C_old': -7.6}}
    return Model(df, sas_specs=sas_specs, solute_parameters=solute_parameters, dt=1, influx='influx (mm/hr)')
                  


# calc error when quickflow is small  as determined by baseflow separation
# step-wise calibration of parameters to capture certain parts of variability
# Function that builds SAS model and returns RMSE
def minimize_me(params):
   model = make_gamma_model_from(params) #***edit which distribution to minimize***
   model.run()
#    obs = model.data_df['ORPB 18O'][isbaseflow].to_numpy() # for baseflow or quickflow
#    pred = (1-model.data_df['bf1_weight'][isquickflow].to_numpy())*(model.data_df['precip 18O --> quickflow (mm/hr)'][isquickflow].to_numpy()) #for quickflow
#    pred = model.data_df['bf1_weight'][isbaseflow].to_numpy()*(model.data_df['precip 18O --> baseflow 1 (mm/hr)'][isbaseflow].to_numpy()) # for baseflow

   obs = model.data_df['ORPB 18O'][issample].to_numpy() # for all data valid with ORPB 18O sample
#    pred = (1-model.data_df['bf1_weight'][issample].to_numpy())*(model.data_df['precip 18O --> quickflow (mm/hr)'][issample].to_numpy()) + model.data_df['bf1_weight'][issample].to_numpy()*(model.data_df['precip 18O --> baseflow 1 (mm/hr)'][issample].to_numpy())
   pred = model.data_df['precip 18O --> discharge (mm/hr)'][issample] # for data not separated into qf and bf
   RMSE = np.sqrt(np.mean((pred-obs)**2))
   print(f'RMSE = {RMSE} for params = {params}')
   return RMSE

# Function that finds best KGE (-inf, 1] where 1 is perfect fit)
def maximize_me(params):
    model = make_Putnam_model_from(params) #***edit which distribution to maximize***
    model.run()
    obs = model.data_df['ORPB 18O'][issample].to_numpy()
    pred = (1-model.data_df['bf1_weight'][issample].to_numpy())*(model.data_df['precip 18O --> quickflow (mm/hr)'][issample].to_numpy()) + model.data_df['bf1_weight'][issample].to_numpy()*(model.data_df['precip 18O --> baseflow 1 (mm/hr)'][issample].to_numpy())
    evaluator = RegressionMetric(obs, pred)
    kge = evaluator.kling_gupta_efficiency()
    print(f'KGE = {kge} for params = {params}')
    return abs(kge-1) #return negative since fmin only minimizes to zero


#%%
# -------------------Set parameters ----------------------------

#--------storage----------
# Now let's supply initial estimates of parameters #opt gamma rmse: 0.5374013865404, b-rmse: 0.53211672085
# S_0 = 5701.46684 # latest g-opt: 5701.46684, b-opt: 6245.18546 (mm) initial storage, can be set to any value
# Try keeping S_0 constant and optimizing other parameters

# df['abs_storage (mm)'] = df['storage (mm)'] + S_0 # (mm)
# S_min = df['storage (mm)'].min() # (mm)
# S_max = df['storage (mm)'].max() # (mm)
# df['wetness'] = (df['storage (mm)'] - S_min) / (S_max - S_min) # catchment wetness

#--------ET storage----------
# S_Tet = 0.02822 #0.58 # normalize to S_max-S_min = 1530.484 (mm)
# LnS_Tet = np.log(S_Tet) # optimize ln(S_Tet) to ensure S_Tet is always positive
et_scale = 40.827 #43.4542481 # (mm) from pmcmc

#--------solute parameters----------
c18O_old = -7.6 #-7.3898 #-7.790043373 #normalize to -7.6 for Putnam model #df['precip 18O'].mean() -1 # latest g-opt: -7.65909, b-opt: -7.720487 (per mil)
# 0.98*-7.6 = -7.44 for starting date 2014-08-01 (Putnam's subset)
#--------distribution parameters--------
a = 1.26 #1.26399 #1.249941962 # laetest g-opt: 0.8293324, b-opt: 0.733789
# b = 0.9923319 #latest b-opt: 0.945078
# a_bf = 1.26 #(df['baseflow 1 (mm/hr)'].mean())**2/(df['baseflow 1 (mm/hr)'].std())**2 # mean^2/std^2 = (df['baseflow 1 (mm/hr)'].mean())**2/(df['baseflow 1 (mm/hr)'].std())**2
# t_bf = 1.48 #1.36-2.23 normalized to storage (df['baseflow 1 (mm/hr)'].std())**2/df['baseflow 1 (mm/hr)'].mean() #variance/mean, or should be mean storage that contains baseflow
# q_max = 1 #normalized to wetness
# lamda = -109.184 #-115.0382 #-106.9838848
# S_c = 51.5895 #49.9052 #50.67889314

# calculate realistic scale and lamda
# S_scale = controls mean travel time and how it responds to storage, lamda = senstivity of transit time to storage, S_c = storage threshold where travel time begins to change ensure positivity by making storage>S_c so S_c<=min(storage)
S_c=df['storage (mm)'].min()-50 #min(df['storage (mm)'])-50 # this is solid, should be less than dynamic storage
lamda=.5 #12.1017 # 10/(df['storage (mm)'].median()-S_c) #median of storage is s_ref
#%%
# sscale=lamda*(scale-sc)
# plt.plot(scale,sscale) #increase monotonically
# plt.plot(data_df['storage (mm)'], sscale) # check that they increase
# hard to tell if sscale is good for discharge, try binning by sscale
# data_df['sscale']=sscale
# df = data_df[['sscale', 'discharge (mm/hr)']].dropna()
# df['Q'] = df['discharge (mm/hr)']
# df['S_bin'] = pd.qcut(df['sscale'], q=10, duplicates='drop')
#binned = (
#     df
#     .groupby('S_bin')
#     .agg(
#         S_scale_mid=('sscale', 'median'),
#         Q_mean=('Q', 'mean'),
#         Q_std=('Q', 'std'),
#         Q_p90=('Q', lambda x: np.percentile(x, 90)),
#         Q_p95=('Q', lambda x: np.percentile(x, 95)),
#         n=('Q', 'size')
#     )
#     .reset_index()
# )
# variance test: should increase monotonically
# plt.figure()
# plt.plot(binned['S_scale_mid'], binned['Q_std'], marker='o')
# plt.xlabel('S_scale')
# plt.ylabel('std(Q)')
# plt.title('Discharge variability vs storage scale')
# plt.show()
# upper-tail test: confirms scale behavior
# plt.figure()
# plt.plot(binned['S_scale_mid'], binned['Q_std'], marker='o')
# plt.xlabel('S_scale')
# plt.ylabel('std(Q)')
# plt.title('Discharge variability vs storage scale')
# plt.show()

#%%
#--------set params_init---------
params_init = c18O_old, a, lamda, S_c, et_scale #***edit for distribution type***
params=params_init

#consider normalizing parameters to obtain better convergence of optimization
# Beta notes: a<1, b=1 young water prefernce, a=1,b<1 old water preference, a=b=1 uniform selection
#Gamma params: s_0=5701.46684 [c18O_old, a, et_storage] RMSE=0.5620217558825428 for params = [-7.60000000e+00  9.09859331e-01  1.58150050e+03]
#Kumaraswamy params: s_0=5701.46684 [c18O_old, a, b, et_storage] RMSE = = 0.5576700088867051 for params = [-7.600000e+00  9.098500e-01  9.923319e-01  1.581500e+03]



#%%
# ---------------------Optimize parameters ----------------------------
# Then we can optimize by fmin
from scipy.optimize import fmin
# params = fmin(minimize_me, params_init, ftol=1e-1, maxiter=100, disp=True)
# params = fmin(maximize_me, params_init, maxiter=100, disp=True)


# Try basinhopping
from scipy.optimize import basinhopping
# params = basinhopping(maximize_me, params_init, niter=25, T=0.01)
params = basinhopping(minimize_me, params_init, niter=25, T=0.01, stepsize=0.1)

#%%
# ------------------Try Monte Carlo sampling----------------------
from tqdm import tqdm
from scores.continuous import nse
import xarray as xr
normalized_storage = df['storage (mm)']-df['storage (mm)'].mean()
R=normalized_storage.max()-normalized_storage.min()
# param bounds
# bounds = {# for gamma model
#     "c18O_old": (-7.819719, -7.5),
#     "a": (0.3, 5.0),
#     "lamda": (0.01, 2), #(10, 13),
#     "S_c": (int(normalized_storage.min()-.75*R), int(normalized_storage.min()-.05*R)), #(-300, -150), #(-330, -310), normlized storage to mean.min()-.75*range to just min-5, scale shrinks to 0 as catchment approaches extreme dryness
#     "et_scale": (5, 100)
# }
bounds = {# for gamma split model
    "c18O_old": (-16.9, -0.445),
    "a_qf": (0.1, 1.0),      # quickflow: young water, shorter travel time
    "t_qf": (0.1, 1.0),      # quickflow: young water, shorter travel time
    "a_bf": (1.0, 5.0),      # baseflow: old water, longer travel time
    "lamda": (0.01, 2.0),
    "S_c": ((df['storage (mm)'].min()-.75*R), (df['storage (mm)'].min()-.05*R)), #(int(normalized_storage.min()-.75*R), int(normalized_storage.min()-.05*R)),
    "et_scale": (5, 100)
}
param_names = list(bounds.keys())

def sample_params(bounds):
    return [np.random.uniform(low, high) for low, high in bounds.values()]

def evaluate_model(params):
    try:
        model=make_gamma_split_model_from(params)
        model.run()

        obs = model.data_df['ORPB 18O'][issample].to_numpy()
        pred = model.data_df['precip 18O --> discharge (mm/hr)'][issample].to_numpy()
        del model #free up memory
        if np.any(np.isnan(pred)) or np.any(np.isinf(pred)):
            print('nan results in pred or obs')
            return np.nan
        RMSE = np.sqrt(np.mean((pred-obs)**2))
        obs_xr = xr.DataArray(obs)
        pred_xr = xr.DataArray(pred)
        NSE = nse(pred_xr, obs_xr)
        return RMSE, NSE.item()
    except Exception as e: # catches unstable parameter combinations
        return np.nan

N = 1000 #samples
results = []
for _ in tqdm(range(N)):
    params = sample_params(bounds)
    result = evaluate_model(params)
    if result is not np.nan and not np.isnan(result[0]):
        RMSE, NSE = result
        results.append((params + [RMSE, NSE]))

cols = param_names + ['RMSE','NSE']
results_df = pd.DataFrame(results, columns=cols)
print(results_df.describe())


#%%
# load results_df from Rockfish
results_df = pd.read_csv('results_df.csv')
results_df.drop(columns='Unnamed: 0', inplace=True) # drop index column if it exists

threshold = results_df['RMSE'].quantile(0.05) # top 5% of samples
thresholdnse = results_df['NSE'].quantile(0.95)
best = results_df[(results_df['RMSE'] <= threshold) & (results_df['NSE']>=thresholdnse)]
print("RMSE threshold: ", threshold)
print(best.describe())

best.hist(bins=20, figsize=(12, 8))
plt.tight_layout()
plt.legend()
plt.show()
#%%
#extract results and visualize
# p = best.loc[best['RMSE']==best['RMSE'].min()]
p = best.loc[best['NSE']==best['NSE'].max()].iloc[0] # take first row if multiple have same max NSE
p = p.to_numpy()
p = np.delete(p, -1) # take off RMSE column for params
p = np.delete(p, -1) # take of NSE column
print(p)
model = make_gamma_split_model_from(p) #***edit which distribution***
model.run()


#%%
# ------------------Build the model --------------------
# Now build a model with parameters
from mesas.sas.model import Model
model = make_gamma_model_from(params) #***edit which distribution***
model.run()

#%%------------------Save model to pickle-------------------
import pickle
pickle.dump(model, open('Putnam_model_n1.01.50.950.2.pkl', 'wb'))
# model = pickle.load(open('beta_model.pkl', 'rb'))
# model = pd.read_pickle('gamma_model_constS0.pkl')



#%%
# ------------------Visualize the output-------------------
# data_df = model.data_df #don't forget to do this
from mesas.utils import vis
fig = plt.figure(figsize=[14,4])

# plt.plot(isbaseflow, model.data_df['ORPB 18O'][isbaseflow],'-', label='Observed 18O baseflow')
# pred = model.data_df['bf1_weight'][isbaseflow].to_numpy()*(model.data_df['precip 18O --> baseflow 1 (mm/hr)'][isbaseflow].to_numpy())
# plt.plot(isbaseflow, pred,'-', label='Predicted 18O baseflow')

plt.plot(model.data_df.index[issample], model.data_df['ORPB 18O'][issample],'.-', color='grey', label='Observed 18O outflow')

# ----for quickflow/baseflow split
# pred = (1-model.data_df['bf1_weight'][issample].to_numpy())*(model.data_df['precip 18O --> quickflow (mm/hr)'][issample].to_numpy()) + model.data_df['bf1_weight'][issample].to_numpy()*(model.data_df['precip 18O --> baseflow 1 (mm/hr)'][issample].to_numpy())
# plt.plot(model.data_df.index[issample], pred,'-', color='cyan', label='Predicted 18O outflow')
# plt.plot(model.data_df.index[issample], model.data_df['precip 18O --> quickflow (mm/hr)'][issample], '.', label='Predicted 18O quickflow')
# plt.plot(model.data_df.index[issample], model.data_df['precip 18O --> baseflow 1 (mm/hr)'][issample], label='Predicted 18O baseflow 1')

# ----for combined discharge
plt.plot(model.data_df.index, model.data_df['precip 18O --> discharge (mm/hr)'], '.-', color='orange', alpha=0.3,label='Predicted 18O outflow')
# plt.axvspan(pd.Timestamp('2014-10-01'), pd.Timestamp('2015-09-30'), color='lightgrey', alpha=0.5)
# plt.axvspan(pd.Timestamp('2016-10-01'), pd.Timestamp('2017-09-30'), color='lightgrey', alpha=0.5)
# plt.axvspan(pd.Timestamp('2018-10-01'), pd.Timestamp('2019-09-30'), color='lightgrey', alpha=0.5)
plt.legend()
plt.title('Isotope outflow at ORPB')
# plt.xlim([pd.Timestamp('2014-01-01'), pd.Timestamp('2015-12-31')])
# plt.xlim([pd.Timestamp('2017-01-01'), pd.Timestamp('2018-12-31')])
# plt.xlim([pd.Timestamp('2016-03-10'), pd.Timestamp('2016-03-25')])
# plt.xlim([pd.Timestamp('2014-08-01'), pd.Timestamp('2016-08-31')])




# ax1 = plt.subplot2grid((1,2), (0,0))
# vis.plot_SAS_cumulative(model, 'discharge (mm/hr)', ax=ax1)
# ax1.set_title('Cumulative SAS storage')

# ax2 = plt.subplot2grid((1,2), (0,1))
# ax2.plot(model.data_df.index, model.data_df['ORPB 18O'],'.', label='Observed 18O outflow')
# ax2.plot(model.data_df.index, model.data_df['precip 18O --> discharge (mm/hr)'], label='Predicted 18O outflow')
# ax2.legend()
# ax2.set_title('Isotope outflow at ORPB')

#%%
# Plot TTD
pq = model.get_pQ('discharge (mm/hr)')
cumTTD = np.cumsum(pq, axis=0) * model.options['dt']
arS = model.get_ST()[:,1]
ars = model.get_sT()[:,1]
time = model.options['dt']*np.arange(model.options['max_age'])
fig = plt.figure(figsize=[11,4])
ax1 = plt.subplot2grid((1,2), (0,0))
ax1.plot(arS, cumTTD)
ax1.set_ylim([0, 1])
ax1.set_xlim(xmin=0)
ax1.plot(ax1.get_xlim(), [1, 1], color='0.1', lw=0.8, ls=':')
ax1.plot(ax1.get_xlim(), [0, 0], color='0.1', lw=0.8, ls=':')
ax1.set_xlabel('$S_T$')
ax1.set_ylabel('$\Omega(S_T)$')
ax1.set_title('Cumulative TTD')

ax2 = plt.subplot2grid((1,2), (0,1))
ax2.plot(arS, pq)
ax2.set_xlabel('$S_T$')
ax2.set_ylabel('$\omega(S_T)$')
ax2.set_title('TTD')
plt.tight_layout()




# %%
