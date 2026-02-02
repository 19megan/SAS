# This script will scale the resolution of isotope timeseries from ORPB_isotope_data.csv

# Date: 10/01/2025

#%%
# ----------------import dataset----------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv('ORPB_isotope_data.csv', index_col=0, parse_dates=[0]) #len 51745
# data = data.loc[pd.Timestamp('2014-01-01'): pd.Timestamp('2014-07-01')]

# clean up and find missing samples
data = data.drop(data[data['Sample Name'].isna()].index)
issample = ~data['Sample Name'].duplicated(keep='last')
data['is_weekly'] = issample
#%% ----------------stream isotopes-------------------------------------------------
# note: if using stream isotopes, need to change the mass flux calculation to use discharge instead of precipitation
# isotopes = data[['ORPB 2H', 'ORPB 18O', 'ORPB 17O']]
# iso = 'ORPB 2H'
# # ---------------plot original data----------------
# # plot fluxes and concentrations - Esther's Figure 4.4
# fig,(ax1,ax2, ax3)=plt.subplots(nrows=3,ncols=1,figsize=[10,9])
# ax1.plot(data.index, data['discharge (mm/hr)'], color='blue', label='discharge (mm/hr)')
# ax1.set_xlabel('Date')
# ax1.set_ylabel('Discharge (mm/hr)', color='blue')
# ax1.tick_params(axis='y', labelcolor='blue')
# ax12 = ax1.twinx()
# ax12.bar(data.index, data['rainfall (mm/hr)'], color='gray', alpha=0.6, width=0.01, label='rainfall (mm/hr)')
# ax12.bar(data.index, data['snowfall SWE (mm/hr)'], color='black', alpha=0.6, width=0.01, label='snowfall SWE (mm)')
# ax12.set_ylabel('(mm/hr)', color='gray')
# ax12.tick_params(axis='y', labelcolor='gray')
# ax12.set_ylim([0,0.75])#max(max(data_df['rainfall (mm/hr)']), max(data_df['snowfall SWE (mm/hr)']))*1.1])
# ax12.invert_yaxis()
# ax1.set_title('ORPB Discharge and Rainfall')
# ax1.legend()
# ax12.legend()

# issample = np.logical_not(np.isnan(data[iso]))
# ax2.plot(data.index[issample], 
#          data[iso][issample],'.', label=f'hourly {iso}')
# # aggregate to weekly Esther Equation 4.1
# isotopes_weekly = isotopes[issample].resample('W').mean()
# rain = data['rainfall (mm/hr)'][issample]
# num = (isotopes[iso][issample]*rain).resample('W').sum()
# den = rain.resample('W').sum()
# isotopes_weekly[iso] = num/den
# ax2.plot(isotopes_weekly.index, isotopes_weekly[iso], '.', label=f'weekly agg {iso}')
# ax2.legend()
# ax2.set_title('ORPB Tracer at stream')
# ax3.plot(data.index,
#          data[f'precip {iso[-2:]}'], label=f'hourly precip {iso[-2:]}')
# # aggregate to weekly Esther Equation 4.1
# pisotopes_weekly = isotopes.resample('W').mean()
# num = (data[f'precip {iso[-2:]}']*rain).resample('W').sum()
# den = rain.resample('W').sum()
# pisotopes_weekly[f'precip {iso[-2:]}'] = num/den
# ax3.plot(pisotopes_weekly.index, pisotopes_weekly[f'precip {iso[-2:]}'], '-', label=f'weekly agg precip {iso[-2:]}')
# ax3.legend()
# ax3.set_title('ORPB Tracer in precipitation')
# plt.tight_layout()
# plt.show()


#%%-----------------precip isotopes-------------------------------------------------
isotopes = data[['precip 2H', 'precip 18O', 'precip 17O']] # we want to scale the input isotopes
iso = 'precip 18O' #'precip 2H'

def weighted_average(df, col='precip 2H', agg_freq='W'):
    weighted_avg = df[col].mul(df['rainfall (mm/hr)']).resample(agg_freq).sum() / df['rainfall (mm/hr)'].resample(agg_freq).sum()
    return weighted_avg

#%%
# ---------------plot original data and aggregations----------------
# plot fluxes and concentrations - Esther's Figure 4.4
fig,(ax1,ax2, ax3)=plt.subplots(nrows=3,ncols=1,figsize=[10,9])
ax1.plot(data.index, data['discharge (mm/hr)'], color='blue', label='discharge (mm/hr)')
ax1.set_xlabel('Date')
ax1.set_ylabel('Discharge (mm/hr)', color='blue')
ax1.tick_params(axis='y', labelcolor='blue')
ax12 = ax1.twinx()
ax12.bar(data.index, data['rainfall (mm/hr)'], color='gray', alpha=0.6, width=0.01, label='rainfall (mm/hr)')
ax12.bar(data.index, data['snowfall SWE (mm/hr)'], color='black', alpha=0.6, width=0.01, label='snowfall SWE (mm)')
ax12.set_ylabel('(mm/hr)', color='gray')
ax12.tick_params(axis='y', labelcolor='gray')
ax12.set_ylim([0,0.75])#max(max(data_df['rainfall (mm/hr)']), max(data_df['snowfall SWE (mm/hr)']))*1.1])
ax12.invert_yaxis()
ax1.set_title('ORPB Discharge and Rainfall')
ax1.legend()
ax12.legend()

isPsample = issample & data[iso].notna()
ax2.plot(data.index[isPsample], data[iso][isPsample], '.', label=f'weekly {iso}')
# aggregate to monthly Esther Equation 4.1
piso_agg = weighted_average(data.loc[isPsample], col=iso, agg_freq='BME')
ax2.plot(piso_agg.index, piso_agg, '.', markersize=13, label=f'bimonthly agg {iso}')
ax2.set_xlabel('Date')
ax2.set_ylabel('Concentration [‰]')
ax2.legend()
ax2.set_title('ORPB Tracer in precipitation')

isSsample = np.logical_not(np.isnan(data[f'ORPB {iso[-2:]}']))
ax3.plot(data.index[isSsample],
         data[f'ORPB {iso[-2:]}'][isSsample], '.', label=f'observed mixed freq ORPB {iso[-3:]}')
# aggregate using Esther Equation 4.1
siso_agg = weighted_average(data[isSsample], col=f'ORPB {iso[-2:]}', agg_freq='W')
ax3.plot(siso_agg.index, siso_agg, '.', markersize=7, label=f'weekly agg ORPB {iso[-3:]}')
ax3.set_xlabel('Date')
ax3.set_ylabel('Concentration [‰]')
ax3.legend()
ax3.set_title('ORPB Tracer at stream')
plt.tight_layout()
plt.show()

#%%
# ---------------down/upsample to bimonthly (BME), monthly (ME), daily (D), etc.----------------
agg_freq = 'W'
pisotopes_agg = isotopes.resample(agg_freq).mean()
pisotopes_obs = isotopes.resample(agg_freq).last() # at weekly cadence
pisotopes_agg[iso] = weighted_average(data, col=iso, agg_freq=agg_freq)
fig, ax = plt.subplots()
ax.scatter(pisotopes_obs[iso], pisotopes_agg[iso], label=f'$\Delta {iso[-2:]}$')
ax.set_xlabel(f'Observed {iso} at {agg_freq} cadence [‰]')
ax.set_ylabel(f'{agg_freq} aggregation {iso}')
ax.plot()
oneone = [min(pisotopes_obs[iso].min(), pisotopes_agg[iso].min()), 
          max(pisotopes_obs[iso].max(), pisotopes_agg[iso].max())]
ax.plot(oneone, oneone, 'k--', alpha=0.75, zorder=0, label='1:1 line')
ax.set_title(f'{iso}: {agg_freq} aggreation vs observations')
ax.legend()
plt.show() # Like Xu Fei Figure 4.5
print('Correlation: ', pisotopes_obs[iso].corr(pisotopes_agg[iso]))


# %%
# ---------------upsample using GP minus seasonal trend----------------

# First, define data with coarse resolution to new resolution with weighted average on the isotope data
c_res = 'D'
precip = data['rainfall (mm/hr)'].resample(c_res).sum()
precip.name = 'rainfall (mm/hr)'
precip = precip.apply(lambda x: round(x/2.54e-3)*2.54e-3)
c_iso = data[[iso, 'Sample Name']].asfreq(c_res) #just reindex no change to values
c_iso['weekly_obs'] = ~c_iso['Sample Name'].duplicated(keep='last')
df = pd.concat([precip, c_iso],axis=1).loc[c_iso.index]
df.ffill(inplace=True)
df['cumP'] = df['rainfall (mm/hr)'].cumsum()
df['cumP'] = df['cumP'].apply(lambda x: round(x/2.54e-3)*2.54e-3)

# define aggregation resolution
agg_res = c_res
#Use Method II: sinusoidal trend
from scipy.optimize import curve_fit

P_n = df.resample(agg_res).sum()['rainfall (mm/hr)']
iso_n = df.resample(agg_res).apply(weighted_average, col=iso, agg_freq=agg_res)
resampled = pd.concat([P_n, iso_n], axis=1)
resampled.columns = ['rainfall (mm/hr)', iso]
resampled.ffill(inplace=True)
resampled.bfill(inplace=True)
t_agg = resampled.index.dayofyear

# this makes a continuous day counter across multiple years
def seasonal_sine(t, A, phi, b):
    return A * np.sin(2 * np.pi * t / 365 - phi) + b

params_init = [5, 1, -7.5] # initial guess
params, cov = curve_fit(seasonal_sine, t_agg, resampled[iso], params_init)
A_fit, phi_fit, b_fit = params
iso_fit = seasonal_sine(t_agg, *params) # This is c_bar in Xu Fei eq. 4.7

# Plot fitted sine
plt.figure(figsize=[10,3])
plt.scatter(resampled.index, resampled[iso], marker='.', label=f'Observed {agg_res}', alpha = 0.7)
plt.plot(resampled.index, iso_fit, 'r--', label='Fitted sine')
plt.xlabel('Date')
plt.ylabel(f'{iso} [‰]')
plt.title(f'{iso}: Seasonal sine fit')
plt.legend()
plt.show() # Like Xu Fei Figure 4.6
# print('Correlation: ', isotopes['ORPB 18O'].corr(iso_fit))

#%%
# take out seasonality from isotope data - Esther Figure 4.6
resampled[f'{iso} deseasoned'] = resampled[iso] - iso_fit

plt.figure(figsize=[15,5])
plt.plot(resampled.index, resampled[f'{iso} deseasoned'], marker='.', label=f'Deseasoned {iso}')
plt.plot(resampled.index, resampled[iso]-resampled[iso].mean(), marker='.', label=f'Standardized {iso}', alpha=0.5)
plt.xlabel('Date')
plt.ylabel(f'{iso} [‰]')
plt.title(f'{iso}: Deseasoned vs standardized')
plt.legend()
plt.show()

print('Correlation: ', (resampled[iso]-resampled[iso].mean()).corr(resampled[f'{iso} deseasoned']))

#%%
# plot mass fluxes of deseasoned and standardized istotope data - Esther Figure 4.6
resampled[f'mass standrd {iso}'] = (resampled[iso]-resampled[iso].mean())*resampled['rainfall (mm/hr)'] #['discharge (mm/hr)'] #discharge for stream isotopes
resampled[f'mass {iso} deseasoned'] = resampled[f'{iso} deseasoned']*resampled['rainfall (mm/hr)'] #['discharge (mm/hr)']

plt.figure(figsize=[10,5])
plt.plot(resampled.index, resampled[f'mass {iso} deseasoned'], marker='.', label=f'Mass flux deseasoned {iso}')
plt.plot(resampled.index, resampled[f'mass standrd {iso}'], marker='.', label=f'Mass flux standardized {iso}', alpha=0.5)
plt.xlabel('Date')
plt.ylabel(f'Mass flux {iso} mm')
plt.title(f'{iso}: Mass flux deseasoned vs standardized')
# plt.ylim([-1000,1000])
plt.legend()
plt.show()

print('Correlation: ', (resampled[f'mass standrd {iso}']).corr(resampled[f'mass {iso} deseasoned']))

# %%
# Now fit GP regressor to M(t) and precip
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel as C, Matern

df_gp = pd.concat([df['rainfall (mm/hr)'],df['cumP'], resampled[f'{iso} deseasoned'], df['weekly_obs']], axis=1)
precip_min, precip_max = df_gp.loc[df_gp['rainfall (mm/hr)'] != 0, 'rainfall (mm/hr)'].min(), df_gp['cumP'].max() #st and et of cumP
df_gp[f'cum mass {iso} deseasoned'] = (df_gp[f'{iso} deseasoned']*df_gp['rainfall (mm/hr)']).cumsum()
continuous_precip = pd.DataFrame(np.arange(0.0, precip_max+precip_min, precip_min), columns=['precip_continuous']) # creates an integral cumP with equal intervals for better GP prediction
df_gp = df_gp.merge(continuous_precip, how='right', left_on='cumP', right_on='precip_continuous')
df_gp[['cumP', f'cum mass {iso} deseasoned']] = df_gp[['cumP', f'cum mass {iso} deseasoned']].interpolate('spline', order=2)
df_gp['weekly_obs'] = df_gp['weekly_obs'].fillna(False)
X = df_gp['precip_continuous'][df_gp['weekly_obs']].values.reshape(-1,1) #convert to 2D array
y = df_gp[f'cum mass {iso} deseasoned'][df_gp['weekly_obs']].values

# Kernel = constant * RBF + noise
kernel = C(1.0, (1e-2, 1e2)) * RBF(length_scale=1.0, length_scale_bounds=(1e-2, 1e2)) \
         + WhiteKernel(noise_level=1e-3, noise_level_bounds=(1e-5, 1e1))

gp = GaussianProcessRegressor(kernel=Matern(), alpha=0.0000001, n_restarts_optimizer=50, normalize_y=True)

# Fit GP to data
gp.fit(X, y)

# Predict on test points
Xtest = df_gp['precip_continuous'].values.reshape(-1, 1)
y_mean, y_std = gp.predict(Xtest, return_std=True)


#%%
# ---------------plot GP predicted deseasoned mass fluxes----------------
plt.figure(figsize=[10,5])
plt.plot(df_gp.index, y_mean, '.-', label='GP predicted mass flux deseasoned')
plt.fill_between(df_gp.index, y_mean - 2*y_std, y_mean + 2*y_std, 
                 alpha=0.2, color="blue", label="Confidence interval")
plt.plot(df_gp.index, df_gp[f'cum mass {iso} deseasoned'], marker='.', zorder=10, label=f'Observed mass flux deseasoned {iso}', alpha=0.2, color='orange')
plt.xlabel('Date')
plt.ylabel(f'Mass flux deseasoned {iso} [‰ mm/hr]')
plt.title(f'{iso}: GP predicted mass flux deseasoned')
plt.legend()
plt.show()

# %%
# ---------------plot GP predicted isotope fluxes----------------
# map predictions back to df.index
time_map = np.round(df['cumP']/precip_min) * precip_min
df_rev = pd.DataFrame(
    {
        "mean": np.round(y_mean, 5),
        "upper": np.round(y_mean+1.96*y_std, 5),
        "lower": np.round(y_mean-1.96*y_std, 5),
        "precip": Xtest.flatten(),
    }
).drop_duplicates()
df_rev = df_rev.merge(pd.DataFrame(time_map), left_on='precip', right_on='cumP', how='right')
df_rev.index = df.index
df_rev['rainfall (mm/hr)'] = df_rev['precip'].diff()
df_rev['rainfall (mm/hr)'][0] = 0.0
df_rev['iso_fit'] = iso_fit
df_rev['mean_c'] = df_rev['mean'] + (df_rev['iso_fit']*df_rev['rainfall (mm/hr)']).cumsum()
df_rev['mean_c'] = df_rev['mean_c'].diff()/df_rev['rainfall (mm/hr)']
df_rev.loc[df_rev['rainfall (mm/hr)']==0, 'mean_c'] = np.nan
df_rev['upper_c'] = df_rev['upper'] + (df_rev['iso_fit']*df_rev['rainfall (mm/hr)']).cumsum()
df_rev['upper_c'] = df_rev['upper_c'].diff()/df_rev['rainfall (mm/hr)']
df_rev.loc[df_rev['rainfall (mm/hr)']==0, 'upper_c'] = np.nan
df_rev['lower_c'] = df_rev['lower'] + (df_rev['iso_fit']*df_rev['rainfall (mm/hr)']).cumsum()
df_rev['lower_c'] = df_rev['lower_c'].diff()/df_rev['rainfall (mm/hr)']
df_rev.loc[df_rev['rainfall (mm/hr)']==0, 'lower_c'] = np.nan

plt.figure(figsize=[15,5])
plt.scatter(df_rev.index, df_rev['mean_c'], label=f'GP predicted {iso} {agg_res}', marker='.', alpha=0.3, zorder=10)
plt.fill_between(df_rev.index, df_rev['lower_c'], df_rev['upper_c'], 
                 alpha=0.1, color="purple", label="Confidence interval")
plt.scatter(df.index[(df['weekly_obs']==True) & (df['rainfall (mm/hr)']>0)], df.loc[(df['weekly_obs']==True) & (df['rainfall (mm/hr)']>0), iso], color='orange', marker='.', s=10**2, label=f'Observed {iso}')
plt.xlabel('Date')
plt.ylabel(f'{iso} concentration [‰]')
plt.title(f'{iso}: GP predicted isotope concentration')
plt.ylim([isotopes[iso].min()-5, isotopes[iso].max()+5])
plt.legend()
plt.show()



#%% #fix this later***
#-------------------Create resampled data with GP predicted isotope values----------------
# assign gp predicted isotope values to resampled dataframe

P_n = data.resample(agg_res).sum()['rainfall (mm/hr)']
S_n = data.resample(agg_res).sum()['snowmelt (mm/hr)']
D_n = data.resample(agg_res).sum()['discharge (mm/hr)']
ET_n = data.resample(agg_res).sum()['ET (mm/hr)']
ORPB18O_n = data.resample(agg_res)['ORPB 18O'].mean()
resampled = pd.concat([P_n, S_n, D_n, ET_n, ORPB18O_n], axis=1)
resampled.columns = [f'rainfall (mm/{agg_res})', f'snowmelt (mm/{agg_res})', f'discharge (mm/{agg_res})', f'ET (mm/{agg_res})', 'ORPB 18O']
resampled.ffill(inplace=True)
resampled.bfill(inplace=True)
resampled[f'{iso}'] = y_mean_iso
resampled[f'GP upper {iso}'] = y_upper
resampled[f'GP lower {iso}'] = y_lower

# %%
# save to csv files

resolution = 'monthly' #bimonthly, weekly, daily, monthly #EDIT HERE depending on agg_res

resolution_dataset_root = '/Users/simon/Desktop/ORPB_resolution_datasets'
import os
if not os.path.exists(resolution_dataset_root):
    os.makedirs(resolution_dataset_root)

resampled.to_csv(f"{resolution_dataset_root}/ORPB_isotope_data_{resolution}.csv")


# %%
