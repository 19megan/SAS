# This script will scale the resolution of isotope timeseries from ORPB_isotope_data.csv

# Date: 10/01/2025

#%%
# ----------------import dataset----------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv('ORPB_isotope_data.csv', index_col=0, parse_dates=[0])
data[['rainfall (mm/hr)', 'precip 2H', 'precip 18O', 'precip 17O', 'ORPB 2H', 'ORPB 18O', 'ORPB 17O']] = data[['rainfall (mm/hr)', 'precip 2H', 'precip 18O', 'precip 17O', 'ORPB 2H', 'ORPB 18O', 'ORPB 17O']].astype(np.float32)
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
iso = 'precip 2H'

# ---------------plot original data----------------
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

issample = np.logical_not(np.isnan(data[iso]))
ax2.plot(data.index, data[iso], label=f'hourly {iso}')
# aggregate to weekly Esther Equation 4.1
pisotopes_weekly = isotopes.resample('W').mean()
rain = data['rainfall (mm/hr)'][issample]
num = (isotopes[iso]*rain).resample('W').sum()
den = rain.resample('W').sum()
pisotopes_weekly[iso] = num/den
ax2.plot(pisotopes_weekly.index, pisotopes_weekly[iso], '-', label=f'weekly agg {iso}')
ax2.legend()
ax2.set_title('ORPB Tracer in precipitation')

isSsample = np.logical_not(np.isnan(data[f'ORPB {iso[-2:]}']))
ax3.plot(data.index[isSsample],
         data[f'ORPB {iso[-2:]}'][isSsample], '.', label=f'hourly ORPB {iso[-3:]}')
# aggregate to weekly Esther Equation 4.1
isotopes_weekly = isotopes[isSsample].resample('W').mean()
rain = data['rainfall (mm/hr)'][isSsample]
num = (data[f'ORPB {iso[-2:]}'][isSsample]*rain).resample('W').sum()
den = rain.resample('W').sum()
isotopes_weekly[f'ORPB {iso[-2:]}'] = num/den #weighted average
ax3.plot(isotopes_weekly.index, isotopes_weekly[f'ORPB {iso[-2:]}'], '.', label=f'weekly agg ORPB {iso[-3:]}')
ax3.legend()
ax3.set_title('ORPB Tracer at stream')
plt.tight_layout()
plt.show()

#%%
# ---------------downsample to weekly----------------
pisotopes_weekly = isotopes.resample('W').mean()
pisotopes_hourly = isotopes.resample('W').last() # at weekly cadence
# pisotopes_weekly[iso] = (data[iso].mul(data['rainfall (mm/hr)']).resample('W').sum()/data['rainfall (mm/hr)'].resample('W').sum())
# pisotopes_hourly[iso] = (data[iso].mul(data['rainfall (mm/hr)']).resample('W').last()/data['rainfall (mm/hr)'].resample('W').last())
fig, ax = plt.subplots()
ax.scatter(pisotopes_hourly[iso], pisotopes_weekly[iso], label=f'$\Delta {iso[-2:]}$')
ax.set_xlabel(f'Hourly {iso} at weekly cadence [‰]')
ax.set_ylabel(f'Weekly aggregation {iso}')
ax.plot()
oneone = [min(pisotopes_hourly[iso].min(), pisotopes_weekly[iso].min()), 
          max(pisotopes_hourly[iso].max(), pisotopes_weekly[iso].max())]
ax.plot(oneone, oneone, 'k--', alpha=0.75, zorder=0, label='1:1 line')
ax.set_title(f'{iso}: Weekly aggreation vs observations')
ax.legend()
plt.show() # Like Xu Fei Figure 4.5
print('Correlation: ', pisotopes_hourly[iso].corr(pisotopes_weekly[iso]))


# %%
# ---------------upsample weekly to hourly----------------

# first agg to make data weekly
rain = data['rainfall (mm/hr)'].resample('W').sum()
rain.name = 'rainfall (mm/hr)'
pisotopes_weekly = isotopes.asfreq('W')
data_weekly = pd.concat([rain, pisotopes_weekly], axis=1)[pisotopes_weekly.index[0]:pisotopes_weekly.index[-1]]
data_weekly.fillna(method='ffill', inplace=True)
#*********change to data_weekly below todo

#%%
#Use Method II: sinusoidal trend
from scipy.optimize import curve_fit

issample = np.logical_not(np.isnan(isotopes_weekly[iso]))
issample_h = np.logical_not(np.isnan(isotopes[iso]))
t_w = isotopes_weekly.index[issample].dayofyear.values
t_w = isotopes_weekly.index[issample].dayofyear# + (isotopes_weekly.index[issample].year - isotopes_weekly.index[issample].year.min())*365
# this makes a continuous day counter across multiple years
def seasonal_sine(t, A, phi, b):
    return A * np.sin(2 * np.pi * t / 365 - phi) + b

params_init = [1, 1, -7.5] # initial guess

params, cov = curve_fit(seasonal_sine, t_w, isotopes_weekly[iso][issample], params_init)
A_fit, phi_fit, b_fit = params
# Now evaluate at hourly timepoints
t_h = isotopes.index.dayofyear# + (isotopes.index.year - isotopes.index.year.min())*365
iso_fit = seasonal_sine(t_h, *params) # This is c_bar in Xu Fei eq. 4.7
iso_fit = iso_fit.astype(np.float32)
plt.figure(figsize=[10,5])
plt.scatter(isotopes_weekly.index[issample], isotopes_weekly[iso][issample], marker='.', label='Observed weekly', alpha = 0.7)
plt.plot(isotopes.index, iso_fit, 'r--', label='Fitted sine')
plt.xlabel('Date')
plt.ylabel(f'{iso} [‰]')
plt.title(f'{iso}: Seasonal sine fit')
plt.legend()
plt.show() # Like Xu Fei Figure 4.6
# print('Correlation: ', isotopes['ORPB 18O'].corr(iso_fit))

#%%
# take out seasonality from isotope data - Esther Figure 4.6
data[f'{iso} deseasoned'] = data[iso] - iso_fit

plt.figure(figsize=[15,5])
plt.plot(data.index, data[f'{iso} deseasoned'], marker='.', label=f'Deseasoned {iso}')
plt.plot(data.index, data[iso]-data[iso].mean(), marker='.', label=f'Original {iso} - mean', alpha=0.5)
plt.xlabel('Date')
plt.ylabel(f'{iso} [‰]')
plt.title(f'{iso}: Deseasoned vs original-mean')
plt.legend()
plt.show()

print('Correlation: ', (data[iso]-data[iso].mean()).corr(data[f'{iso} deseasoned']))

#%%
# plot mass fluxes of deseasoned and standardized istotope data - Esther Figure 4.6
data[f'mass standrd {iso}'] = (data[iso]-data[iso].mean())*data['rainfall (mm/hr)'] #['discharge (mm/hr)'] #discharge for stream isotopes
data[f'mass {iso} deseasoned'] = data[f'{iso} deseasoned']*data['rainfall (mm/hr)'] #['discharge (mm/hr)']
data[f'mass standrd {iso}'] = data[f'mass standrd {iso}'].astype(np.float32)
data[f'mass {iso} deseasoned'] = data[f'mass {iso} deseasoned'].astype(np.float32)
plt.figure(figsize=[10,5])
plt.plot(data.index, data[f'mass {iso} deseasoned'], marker='.', label=f'Mass flux deseasoned {iso}')
plt.plot(data.index, data[f'mass standrd {iso}'], marker='.', label=f'Mass flux standardized {iso}', alpha=0.5)
plt.xlabel('Date')
plt.ylabel(f'Mass flux {iso} [‰ mm/hr]')
plt.title(f'{iso}: Mass flux deseasoned vs standardized')
# plt.ylim([-5,5]) #only for stream isotopes
plt.legend()
plt.show()

print('Correlation: ', (data[f'mass standrd {iso}']).corr(data[f'mass {iso} deseasoned']))

#%%
# ---------------Find M(t)----------------
# Now find true M(t) process with unknown discretized GP:
# data['precip (mm/hr)'] = data[['rainfall (mm/hr)','snowmelt (mm/hr)']].sum(axis=1)
# M = np.zeros(len(data[issample]))
# dM = np.zeros_like(M)
# for i in range(1, len(data[issample])):
#     if np.isnan(data[iso][issample].iloc[i]):
#         print('Not sample')
#         break
#     for k in range(1, i):
#         M[i] += (data[iso][issample].iloc[k] - iso_fit[k]) * (data['precip (mm/hr)'][issample].iloc[k] - data['precip (mm/hr)'][issample].iloc[k-1])
#     dM[i] = (data[iso][issample].iloc[i] - iso_fit[i]) * (data['precip (mm/hr)'][issample].iloc[i] - data['precip (mm/hr)'][issample].iloc[i-1])
    

# plt.figure(figsize=[10,5])
# plt.plot(data[issample].index, M, label='M(t)')
# plt.legend()

# %%
# Now fit GP regressor to M(t) and precip
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel as C

# Training data #I have to concat it to 1000 length to avoid memory issues
X = data['rainfall (mm/hr)'][issample_h][20000:20200].values.reshape(-1, 1) #convert to 2D array
y = data[f'mass {iso} deseasoned'][issample_h][20000:20200].values

# Kernel = constant * RBF + noise
kernel = C(1.0, (1e-2, 1e2)) * RBF(length_scale=1.0, length_scale_bounds=(1e-2, 1e2)) \
         + WhiteKernel(noise_level=1e-3, noise_level_bounds=(1e-5, 1e1))

gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=10)

# Fit GP to data
gp.fit(X, y)

# Predict on test points
Xtest = data['rainfall (mm/hr)'].values.reshape(-1, 1)
y_mean, y_std = gp.predict(Xtest, return_std=True)

# # Plot
# plt.figure()
# plt.plot(X, y, "kx", label="Training data")
# plt.plot(Xtest, y_mean, "b", label="Mean prediction")
# plt.fill_between(Xtest.ravel(), y_mean - 2*y_std, y_mean + 2*y_std, 
#                  alpha=0.2, color="blue", label="Confidence interval")
# plt.legend()
# plt.show()

# print("Optimized kernel:", gp.kernel_)

#%%
# ---------------plot GP predicted deseasoned mass fluxes----------------
plt.figure(figsize=[10,5])
plt.plot(data.index, y_mean, label='GP predicted mass flux deseasoned')
plt.fill_between(data.index, y_mean - 2*y_std, y_mean + 2*y_std, 
                 alpha=0.2, color="blue", label="Confidence interval")
plt.plot(data.index, data[f'mass {iso} deseasoned'], marker='.', zorder=-1, label=f'Observed mass flux deseasoned {iso}', alpha=0.2, color='orange')
plt.xlabel('Date')
plt.ylabel(f'Mass flux deseasoned {iso} [‰ mm/hr]')
plt.title(f'{iso}: GP predicted mass flux deseasoned')
plt.legend()
plt.show()

# %%
# ---------------plot GP predicted isotope fluxes----------------
y_mean_iso = y_mean/data['rainfall (mm/hr)'] + iso_fit #convert mass flux back to isotope concentration and re-add seasonality -- 'discharge (mm/hr)' for stream isotopes
y_std_iso = y_std/data['rainfall (mm/hr)']
plt.figure(figsize=[15,5])
plt.scatter(data.index, y_mean_iso, label='GP predicted precip hourly', marker='.',alpha=0.2)
plt.fill_between(data.index, y_mean_iso - 2*y_std_iso, y_mean_iso + 2*y_std_iso, 
                 alpha=0.1, color="purple", label="Confidence interval")
plt.scatter(isotopes_weekly.index[issample], isotopes_weekly[iso][issample], color='orange', marker='.', label='Observed weekly', alpha = 0.7)
plt.xlabel('Date')
plt.ylabel(f'{iso} concentration [‰]')
plt.title(f'{iso}: GP predicted isotope concentration')
plt.ylim([isotopes[iso].min()-5, isotopes[iso].max()+5])
plt.legend()
plt.show()

# assign gp predicted isotope values to dataframe
data[f'GP predicted {iso}'] = y_mean_iso
# %%

