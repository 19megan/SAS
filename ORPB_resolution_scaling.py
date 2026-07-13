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
data = data.drop(data[data['Sample Name'].isna()].index) #drop last week of data with no isotope samples and no rain - still keeps h res
issample = ~data.loc[data['rainfall (mm/hr)']>0, 'Sample Name'].duplicated(keep='last')
data['is_weekly'] = issample
data.loc[data['is_weekly'].isna(), 'is_weekly'] = False #when rain is 0, set to False
data['cumP'] = data['rainfall (mm/hr)'].cumsum()
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
iso = 'precip 18O' #'precip 18O' #'precip 2H'
if iso == 'precip 2H': #number of digits needed for isotope type
    ison = 2
elif iso == 'precip 18O':
    ison = 3
else:
    ison = 3

def weighted_average(df, col='precip 2H', agg_freq='W'):
    weighted_avg = df[col].mul(df['rainfall (mm/hr)']).resample(agg_freq).sum(min_count=1) / df['rainfall (mm/hr)'].resample(agg_freq).sum()
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

isSsample = np.logical_not(np.isnan(data[f'ORPB {iso[-ison:]}']))
ax3.plot(data.index[isSsample],
         data[f'ORPB {iso[-ison:]}'][isSsample], '.', label=f'observed mixed freq ORPB {iso[-ison:]}')
# aggregate using Esther Equation 4.1
siso_agg = weighted_average(data[isSsample], col=f'ORPB {iso[-ison:]}', agg_freq='W')
ax3.plot(siso_agg.index, siso_agg, '.', markersize=7, label=f'weekly agg ORPB {iso[-ison:]}')
ax3.set_xlabel('Date')
ax3.set_ylabel('Concentration [‰]')
ax3.legend()
ax3.set_title('ORPB Tracer at stream')
plt.tight_layout()
plt.show()

#%%
# ---------------down/upsample to bimonthly (BME), monthly (ME), daily (D), etc.----------------
agg_freq = 'W' #hourly (h), daily (D), weekly (W), biweekly (2W), monthly (ME)
pisotopes_agg = isotopes.resample(agg_freq).mean()
pisotopes_obs = isotopes.resample(agg_freq).last() # at weekly cadence
pisotopes_agg[iso] = weighted_average(data, col=iso, agg_freq=agg_freq)
fig, ax = plt.subplots()
ax.scatter(pisotopes_obs[iso], pisotopes_agg[iso], label=f'$\Delta {iso[-ison:]}$')
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
# ---------------upsample using GP minus trend----------------

# First, define data with coarse resolution to new resolution with weighted average on the isotope data
c_res = 'ME' #hourly (h), daily (D), weekly (W), biweekly (2W), monthly (ME) #EDIT HERE depending on agg_res
precip = data['rainfall (mm/hr)'].resample(c_res).sum()
precip.name = 'rainfall (mm/hr)'
precip = precip.apply(lambda x: round(x/2.54e-3)*2.54e-3)
c_iso = weighted_average(data, col=iso, agg_freq=c_res) #better to use weighted average to aggregate before downsampling
c_iso = pd.concat([c_iso,data['Sample Name'].asfreq(c_res)],axis=1) #just reindex no change to values
c_iso.columns = [iso, 'Sample Name']
df = pd.concat([precip, c_iso],axis=1).loc[c_iso.index]
df['weekly_obs'] = ~df.loc[df['rainfall (mm/hr)']>0, 'Sample Name'].duplicated(keep='last') # make sure weekly values are only for when rain was actually observed (1066 obs values vs only 235 when observed without rain restriction)
df.loc[df['weekly_obs'].isna(), 'weekly_obs'] = False #set non-masked values to False i.e. where no rain observed
df.loc[df['weekly_obs']==False, iso] = np.nan # set non-weekly observed isotope values to nan so they don't influence GP fit
df.bfill(inplace=True)
df.ffill(inplace=True) #to get last 23 rows of data that have zero rainfall
df['cumP'] = df['rainfall (mm/hr)'].cumsum()
df['cumP'] = df['cumP'].apply(lambda x: round(x/2.54e-3)*2.54e-3)

# define aggregation resolution
agg_res = 'h'

P_n = df.resample(agg_res).sum()['rainfall (mm/hr)']
iso_n = df[[iso, 'Sample Name']].asfreq(agg_res) #just reindex to finer res no change to values, does the same thing as just assigning the coarse res column to the finer res df
resampled = data.copy()
resampled['is_weekly'] = df['weekly_obs'] #makes mask match coarse data
resampled.loc[resampled['is_weekly'].isna(), 'is_weekly'] = False
resampled[iso] = iso_n[iso]
resampled.bfill(inplace=True)
resampled.ffill(inplace=True)
t_agg = resampled.index.dayofyear

#note: after checking data[iso] and resampled[iso], they are very different post 2019-06-10 because there is no rain, but data[iso] still somehow has observations.

#%%
# after this point, created backfilled interpolation, otherwise skip to trend removal for GP interpolation
# resampled['mean_c'] = resampled[iso]
# data = data.join(resampled['mean_c'])
# data['mean_c'] = data['mean_c'].bfill().ffill()
# resolution='monthly' # == c_res to indicate the resolution of the coarse iso
# resolution_dataset_root = '/Users/simon/Desktop/ORPB_resolution_datasets'
# data.to_csv(f"{resolution_dataset_root}/ORPB_isotope_data_bfill_{resolution}_{iso}.csv")

#%%
#----CHOOSE METHOD FOR TREND REMOVAL: Method I: isoMAP based trend, Method II: sinusoidal trend

#Use Method I: isoMAP based trend
isomap = pd.read_csv('isomap.csv')
isomap = isomap[['H2_PB', 'O18_PB']]
isomap.columns = ['2H', '18O']

st, et = resampled.index[0], resampled.index[-1]
delta_year = et.year - st.year + 1
st_month, et_month = st.month, et.month

# label isomap index by creating fake month
mm = pd.DataFrame(
    {"2H": np.tile(isomap['2H'].values, delta_year)[st_month-1:et_month-12],
     "18O": np.tile(isomap['18O'].values, delta_year)[st_month-1:et_month-12]
     }
)
mm.index = pd.to_datetime(pd.date_range(st, et + pd.DateOffset(days=29), freq='ME').strftime("%Y-%m-%d %H:00:00"))

# create compatible index for monthly mean for calculation
temp = pd.concat([resampled[iso], mm],axis=1)
temp.columns = [iso, '2H_mm', '18O_mm']
temp.bfill(inplace=True)
temp.dropna(inplace=True)
temp = temp.loc[resampled[iso].index]

trend = temp[['2H_mm', '18O_mm']]
trend.columns = ['2H', '18O']

# define iso_fit trend to subtract from iso
iso_fit = temp[f'{iso[-ison:]}_mm']

#%%
#Use Method II: sinusoidal trend
from scipy.optimize import curve_fit

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
from sklearn.preprocessing import StandardScaler
gp_alpha = 1e-6 #change to smallest value that doesn't give ABNORMAL and also LOO coverage and mean Z in next block

df_gp = pd.concat([data['rainfall (mm/hr)'],data['cumP'], resampled[f'{iso} deseasoned'], resampled['is_weekly']], axis=1)
df_gp['is_weekly'] = df_gp['is_weekly'].fillna(False).infer_objects(copy=False).astype(bool) #infer_objects forces the type without silently downcasting that pandas will depreciate
precip_min, precip_max = df_gp.loc[df_gp['rainfall (mm/hr)'] != 0, 'rainfall (mm/hr)'].min(), df_gp['cumP'].max() #st and et of cumP
df_gp[f'cum mass {iso} deseasoned'] = (df_gp[f'{iso} deseasoned']*df_gp['rainfall (mm/hr)']).cumsum()
continuous_precip = pd.DataFrame(np.arange(0.0, precip_max+precip_min, precip_min), columns=['precip_continuous']) # creates an integral cumP with equal intervals for better GP prediction
# df_gp = df_gp.merge(continuous_precip, how='right', left_on='cumP', right_on='precip_continuous')
# df_gp[['cumP', f'cum mass {iso} deseasoned']] = df_gp[['cumP', f'cum mass {iso} deseasoned']].interpolate('spline', order=2)
# X = df_gp['precip_continuous'][df_gp['weekly_obs']].values.reshape(-1,1) #convert to 2D array

#---------- if merge cuts out all obs_weekly rows:
# 1. FORCE continuous_precip to be a 1D numpy array
continuous_precip = np.asarray(continuous_precip).astype(float).ravel()

# 2. FORCE observation cumP to be a 1D numpy array
mask = np.asarray(df_gp['is_weekly']).astype(bool)
cumP_obs = np.asarray(df_gp.loc[mask, 'cumP']).astype(float).ravel()

# 3. Map observations to grid
idx = np.searchsorted(continuous_precip, cumP_obs)
idx = np.clip(idx, 0, continuous_precip.size - 1)

# 4. Build X and y (both numpy arrays)
X = continuous_precip.take(idx).reshape(-1, 1)
x_scaler = StandardScaler()
Xs = x_scaler.fit_transform(X) #scale the input
y = np.asarray(df_gp.loc[mask, f'cum mass {iso} deseasoned']).astype(float).ravel()

#%%
# --- OLD: per-point measurement alpha for use with RBF + learned WhiteKernel. Kept for reference.
# ---------- per-observation measurement-error variance (alpha array) ----------
# LGR 1-SD analytical precision for 18O (per mil).
# sigma_18O = 0.25  # ‰

# # Each weekly sample's measured value is forward-filled (line ~171) across a block
# # of hourly rows until the next sample, so the SAME measurement error applies to
# # every hour in that block. Block j therefore contributes (sigma * W_j) to the
# # cumulative mass, where W_j is the total rainfall its value is applied to.
# raw_iso = resampled[iso].reindex(df_gp.index)                 # ffilled, piecewise-constant
# block_id = (raw_iso != raw_iso.shift()).cumsum()              # one id per sample's block

# # rainfall carried by each block (ignore rows before the first sample, where iso is NaN)
# rain_valid = df_gp['rainfall (mm/hr)'].where(raw_iso.notna(), 0.0)
# block_W = rain_valid.groupby(block_id).sum()                  # W_j per block [mm]

# # independent per-block errors accumulate: Var(cum mass at k) = sigma^2 * sum_{j<=k} W_j^2
# cum_var_by_block = (sigma_18O ** 2) * (block_W ** 2).cumsum()  # [‰^2 mm^2]
# var_per_row = block_id.map(cum_var_by_block).to_numpy()
# alpha_meas = var_per_row[mask]                               # aligned to y

# # normalize_y=True rescales y by its std (ddof=0), so alpha must live in that space
# alpha_norm = alpha_meas / np.var(y)
# alpha_norm = np.clip(alpha_norm, 1e-10, None)               # keep strictly positive
# # -----------------------------------------------------------------------------


#%%
# --- Esther's reconstruction: plain Matern (nu=1.5), fixed tiny jitter, NO WhiteKernel. ---
# StandardScaler is kept only for numerical stability; Matern default length-scale bounds
# (1e-5, 1e5) are wide enough that scaling does not distort the fit vs her unscaled setup.
kernel = Matern()
# alpha regularizes the (near-singular) kernel matrix. Esther used 1e-7, but her 12H data was
# better conditioned; your closely-spaced hourly cumP points need more jitter or L-BFGS fails
# (ABNORMAL). 1e-3 ~ the noise the WhiteKernel previously found. Lower it toward 1e-5 if it
# still converges cleanly (smaller = closer to Esther); raise it if ABNORMAL persists.
#"jitter raised from Esther's 1e-7 to 1e-6 for numerical stability on the higher-resolution hourly data; selected by LOO calibration."
gp = GaussianProcessRegressor(kernel=kernel, alpha=gp_alpha, n_restarts_optimizer=50, normalize_y=True)

# --- OLD (calibrated) setup: RBF + learned WhiteKernel + per-point measurement alpha. Kept for reference.
# kernel = C(1.0, (1e-2, 1e2)) * RBF(length_scale=1.0, length_scale_bounds=(3e-1, 1e2)) \
#          + WhiteKernel(noise_level=1e-1, noise_level_bounds=(1e-3, 1e2))
# gp = GaussianProcessRegressor(kernel=kernel, alpha=alpha_norm, n_restarts_optimizer=50, normalize_y=True)

# Fit GP to data
gp.fit(Xs, y)
print('fitted kernel:', gp.kernel_)

#%% ---------------- validate uncertainty: leave-one-out coverage ----------------
# Refit with the ALREADY-optimized kernel (optimizer=None) dropping each point, then
# check what fraction of held-out observations land inside their 95% predictive band.
# ~95% coverage & mean|z|~0.8 => uncertainty is calibrated; far below => band too narrow.
n = len(y)
loo_mean = np.empty(n)
loo_std  = np.empty(n)
for i in range(n):
    keep = np.arange(n) != i
    g = GaussianProcessRegressor(kernel=gp.kernel_, optimizer=None,
                                 alpha=gp_alpha, normalize_y=True)  # match the fixed-jitter fit
    g.fit(Xs[keep], y[keep])
    m, s = g.predict(Xs[i].reshape(1, -1), return_std=True)
    loo_mean[i], loo_std[i] = m[0], s[0]

# predictive std for a held-out OBSERVATION: latent + WhiteKernel(repr) + its own measurement var
sigma_pred = np.sqrt(loo_std**2)# + alpha_meas)          # alpha_meas is in y-units (‰·mm)^2
z = (y - loo_mean) / sigma_pred
print(f'95% LOO coverage: {np.mean(np.abs(z) <= 1.96):.1%}  (target ~95%)') # >>95% bands are overcovering, <<95% (40-70) bands are too narrow
print(f'mean |z|:         {np.mean(np.abs(z)):.2f}  (target ~0.80 if calibrated)') #average standardized residual. <<0.8 band too wide, >>1 residuals are bigger than band claims 

#%%
# Predict on test points
# Xtest = df_gp['precip_continuous'].values.reshape(-1, 1)
# --- Direct prediction at each timestamp's cumP (no giant grid / interpolation needed) ---
# We only ever use predictions at data.index, so predict there directly. K_trans is
# (n, n_train) in memory (~cheap), avoiding the memory blow-up from predicting on the
# full, millions-long continuous_precip. precip is the EXACT cumP (no rounding/snapping).
X_pred = x_scaler.transform(data['cumP'].to_numpy().reshape(-1, 1))
y_mean, y_std = gp.predict(X_pred, return_std=True)

df_rev = pd.DataFrame(
    {
        "mean":  y_mean,
        "upper": y_mean + 1.96 * y_std,
        "lower": y_mean - 1.96 * y_std,
        "precip": data['cumP'].to_numpy(),
    },
    index=data.index,
)

# --- OLD: coarse-grid predict + cubic interpolation + searchsorted map-back. Needed ONLY
# --- for the MC latent-sampling band (return_cov cannot be formed at ~44k points). Restore
# --- this block (defines Xtests / grid / idx) if you switch back to the rigorous MC band.
# Xtest = continuous_precip.reshape(-1, 1)
# Xtests = x_scaler.transform(Xtest)
# from scipy.interpolate import interp1d
# grid = np.linspace(Xtests.min(), Xtests.max(), 2000).reshape(-1, 1)
# y_mean_grid, y_std_grid = gp.predict(grid, return_std=True)
# mean_fn = interp1d(grid.ravel(), y_mean_grid, kind='cubic')
# std_fn  = interp1d(grid.ravel(), y_std_grid,  kind='cubic')
# y_mean = mean_fn(Xtests.ravel()); y_std = std_fn(Xtests.ravel())
# cumP = data['cumP'].to_numpy()
# idx = np.searchsorted(continuous_precip, cumP)
# idx = np.clip(idx, 1, len(continuous_precip) - 1)
# idx -= (cumP - continuous_precip[idx - 1]) < (continuous_precip[idx] - cumP)
# df_rev = pd.DataFrame(
#     {"mean": y_mean[idx], "upper": y_mean[idx] + 1.96 * y_std[idx],
#      "lower": y_mean[idx] - 1.96 * y_std[idx], "precip": continuous_precip[idx]},
#     index=data.index,
# )

#%%
# ---------------plot GP predicted deseasoned mass fluxes----------------
plt.figure(figsize=[10,5])
plt.plot(df_rev.index, df_rev['mean'], '.-', label='GP predicted mass flux deseasoned')
plt.fill_between(df_rev.index, df_rev['lower'], df_rev['upper'], 
                 alpha=0.2, color="blue", label="Confidence interval")
plt.plot(df_gp.index, df_gp[f'cum mass {iso} deseasoned'], marker='.', zorder=10, label=f'Observed mass flux deseasoned {iso}', alpha=0.2, color='orange')
plt.xlabel('Date')
plt.ylabel(f'Mass flux deseasoned {iso} [‰ mm/hr]')
plt.title(f'{iso}: GP predicted mass flux deseasoned')
plt.legend()
plt.show()

# %%
# ---------------plot GP predicted isotope fluxes----------------

df_rev[f'rainfall (mm/{agg_res})'] = df_rev['precip'].diff()
df_rev.loc[df_rev.index[0], f'rainfall (mm/{agg_res})'] = 0.0
df_rev['iso_fit'] = iso_fit
# --- Esther concentration band: differenced the cumulative-mass ENVELOPE, which drops the
# --- covariance between adjacent points and collapses the band.
df_rev['mean_c'] = df_rev['mean'] + (df_rev['iso_fit']*df_rev[f'rainfall (mm/{agg_res})']).cumsum()
df_rev['mean_c'] = df_rev['mean_c'].diff()/(df_rev[f'rainfall (mm/{agg_res})']+1e-6)
df_rev.loc[df_rev[f'rainfall (mm/{agg_res})']==0, 'mean_c'] = np.nan
df_rev['upper_c'] = df_rev['upper'] + (df_rev['iso_fit']*df_rev[f'rainfall (mm/{agg_res})']).cumsum()
df_rev['upper_c'] = df_rev['upper_c'].diff()/(df_rev[f'rainfall (mm/{agg_res})']+1e-6)
df_rev.loc[df_rev[f'rainfall (mm/{agg_res})']==0, 'upper_c'] = np.nan
df_rev['lower_c'] = df_rev['lower'] + (df_rev['iso_fit']*df_rev[f'rainfall (mm/{agg_res})']).cumsum()
df_rev['lower_c'] = df_rev['lower_c'].diff()/(df_rev[f'rainfall (mm/{agg_res})']+1e-6)
df_rev.loc[df_rev[f'rainfall (mm/{agg_res})']==0, 'lower_c'] = np.nan

# --- OLD: if using RBF kernel ->concentration band by Monte-Carlo through the derivative ---
# Draw posterior samples of the cumulative-mass curve and differentiate EACH sample,
# then take percentiles. This propagates the full correlated posterior correctly.
# mean_c stays the deterministic GP mean; only the band comes from the samples.
# df_rev['mean_c'] = df_rev['mean'] + (df_rev['iso_fit']*df_rev[f'rainfall (mm/{agg_res})']).cumsum()
# df_rev['mean_c'] = df_rev['mean_c'].diff()/(df_rev[f'rainfall (mm/{agg_res})']+1e-6)
# df_rev.loc[df_rev[f'rainfall (mm/{agg_res})']==0, 'mean_c'] = np.nan

# N = 300
# # Sample the LATENT cumulative-mass function (exclude WhiteKernel noise): predict the full
# # posterior covariance, subtract the observation-noise variance off the diagonal, then draw.
# # Latent samples are smooth, so differentiating them does not amplify iid noise into an
# # exploding band (which is what gp.sample_y did -- it includes the WhiteKernel).
# mean_g, cov_g = gp.predict(grid, return_cov=True)                   # y-units, incl. WhiteKernel on diag
# noise_var_y = gp.kernel_.k2.noise_level * np.var(y)                 # WhiteKernel var, back to y-units
# cov_latent = cov_g - noise_var_y * np.eye(len(mean_g))             # remove observation noise
# cov_latent = (cov_latent + cov_latent.T) / 2 + 1e-8 * np.eye(len(mean_g))  # symmetrize + jitter for PSD
# rng = np.random.default_rng(0)
# M_grid_samples = rng.multivariate_normal(mean_g, cov_latent, size=N, method='eigh').T  # (len(grid), N)
# samp_fn = interp1d(grid.ravel(), M_grid_samples, axis=0, kind='cubic', fill_value='extrapolate')
# M_rev = samp_fn(Xtests.ravel())[idx, :]                             # (len(df_rev), N) at df_rev rows

# dP = df_rev[f'rainfall (mm/{agg_res})'].to_numpy()                  # ΔcumP per row
# iso_fit_arr = df_rev['iso_fit'].to_numpy()
# # concentration = d(deseasoned cum mass)/dP + seasonal fit. Differencing the seasonal
# # cumsum reduces to just adding iso_fit back, so no cumsum is needed -- which also avoids
# # np.cumsum propagating a single NaN across the whole column (the old pandas .cumsum skipped NaNs).
# C_samples = np.diff(M_rev, axis=0) / (dP[1:, None] + 1e-6) + iso_fit_arr[1:, None]
# C_samples = np.vstack([np.full((1, N), np.nan), C_samples])        # align, first row NaN

# print('M_rev all finite?', np.isfinite(M_rev).all(), '| iso_fit NaNs:', np.isnan(iso_fit_arr).sum())
# df_rev['lower_c'] = np.nanpercentile(C_samples, 2.5, axis=1)
# df_rev['upper_c'] = np.nanpercentile(C_samples, 97.5, axis=1)
# df_rev.loc[df_rev[f'rainfall (mm/{agg_res})']==0, ['lower_c', 'upper_c']] = np.nan

plt.figure(figsize=[15,5])
plt.scatter(df_rev.index[df_rev[f'rainfall (mm/{agg_res})']>0], df_rev.loc[df_rev[f'rainfall (mm/{agg_res})']>0, 'mean_c'], color='darkblue', s=10**2, label=f'GP predicted {iso} {agg_res}', marker='.', alpha=0.3, zorder=10)
plt.fill_between(df_rev.index[df_rev[f'rainfall (mm/{agg_res})']>0], df_rev.loc[df_rev[f'rainfall (mm/{agg_res})']>0, 'lower_c'], df_rev.loc[df_rev[f'rainfall (mm/{agg_res})']>0, 'upper_c'], 
                 alpha=0.1, color="purple", label="Confidence interval")
plt.scatter(df.index[(df['weekly_obs']==True) & (df['rainfall (mm/hr)']>0)], df.loc[(df['weekly_obs']==True) & (df['rainfall (mm/hr)']>0), iso], color='red', marker='.', s=10**2, label=f'Observed {iso} {c_res}') #rainfall mask technically isn't needed since weekly_obs only takes rain days
plt.xlabel('Date')
plt.ylabel(f'{iso} concentration [‰]')
plt.title(f'{iso}: GP predicted isotope concentration')
plt.ylim([isotopes[iso].min()-5, isotopes[iso].max()+5])
plt.legend()
plt.show()



#%% #fix this later***
#-------------------Create resampled data with GP predicted isotope values----------------

# S_n = data.resample(agg_res).sum()['snowmelt (mm/hr)']
# D_n = data.resample(agg_res).sum()['discharge (mm/hr)']
# ET_n = data.resample(agg_res).sum()['ET (mm/hr)']
# ORPB18O_n = data.resample(agg_res)['ORPB 18O'].mean()
# df_rev[[f'snowmelt (mm/{agg_res})', f'snowfall SWE (mm/{agg_res})', f'discharge (mm/{agg_res})', 'storage (mm)', f'ET (mm/{agg_res})', f'baseflow 1 (mm/{agg_res})', 'ORPB 18O', 'ORPB 2H', 'is_weekly']] = data.resample(agg_res).agg({'snowmelt (mm/hr)': 'sum', 'snowfall SWE (mm/hr)': 'sum', 'discharge (mm/hr)': 'sum',  'storage (mm)': 'mean', 'ET (mm/hr)': 'sum', 'baseflow 1 (mm/hr)': 'sum', 'ORPB 18O': 'mean', 'ORPB 2H': 'mean', 'is_weekly': 'any'})
# This is an issue, we want the same dataframe, data, just add the new columns from df_rev to it

data = data.join(df_rev[['mean_c', 'lower_c', 'upper_c']])

# %%
# save to csv files

resolution = 'hourly' #bimonthly, weekly, daily, monthly #EDIT HERE depending on agg_res

resolution_dataset_root = '/Users/simon/Desktop/ORPB_resolution_datasets'
import os
if not os.path.exists(resolution_dataset_root):
    os.makedirs(resolution_dataset_root)

# df_rev.to_csv(f"{resolution_dataset_root}/ORPB_isotope_data_isoMAP_{iso}_{resolution}.csv")
# df_rev.to_csv(f"{resolution_dataset_root}/ORPB_isotope_data_sinusoid_{iso}_{resolution}.csv")

data.to_csv(f"{resolution_dataset_root}/ORPB_isotope_data_bfill_{resolution}_{iso}.csv")
# %%



# additonal visualization
# note: rerun agg_res and c_res df creation block
resolution = 'weekly'
resolution_dataset_root = '/Users/simon/Desktop/ORPB_resolution_datasets'
df_rev = pd.read_csv(f"{resolution_dataset_root}/ORPB_isotope_data_isoMAP_{resolution}_{iso}.csv", index_col=0, parse_dates=[0])

plt.figure(figsize=[15,5])
plt.scatter(df_rev.index[df_rev[f'rainfall (mm/hr)']>0], df_rev.loc[df_rev[f'rainfall (mm/hr)']>0, 'mean_c'], color='darkblue', s=10**2, label=f'GP predicted {iso} {agg_res}', marker='.', alpha=0.3, zorder=10)
plt.fill_between(df_rev.index[df_rev[f'rainfall (mm/hr)']>0], df_rev.loc[df_rev[f'rainfall (mm/hr)']>0, 'lower_c'], df_rev.loc[df_rev[f'rainfall (mm/hr)']>0, 'upper_c'], 
                 alpha=0.1, color="purple", label="Confidence interval")
plt.scatter(df.index[(df['weekly_obs']==True) & (df['rainfall (mm/hr)']>0)], df.loc[(df['weekly_obs']==True) & (df['rainfall (mm/hr)']>0), iso], color='red', marker='.', s=10**2, label=f'Observed {iso} {c_res}') #rainfall mask technically isn't needed since weekly_obs only takes rain days
plt.xlabel('Date')
plt.ylabel(f'{iso} concentration [‰]')
plt.title(f'{iso}: GP predicted isotope concentration')
plt.xlim([pd.Timestamp('2018-01-01'), pd.Timestamp('2018-03-01')])
plt.legend()
plt.show()
#%%
#precip and discharge plot for reference
fig,ax = plt.subplots(nrows=1,ncols=1,figsize=[12,4])
ax.plot(data.index, data['discharge (mm/hr)'], color='blue', label='discharge (mm/hr)')
ax.set_xlabel('Date')
ax.set_ylabel('Discharge (mm/hr)', color='blue')
ax.tick_params(axis='y', labelcolor='blue')
ax1 = ax.twinx()
ax1.bar(data.index, data['rainfall (mm/hr)'], color='gray', alpha=0.6, width=0.01, label='rainfall (mm/hr)')
ax1.bar(data.index, data['snowfall SWE (mm/hr)'], color='black', alpha=0.6, width=0.01, label='snowfall SWE (mm)')
ax1.set_ylabel('(mm/hr)', color='gray')
ax1.tick_params(axis='y', labelcolor='gray')
ax1.set_ylim([0,0.75])#max(max(data_df['rainfall (mm/hr)']), max(data_df['snowfall SWE (mm/hr)']))*1.1])
ax1.invert_yaxis()
ax.set_title('ORPB Discharge and Rainfall')
ax.legend()
ax1.legend()
ax.set_xlim([pd.Timestamp('2014-08-01'), pd.Timestamp('2015-08-01')])
ax1.set_xlim([pd.Timestamp('2014-08-01'), pd.Timestamp('2015-08-01')])
plt.tight_layout()
# %%



# plot with Olin hall daily data
  

df_oh = pd.read_csv("olin_hall_precip.csv", index_col=1, parse_dates=True)
starttime = max(df_oh.index[0], df.index[0])
endtime = min(df_oh.index[-1], df.index[-1])

plt.figure(figsize=(16, 7))
plt.plot(
    df_oh[f"O18_daily"][df_oh["daily obs"] == True],
    ".",
    color="orange",
    ms=10,
    label="Olin Hall daily observation",
)
plt.plot(
    df_rev["mean_c"], ".", color="red", label=f"Predicted daily from {resolution} resolution"
)
plt.legend(ncol=5, frameon=False)
plt.title(r"$^{18}O$", fontsize=16)
plt.xlim([starttime, endtime])
# fig.supylabel("Delta isotope concentration (‰)", fontsize=16)
plt.tight_layout()
# %%
