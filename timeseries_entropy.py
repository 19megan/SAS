# This script will calculated the entropy of an isotope timeseries.
# Date: 02/06/2026

#%%
import numpy as np
import pandas as pd
from astropy.stats import knuth_bin_width # maximizes posterior probability of the histogram - closely related to Shannon entropy
import matplotlib.pyplot as plt

resolution_dataset_root = '/Users/simon/Desktop/ORPB_resolution_datasets'

iso = 'precip 18O'
resolution = ['daily', 'weekly', 'bimonthly', 'monthly']
df = ['df1', 'df2', 'df3', 'df4']
data = pd.read_csv('ORPB_isotope_data.csv', index_col=0, parse_dates=[0]) #len 51745

for i in range(len(resolution)):
    globals()[df[i]] = pd.read_csv(f"{resolution_dataset_root}/ORPB_isotope_data_{iso}_{resolution[i]}.csv", index_col=0, parse_dates=[0])


# clean up and find missing samples
data = data.drop(data[data['rainfall (mm/hr)']<=0].index) #for precip we don't need non rain days
data = data.drop(data[data['Sample Name'].isna()].index)
issample = ~data['Sample Name'].duplicated(keep='last')
data['is_weekly'] = issample

x = data.loc[data['is_weekly']==True, 'precip 18O'].dropna().values
# try using GP predicted values from resolution scaling script

# %%
bin_width = knuth_bin_width(x)
bins = np.arange(x.min(), x.max() + bin_width, bin_width)
plt.hist(x, bins=bins, density=True)
plt.xlabel(r'$\delta^{18}O$ (‰)')
plt.ylabel('Probability Density')
plt.title("Knuth's Rule Histogram (Entropy-based)")
plt.show()

#%%
counts, bin_edges = np.histogram(x, bins=bins)
probabilities = counts / counts.sum()
entropy = -np.sum(probabilities * np.log(probabilities + 1e-10)) # add small value to avoid log(0)
print(f'Entropy of the weekly precip 18O timeseries: {entropy:.4f}')

H_max = np.log(len(probabilities))
print(f'Maximum possible entropy (uniform distribution): {H_max:.4f}')
H = -(probabilities * np.log(probabilities + 1e-10))
print(f'Entropy contributions per bin: {H}')



#%%
# Issue with bin width for unique concentration timeseries is that it will be different when comparing other timeseries
# determine bin width once for pooled data:

# Pool all isotope data
x_all = np.concatenate([
    df1['mean_c'].dropna().values,
    df2['mean_c'].dropna().values,
    df3['mean_c'].dropna().values,
    df4['mean_c'].dropna().values
])

# Get global optimal bin width
from astropy.stats import knuth_bin_width
bw = knuth_bin_width(x_all)

bins = np.arange(x_all.min(), x_all.max() + bw, bw)

# then for each series:
for i in range(len(resolution)):
    x_i = globals()[df[i]]['mean_c'].dropna().values
    counts, _ = np.histogram(x_i, bins=bins)
    p = counts / counts.sum()
    entropy_i = -np.sum(p[p > 0] * np.log(p[p > 0]))
    print(f'Entropy of {resolution[i]} series: {entropy_i:.4f}')


# otherwise, could keep optimal bin width and use normralized entropy
# H_norm = H / np.log(N_bins) # rescales entropy to [0,1]

#could use both and show robustness

# %%

