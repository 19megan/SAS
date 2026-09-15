# This script will explore some of the different SAS shapes
# and functions used in the SAS studies in my meta-analysis

# Date: 06/08/2026

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gamma, beta
from mesas.sas.model import Model
import pandas as pd


dates = pd.date_range(start='2020-01-01', end='2020-12-31', freq='D')
df = pd.DataFrame(data ={'date': dates}, index=dates)
df_harman = pd.read_csv('LowerHafrenMESAS_data.csv')

df_ORPB = pd.read_csv('ORPB_isotope_data.csv', index_col=0, parse_dates=[0])
df_ORPB = df_ORPB.loc[pd.Timestamp('2014-01-01'): pd.Timestamp('2014-12-31')]
df_ORPB['ws'] = (df_ORPB['storage (mm)']-df_ORPB['storage (mm)'].min())/(df_ORPB['storage (mm)'].max()-df_ORPB['storage (mm)'].min()) #normalize storage to get ws between 0 and 1

df_borriero = pd.read_table('Selke_hydroclim_data.txt', index_col=0, parse_dates=[0]) #daily
iso_p = pd.read_table('Selke_d18O_P_raw.txt', index_col=0, parse_dates=[0]) #monthly
iso_q = pd.read_table('Selke_d18O_Q.txt', index_col=0, parse_dates=[0]) #monthly
df_borriero = df_borriero.join(iso_p[['d18O_P raw']])
df_borriero = df_borriero.join(iso_q[['d18O_Q']])
S0 = 1778
k1 = 0.675
k2 = 1.165
df_borriero['S'] = df_borriero['J [mm/d]']-df_borriero['ET [mm/d]']-df_borriero['Q [mm/d]']+S0
df_borriero['w'] = (df_borriero['S']-df_borriero['S'].min())/(df_borriero['S'].max()-df_borriero['S'].min())

# df_br = pd.read_excel('Weierbach_rainfall_Holtz_2009-2019.xlsx', index_col=0, parse_dates=[0], header=3) #start: 2009-01-01
# df_br['rainfall (mm)'] = df_br['rainfall (mm)'].apply(pd.to_numeric, errors = 'coerce')
# df_br = df_br.resample('4h').sum() # model run at 4h timestep
# stream = pd.read_excel('Weierbach_stream discharge_2009-2019.xlsx', index_col=0, parse_dates=[0], header=3) #start: 2009-09-01
# stream['Q (m3/s)'] = stream['Q (m3/s)'].apply(pd.to_numeric, errors = 'coerce')
# stream['Q (m3/4h)'] = stream['Q (m3/s)']*60*15 # convert to m3/15min
# stream = stream.resample('4h').sum()
# df_br = df_br.join(stream[['Q (m3/4h)']])
# iso_p = pd.read_excel('Weierbach_OH_rainfall_2009-2019.xlsx', index_col=3, parse_dates=[3], header=3) #fortnightly start: 2009-12-04
# df_br = df_br.join(iso_p[['d18O (permil)', 'd2H (permil)']])
# iso_q = pd.read_excel('Weierbach_OH_streamwater_2009-2019.xlsx', index_col=2, parse_dates=[2]) # start: 2009-09-21
# iso_q = iso_q[(iso_q['sample_type']=='streamwater') & (iso_q['sampling_location']=='SW1')]
# iso_q[['Q d18O (permil)', 'Q d2H (permil)']] = iso_q[['d18O (permil)', 'd2H (permil)']]
# df_br = df_br.join(iso_q[['Q d18O (permil)', 'Q d2H (permil)']])
# iso_q3H = pd.read_excel('Weierbach_tritium_2011-2017.xlsx', index_col=3, parse_dates=[3], header=3)
# iso_q3H.index = iso_q3H.index.normalize() # normalize to 00:00:00 # only 27 samples start: 2011-06-10
# iso_q3H['Q 3H (TU)'] = iso_q3H['3H (TU)']
# df_br = df_br.join(iso_q3H[['Q 3H (TU)']])
# # start where first isotope precip samples start (2009-12-04)
# df_br = df_br.loc[pd.Timestamp('2009-12-04'):]

df_sprenger_cal = pd.read_csv('CanVilla_Calibration2015_2017_model_Input.csv', sep=';', index_col=0, parse_dates=[0]) #hourly
xylem_ET = pd.read_csv('CanVilla_dO18_Xylem_ET.csv', sep=';', index_col=0, parse_dates=[0]) #monthly for only 2015
df_sprenger_cal = df_sprenger_cal.join(xylem_ET[['d18O_ET']])
df_sprenger_cal.loc[df_sprenger_cal['measC_Q [-]']==-999, 'measC_Q [-]'] = np.nan #fill in missing values for calibration
df_sprenger_val = pd.read_csv('CanVilla_Validation2011_2013_model_Input.csv', sep=';', index_col=0, parse_dates=[0]) #hourly
df_sprenger_val.loc[df_sprenger_val['measC_Q [-]']==-999, 'measC_Q [-]'] = np.nan #fill in missing values for validation
df_sprenger_cal['S'] = df_sprenger_cal['J [mm/h]']-df_sprenger_cal['ET [mm/h]']-df_sprenger_cal['Q [mm/h]']+389.91

inflow = pd.read_csv('DryCreek_isotopes_TTD_precip.csv', index_col=0, parse_dates=[0]) #start: 2015-12-01 00:00:00 hourly
stream = pd.read_csv('DryCreek_TTD_runoff.csv', index_col=0, parse_dates=[0]) #start: 2016-10-01 00:00:00 hourly
df_l = inflow.join(stream[['Q [mm/h]', 'measC_Q']])
df_l = df_l.resample('4h').sum() # model run at 4h timestep
df_l['Cin [-]'] = inflow['Cin [-]'].resample('4h').mean()
df_l['measC_Q'] = stream['measC_Q'].resample('4h').mean()
et = pd.read_csv('DryCreek_PET.csv', index_col=0, parse_dates=[0]) #start: 2016-10-01 00:00:00 hourly
et = et.resample('4h').sum()
df_l = df_l.join(et[['ET [mm/h]']])
df_l[['J [mm/4h]', 'Q [mm/4h]', 'ET [mm/4h]']] = df_l[['J [mm/h]', 'Q [mm/h]', 'ET [mm/h]']]
df_l['Cin [-]'] = df_l['Cin [-]'].bfill().ffill()
df_l['S'] = df_l['J [mm/4h]']-df_l['ET [mm/4h]']-df_l['Q [mm/4h]']+89.23
df_l['wi [-]'] = np.log(df_l['Q [mm/4h]']).values/np.log(df_l['Q [mm/4h]']).max()
logfactor_Q = 25.11
k1 = 0.935
k2 = 1.311
df_l['k'] = k1+(k2-k1)/np.log(logfactor_Q)*np.log(logfactor_Q-(logfactor_Q-1)*df_l['wi [-]'])
df_l_val = df_l.loc[pd.Timestamp('2019-12-02 08:00:00'):pd.Timestamp('2020-12-31 23:00:00')] #calibration period: 2019 water year
df_l = df_l.loc[pd.Timestamp('2018-10-01 00:00:00'):pd.Timestamp('2019-07-19 08:00:00')] #calibration period: 2019 water year
# %%

def make_LowerHafren_model_from(params):
    alpha, S_scale, et_scale = params
    scale = alpha
    shape = S_scale #use last value of S_scale as shape parameter
    # Generate x values based on 0.1st and 99.9th percentile of dist
    x=np.linspace(gamma.ppf(0.0001, shape, scale=scale),
                  gamma.ppf(0.999, shape, scale=scale), 100)
    y=gamma.pdf(x, shape, scale=scale)

    # plot the distribution
    plt.figure(figsize=(8, 4))
    plt.plot(x,y, 'r-', lw=2, label=f'TV gamma PDF (shape={shape}, scale={scale})')
    plt.title('Lower Hafren best Q SAS function')
    plt.xlabel('Age (hours)')
    plt.ylabel('Probability Density')
    plt.legend()
    plt.show()

def power_law(x, k):
    return k * x ** (k-1)
def TV_power_law(x, ws, k1, k2):
    k = k1 + (1-ws) * (k2-k1)
    return k * x ** (k-1)

def make_Selke_model_from(params):
    k1, k2, ket, S0 = params
    # Generate x values based on 0.1st and 99.9th percentile of dist
    x=np.linspace(0, 2000, 100) #df_borriero['S']
    y=TV_power_law(x, df_borriero['w'].iloc[-1], k1, k2)

    # plot the distribution
    plt.figure(figsize=(8, 4))
    plt.plot(x,y, 'r-', lw=2, label=f'TV power law PDF (k1={k1}, k2={k2})')
    plt.title('Selke best Q SAS function')
    plt.xlabel('Age (days)')
    plt.ylabel('Probability Density')
    plt.legend()
    plt.show()

def make_CanVilla_model_from(params):
    k1, k2, ket, S0 = params
    # Generate x values based on 0.1st and 99.9th percentile of dist
    x=df_sprenger_cal['S'] #np.linspace(0, 500, 100)  #df_sprenger_cal['S']
    y=TV_power_law(x, df_sprenger_cal['wi [-]'].iloc[-1], k1, k2)

    # plot the distribution
    plt.figure(figsize=(8, 4))
    plt.plot(x,y, 'r-', lw=2, label=f'TV power law PDF (k1={k1}, k2={k2})')
    plt.title('Can Villa best Q SAS function')
    plt.xlabel('Age (hours)')
    plt.ylabel('Probability Density')
    plt.legend()
    plt.show()

def make_DryCreek_model_from(params):
    k1, k2, ket, S0 = params
    # Generate x values based on 0.1st and 99.9th percentile of dist
    x=np.linspace(0, 500, 100)  #df_l['S']
    y=power_law(x, df_l['k'].iloc[-1])

    # plot the distribution
    plt.figure(figsize=(8, 4))
    plt.plot(x,y, 'r-', lw=2, label=f'TV power law PDF (k1={k1}, k2={k2})')
    plt.title('Dry Creek best Q SAS function')
    plt.xlabel('Age (4hours)')
    plt.ylabel('Probability Density')
    plt.legend()
    plt.show()

def make_BruntlandBurn_model_from(params):
    kQ1, kQ2, kET, S0, alpha = params
    # Generate x values based on 0.1st and 99.9th percentile of dist
    x=np.linspace(0, 3000, 100)
    y=TV_power_law(x, df_ORPB['ws'].iloc[-1], kQ1, kQ2)

    # plot the distribution
    plt.figure(figsize=(8, 4))
    plt.plot(x,y, 'r-', lw=2, label=f'TV power law PDF (kQ1={kQ1}, kQ2={kQ2})')
    plt.title('Bruntland-Burn best Q SAS function')
    plt.xlabel('(some) Age (hours)')
    plt.ylabel('Probability Density')
    plt.legend()
    plt.show()

import matplotlib.cm as cm

def make_ProvidenceCreek_model_from(params):
    aQmin, bQmin, S_A = params
    Q98 = df_ORPB['discharge (mm/hr)'].quantile(0.98)
    Q2 = df_ORPB['discharge (mm/hr)'].quantile(0.02)
    i = aQmin + (2-bQmin-aQmin)*np.sqrt(Q98)-np.sqrt(df_ORPB['discharge (mm/hr)'])/(np.sqrt(Q98)-np.sqrt(Q2))

    cmap = plt.get_cmap('viridis')
    colors = [cmap(i) for i in np.linspace(0,1,len(df_ORPB))]
    fig, ax = plt.subplots(figsize=[8,4])

    # plot over time
    for j in range(len(df_ORPB)):
        aQmin = min(i.iloc[j],1)
    # Generate x values based on 0.1st and 99.9th percentile of dist
        x=np.linspace(beta.ppf(0.0001, aQmin, bQmin),
                  beta.ppf(0.999, aQmin, bQmin), 100)
        y=beta.pdf(x, aQmin, bQmin)
        ax.plot(x,y, 'r-', lw=2, label=f'TV beta PDF (aQmin={aQmin}, bQmin={bQmin})')
    ax.set_title('Providence Creek best Q SAS function')
    ax.set_xlabel('Age (hours)')
    ax.set_ylabel('Probability Density')
    # ax.legend()
    plt.show()
    sm = cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, 
                                                         vmax=len(df_ORPB)-1))
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label='Time Index', fraction=0.046, pad=0.04)


# %%

make_LowerHafren_model_from([0.69, df_harman['S_scale'].iloc[-1], 398]) #params: alpha, S_scale, et_scale
#%%
make_Selke_model_from([0.675, 1.165, 0.3, 1778]) #params: k1, k2, ket, S0
# %%
make_CanVilla_model_from([0.28, 1.26, 3.26, 389.91]) #params: k1, k2, ket, S0
# %%
make_DryCreek_model_from([0.935, 1.311, 0.129, 89.23]) #params: k1, k2, ket, S0
# %%
make_BruntlandBurn_model_from([0.36, 0.8, 0.9, 2400, 0.991]) #params: kQ1, kQ2, kET, S0, alpha
# %%
make_ProvidenceCreek_model_from([0.7, 0.3, 3]) #params: aQmin, bQmin, S_A
# %%
