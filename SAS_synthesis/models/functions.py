# This has all functions needed for SAS_synthesis_script.py
# Date created: 09/14/2026

import pandas as pd
import numpy as np
from mesas.sas.model import Model
from SAS_synthesis.models.SAS_models import make_SAS_model
from tqdm import tqdm
from scores.continuous import nse
import xarray as xr

def load_data(location, data_file_path):
    """
    Load data from a CSV file.

    Parameters:
    location (str): The name of the catchment
    data_file_path (str): The path to the data files

    Returns:
    df (pd.DataFrame): The loaded data as a pandas DataFrame.
    spinup (pd.DataFrame): The spin-up data as a pandas DataFrame.
    df_validation (pd.DataFrame): The validation data as a pandas DataFrame.
    issample (list): list of booleans indicating if the sample is from the catchment or not.
    influx (list): name of influx colum.
    et (list): name of evapotranspiration column.
    discharge (list): name of discharge column.
    iso_out (list): name of outflow isotope column.
    iso_in (list): name of inflow isotope column.
    age_unit (str): The unit of age used in the data.
    """

    if location == 'Lower Hafren':
        df_harman = pd.read_csv(f'{data_file_path}/LowerHafrenMESAS_data.csv', index_col=1, parse_dates=[1])
        df_harman = df_harman.loc[pd.Timestamp('1999-01-01'): pd.Timestamp('2008-12-31')]
        issample = df_harman['Q Cl mg/l'].notna()
        influx = 'J'
        et = 'ET'
        discharge = 'Q'
        iso_out = 'Q Cl mg/l'
        iso_in = 'Cl mg/l'
        age_unit = 'day'
        return df_harman, df_harman, df_harman, issample, influx, et, discharge, iso_out, iso_in, age_unit

    elif location == 'Bruntland Burn':

        return np.nan


    elif location == 'Providence Creek':

        return np.nan
    

    elif location == 'Weierbach (2019)':
        df_r = pd.read_excel(f'{data_file_path}/Weierbach_rainfall_Holtz_2009-2019.xlsx', index_col=0, parse_dates=[0], header=3) #start: 2009-01-01
        df_r['rainfall (mm)'] = df_r['rainfall (mm)'].apply(pd.to_numeric, errors = 'coerce')
        df_r = df_r.resample('4h').sum() # model run at 4h timestep
        stream = pd.read_excel(f'{data_file_path}/Weierbach_stream discharge_2009-2019.xlsx', index_col=0, parse_dates=[0], header=3) #start: 2009-09-01
        stream['Q (m3/s)'] = stream['Q (m3/s)'].apply(pd.to_numeric, errors = 'coerce')
        stream['Q (m3/4h)'] = stream['Q (m3/s)']*60*15 # convert to m3/15min
        stream = stream.resample('4h').sum()
        df_r = df_r.join(stream[['Q (m3/4h)']])
        iso_p = pd.read_excel(f'{data_file_path}/Weierbach_OH_rainfall_2009-2019.xlsx', index_col=3, parse_dates=[3], header=3) #fortnightly start: 2009-12-04
        df_r = df_r.join(iso_p[['d18O (permil)', 'd2H (permil)']])
        iso_q = pd.read_excel(f'{data_file_path}/Weierbach_OH_streamwater_2009-2019.xlsx', index_col=2, parse_dates=[2]) # start: 2009-09-21
        iso_q = iso_q[(iso_q['sample_type']=='streamwater') & (iso_q['sampling_location']=='SW1')]
        iso_q[['Q d18O (permil)', 'Q d2H (permil)']] = iso_q[['d18O (permil)', 'd2H (permil)']]
        df_r = df_r.join(iso_q[['Q d18O (permil)', 'Q d2H (permil)']])

        # start where first isotope precip samples start (2009-12-04)
        df_r = df_r.loc[pd.Timestamp('2009-12-04'):]

        influx = 'rainfall (mm)'
        pet = 'PET'
        discharge = 'Q (m3/4h)'
        iso_out = 'Q d2H (permil)'
        iso_in = 'd2H (permil)'
        age_unit = '4h'

        df_r[iso_in] = df_r[iso_in].bfill().ffill()
        
        Sref = 2426 #mm
        Sroot=Sref-150 #total storage threshold
        n=20
        #*****************need to fix and get PET data
        df_r['ET'] = 0
        df_r['S'].iloc[0] = Sref
        for i in range(len(df_r)-1):
            df_r['ET'].iloc[i] = df_r[pet].iloc[i]*np.tanh((df_r['S'].iloc[i]/Sroot)**n)
            df_r['S'].iloc[i+1] = df_r[influx].iloc[i]-df_r['ET'].iloc[i]-df_r[discharge].iloc[i]+Sref
        df_r['ET'].iloc[len(df_r)] = df_r[pet].iloc[len(df_r)]*np.tanh((df_r['S'].iloc[len(df_r)]/Sroot)**n)
        # ET = df_r[pet]*np.tanh((df_r['S']/Sroot)**n)

        # model params
        lamda1s = 0.11
        f0 = 0.1
        Smin = df_r['S'].min()
        Sth = 105
        dSth = 3.97 #storage variation threshold for flashy events
        m = 1000 #fixed for threshold behavior
        df_r['f'] = f0*(1-np.tanh((df_r['S']/(Smin+Sth))**m))
        df_r['dS'] = df_r[influx] - df_r[discharge] - df_r[et]
        df_r['dSbar'] = 0
        for i in range(len(df_r)):
            df_r['dSbar'].iloc[i] = np.max(1/3*np.sum([df_r['dS'].iloc[i], df_r['dS'].iloc[i-1], df_r['dS'].iloc[i-2]]), 0)
        df_r['g'] = 1-np.exp(-df_r['dSbar']/dSth)
        df_r['lamda1'] = lamda1s * (df_r['f']+(1-df_r['f'])*df_r['g'])
        df_r['lamda2'] = 0.32 #constant
        df_r['lamda3'] = np.ones(len(df_r))-df_r['lamda2']-df_r['lamda1']

        # "The input data we used for the spin-up corresponds to the input data from October 2010 to October 2015 that we looped back over periods of 5 years."
        spinup = pd.concat([df_r.loc[pd.Timestamp('2010-10-01'): pd.Timestamp('2015-09-31 20:00:00')]]*2, ignore_index=True) #should be 100 yr spinup, but too much memory here
        newd = pd.date_range(start='2005-10-01', end='2015-09-31', freq='4h')
        spinup.index=newd

        df_r = df_r.loc[pd.Timestamp('2015-10-01'):pd.Timestamp('2017-09-31 20:00:00')]
        issample = df_r['Q d2H (permil)'].notna()

        return df_r, spinup, df_r, issample, influx, pet, discharge, iso_out, iso_in, age_unit
    
    elif location == 'Corin':

        return np.nan


    elif location == 'Weierbach (2021)':
        df_br = pd.read_excel(f'{data_file_path}/Weierbach_rainfall_Holtz_2009-2019.xlsx', index_col=0, parse_dates=[0], header=3) #start: 2009-01-01
        df_br['rainfall (mm)'] = df_br['rainfall (mm)'].apply(pd.to_numeric, errors = 'coerce')
        df_br = df_br.resample('4h').sum() # model run at 4h timestep
        stream = pd.read_excel(f'{data_file_path}/Weierbach_stream discharge_2009-2019.xlsx', index_col=0, parse_dates=[0], header=3) #start: 2009-09-01
        stream['Q (m3/s)'] = stream['Q (m3/s)'].apply(pd.to_numeric, errors = 'coerce')
        stream['Q (m3/4h)'] = stream['Q (m3/s)']*60*15 # convert to m3/15min
        stream = stream.resample('4h').sum()
        df_br = df_br.join(stream[['Q (m3/4h)']])
        iso_p = pd.read_excel(f'{data_file_path}/Weierbach_OH_rainfall_2009-2019.xlsx', index_col=3, parse_dates=[3], header=3) #fortnightly start: 2009-12-04
        df_br = df_br.join(iso_p[['d18O (permil)', 'd2H (permil)']])
        iso_q = pd.read_excel(f'{data_file_path}/Weierbach_OH_streamwater_2009-2019.xlsx', index_col=2, parse_dates=[2]) # start: 2009-09-21
        iso_q = iso_q[(iso_q['sample_type']=='streamwater') & (iso_q['sampling_location']=='SW1')]
        iso_q[['Q d18O (permil)', 'Q d2H (permil)']] = iso_q[['d18O (permil)', 'd2H (permil)']]
        df_br = df_br.join(iso_q[['Q d18O (permil)', 'Q d2H (permil)']])
        iso_q3H = pd.read_excel(f'{data_file_path}/Weierbach_tritium_2011-2017.xlsx', index_col=3, parse_dates=[3], header=3)
        iso_q3H.index = iso_q3H.index.normalize() # normalize to 00:00:00 # only 27 samples start: 2011-06-10
        iso_q3H['Q 3H (TU)'] = iso_q3H['3H (TU)']
        df_br = df_br.join(iso_q3H[['Q 3H (TU)']])

        # start where first isotope precip samples start (2009-12-04)
        df_br = df_br.loc[pd.Timestamp('2009-12-04'):]
        issample = df_br['Q d2H (permil)'].notna()
        spinup = df_br.loc[pd.Timestamp('2009-12-04'): pd.Timestamp('2015-09-31 20:00:00')]
        
        influx = 'rainfall (mm)'
        pet = 'PET'
        discharge = 'Q (m3/4h)'
        iso_out = 'Q d2H (permil)'
        iso_in = 'd2H (permil)'
        age_unit = '4h'

        df_br[iso_in] = df_br[iso_in].bfill().ffill()


        Sref = 2000 #mm
        Sroot=Sref-150 #total storage threshold
        n=20
        #*****************need to fix and get PET data
        df_br['S'] = df_br[influx]-df_br[et]-df_br[discharge]+Sref
        ET = df_br[pet]*np.tanh((df_br['S']/Sroot)**n)
        w = (df_br['S']-df_br['S'].min())/(df_br['S'].max()-df_br['S'].min())
        df_br['k'] = k1+(1-w)*(k2-k1)
        # add spin-up years 2008-2012
        spinup = pd.concat([df_benettin.loc['2013-02-04':'2015-05-12']]*2, ignore_index=True)
        newd = pd.date_range(start='2008-07-24', end='2013-02-03', freq='D')
        spinup.index=newd

        return df_br, spinup, df_br, issample, influx, pet, discharge, iso_out, iso_in, age_unit


    elif location == 'Chenqi':
        df_z = pd.read_excel(f'{data_file_path}/data of Chenqi in China - karst catchment-by Zhicai Zhang.xlsx', index_col=0, parse_dates=[0]) #start: 2016-08-01 - 2019-08-30
        df_z['outlet-D (‰)'] = df_z['outlet-D (‰)'].apply(pd.to_numeric, errors = 'coerce')
        rain_D = df_z[['Unnamed: 10', 'rain-D (‰)']].copy() #rain-D is on a separate date column 'Unnamed: 10'
        rain_D.index = rain_D['Unnamed: 10']
        rain_D = rain_D.dropna()
        df_z['rain-D (‰)'] = rain_D['rain-D (‰)']
        df_z['Q (m3/D)'] = df_z['outlet-Q (*10-3  m3/s)']*1e-3*60*60*24 #convert to m3/day
        # fill gaps in the precipitation isotopes: use the adjacent days (mean if both
        # the day before and the day after are available), otherwise the study period mean
        D_mean = -60.4 # value reported in the paper (arithmetic mean of the 317 obs here is -58.6)
        obs = df_z['rain-D (‰)'].copy() # neighbours are taken from observations only, never from filled values
        adjacent = pd.concat([obs.shift(1), obs.shift(-1)], axis=1).mean(axis=1) # daily gap-free index, so shift = day before/after
        rainy = df_z['P(mm)']>0
        df_z['rain-D (‰)'] = obs.fillna(adjacent.where(rainy)) # adjacent-day fill on days with precipitation
        df_z['rain-D (‰)'] = df_z['rain-D (‰)'].fillna(D_mean) # no adjacent data -> mean (dry days get it too, J=0 so it is unused)
        # model run at daily timesteps


        influx = 'P(mm)'
        pet = 'PET (mm)'
        et = 'ET (mm/day)'
        dischargeflux = 'Q (m3/D)'
        discharge = 'Q (mm/day)'
        iso_out = 'outlet-D (‰)'
        iso_in = 'rain-D (‰)'
        age_unit = 'D'

        S0 = 499
        k1 = 0.79
        k2 = 1.31
        ket = 0.42
        C_old = D_mean
        alpha = 0.92 #adjustment factor for actual drainage area
        area = 1.25*1e6 #km2->m2
        f = 0.005 #fractionation for ET

        df_z[discharge] = df_z[dischargeflux]/(alpha*area)*1000 #convert to mm/day
        # ET(t) = beta*PET(t), with beta from the water balance closed over the entire study period:
        # sum(P) = sum(ET) + sum(Q) + dS, and dS ~ 0 over the 3 years -> beta = (sum(P)-sum(Q))/sum(PET)
        beta = (df_z[influx].sum()-df_z[discharge].sum())/df_z[pet].sum() #0.580 (0.584 and 0.617 for the two full water years)
        df_z[et] = beta*df_z[pet]
        df_z['S'] = df_z[influx]-df_z[et]-df_z[discharge]+S0
        df_z['w'] = (df_z[discharge]-df_z[discharge].min())/(df_z[discharge].max()-df_z[discharge].min())
        df_z['k'] = k1 + (1-df_z['w'])*(k2-k1)

        # spin-up period: calibration period 08-01-2016 - 07-31-2018 looped n_loops times.
        # the IC memory takes ~2yrs to wash out, so one pass from empty storage is not enough to equilibrate
        n_loops = 5
        spinup = pd.concat([df_z.loc['2016-08-01':'2018-07-31']]*n_loops, ignore_index=True)
        newd = pd.date_range(end='2018-07-31', periods=len(spinup), freq='D')
        spinup.index=newd
        # validation period: 08-01-2018 - 07-31-2019
        df_z_val = df_z.loc[pd.Timestamp('2018-08-01'):pd.Timestamp('2019-08-01')]
        df_z = df_z.loc[pd.Timestamp('2016-08-01'):pd.Timestamp('2018-07-31')] #calibration period
        issample = (df_z[iso_out].notna()) & (df_z[discharge]>0)

        return df_z, spinup, df_z_val, issample, influx, et, discharge, iso_out, iso_in, age_unit


    elif location == 'Can Villa':
        df_sprenger_cal = pd.read_csv(f'{data_file_path}/CanVilla_Calibration2015_2017_model_Input.csv', sep=';', index_col=0, parse_dates=[0]) #hourly
        xylem_ET = pd.read_csv(f'{data_file_path}/CanVilla_dO18_Xylem_ET.csv', sep=';', index_col=0, parse_dates=[0]) #monthly for only 2015
        df_sprenger_cal = df_sprenger_cal.join(xylem_ET[['d18O_ET']])
        df_sprenger_cal.loc[df_sprenger_cal['measC_Q [-]']==-999, 'measC_Q [-]'] = np.nan #fill in missing values for calibration
        df_sprenger_val = pd.read_csv(f'{data_file_path}/CanVilla_Validation2011_2013_model_Input.csv', sep=';', index_col=0, parse_dates=[0]) #hourly
        df_sprenger_val.loc[df_sprenger_val['measC_Q [-]']==-999, 'measC_Q [-]'] = np.nan #fill in missing values for validation
        issample = (df_sprenger_cal['measC_Q [-]'].notna()) & (df_sprenger_cal['Q [mm/h]']>0)
        influx = 'J [mm/h]'
        et = 'ET [mm/h]'
        discharge = 'Q [mm/h]'
        iso_out = 'measC_Q [-]'
        iso_in = 'Cin [-]'
        location = 'Can Villa'
        age_unit = 'hour'

        S0 = 389.91 #mm
        k1 = 0.28
        k2 = 1.26
        df_sprenger_cal['S'] = df_sprenger_cal[influx]-df_sprenger_cal[et]-df_sprenger_cal[discharge]+S0
        df_sprenger_cal['k'] = k1+(1-df_sprenger_cal['wi [-]'])*(k2-k1)
        # spin-up 2 years using first year of calibration data (2015-2017)
        spinup = pd.concat([df_sprenger_cal.loc['2015-05-01 00:00:00':'2015-12-31 23:00:00']]*2, ignore_index=True)
        newd = pd.date_range(start='2013-12-27 00:00:00', end='2015-04-30 23:00:00', freq='h')
        spinup.index=newd
        return df_sprenger_cal, spinup, df_sprenger_val, issample, influx, et, discharge, iso_out, iso_in, age_unit


    elif location == 'Dry Creek':
        inflow = pd.read_csv(f'{data_file_path}/DryCreek_isotopes_TTD_precip.csv', index_col=0, parse_dates=[0]) #start: 2015-12-01 00:00:00 hourly
        stream = pd.read_csv(f'{data_file_path}/DryCreek_TTD_runoff.csv', index_col=0, parse_dates=[0]) #start: 2016-10-01 00:00:00 hourly
        df_l = inflow.join(stream[['Q [mm/h]', 'measC_Q']])
        df_l = df_l.resample('4h').sum() # model run at 4h timestep
        df_l['Cin [-]'] = inflow['Cin [-]'].resample('4h').mean()
        df_l['measC_Q'] = stream['measC_Q'].resample('4h').mean()
        et = pd.read_csv(f'{data_file_path}/DryCreek_PET.csv', index_col=0, parse_dates=[0]) #start: 2016-10-01 00:00:00 hourly
        et = et.resample('4h').sum()
        df_l = df_l.join(et[['ET [mm/h]']])
        df_l[['J [mm/4h]', 'Q [mm/4h]', 'ET [mm/4h]']] = df_l[['J [mm/h]', 'Q [mm/h]', 'ET [mm/h]']]
        df_l['Cin [-]'] = df_l['Cin [-]'].bfill().ffill()

        influx = 'J [mm/4h]'
        et = 'ET [mm/4h]'
        discharge = 'Q [mm/4h]'
        iso_out = 'measC_Q'
        iso_in = 'Cin [-]'
        age_unit = '4h'

        S0 = 89.23 #168.18
        k1 = 0.935 #0.45
        k2 = 1.311 #1.04
        logfactor_Q = 25.11 #44.6
        ket = 0.129 #0.94
        C_old = -146.78 #-93.28
        df_l['S'] = df_l[influx]-df_l[et]-df_l[discharge]+S0
        df_l['wi [-]'] = np.log(df_l['Q [mm/4h]']).values/np.log(df_l['Q [mm/4h]']).max()
        # df_l['k'] = k1+(k2-k1)*np.log((1-logfactor_Q)*df_l['wi [-]']) #in the paper
        df_l['k'] = k1+(k2-k1)/np.log(logfactor_Q)*np.log(logfactor_Q-(logfactor_Q-1)*df_l['wi [-]']) #in their code

        # add spin-up years 2008-2012
        spinup = pd.concat([df_l.loc['2017-10-01 00:00:00':'2018-10-01 00:00:00']]*2, ignore_index=True) #repeated 10x for IC on age-rank storage
        newd = pd.date_range(start='2016-09-30 20:00:00', end='2018-10-01 00:00:00', freq='4h')
        spinup.index=newd

        df_l_val = df_l.loc[pd.Timestamp('2019-12-02 08:00:00'):pd.Timestamp('2020-12-31 23:00:00')] #calibration period: 2019 water year
        df_l = df_l.loc[pd.Timestamp('2018-10-01 00:00:00'):pd.Timestamp('2019-07-19 08:00:00')] #calibration period: 2019 water year
        issample = (df_l['measC_Q'].notna()) & (df_l['Q [mm/4h]']>0)

        return df_l, spinup, df_l_val, issample, influx, et, discharge, iso_out, iso_in, age_unit


    elif location == 'Selke':
        df_borriero = pd.read_table(f'{data_file_path}/Selke_hydroclim_data.txt', index_col=0, parse_dates=[0]) #daily
        iso_p = pd.read_table(f'{data_file_path}/Selke_d18O_P_raw.txt', index_col=0, parse_dates=[0]) #monthly
        iso_q = pd.read_table(f'{data_file_path}/Selke_d18O_Q.txt', index_col=0, parse_dates=[0]) #monthly
        df_borriero = df_borriero.join(iso_p[['d18O_P raw']])
        df_borriero = df_borriero.join(iso_q[['d18O_Q']])
        issample = df_borriero['d18O_Q'].notna()
        influx = 'J [mm/d]'
        et = 'ET [mm/d]'
        discharge = 'Q [mm/d]'
        iso_out = 'd18O_Q'
        iso_in = 'd18O_P raw'
        age_unit = 'day'

        #monthly interpolation to daily with step
        df_borriero[iso_in] = df_borriero[iso_in].bfill().ffill()
        S0 = 1778
        k1 = 0.675
        k2 = 1.165
        df_borriero['S'] = df_borriero[influx]-df_borriero[et]-df_borriero[discharge]+S0
        w = (df_borriero['S']-df_borriero['S'].min())/(df_borriero['S'].max()-df_borriero['S'].min())
        df_borriero['k'] = k1+(1-w)*(k2-k1)
        # add spin-up years 2008-2012
        spinup = pd.concat([df_borriero.loc['2013-02-04':'2015-05-12']]*2, ignore_index=True)
        newd = pd.date_range(start='2008-07-24', end='2013-02-03', freq='D')
        spinup.index=newd

        return df_borriero, spinup, df_borriero, issample, influx, et, discharge, iso_out, iso_in, age_unit


    elif location == 'Hemuqiao':
        return np.nan
    

    else:
        raise ValueError(f"Unknown location: {location}. Please provide a valid location.")





def make_prior_bounds(location):
    """
    Create prior bounds for model parameters based on the specified location.

    Parameters:
    location (str): The name of the catchment.

    Returns:
    bounds (dict): A dictionary containing parameter names as keys and their corresponding (low, high) bounds as values.
    """
    if location == 'Chenqi':
        bounds = {
            "k1": (0.05, 10), #(low, high) for 95% CI
            "k2": (0.05, 10),
            "ket": (0.05, 10),
            "S0": (50, 6000),
            "f": (0, 0.1),
            "alpha": (0.15, 1.5),
            "C_old": (-146.71, 29.02),
        }
        return bounds
    elif location == 'Dry Creek':
        bounds = { # for lapides/DryCreek model
            "S0": (60, 270), #(low, high) for 95% CI
            "k1": (0., 1),
            "k2": (0.75, 1.25),
            "logfactor_Q": (0, 100),
            "ket": (0., 2),
            "C_old": (-200, 0),
        }
        return bounds
    elif location == 'Selke':
        bounds = { # for benettin/Selke model
            "S0": (681, 2875), #(low, high) for 95% CI
            "k1": (0.4, 0.95),
            "k2": (0.53, 1.8),
            "ket": (0.2, 1.95),
        }
        return bounds
    else:
        raise ValueError(f"Unknown location: {location}. Please provide a valid location.")



def make_model(location, data_file_path, spin_up=False, MC=False, validate=False):
    """
    Create a model instance based on the specified location.

    Parameters:
    location (str): The name of the catchment.
    data_file_path (str): The path to the data files.
    spinup (bool): Whether to use spin-up data or not.
    MC (bool): Whether to use Monte Carlo simulation or not.

    Returns:
    model: An instance of the Model class for the specified location.
    """
    spin_done = False
    df, spinup, df_val, issample, influx, et, discharge, iso_out, iso_in, age_unit = load_data(location, data_file_path)
    if spin_up == True:
        if MC == True:
            # Monte Carlo simulation logic here

            bounds = make_prior_bounds(location)
            param_names = list(bounds.keys())

            def sample_params(bounds):
                return [np.random.uniform(low, high) for low, high in bounds.values()]

            def evaluate_model(params):
                try:
                    # Run spin-up model to get sT_init and mT_init
                    spin_done = False
                    model = make_SAS_model(location, df, spinup, spin_done, influx, params=params)
                    model.run()
                    if len(model.get_mT(iso_in)[:,-1]) < len(df): #CanVilla
                        df['mT_init'] = pd.concat([pd.Series(model.get_mT(iso_in)[:,-1]), pd.Series(np.zeros(len(df)-len(model.get_mT(iso_in)[:,-1])))], ignore_index=True).values
                        sT_init = pd.concat([pd.Series(model.get_sT()[:,-1]), pd.Series(np.zeros(len(df)-len(model.get_sT()[:,-1])))], ignore_index=True).values
                    else:
                        df['mT_init'] = model.get_mT(iso_in)[:len(df),-1]
                        sT_init = model.get_sT()[:len(df),-1]
                    spin_done = True
                    model=make_SAS_model(location, df, spinup, spin_done, influx, sT_init, params=params)
                    model.run()

                    obs = model.data_df[iso_out][issample].to_numpy()
                    pred = model.data_df[f'{iso_in} --> {discharge}'][issample].to_numpy()
                    del model #free up memory
                    RMSE = np.sqrt(np.mean((pred-obs)**2))
                    obs_xr = xr.DataArray(obs)
                    pred_xr = xr.DataArray(pred)
                    NSE = nse(pred_xr, obs_xr)
                    return RMSE, NSE.item()
                except Exception as e: # catches unstable parameter combinations
                    return np.nan

            N = 100 #samples
            results = []
            for _ in tqdm(range(N)):
                params = sample_params(bounds)
                result = evaluate_model(params)
                if result is not np.nan and not np.isnan(result[0]):
                    RMSE, NSE = result
                    results.append((params + [RMSE, NSE]))
            
            cols = param_names + ['RMSE','NSE']
            results_df = pd.DataFrame(results, columns=cols)
            p = results_df.loc[results_df['NSE']==results_df['NSE'].max()] #[251.30273735,   0.71963209,   1.20491725,  27.47813886, 1.2644441 , -46.05504647]
            p = p.to_numpy()
            p = np.delete(p, -1) # take off RMSE column for params
            p = np.delete(p, -1) # take of NSE column
            spin_done = False
            if validate == True:
                model = make_SAS_model(location, df_val, spinup, spin_done, influx, params=p)
                model.run()
                if len(model.get_mT(iso_in)[:,-1]) < len(df_val): #df_val
                    df_val['mT_init'] = pd.concat([pd.Series(model.get_mT(iso_in)[:,-1]), pd.Series(np.zeros(len(df_val)-len(model.get_mT(iso_in)[:,-1])))], ignore_index=True).values
                    sT_init = pd.concat([pd.Series(model.get_sT()[:,-1]), pd.Series(np.zeros(len(df_val)-len(model.get_sT()[:,-1])))], ignore_index=True).values
                else:
                    df_val['mT_init'] = model.get_mT(iso_in)[:len(df_val),-1]
                    sT_init = model.get_sT()[:len(df_val),-1]
                spin_done = True
                model=make_SAS_model(location, df_val, spinup, spin_done, influx, sT_init, params=p)

            else:
                model = make_SAS_model(location, df, spinup, spin_done, influx, params=p)
                model.run()
                if len(model.get_mT(iso_in)[:,-1]) < len(df): #df_val
                    df['mT_init'] = pd.concat([pd.Series(model.get_mT(iso_in)[:,-1]), pd.Series(np.zeros(len(df)-len(model.get_mT(iso_in)[:,-1])))], ignore_index=True).values
                    sT_init = pd.concat([pd.Series(model.get_sT()[:,-1]), pd.Series(np.zeros(len(df)-len(model.get_sT()[:,-1])))], ignore_index=True).values
                else:
                    df['mT_init'] = model.get_mT(iso_in)[:len(df),-1]
                    sT_init = model.get_sT()[:len(df),-1]
                spin_done = True
                model=make_SAS_model(location, df, spinup, spin_done, influx, sT_init, params=p)

        else: #No MC
            # Run spin-up model to get sT_init and mT_init
            model = Model(data_df=spinup, config=f'{data_file_path}/../models/{location}_config_spinup.json', influx=influx)
            model.run()
            if len(model.get_mT(iso_in)[:,-1]) < len(df): #CanVilla
                df['mT_init'] = pd.concat([pd.Series(model.get_mT(iso_in)[:,-1]), pd.Series(np.zeros(len(df)-len(model.get_mT(iso_in)[:,-1])))], ignore_index=True).values
                sT_init = pd.concat([pd.Series(model.get_sT()[:,-1]), pd.Series(np.zeros(len(df)-len(model.get_sT()[:,-1])))], ignore_index=True).values
            else:
                df['mT_init'] = model.get_mT(iso_in)[:len(df),-1]
                sT_init = model.get_sT()[:len(df),-1]
            spin_done = True
            model = make_SAS_model(location, df, spinup, spin_done, influx, sT_init)

    else: #No spin-up
        model = Model(data_df=df, config=f'{data_file_path}/../models/{location}_config.json', influx=influx)

    return model
        