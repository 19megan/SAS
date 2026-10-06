# This has all functions needed for SAS_synthesis_script.py
# Date created: 09/14/2026

import pandas as pd
import numpy as np
from mesas.sas.model import Model
from SAS_synthesis.models.SAS_models import make_SAS_model, weierbach_hydrology, WEIERBACH_PARAMS
from tqdm import tqdm
from scores.continuous import nse
import xarray as xr


def _bb_storage(P, PET, Q, tol=1e-8, max_iter=50):
    """
    Bruntland Burn relative storage, following Benettin et al. (2017) section 3.2:
    dS(t)/dt = P(t) - lf(t)*PET(t) - Q(t), with lf(t) = min{1, S(t)/(2*std(S))}.

    lf depends on S and on std(S) over the whole record, so the balance is stepped forward with
    std(S) held fixed and std(S) is then updated, until it converges. The returned S is offset so
    its minimum is zero (lf must stay in [0,1], which fixes the otherwise free datum).

    Parameters:
    P, PET, Q (np.ndarray): influx, potential ET and discharge, same units per timestep.

    Returns:
    S (np.ndarray): relative storage, minimum 0.
    ET (np.ndarray): actual evapotranspiration, lf*PET.
    """
    n = len(P)
    sd = np.std(np.cumsum(P - PET - Q)) # first guess, lf = 1 throughout
    for _ in range(max_iter):
        S = np.empty(n)
        ET = np.empty(n)
        S[0] = 0.0
        for i in range(n):
            lf = min(1.0, max(0.0, S[i]/(2*sd)))
            ET[i] = lf*PET[i]
            if i+1 < n:
                S[i+1] = S[i] + P[i] - ET[i] - Q[i]
        S = S - S.min()
        sd_new = np.std(S)
        if abs(sd_new-sd) < tol:
            break
        sd = sd_new
    return S, ET


def _weierbach_data(data_file_path):
    """
    Weierbach input data at the 4 h model timestep, following Rodriguez and Klaus (2019) section 2.5
    and Rodriguez et al. (2021) section 2.2.

    Returns:
    df (pd.DataFrame): October 2010 - October 2017. 'main' marks the period of interest
        (October 2015 - October 2017), 'record' the daily steps of it that the model records.
    spinup (pd.DataFrame): the October 2010 - October 2015 cycle looped back twice before 2010.
    """
    start, main_start, end = pd.Timestamp('2010-10-01'), pd.Timestamp('2015-10-01'), pd.Timestamp('2017-09-30 20:00')
    area = 0.42e6 # m2, the 42 ha catchment above the SW1 weir
    grid = pd.date_range(start, end, freq='4h')

    def read_15min(file, col):
        s = pd.read_excel(f'{data_file_path}/{file}', index_col=0, parse_dates=[0], header=3)[col]
        return pd.to_numeric(s, errors='coerce') # 'no data' -> NaN

    # precipitation (mm): 15 min tipping bucket, 10 min from 2016. Depths, so they sum to 4 h
    # directly. The few 'no data' stamps count as 0.
    P = read_15min('Weierbach_rainfall_Holtz_2009-2019.xlsx', 'rainfall (mm)')
    J = P[start:end+pd.Timedelta('4h')].resample('4h').sum().reindex(grid)

    # discharge at SW1 (m3/s, 15 min) -> mm per 15 min over the catchment area -> mm/4h.
    # Two gaps: 38 days in 2012 (spin-up cycle) and 2017-01-01 - 2017-01-30 (main period; none of
    # the four gauges cover it). Interpolated linearly in time so they are not summed as zero flow.
    Q = read_15min('Weierbach_stream discharge_2009-2019.xlsx', 'Q (m3/s)')
    Q = Q[start-pd.Timedelta('30D'):end+pd.Timedelta('30D')]
    gap = Q.isna().resample('4h').max().reindex(grid)
    Q = Q.interpolate('time', limit_area='inside')
    Q = (Q*900/area*1000)[start:end+pd.Timedelta('4h')].resample('4h').sum().reindex(grid)

    # potential ET: FAO Penman-Monteith reference ET0 at 15 min in mm/h (Glaser et al., 2016),
    # stamped with MATLAB datenums (days since year 0; 719529 = 1970-01-01). The datenums are
    # stored to 4 decimals (~9 s), so round to the 15 min grid. mm/h * 0.25 h = mm per 15 min.
    et0 = pd.read_table(f'{data_file_path}/ET0_15min_mmh_oct10-jan18_Weierbach.txt', header=None, names=['datenum', 'ET0 [mm/h]'])
    et0.index = pd.to_datetime(et0['datenum']-719529, unit='D').dt.round('15min')
    PET = (et0['ET0 [mm/h]']*0.25)[start:end+pd.Timedelta('4h')].resample('4h').sum().reindex(grid)

    df = pd.DataFrame({'J [mm/4h]': J, 'Q [mm/4h]': Q, 'PET [mm/4h]': PET, 'Q gapfilled': gap}, index=grid)

    # d2H in precipitation: bulk samples dated at collection. "The time series of tracer in
    # precipitation was interpolated between two consecutive samples (e.g., A and B) as being equal
    # to the value of the next sample (i.e., B)" -> backfill. NB the papers also inserted the
    # sequential rainfall samples (~23 h) for 2015-2017; those are not in this file.
    iso_p = pd.read_excel(f'{data_file_path}/Weierbach_OH_rainfall_2009-2019.xlsx', header=3)
    iso_p = pd.Series(pd.to_numeric(iso_p['d2H (permil)'], errors='coerce').values,
                      index=pd.to_datetime(iso_p['sampling_date'])).dropna()
    df['Cin d2H'] = iso_p.reindex(grid.union(iso_p.index)).bfill().reindex(grid)

    # d2H in the stream: grab samples at SW1 (dates only, placed at 00:00 of the sampling day)
    iso_q = pd.read_excel(f'{data_file_path}/Weierbach_OH_streamwater_2009-2019.xlsx')
    iso_q = iso_q[(iso_q['sample_type']=='streamwater') & (iso_q['sampling_location']=='SW1')]
    iso_q = pd.Series(pd.to_numeric(iso_q['d2H (permil)'], errors='coerce').values,
                      index=pd.to_datetime(iso_q.iloc[:, 2])).dropna()
    df['measC_Q d2H'] = iso_q.groupby(level=0).mean().reindex(grid)

    # tritium in the stream (Rodriguez et al., 2021 used 24 of these), placed in their 4 h step
    iso_3H = pd.read_excel(f'{data_file_path}/Weierbach_tritium_2011-2017.xlsx', header=3)
    iso_3H = pd.Series(iso_3H['3H (TU)'].values, index=pd.to_datetime(iso_3H['sampling_date']).dt.floor('4h'))
    df['measC_Q 3H'] = iso_3H.groupby(level=0).mean().reindex(grid)

    df['main'] = df.index >= main_start
    # record the state once a day (last step of the day) over the main period only; the recorded
    # arrays are (max_age x recorded steps), so this keeps them ~100 MB instead of ~2 GB
    df['record'] = df['main'] & (df.index.hour == 20)

    # "The input data we used for the spin-up corresponds to the input data from October 2010 to
    # October 2015 that we looped back over periods of 5 years." The water balance (S, ET) is run
    # over the papers' full 100 years in weierbach_hydrology. The transport spin-up only has to be
    # longer than the main run's age window (7 years): water younger than T depends only on the
    # last T of forcing when the SAS functions are functions of ST.
    n_loops = 2
    spinup = pd.concat([df.loc[~df['main']]]*n_loops, ignore_index=True)
    spinup.index = pd.date_range(end=start-pd.Timedelta('4h'), periods=len(spinup), freq='4h')
    spinup['main'] = False
    spinup['record'] = False

    return df, spinup


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
        # Benettin et al. (2017), WRR 53, doi:10.1002/2016WR020117. Model run at hourly timesteps
        # (the paper's marginal TTD averages "more than 26,000" curves = the 26,280 hourly steps
        # of 2011-06-01 to 2014-05-31).
        # NB the compiled csv columns are named *_mmd but hold PER-HOUR depths: daily sums of
        # etpot_mmd reproduce 'ETP (mm d-1)' in the daily xlsx to 3 decimals. Do not divide by 24.
        df_benettin = pd.read_csv(f'{data_file_path}/BruntlandBurn_data_compiled_2011-2014.csv',
                                  index_col=0, parse_dates=[0], date_format='%d-%b-%Y %H:%M:%S') #hourly 2011-06-02 - 2014-09-30
        df_benettin.columns = ['J [mm/h]', 'Q [mm/h]', 'PET [mm/h]'] # P, Q, PET from Soulsby et al. (2015)

        influx = 'J [mm/h]'
        pet = 'PET [mm/h]'
        et = 'ET [mm/h]'
        discharge = 'Q [mm/h]'
        iso_out = 'measC_Q d2H'
        iso_in = 'Cin d2H'
        age_unit = 'hour'

        # patch the 204 missing discharge hours (4 gaps, the longest 196 h over 2012-05-30 - 2012-06-07)
        # with the daily record, which is gap-free over the same period
        daily = pd.read_excel(f'{data_file_path}/BruntlandBurn_ORIGINAL_BBdaily-all_2011_2014.xlsx', index_col=0, parse_dates=[0])
        Q_from_daily = (daily['BB_Q (mm/d)']/24).reindex(df_benettin.index, method='ffill')
        gap = df_benettin[discharge].isna()
        df_benettin.loc[gap, discharge] = Q_from_daily[gap]

        # measured deuterium. Use the EcH2Oiso workbook, not the P_d2Hmodel/Q_d2H columns of the daily
        # xlsx: this sheet reproduces the flow-weighted mean precip d2H (-59.55) and mean stream d2H
        # (-57.55) the paper reports (-59.8 and -57.4), the daily xlsx columns are a modelled variant.
        iso = pd.read_excel(f'{data_file_path}/BruntlandBurn_Isotopes_usedfor_EcH2Oiso.xlsx', sheet_name='Precip-Stream', header=[0,1])
        iso.columns = ['date', 'P d2H', 'P d18O', 'Q d2H', 'Q d18O'] # two-row header: (Precip|Stream) x (dD|d18O)
        iso = iso.set_index(pd.to_datetime(iso.pop('date'))) #daily 2011-06-01 - 2016-09-19

        day = df_benettin.index.normalize()
        # precipitation samples are cumulative over 24 h -> the daily value applies to all 24 hours
        # (every rainy hour in the record gets a value, so no filling is needed)
        df_benettin[iso_in] = iso['P d2H'].reindex(day).to_numpy()
        # stream samples are instantaneous grabs at 9 A.M. -> only the 09:00 hour is an observation.
        # Broadcasting these over the day instead would inflate the sample count 24x and make the
        # NSE incomparable with the other catchments.
        df_benettin[iso_out] = np.where(df_benettin.index.hour==9, iso['Q d2H'].reindex(day).to_numpy(), np.nan)

        # relative storage from the water balance dS/dt = P - lf*PET - Q, with the PET limiting
        # factor lf = min(1, S/(2*std(S))). S appears on both sides, so step the balance forward and
        # iterate on std(S) until it stops moving; S is then offset so its minimum is zero.
        df_benettin['S_rel'], df_benettin[et] = _bb_storage(df_benettin[influx].to_numpy(),
                                                            df_benettin[pet].to_numpy(),
                                                            df_benettin[discharge].to_numpy())
        # check from the paper: "storage variations computed with this simplified methodology show no
        # trend in the observed 3 year period" -> 13.7 mm/yr drift across a 246 mm range. ET/PET = 0.68

        # spin-up: 8 years repeating the hydrologic data of the first year of measurements
        n_loops = 8
        spinup = pd.concat([df_benettin.loc['2011-06-01':'2012-05-31']]*n_loops, ignore_index=True)
        spinup.index = pd.date_range(end='2011-06-01 23:00:00', periods=len(spinup), freq='h')

        # calibration period: 1 June 2011 to 1 June 2014 (1095 days). The record runs 4 months longer
        # (to 2014-09-30) and the isotopes to 2016, so there is spare data if you want to extend it.
        df_benettin_val = df_benettin.loc[pd.Timestamp('2014-06-01'):]
        df_benettin = df_benettin.loc[pd.Timestamp('2011-06-01'):pd.Timestamp('2014-05-31 23:00:00')]
        issample = (df_benettin[iso_out].notna()) & (df_benettin[discharge]>0) #992 stream samples

        return df_benettin, spinup, df_benettin_val, issample, influx, et, discharge, iso_out, iso_in, age_unit



    elif location == 'Providence Creek':

        return np.nan
    

    elif location in ('Weierbach (2019)', 'Weierbach (2021)'):
        # Rodriguez and Klaus (2019), WRR 55, doi:10.1029/2019WR024973, and Rodriguez et al. (2021),
        # HESS 25, doi:10.5194/hess-25-401-2021. Same data, timestep (4 h), periods and ET/storage
        # model; they differ in Sref (calibrated in 2019, fixed at 2000 mm in 2021) and in the
        # calibration (2021 adds tritium), which is handled in make_SAS_model.
        df_w, spinup = _weierbach_data(data_file_path)
        influx = 'J [mm/4h]'
        et = 'ET [mm/4h]'
        discharge = 'Q [mm/4h]'
        iso_out = 'measC_Q d2H'
        iso_in = 'Cin d2H'
        age_unit = '4h'

        # S, ET and the composite SAS weights for the default (paper) parameters. make_SAS_model
        # recomputes them when the parameters are sampled (ET depends on Sref through Sroot).
        p = WEIERBACH_PARAMS[location]
        weierbach_hydrology(spinup, df_w, p['Sref'], p['Sth'], p['dSth'], p['f0'], p['lamda1s'], p['lamda2'])

        # performance is evaluated over October 2015 - October 2017 only. df_w also carries the
        # October 2010 - October 2015 cycle in front of it so the main run can track water up to
        # 7 years old (mesas caps max_age at the run length; older water gets C_old).
        issample = df_w['main'] & df_w[iso_out].notna() & (df_w[discharge]>0)

        return df_w, spinup, df_w, issample, influx, et, discharge, iso_out, iso_in, age_unit

    elif location == 'Corin':

        return np.nan


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
        # typo in the source file: one sample is dated 2015-12-04, which sits in the December slot of
        # an otherwise strictly monthly sequence (2014-11-06 -> 2015-01-14) but falls 7 months past
        # the end of the hydroclimatic record, so the join below would silently drop it and calibrate
        # on 26 of the 27 samples. Corrected here rather than in the data file.
        iso_q = iso_q.rename(index={pd.Timestamp('2015-12-04'): pd.Timestamp('2014-12-04')}).sort_index()
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
        df_borriero[iso_in] = df_borriero[iso_in].ffill().bfill()
        # storage variations, paper equation 1: S(t) = S0 + V(t). V is a running sum of the fluxes,
        # not the per-timestep difference (which spans only -5.8 to 36.8 mm against -149.2 to
        # 24.2 mm accumulated, and would give PLTV an essentially unrelated wetness series).
        df_borriero['V'] = (df_borriero[influx]-df_borriero[et]-df_borriero[discharge]).cumsum()

        # fallback parameters, used when the model is built without the Monte Carlo search
        S0 = 1778
        k1 = 0.675
        k2 = 1.165
        df_borriero['S'] = S0 + df_borriero['V']
        w = (df_borriero['S']-df_borriero['S'].min())/(df_borriero['S'].max()-df_borriero['S'].min())
        df_borriero['k'] = k1+(1-w)*(k2-k1)

        # warm-up: "A 5-year warm-up period (i.e. repetition of the input data) from February 2008
        # to January 2013". The 828 d record is tiled, so V restarts each cycle instead of
        # accumulating its -67.5 mm net imbalance once per cycle (-338 mm over the five years).
        spin_end, spin_start = pd.Timestamp('2013-02-03'), pd.Timestamp('2008-02-01')
        n_loops = int(np.ceil(((spin_end-spin_start).days+1)/len(df_borriero)))
        spinup = pd.concat([df_borriero]*n_loops, ignore_index=True)
        spinup.index = pd.date_range(end=spin_end, periods=len(spinup), freq='D')
        spinup = spinup.loc[spin_start:]

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
    elif location == 'Bruntland Burn':
        bounds = { # Benettin et al. (2017) Table 1, Model 2 (time-variant Q SAS)
            "S0": (500, 4000), #(low, high); the paper's "average total storage" Stot
            "k1": (0.2, 3),    #kQ1, the wet-state exponent
            "k2": (0.2, 3),    #kQ2, the dry-state exponent
            "ket": (0.2, 3),
            "alpha": (0.95, 1.00), #evaporative fractionation factor for ET
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
            # Borriero et al. (2023) Table 2 search ranges: S0 300-3000 mm, every SAS parameter
            # 0.1-2. The commented values are the behavioral ranges instead (Table S1); the text
            # reports S0 acceptable anywhere in 335-2895 mm and kET unidentifiable over the whole
            # search range, so expect S0 and ket to stay flat however many samples are drawn.
            "S0": (300, 3000), #(618, 2875), #(low, high) for 95% CI
            "k1": (0.2, 2.0), #(0.4, 0.95),
            "k2": (0.2, 2.0), #(0.53, 1.8),
            "ket": (0.2, 2.0), #(0.2, 1.95),
        }
        return bounds
    elif location == 'Weierbach (2019)':
        bounds = { # Rodriguez and Klaus (2019) Table 2 initial ranges. Keys must stay in the order of WEIERBACH_PARAM_NAMES
            "Sref": (1500, 2500),
            "Sth": (80, 180),
            "f0": (0, 1),
            "dSth": (0.1, 5),
            "lamda1s": (0, 0.5), #lamda1s + lamda2 > 1 gives negative lamda3; make_SAS_model rejects those
            "Su": (1, 10),
            "lamda2": (0, 1),
            "mu2": (300, 700),
            "theta2": (0, 200),
            "mu3": (700, 1750),
            "theta3": (0, 200),
            "muET": (300, 1100),
            "thetaET": (0, 200),
        }
        return bounds
    elif location == 'Weierbach (2021)':
        bounds = { # Rodriguez et al. (2021) Table 1 initial ranges; Sref fixed at 2000 mm. Keys in the order of WEIERBACH_PARAM_NAMES
            "Sth": (20, 200),
            "dSth": (0.1, 20),
            "Su": (1, 50),
            "f0": (0, 1),
            "lamda1s_frac": (0, 1), #"lamda1s is uniformly sampled between 0 and 1-lamda2", so this is lamda1s/(1-lamda2)
            "lamda2": (0, 1),
            "mu2": (0, 1600),
            "theta2": (0, 100),
            "mu3": (0, 1600),
            "theta3": (0, 100),
            "muET": (0, 1600),
            "thetaET": (0, 100),
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
            p = results_df.loc[results_df['NSE']==results_df['NSE'].max()]
            p = p.to_numpy()
            p = np.delete(p, -1) # take off RMSE column for params
            p = np.delete(p, -1) # take of NSE column
            print('Best params:')
            for i in range(len(param_names)):
                print(f'{param_names[i]}: {p[i]}')
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
            if location.startswith('Weierbach'):
                # the composite SAS weights are data_df columns built from the paper parameters
                model = make_SAS_model(location, df, spinup, spin_done, influx)
            else:
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


def save_model_run(path, model, issample, **info):
    """
    Pickle a model that has been run so it can be reloaded without rerunning the spin-up/MC.

    The whole Model is saved, so the reloaded one keeps data_df (with the predicted
    '<iso_in> --> <discharge>' columns), options, the recorded-step index and the result arrays,
    i.e. get_pQ, get_sT, get_mT, ... all work as before. issample is saved with it because some
    catchments overwrite the one from load_data after the run.

    Parameters:
    path (str or Path): The .pkl file to write.
    model (Model): A model that has been run.
    issample (pd.Series): Boolean mask of the observed samples.
    **info: Anything else to keep with the run (e.g. spin_up, MC), returned on load.
    """
    import pickle
    from pathlib import Path
    from datetime import datetime
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    run = {'model': model, 'issample': issample, 'saved': datetime.now().isoformat(timespec='seconds'), **info}
    with open(path, 'wb') as f:
        pickle.dump(run, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f'Saved model run to {path}')


def load_model_run(path):
    """
    Reload a model run saved with save_model_run.

    Returns:
    dict: 'model', 'issample', 'saved' and any extra info given when saving.
    """
    import pickle
    with open(path, 'rb') as f:
        run = pickle.load(f)
    print(f"Loaded model run from {path} (saved {run['saved']})")
    return run
        