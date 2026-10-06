# Contains all the SAS models with necessary mT_init as a string able to be passed
# Date created: 09/14/2026

from mesas.sas.model import Model
import numpy as np
import pandas as pd
from math import tanh


# Weierbach composite SAS model parameters
WEIERBACH_PARAM_NAMES = { # order of the sampled parameter vector (make_prior_bounds keys)
    'Weierbach (2019)': ['Sref', 'Sth', 'f0', 'dSth', 'lamda1s', 'Su', 'lamda2', 'mu2', 'theta2', 'mu3', 'theta3', 'muET', 'thetaET'],
    'Weierbach (2021)': ['Sth', 'dSth', 'Su', 'f0', 'lamda1s_frac', 'lamda2', 'mu2', 'theta2', 'mu3', 'theta3', 'muET', 'thetaET'],
}
_WEIERBACH_2019 = dict(Sref=2426, Sth=105, f0=0.1, dSth=3.97, lamda1s=0.11, Su=4.1, lamda2=0.32,
                       mu2=575, theta2=34, mu3=1020, theta3=87, muET=305, thetaET=77)
WEIERBACH_PARAMS = { # defaults, used without Monte Carlo
    # Rodriguez and Klaus (2019) Table 2, calibrated values (the third best set, NSE = 0.255)
    'Weierbach (2019)': _WEIERBACH_2019,
    # Rodriguez et al. (2021) fix Sref = 2000 mm and report only behavioral ranges, no single best
    # set, so the 2019 values stand in for the SAS parameters until this is calibrated
    'Weierbach (2021)': {**_WEIERBACH_2019, 'Sref': 2000},
}


def weierbach_hydrology(spinup, df, Sref, Sth, dSth, f0, lamda1s, lamda2, n_cycles=20, n=20, m=1000):
    """
    Storage, actual ET and composite SAS weights for the Weierbach (Rodriguez and Klaus, 2019,
    equations 4 and 13-19; Rodriguez et al., 2021, equations A5-A9). Writes the columns 'S',
    'ET [mm/4h]', 'lamda1', 'lamda2' and 'lamda3' into spinup and df in place.

    S(t) = Sref + integral(J - Q - ET), starting at the beginning of the 100 year spin-up, with
    ET(t) = PET(t)*tanh((S(t)/Sroot)^n) and Sroot = Sref - 150 mm in both papers (2019 section 2.7,
    2021 eq. A9; the "150" of the 2021 Table 1 is the water accessible by ET, Sref - Sroot). ET
    depends on S, so the balance is stepped forward. The October 2010 - October 2015 cycle is
    looped n_cycles-1 times before df (whose first 5 years are that cycle again), i.e. 100 years in
    total for n_cycles = 20. S settles where the reduced ET closes the cycle's water balance.

    Parameters:
    spinup, df (pd.DataFrame): from _weierbach_data, spinup being whole loops of the cycle.
    Sref (float): initial storage (mm).
    Sth, dSth, f0, lamda1s, lamda2 (float): parameters of the time-variant weights.
    """
    J, Q, PET = 'J [mm/4h]', 'Q [mm/4h]', 'PET [mm/4h]'
    cycle = df.loc[~df['main']]
    assert len(spinup) % len(cycle) == 0 and len(spinup) <= (n_cycles-1)*len(cycle)

    def chain(col):
        return np.concatenate([np.tile(cycle[col].to_numpy(), n_cycles-1), df[col].to_numpy()])
    j, q, pet = chain(J), chain(Q), chain(PET)

    Sroot = Sref-150
    S = np.empty(len(j))
    ET = np.empty(len(j))
    s = Sref
    for i in range(len(j)):
        S[i] = s
        ET[i] = pet[i]*tanh((s/Sroot)**n)
        s = s + j[i] - q[i] - ET[i]

    # storage variations over the moving window dt* = 2 dt, i.e. this step and the two before (eq. 16-17)
    dS = j - q - ET
    dSbar = np.maximum(pd.Series(dS).rolling(3, min_periods=1).mean().to_numpy(), 0)
    # Smin over the period of interest, October 2015 - October 2017 (cf. Figure 5 of the 2019 paper)
    main = np.zeros(len(j), dtype=bool)
    main[len(j)-len(df):] = df['main'].to_numpy()
    Smin = S[main].min()
    with np.errstate(over='ignore'): # (S/(Smin+Sth))^1000 overflows to inf when S is well above, tanh -> 1
        f = f0*(1-np.tanh((S/(Smin+Sth))**m)) # eq. 14
    g = 1-np.exp(-dSbar/dSth) # eq. 15
    lamda1 = lamda1s*(f+(1-f)*g) # eq. 13

    for frame, seg in ((df, slice(len(j)-len(df), None)), (spinup, slice(len(j)-len(df)-len(spinup), len(j)-len(df)))):
        frame['S'] = S[seg]
        frame['ET [mm/4h]'] = ET[seg]
        frame['lamda1'] = lamda1[seg]
        frame['lamda2'] = lamda2
        frame['lamda3'] = 1-lamda2-lamda1[seg]


def make_SAS_model(location, df, spinup, spin_done, influx, sT_init=None, params=None):
    """
    Create a SAS model for a given location and data frame.
    
    Parameters:
    location (str): The location for which the model is created.
    df (pd.DataFrame): The data frame containing the data.
    spinup (pd.DataFrame): The data frame containing the spin-up data.
    spin_done (bool): Whether the spin-up has been completed.
    influx (str): The name of the influx column in the data frame.
    sT_init (pd.Series): The initial storage state for the model.
    params (list): A list of parameter values for the model.
    
    Returns:
    Model: An instance of the Model class.
    
    """
    if location in ("Weierbach (2019)", "Weierbach (2021)"):
        # Composite Q SAS function (2019 eq. 9-12, 2021 eq. 5 and A1-A4): a uniform on [0, Su] for
        # young event water plus two gammas for older water, weighted by the data_df columns
        # 'lamda1' (time-variant), 'lamda2' (constant) and 'lamda3' = 1 - lamda1 - lamda2 (mesas
        # takes a component's weight from the column with its name). ET SAS is a single gamma.
        # The gammas are given by mean mu and scale theta, so the shape is a = mu/theta.
        # Both papers use deuterium; 2021 also calibrates on tritium, which needs 3H in
        # precipitation (GNIP Trier, IAEA WISER) that is not in the data folder yet.
        def make_Weierbach_model_from(params, data_df):
            p = dict(WEIERBACH_PARAMS[location])
            if params is not None:
                p.update(zip(WEIERBACH_PARAM_NAMES[location], params))
            if 'lamda1s_frac' in p: # 2021: lamda1s sampled in [0, 1-lamda2] so the weights stay >= 0
                p['lamda1s'] = p.pop('lamda1s_frac')*(1-p['lamda2'])
            if p['lamda1s']+p['lamda2'] > 1:
                raise ValueError('lamda1s + lamda2 > 1 gives a negative lamda3')
            # S, ET and the weights depend on Sref and the weight parameters, so rebuild them
            weierbach_hydrology(spinup, df, p['Sref'], p['Sth'], p['dSth'], p['f0'], p['lamda1s'], p['lamda2'])
            sas_specs = {
                "Q [mm/4h]":{
                    "lamda1":{"func": "kumaraswamy",
                        "args": {"loc": 0.0, "scale": p['Su'], "a": 1.0, "b": 1.0}}, # uniform on [0, Su]
                    "lamda2":{"func": "gamma",
                        "args": {"loc": 0.0, "scale": p['theta2'], "a": p['mu2']/p['theta2']}},
                    "lamda3":{"func": "gamma",
                        "args": {"loc": 0.0, "scale": p['theta3'], "a": p['mu3']/p['theta3']}}},
                "ET [mm/4h]":{
                    "ET SAS function":{"func": "gamma",
                        "args": {"loc": 0.0, "scale": p['thetaET'], "a": p['muET']/p['thetaET']}}}}
            # C_old: the 2019 initial storage d2H, Cref = -50.6 permil (flux-weighted mean stream d2H
            # 2015-2017). Here it is the d2H of water older than the tracked ages.
            solute_parameters = {'Cin d2H': {'C_old': -50.6, 'observations': 'measC_Q d2H'}}
            if spin_done:
                solute_parameters['Cin d2H']['mT_init'] = 'mT_init'
                # record daily over October 2015 - October 2017 only ('record' column); max_age is
                # len(data_df) = 7 years (from len(sT_init))
                return Model(data_df=data_df, sas_specs=sas_specs, solute_parameters=solute_parameters,
                             verbose=False, n_substeps=1, influx=influx, sT_init=sT_init,
                             record_state='record', record_arrays=('sT', 'pQ'), record_dtype='float32')
            return Model(data_df=data_df, sas_specs=sas_specs, solute_parameters=solute_parameters,
                         verbose=False, n_substeps=1, influx=influx, record_state=False)

        return make_Weierbach_model_from(params, df if spin_done else spinup)

    
    elif location == "Bruntland Burn":
        # Benettin et al. (2017), Model 2. The SAS functions are power laws of the normalized rank
        # storage, Omega(Ps) = Ps^k (equation 5), which is kumaraswamy with b = 1 and scale = S.
        # kQ is time-variant: ws = (S-Smin)/(Smax-Smin), kQ(t) = kQ1 + (1-ws(t))*(kQ2-kQ1), so kQ1
        # is the wet-state exponent and kQ2 the dry-state one. kET is constant.
        # Reference (best-performing) parameters: kQ1 = 0.36, kQ2 = 0.80, kET = 0.90, S0 = 2400,
        # alpha = 0.991. The initial storage isotopic content is fixed at -58 permil, not calibrated.
        def make_BruntlandBurn_model_from(params, data_df):
            S0, k1, k2, ket, alpha = params
            # S_rel is the relative storage from the water balance and does not depend on the
            # calibrated parameters, so only the offset is applied here. S0 is treated as the
            # AVERAGE total storage (the symbol used in the paper's Table 1); it is also referred
            # to as the initial storage in the text, in which case use "+ S0" alone.
            data_df['S'] = data_df['S_rel'] - data_df['S_rel'].mean() + S0
            ws = (data_df['S']-data_df['S'].min())/(data_df['S'].max()-data_df['S'].min())
            data_df['k'] = k1 + (1-ws)*(k2-k1)
            sas_specs = {
                "Q [mm/h]":{
                    "Q SAS function":{
                        "func": "kumaraswamy",
                        "args": {"loc": 0.0, "scale": "S", "a": "k", "b": 1.0}}},
                "ET [mm/h]":{
                    "ET SAS function":{
                        "func": "kumaraswamy",
                        "args": {"loc": 0.0, "scale": "S", "a": ket, "b": 1.0}}}}
            # evaporative fractionation: ET leaves depleted by alpha relative to storage (Appendix B)
            solute_parameters = {'Cin d2H': {'C_old': -58.0, 'observations': 'measC_Q d2H',
                                             'alpha': {'Q [mm/h]': 1.0, 'ET [mm/h]': alpha}}}
            # The recorded state is (max_age x timesteps), so at an hourly timestep record_state=True
            # would allocate 26280^2 cells per array (>11 GB for pQ alone, and the 8 yr spin-up is
            # 70080 steps = 73 GB). Neither is needed:
            #  - spin-up: record_state=False keeps only the final timestep, which is all that is read
            #    back for sT_init/mT_init.
            #  - main run: record daily (record_every=24) in float32. The isotope prediction C_Q is
            #    solved at every hourly step regardless of what state is recorded, so NSE/RMSE are
            #    unaffected; only the TTD snapshots are thinned to one per day.
            if spin_done:
                solute_parameters['Cin d2H']['mT_init'] = 'mT_init'
                return Model(data_df=data_df, sas_specs=sas_specs, solute_parameters=solute_parameters,
                             verbose=False, n_substeps=1, influx=influx, sT_init=sT_init,
                             record_every=24, record_dtype='float32')
            return Model(data_df=data_df, sas_specs=sas_specs, solute_parameters=solute_parameters,
                         verbose=False, n_substeps=1, influx=influx, record_state=False)

        return make_BruntlandBurn_model_from(params, df if spin_done else spinup)


    elif location == "Chenqi":
        if spin_done:
            def make_Chenqi_model_from(params):
                k1, k2, ket, S0, f, alpha, C_old = params
                df['Q (mm/day)'] = df['Q (m3/D)']/(alpha*(1.25*1e6))*1000 #convert to mm/day
                beta = (df['P(mm)'].sum()-df['Q (mm/day)'].sum())/df['PET (mm)'].sum() #0.580 (0.584 and 0.617 for the two full water years)
                df['ET (mm/day)'] = beta*df['PET (mm)']
                df['S'] = df['P(mm)']-df['ET (mm/day)']-df['Q (mm/day)']+S0    
                df['w'] = (df['Q (mm/day)']-df['Q (mm/day)'].min())/(df['Q (mm/day)'].max()-df['Q (mm/day)'].min())
                df['k'] = k1 + (1-df['w'])*(k2-k1)

                sas_specs = {
                    "Q (mm/day)":{
                        "Q SAS function":{
                            "func": "kumaraswamy",
                            "args": {"loc": 0.0,"scale": "S", "a": "k", "b": 1.0}}},
                    "ET (mm/day)":{
                        "ET SAS function":{
                            "func": "kumaraswamy",
                            "args": {'loc': 0.0, 'scale': "S", 'a': ket, 'b': 1.0}}}}
                solute_parameters = {'rain-D (‰)': {'C_old': C_old, 'observations': 'outlet-D (‰)', 'alpha': {'Q (mm/day)': 1, 'ET (mm/day)': 1-f}, 'mT_init': 'mT_init'}}
                return Model(data_df=df, sas_specs=sas_specs, solute_parameters=solute_parameters, verbose=False, n_substeps=1, record_state=True, influx=influx, sT_init=sT_init)
            return make_Chenqi_model_from(params)
        else:
            def make_Chenqi_spinup_model_from(params):
                k1, k2, ket, S0, f, alpha, C_old = params
                spinup['Q (mm/day)'] = spinup['Q (m3/D)']/(alpha*(1.25*1e6))*1000 #convert to mm/day
                beta = (spinup[influx].sum()-spinup['Q (mm/day)'].sum())/spinup['PET (mm)'].sum() #0.580 (0.584 and 0.617 for the two full water years)
                spinup['ET (mm/day)'] = beta*spinup['PET (mm)']
                spinup['S'] = spinup[influx]-spinup['ET (mm/day)']-spinup['Q (mm/day)']+S0    
                spinup['w'] = (spinup['Q (mm/day)']-spinup['Q (mm/day)'].min())/(spinup['Q (mm/day)'].max()-spinup['Q (mm/day)'].min())
                spinup['k'] = k1 + (1-spinup['w'])*(k2-k1)

                sas_specs = {
                    "Q (mm/day)":{
                        "Q SAS function":{
                            "func": "kumaraswamy",
                            "args": {"loc": 0.0,"scale": "S", "a": "k", "b": 1.0}}},
                    "ET (mm/day)":{
                        "ET SAS function":{
                            "func": "kumaraswamy",
                            "args": {'loc': 0.0, 'scale': "S", 'a': ket, 'b': 1.0}}}}
                solute_parameters = {'rain-D (‰)': {'C_old': C_old, 'observations': 'outlet-D (‰)', 'alpha': {'Q (mm/day)': 1, 'ET (mm/day)': 1-f}}}
                return Model(data_df=spinup, sas_specs=sas_specs, solute_parameters=solute_parameters, verbose=False, n_substeps=1, record_state=True, influx=influx)
            return make_Chenqi_spinup_model_from(params) 

    elif location == "CanVilla":
        def make_CanVilla_model_from(params):
            S0, k1, k2 = params
            df['S'] = df[influx]-df['ET [mm/h]']-df['Q [mm/h]']+S0
            df['k'] = k1+(1-df['wi [-]'])*(k2-k1)
            sas_specs = {
                "Q [mm/h]":{
                    "Q SAS function":{
                        "func": "kumaraswamy",
                        "args": {"loc": 0.0,"scale": "S", "a": "k", "b": 1.0}}},
                "ET [mm/h]":{
                    "ET SAS function":{
                        "func": "kumaraswamy",
                        "args": {'loc': 0.0, 'scale': "S", 'a': 3.26, 'b': 1.0}}}}
            solute_parameters = {'Cin [-]': {'C_old': -7.7, 'observations': 'measC_Q [-]', 'mT_init': 'mT_init'}}
            return Model(data_df=df, sas_specs=sas_specs, solute_parameters=solute_parameters, verbose=False, n_substeps=1, record_state=True, influx=influx, sT_init=sT_init)
        return make_CanVilla_model_from([389.91, 0.28, 1.26])

    elif location == "Dry Creek":
        if spin_done:
            def make_DryCreek_model_from(params):
                S0, k1, k2, logfactor_Q, ket, C_old = params
                df['S'] = df[influx]-df["ET [mm/4h]"]-df['Q [mm/4h]']+S0    
                df['wi [-]'] = np.log(df['Q [mm/4h]']).values/np.log(df['Q [mm/4h]']).max()
                # # df_l['k'] = k1+(k2-k1)*np.log((1-logfactor_Q)*df_l['wi [-]']) #in the paper
                df['k'] = k1+(k2-k1)/np.log(logfactor_Q)*np.log(logfactor_Q-(logfactor_Q-1)*df['wi [-]']) #in their code

                sas_specs = {
                    "Q [mm/4h]":{
                        "Q SAS function":{
                            "func": "kumaraswamy",
                            "args": {"loc": 0.0,"scale": "S", "a": "k", "b": 1.0}}},
                    "ET [mm/4h]":{
                        "ET SAS function":{
                            "func": "kumaraswamy",
                            "args": {'loc': 0.0, 'scale': "S", 'a': ket, 'b': 1.0}}}}
                solute_parameters = {'Cin [-]': {'C_old': C_old, 'observations': 'measC_Q', 'mT_init': 'mT_init'}}
                return Model(data_df=df, sas_specs=sas_specs, solute_parameters=solute_parameters, verbose=False, n_substeps=1, record_state=True, influx=influx, sT_init=sT_init)
            return make_DryCreek_model_from(params)
        else:
            def make_DryCreek_spinup_model_from(params):
                S0, k1, k2, logfactor_Q, ket, C_old = params
                spinup['S'] = spinup[influx]-spinup["ET [mm/4h]"]-spinup['Q [mm/4h]']+S0    
                spinup['wi [-]'] = np.log(spinup['Q [mm/4h]']).values/np.log(spinup['Q [mm/4h]']).max()
                # # df_l['k'] = k1+(k2-k1)*np.log((1-logfactor_Q)*df_l['wi [-]']) #in the paper
                spinup['k'] = k1+(k2-k1)/np.log(logfactor_Q)*np.log(logfactor_Q-(logfactor_Q-1)*spinup['wi [-]']) #in their code

                sas_specs = {
                    "Q [mm/4h]":{
                        "Q SAS function":{
                            "func": "kumaraswamy",
                            "args": {"loc": 0.0,"scale": "S", "a": "k", "b": 1.0}}},
                    "ET [mm/4h]":{
                        "ET SAS function":{
                            "func": "kumaraswamy",
                            "args": {'loc': 0.0, 'scale': "S", 'a': ket, 'b': 1.0}}}}
                solute_parameters = {'Cin [-]': {'C_old': C_old, 'observations': 'measC_Q'}}
                return Model(data_df=spinup, sas_specs=sas_specs, solute_parameters=solute_parameters, verbose=False, n_substeps=1, record_state=True, influx=influx)
            return make_DryCreek_spinup_model_from(params)

    elif location == "Selke":
        # Borriero et al. (2023), setup e: power law time-variant (PLTV) SAS for Q, which is
        # kumaraswamy with b = 1 and scale = S, and kQ(t) linear in the catchment wetness
        # wi = (S-Smin)/(Smax-Smin). C_old = -9.2 is the paper's initial storage composition, set
        # to the mean measured d18O_Q. No evaporative fractionation (paper equation 2 has none).
        def make_Selke_model_from(params, data_df):
            S0, k1, k2, ket = params
            # paper equation 1: S(t) = S0 + V(t), with the accumulated V built in load_data
            data_df['S'] = S0 + data_df['V']
            w = (data_df['S']-data_df['S'].min())/(data_df['S'].max()-data_df['S'].min())
            data_df['k'] = k1+(1-w)*(k2-k1)
            sas_specs = {
                "Q [mm/d]":{
                    "Q SAS function":{
                        "func": "kumaraswamy",
                        "args": {"loc": 0.0,"scale": "S", "a": "k", "b": 1.0}}},
                "ET [mm/d]":{
                    "ET SAS function":{
                        "func": "kumaraswamy",
                        "args": {'loc': 0.0, 'scale': "S", 'a': ket, 'b': 1.0}}}}
            solute_parameters = {'d18O_P raw': {'C_old': -9.2, 'observations': 'd18O_Q'}}
            if spin_done:
                solute_parameters['d18O_P raw']['mT_init'] = 'mT_init'
                return Model(data_df=data_df, sas_specs=sas_specs, solute_parameters=solute_parameters,
                             verbose=False, n_substeps=1, record_state=True, influx=influx, sT_init=sT_init)
            return Model(data_df=data_df, sas_specs=sas_specs, solute_parameters=solute_parameters,
                         verbose=False, n_substeps=1, record_state=True, influx=influx)

        return make_Selke_model_from(params, df if spin_done else spinup)

