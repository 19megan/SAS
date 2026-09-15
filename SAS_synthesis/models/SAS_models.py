# Contains all the SAS models with necessary mT_init as a string able to be passed
# Date created: 09/14/2026

from mesas.sas.model import Model
import numpy as np




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
    if location == "Weierbach (2019)":
        def make_Weierbach2019_model_from(params):
            _ = params
            sas_specs={
                "Q [m3/4h]":{
                    "lamda1":{"func": "kumaraswamy",
                        "args": {"loc": 0.0,"scale": 4.1,"a": 1.0,"b": 1.0}
                    },
                    "lamda2":{"func": "gamma",
                        "args": {"a": 16.9117647058824,"scale": 34,"loc": 0}
                    },
                    "lamda3":{"func": "gamma",
                        "args": {"a": 11.7241379310345,"scale": 87,"loc": 0}
                    }},
                "ET [m3/4h]":{
                    "ET SAS function":{"func": "gamma",
                        "args": {"loc": 0.0,"scale": 77,"a": 3.96103896103896}
                    }}}
            solute_parameters={"d2H (permil)":{"C_old": -50.6, "observations": "Q d2H (permil)", 'mT_init': 'mT_init'}}
            return Model(data_df=df, sas_specs=sas_specs, solute_parameters=solute_parameters, verbose=False, n_substeps=1, record_state=True, influx=influx, sT_init=sT_init)

        return make_Weierbach2019_model_from(df)


    elif location == "Weierbach (2021)":
        def make_Weierbach2021_model_from(params):
            S0, k1, k2, ket = params
            df_benettin['S'] = df_benettin[influx]-df_benettin[et]-df_benettin[discharge]+S0
            w = (df_benettin['S']-df_benettin['S'].min())/(df_benettin['S'].max()-df_benettin['S'].min())
            df_benettin['k'] = k1+(1-w)*(k2-k1)
            sas_specs = {
                "Q [mm/d]":{
                    "Q SAS function":{
                        "func": "kumaraswamy",
                        "args": {"loc": 0.0,"scale": "S", "a": "k", "b": 1.0}}},
                "ET [mm/d]":{
                    "ET SAS function":{
                        "func": "kumaraswamy",
                        "args": {'loc': 0.0, 'scale': "S", 'a': ket, 'b': 1.0}}}}
            solute_parameters = {'d18O_P raw': {'C_old': -9.2, 'observations': 'd18O_Q', 'mT_init': 'mT_init'}}
            return Model(data_df=df_benettin, sas_specs=sas_specs, solute_parameters=solute_parameters, verbose=False, n_substeps=1, record_state=True, influx=influx, sT_init=sT_init)
        return make_Weierbach2021_model_from(S0=0.0, k1=0.1, k2=0.5, ket=0.3) #wrong. update these

    
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
        if spin_done:
            def make_Selke_model_from(params):
                S0, k1, k2, ket = params
                df['S'] = df[influx]-df["ET [mm/d]"]-df['Q [mm/d]']+S0
                w = (df['S']-df['S'].min())/(df['S'].max()-df['S'].min())
                df['k'] = k1+(1-w)*(k2-k1)
                sas_specs = {
                    "Q [mm/d]":{
                        "Q SAS function":{
                            "func": "kumaraswamy",
                            "args": {"loc": 0.0,"scale": "S", "a": "k", "b": 1.0}}},
                    "ET [mm/d]":{
                        "ET SAS function":{
                            "func": "kumaraswamy",
                            "args": {'loc': 0.0, 'scale': "S", 'a': ket, 'b': 1.0}}}}
                solute_parameters = {'d18O_P raw': {'C_old': -9.2, 'observations': 'd18O_Q', 'mT_init': 'mT_init'}}
                return Model(data_df=df, sas_specs=sas_specs, solute_parameters=solute_parameters, verbose=False, n_substeps=1, record_state=True, influx=influx, sT_init=sT_init)
            return make_Selke_model_from(params)
        else:
            def make_Selke_spinup_model_from(params):
                S0, k1, k2, ket = params
                spinup['S'] = spinup[influx]-spinup["ET [mm/d]"]-spinup['Q [mm/d]']+S0
                w = (spinup['S']-spinup['S'].min())/(spinup['S'].max()-spinup['S'].min())
                spinup['k'] = k1+(1-w)*(k2-k1)
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
                return Model(data_df=spinup, sas_specs=sas_specs, solute_parameters=solute_parameters, verbose=False, n_substeps=1, record_state=True, influx=influx)
            return make_Selke_spinup_model_from(params)

