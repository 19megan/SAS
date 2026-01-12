from pcraster import *
from pcraster.framework import *

def iniParameters(self):
    
    ##########################################################################
    # P A R A M E T E R S                                                    #
    ##########################################################################    
    
    ################################################################## MAPS ##
    self.dem        = readmap(self.input_folder+"/maps/dem_new.map")          # DEM
    self.LAI        = scalar(0)                                               # Leaf Area Index, i.e. the output fro this example correnponds to the "open" scenario, without canopy 

    ################################################################ SCALAR ##
    self.sfCorr   = scalar(0.15)                                              # ADDITIONAL correction factor for snowfall on top of wind correction
    self.ttlow    = scalar(-2)                                                # temperature [C] below which all precip as snow
    self.tthigh   = scalar(2)                                                 # temperature [C] above which all precip as rain
    self.RDthres  = scalar(0.001)                                             # amount of precipitation [m] above which a day is considered cloudy (affects long wave inputs)
    self.albGround= scalar(0.21)                                              # albedo for snowfree areas
    self.deplOffset = scalar(3.5)                                             # offset of snowmelt concentration, equilibrium frationaltion is 3.5 permiull
    self.retcap     = scalar(0.05)                                            # maximum water rentention capacity of the snowpack
    self.ccov       = scalar(0)                                               # canopy cover
    self.IsnowMax   = scalar(0.0001) #self.LAI * 4.4                          # maximum interception capacity [kg /m2], For numerical stability here set to 0.0001, when LAI >0 according to Liston et al 2006. 
    self.IunDD      = scalar(5)                                               # unloading factor for interception storage kg/m2/d/C (Liston et al 2006)
    self.Isroo      = scalar(150)                                             # snow density in snow interception storage [kg/m3]
    metAlt	  = scalar(283)		                                      # elevation of the meteorological measurements
    cAtten        = scalar(1)                                                 # attenuation coeffient for canopy longwave radiation (Montehit and Unsworth 2013 fig 8.2)
    self.altDiff    = self.dem - metAlt                                       # altitude difference between model cell and meteorological station
    self.albPow     = scalar(1)                                               # parameter lowering albedo for aging snow 
    self.Efrac      = scalar(15)                                              # offset tha evaporative fractination is to introduce to the snowpack
    self.Tgrad	    = scalar(-0.006)                                          # temperature lapse rate [C/m]
    
    
       
    ##########################################################################
    # I N I T I A L    C O N D I T I O N S                                   #
    ##########################################################################      

    ################################################### FLUXES AND STORAGES ##
    self.Pmelt=scalar(0)                                                      # water output from the snowpack
    self.liquidP=scalar(0)                                                    # liquid water in the snowpack
    self.Iunload    = scalar(0)                                               # unloading from snow interception storage [kg/m2]
    self.Isnow      = scalar(0)                                               # snow canopy interception storage [kg/m2]
    self.dNoSnow    = scalar(1)                                               # days since snowfall
    self.dSnow      = scalar(0)                                               # days with snow on the ground
    self.dMelt      = scalar(0)                                               # number of snowmelt days
    self.tau0       = exp(-cAtten * self.LAI)                                 # fraction of shortwave radiation transmitted by canopy (Wigmosta 1994)
    self.ST         = scalar(-4)                                              # snow temperature
    self.SWE        = scalar(0.2)                                             # snow water equivalent [m]
    self.Wice       = scalar(0.2)                                             # ice in snowpack
    self.Wliq       = scalar(0)                                               # liquid water in snowpack
    self.Qs         = scalar(0)                                               # sensible heat
    self.Ql         = scalar(0)                                               # latent heat
    self.Qp         = scalar(0)                                               # heat in precip
    self.Qm         = scalar(0)                                               # heat from melt/refreezing
    self.M          = scalar(0)                                               # amount of melt
    self.radRns     = scalar(0)                                               # net radiation
    
    ###################################################### ISOTOPES AND AGE ##
    self.concsn     = scalar(-25)                                             # initial snowpack isotope concentration
    self.concIst    = scalar(-25)                                             # initial snow interception storate isotope composition
    self.agesnow    = scalar(1)                                               # age of snow
    self.ageIst     = scalar(0)                                               # age of interception storage
   
    
    ############################################################################
    # O U T P U T   W R I T I N G                                              #
    ############################################################################
    
    if self.reportTimeseries:
        
        # South facing slope (SF)
        self.tssSWE_SF=TimeoutputTimeseries("Output/SWE_SF", self,self.input_folder+"/maps/output_SF.map",noHeader=True)        # SWE [mm]
        self.tssIsnow_SF=TimeoutputTimeseries("Output/Isnow_SF", self,self.input_folder+"/maps/output_SF.map",noHeader=True)    # interception storage [mm]
        self.tssPmelt_SF=TimeoutputTimeseries("Output/Pmelt_SF", self,self.input_folder+"/maps/output_SF.map",noHeader=True)    # meltwater
        self.tssIunlm_SF=TimeoutputTimeseries("Output/Iunlm_SF", self,self.input_folder+"/maps/output_SF.map",noHeader=True)    # unloading from interception storage
        self.tssconcsn_SF=TimeoutputTimeseries("Output/concsn_SF", self,self.input_folder+"/maps/output_SF.map",noHeader=True)  # snow concentration
        self.tssconcmelt_SF=TimeoutputTimeseries("Output/concmelt_SF", self,self.input_folder+"/maps/output_SF.map",noHeader=True)  # meltwater concentration
        
        # North facing slope
        self.tssSWE_NF=TimeoutputTimeseries("Output/SWE_NF", self,self.input_folder+"/maps/output_NF.map",noHeader=True)
        self.tssIsnow_NF=TimeoutputTimeseries("Output/Isnow_NF", self,self.input_folder+"/maps/output_NF.map",noHeader=True)    # interception storage [mm]
        self.tssPmelt_NF=TimeoutputTimeseries("Output/Pmelt_NF", self,self.input_folder+"/maps/output_NF.map",noHeader=True)    # meltwater
        self.tssIunlm_NF=TimeoutputTimeseries("Output/Iunlm_NF", self,self.input_folder+"/maps/output_NF.map",noHeader=True)    # unloading from interception storage
        self.tssconcsn_NF=TimeoutputTimeseries("Output/concsn_NF", self,self.input_folder+"/maps/output_NF.map",noHeader=True)  # snow concentration
        self.tssconcmelt_NF=TimeoutputTimeseries("Output/concmelt_NF", self,self.input_folder+"/maps/output_NF.map",noHeader=True)  # meltwater concentration
        
        # Catchment averages
        self.tssPmtot=TimeoutputTimeseries("Output/Pmtot", self,self.input_folder+"/maps/output_NF.map",noHeader=True)          # meltwater flux
        self.tssmIsotot=TimeoutputTimeseries("Output/mIsotot", self,self.input_folder+"/maps/output_NF.map",noHeader=True)      # meltwater concentration
        self.tssSWEtot=TimeoutputTimeseries("Output/SWEtot", self,self.input_folder+"/maps/output_NF.map",noHeader=True)        # SWE
        self.tsssnIsotot=TimeoutputTimeseries("Output/snIsotot", self,self.input_folder+"/maps/output_NF.map",noHeader=True)    # snow concentration


###################################