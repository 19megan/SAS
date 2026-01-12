from pcraster import *
from pcraster.framework import *
import os

def reportAll(self):
    
    if self.reportSpatial:
        self.report(self.SWE,"SWE")
                
        # flow variables
        self.report(self.Wliqout*1000,"Wliqout")                                  
        self.report(-1000*self.Eice,"Eice")
        
        #concentrations/ages
        self.report(self.concmelt,"concmelt")
        self.report(self.concsn,"concsn")

        
    else:
        print('no spatial output')

    self.NoCells = scalar(338)    
    
    if self.reportTimeseries:
        
        # South Facing
        self.tssSWE_SF.sample(self.SWE*1000)
        self.tssIsnow_SF.sample(1000*self.Isnow/self.Isroo)
        self.tssPmelt_SF.sample(self.Pmelt)
        self.tssIunlm_SF.sample(self.Iunlm*1000)
        self.tssconcsn_SF.sample(self.concsn)
        self.tssconcmelt_SF.sample(self.concmelt)
        
        # North facing
        self.tssSWE_NF.sample(self.SWE*1000)
        self.tssIsnow_NF.sample(1000*self.Isnow/self.Isroo)
        self.tssPmelt_NF.sample(self.Pmelt)
        self.tssIunlm_NF.sample(self.Iunlm*1000)
        self.tssconcsn_NF.sample(self.concsn)
        self.tssconcmelt_NF.sample(self.concmelt)
        
        # CATCHMENT AVERAGED fluxes and storages
        
        # cathcment averaged meltwater water flux    
        self.Pm_tot = maptotal(self.Pmelt)/self.NoCells
        self.tssPmtot.sample(self.Pm_tot)                                       
        
        # catchment averaged concentration of snowmelt output
        self.mIso_tot = ifthenelse(self.Pm_tot == 0,                          
                                 0, 
                                 maptotal(self.Pmelt * self.concmelt) / 
                                         (maptotal(self.Pmelt)+0.00001))
        self.tssmIsotot.sample(self.mIso_tot)
        
        # catchment averaged SWE
        self.SWE_tot = maptotal(self.SWE*1000)/self.NoCells
        self.tssSWEtot.sample(self.SWE_tot)                                  
        
        # catchment averaged concentration of snowpack
        self.snIso_tot = ifthenelse(self.SWE_tot == 0,                          
                                 0, 
                                 maptotal(self.SWE * self.concsn) / 
                                         (maptotal(self.SWE)+0.000001))
        self.tsssnIsotot.sample(self.snIso_tot)
        


    else:
        print('no timeseries output')
 
