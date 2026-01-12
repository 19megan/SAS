from pcraster import *
from pcraster.framework import *
import waterBalanceCheck
# energy balance snow module in the STARR model


def update(self, P, T, Pconc):                        	                                                
    # specify snow relevant constants
    cSB             = scalar(4.89E-9)                                           # Stefan-Boltzmann constant {MJ/day*m2*K^4)
    cair            = scalar(1.29E-3)                                           # heat capacity of air [MJ/m3*C]
    zu              = scalar(2)                                                 # measurement height of climate variables (Walter et al 2005)
    ds              = scalar(0)                                                 # zero-plane dispalcement for snow [m] (Walter et al 2005)
    zms             = scalar(0.001)                                             # momentum roughness for snow [m] (Walter et al 2005)
    zhs             = scalar(0.0002)                                            # heat and vapour roughness parameter for snow [m] (Walter et al 2005)
    k               = scalar(0.41)                                              # von Karman's constant
    Rv              = scalar(4.63E-3)                                           # gas constant for water vapour
    Rt              = scalar(0.4615)                                            # thermodynamic vapour constant [kJ/kg*K]
    lamv            = scalar(2.800)                                             # latent heat of vaporization [MJ/kg]
    lamf            = scalar(0.333)                                             # latent heat of fusion [MJ/kg]
    rooW            = scalar(1000)                                              # density of water [kg/m3]
    Cw              = scalar(4.2E-3)                                            # specific heat capacity of water [MJ/kg*C]
    Ci              = scalar(2.03E-3)                                           # specific heat capacity of ice [MJ/kg*C]                               

    ## RENAME and convert units of some of the INPUTS
    self.AT 	    = T;                 		                        # comform to naming conventions elsewhere in STARR  
    self.P  	    = P*0.001;		                                        # convert from [mm/d] to [m/d]
    self.WS         = self.WSin * (1-(0.8*self.ccov))                           # reduction of windspeed due to vegetation (Tarboton and Luke 1996)
    self.Pconc      = Pconc                                                     # precipitation isotope timeseries
                    
    ##################################################### RAINFALL SEPARATION ##
    ############################################################################
    self.Pq         = ifthenelse(self.AT <= self.ttlow,1,                       # calculate the thermal quality of snow (ice or liquid)
                                 ifthenelse(pcrand(self.AT > self.ttlow,self.AT < self.tthigh), # with a range where precip is a mix of rain snow defined by temperature threshold parameters parameters
                                            (self.tthigh - self.AT)/(self.tthigh-self.ttlow),
                                            0))
    self.Pliq       = self.P * (1-self.Pq)                                      # all precip that isnt snow is liquid 
    self.liquidP    = self.Pliq * 1000                                          # COMFORM to output to other parts of STARR and convert to [mm/d]                                
    self.cSnow      =  1 / (ifthenelse(self.AT<0,                               # correct to undercatch due to wind (Yang et al 1998)
                                        max(exp(4.606 - 0.036 * self.WSin**1.75),38.6), # correction for cold air temperatures    
                                        max(101.04 - 5.62 * self.WSin,64.5))/100) # correction for range 0-3 C
    self.cSnow      = ifthenelse(self.cSnow < 1.5,self.cSnow+self.sfCorr,self.cSnow) # additional correction coeffient sfCorr that can be used in snow model calibration during low wind speeds      
    
    self.Pice       = self.P * self.Pq * self.cSnow                             # precip as ice
    self.Ptot       = self.Pice+self.Pliq                                       # total corrected precip
       
                
    ################################################# SNOWPACK ENERGY BALANCE ##
    ############################################################################
        
    ## RADIATION components	
    # albedo
    self.dNoSnow    = ifthenelse(self.Pice > 0.001,1,min(30,(self.dNoSnow + 1)))# calculate days since previous snowstorm and set that snowpack doesn't get any "older" after 30 days. Snowfall lower than 5mm do not reset the snow age, the 1mm here could be calibrated                        
    self.alb        = ifthenelse(self.SWE >0,                                   # estimate snow albedo a a function of days since snowstorm (Wigmosta 1994)
                        ifthenelse(self.dSnow>100,(0.94**self.dNoSnow**0.58)**self.albPow, (0.94**self.dNoSnow**0.58)), # note the ^albPow for old snow is own modification to allow more rapid albedo decay for old snowpacks
                        self.albGround)                                         # for snowfree areas specify constant             

    # incoming SHORTWAVE radiation 
    self.radRss     = self.radRs * (1-self.alb) * (self.tau0 * self.ccov + (1-self.ccov)) # [MJ*d/m2] passing through canopy and transmitted by canopy (Wigmosta et al 1996)
    
    # net LONGWAVE radiation in the snowpack amitted by atmosphere, overstorey and lost by snowpack  (Wigmosta et al 1996)
    self.emAir      = ifthenelse(self.P>self.RDthres,                                  # atmosphere emissivity, different for cloudy and clear days (Walter et al 2005/Campbell and Norman 2004)
                                 (0.72+0.005*self.AT)*(1-0.84)+0.84,            # Walter et al use a value of 5mm to tell apart cloudy and sunny days, we use 1 mm because better fits for Krycklan...
                                  0.72+0.005*self.AT) 
    
    self.Ld         = self.emAir * cSB * (273.15 + self.AT)**4                  # atmospheric longwave radiation
    self.L0         = cSB * (273.15 + self.AT)**4                               # longwave emitted by overstorey, assuming emissivity of unity
    self.Ls         = 0.97 * cSB * (273.15 + self.ST)**4                        # longwave emitted by snow, emissivity 0.97 from (Walter et al 2005)
    self.radLs      = (self.L0 * self.ccov + (self.Ld * (1-self.ccov))-self.Ls) # Net longwave radiation [MJ/d*m2]
    
    # net TOTAL radation on the SNOWPACK
    self.radRns     = self.radRss + self.radLs

    # ADVECTION from PRECIPITATION
    self.Qp         = rooW * Cw * self.AT * (self.Pliq + 0.5 * self.Pice)       # Heat from rain, both liquid and solid (Wigmosta et al. 1994), conversion from mm to m and kj to MJ                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
    
    # SENSIBLE HEAT exchange in the snowpack ...snow temperature from previous timestep is takes as input
    self.ras        = ((ln((zu-ds+zms)/zms) * ln((zu-ds+zhs)/zhs)) / (k**2 * self.WS)) / 86400 # resistance to heat transfer (Walter et al 2005)
    self.Qs         = cair * (self.AT - self.ST) / self.ras                     # sensible heat transfer by turbulent convection [MJ/d*m2]    
        
    # HEAT from convective VAPOUR EXCHANGE (evaporation and condensation)  ...snow temperature from previous timestep is takes as input
    self.pVap	    = 0.6108* exp((17.27*self.AT)/(self.AT+237.3))*10 * self.RH/100	# saturation vapour pressure in a given air temperature (Allen et al 2000), converted to mbbar and scaled to actual with relative humidity data 
    self.rooA       = (self.pVap / ((self.AT+273.15)*Rv)) / 1000                # vapour density of air (Dingman 1993, eq D-7a) converted to [kg/m3]   
    self.rooSA      = exp((16.78*self.ST-116.8) / (self.ST + 273.3)) * (1/((273.15+self.ST)*Rt)) # vapour density at the snow surface (Walter et al 2005)
    self.Ql         = lamv * ((self.rooA - self.rooSA)/self.ras)                # latent heat flux [MJ/d*m2] (Walter et al 2005)
   
    
    ############################################# SNOW INTERCEPTION BY CANOPY ##  
    ############################################################################  
    
    # needs to be done here, because interception sublimation is assumed equal to snowpack sublimation (=latent heat exchange). The overall concept for this from Liston et al (2006)

    self.Isnow_o    = self.Isnow                                                # save these for later
    self.concIst_o  = self.concIst                                                
    
    # amount of interception storage   
    self.rooSP      = 150                                                       # snow density in snowfall assumed constant, could be estimated from min(250,67.9 + 52.3*exp(self.AT/2.6)) but storage density and water balance appear to get tricky if interception storage density varies over time              
    self.Pkgm2      = self.Pice * self.rooSP                                    # incoming snow to units kg / m2
                                                                    
    self.Isnow      = min(self.IsnowMax,                                        # interception storage for snow, cannot exceed the maximum interception storage [kg/m2]
                        self.Isnow_o + 0.7*(self.IsnowMax - self.Isnow_o) * (1 - exp(-self.Pkgm2/self.IsnowMax)))     # calculate the addition to snow int storage from rainfall Liston et al (2006)
    self.Isnow      = ifthenelse(defined(self.Isnow),self.Isnow,0)              # replace the mv aroused by pontential LAI=0 values with 0
    self.Isnow      = ifthenelse(self.LAI==0,0,self.Isnow)
    
    
    self.PIthru     = ifthenelse(self.Pkgm2>0,
                                (self.Pkgm2 - max(0,self.Isnow - self.Isnow_o)), 0)     # the snow that does not accumulate to storage, passes through
    self.PIthru     = ifthenelse(self.LAI==0, self.Pice,self.PIthru/self.rooSP) # adjust for LAI and convert to [m]               
                        
    
    # calculate evaporation from interception storage, assuming similar latent heat exchange as in the snowpack
    self.rooSIA     = exp((16.78*min(0,self.AT)-116.8) / (min(0,self.AT) + 273.3)) * (1/((273.15+min(0,self.AT))*Rt)) # vapour density at the snow surface (Walter et al 2005), take intercepted snow as AT, when AT negative
    self.rasSn      = ((ln((zu-ds+zms)/zms) * ln((zu-ds+zhs)/zhs)) / (k**2 * self.WSin)) / 86400 # resistance to heat transfer (Walter et al 2005)
    self.Qlsnow     = lamv * ((self.rooA - self.rooSIA)/self.rasSn)             # latent heat available for sublimation from interception storage
    self.sIntE      = ifthenelse(self.Isnow>0,                                  # sublimation from interception storage [m], cannot exceed storage or happed without storage
                                min(max(0,-self.Qlsnow/(rooW*lamv)),self.Isnow/self.Isroo),
                                0) 
    
    # calculate unloading of the snowpack due melt etc
    self.Iunload    = ifthenelse(self.Isnow>0,
                                min(self.Isnow,max(0,self.IunDD * self.AT)),   # unloading rate with temperature index kg/m2 (Liston et al 2006), unloaded water goes to ground snowpack, where it can actually be melted 
                                0)                                             
    self.Iunload    = min(max(0,self.Iunload - self.sIntE*self.Isroo),self.Iunload)  # allow ET to take place first
    self.Iunlm      = ifthenelse(self.LAI==0,0,self.Iunload/self.Isroo)         # convert unloading to [m] in order to feed to ground snowpack  

    self.Isnow      = max(0,self.Isnow - self.Iunload - self.sIntE*self.Isroo)  # update the Interception storage to account for sublimation and unloading [kg/m2]
    
    # ISOTOPIC COMPOSITION of interception storage
    
    self.concIst    = ifthenelse(self.Pkgm2 > 0,                                # if there is snowfall
                                ifthenelse(self.Isnow/self.Isroo>0.005,        # for numerical stability ignore ET fractionation for low storages (<5mm)
                                            (self.concIst * self.Isnow + self.Pkgm2 * self.Pconc - (self.concIst - self.Efrac) * self.sIntE*self.Isroo) 
                                            / (self.Isnow + self.Pkgm2 - self.sIntE*self.Isroo+0.0001),
                                            (self.concIst * self.Isnow + self.Pkgm2 * self.Pconc) 
                                            / (self.Isnow + self.Pkgm2+0.0001)),
                                ifthenelse(self.Isnow==0,                      # no snow accumulation 
                                            0,                                 # and no interception storage, concentration of 0
                                            ifthenelse(self.Isnow/self.Isroo>0.005,# in an existing interception storage, allow fractionation if storage above a limit
                                                        (self.concIst * self.Isnow - (self.concIst - self.Efrac) * self.sIntE*self.Isroo)
                                                        / (self.Isnow-self.sIntE*self.Isroo+0.0001),
                                                        self.concIst)))                                   
    # AGE of the interception storage
    if self.doAges:
        self.ageIst_o = self.ageIst
        
        self.ageIst = ifthenelse(self.Pkgm2>0,                                  # if there is snow accumulation              
                                ((self.ageIst + 1) * self.Isnow + max(0,self.Isnow - self.Isnow_o))  # age updated to incoming snow that has an age of one
                                / (self.Isnow + max(0,self.Isnow - self.Isnow_o)+0.0001),
                                ifthenelse(self.Isnow==0,                     # no snow accumulation and no snowpack
                                            1,                                 # age set to 0
                                            self.ageIst+1))                    # otherwise interception storage ages one day
        self.ageIst = ifthenelse(defined(self.ageIst),self.ageIst,0)            # set missing values to 0 for numerical stability        
    

    ##################### back to ground snowpack...                                            
                                                                                                                                        
    # sum of ENERGY INPUT/OUTPUT which will results in melting/refreezing and heating/cooling the snowpack
    self.Esum       = self.radRns + self.Qp + self.Qs + self.Ql                 # positive fluxes add energy to the snowpack and negative remove energy from snowpack     
        
    # the portion of the sum energy diverted to snowmelt 
    self.Qm         = ifthenelse(self.SWE==0,
                                 0,                                             # if there is no snow, there can be no melt or refreezing          
                                 ifthenelse(self.Esum > 0,                      # determine if the snowpack is melting or freezing
                                            ifthenelse(self.ST<0,               # MELTING: is there cold content to melt first?
                                                        max(0,self.Esum - (0 - self.ST)*(1/(Ci*(self.SWE+0.00001)*rooW))), # what is available for melt after heating the snowpack. If more cold content than heat, nothing left for melt (add 0.00001 for numerical stability)
                                                        self.Esum),             # snowpack isothermal, all energy is diverted to snowmelt !
                                            ifthenelse(self.Wliq==0,            # FREEZING: is there liquid water to refreeze?
                                                        0,                      # no liquid water, no energy is wasted on freezing water
                                                        max(-rooW*lamf*self.Wliq,self.Esum)))) # either there is energy to refreeze everything or just a fraction is refrozen  

    self.Qc         = self.Esum - self.Qm                                       # snowpack cold content change                                                                                 
        
    ################################################### SNOWPACK MASS BALANCE ##   
    ############################################################################   
    
    # mass balance formulation concept from (Wigmosta 1994) except that sublimation/deposition takes place from ice phase
    
    self.Wliq_o     = self.Wliq                                                 # store values from previous timestep
    self.Wice_o     = self.Wice
    self.SWE_o      = self.SWE

    self.Wice       = max(0,self.Wice+self.PIthru+self.Iunlm)                   # update the snow ice content after adding water from throughfall or unloading                
    self.M          = ifthenelse(self.Qm<0,                                     # decide if there is refreezing or melt taking place based on the Qm 
                                 max(self.Qm/(rooW*lamf),-self.Wliq),           # rate of refreezing [m/d], limited by the liquid water storage available for it (minus stands for negavite energy
                                 min(self.Qm/(rooW*lamf),self.Wice))            # rate of melt, limited by the availability of ice in snowpack
    
    self.Wice       = max(0,self.Wice-self.M)                                   # update the snow ice content after melt/refreeze                                                                  
   
    # calculate and limit evaporation/deposition
    self.Eice       = ifthenelse(self.Wice>0,max(self.Ql/(rooW*lamv),-self.Wice),0) # sublimation/deposition to the solid ice phase, assuming there is ice left !! no evaporation from liquid phase...!                                       
           
    # update the mass balance of the ice phase
    self.Wice       = max(0,self.Wice + self.Eice)                              # update snow ice content after sublimation, cannot go negative
   
    # mass balance for liquid phase
    self.Wliq       = max(0,min(self.retcap*self.Wice,                          # liquid water in the snowpack [m], cannot exceed water retention capacity or go negative
                        self.Wliq + self.Pliq + self.M))                        # water is added via rain and added/removed via melting/freezing
    self.Wliqout    = max(0,(self.Wliq_o + self.Pliq + self.M)-self.retcap*self.Wice) # water flow out of the snowpack, a certain volume retained    
    self.Pmelt      = self.Wliqout * 1000                                       # COMFORM to starr variable naming and convert to [mm/d]
    
    # total water content in snowpack       
    self.SWE        = self.Wice + self.Wliq                                     # total snow water equivalent a sum of liquid and ice fraction
    self.deltaSWE   = self.SWE-self.SWE_o                                       # save the change in SWE 
    
    # update snow temperature according to energy and mass balance                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
    self.ST_o       = self.ST
    # constrain snow temperature
    self.ST         = ifthenelse(self.SWE<0.05,self.AT,                         # with shallow snow depths snow temp equals air temp
                                    self.ST_o + (self.Qc)/(Ci*self.SWE*rooW))   # excess energy from melting 
    self.ST         = ifthenelse(self.Qm>0,0,self.ST)                           # during melt ST = 0 
    self.ST         = ifthenelse(self.ST<-4,                                    # cannot go too cold, miniimise to AT. This breaks the energy conservation so this could be improved. Perhaps introduce soil energy store where excess energy/cold content is diverted  
                        max(self.ST,self.AT),                                   # if snow temperature tries to go below -4, set lower limit to air temperature
                        self.ST)                                                # for moderately cold snow (>-4) leave at simulated ST
                                                                                         #                                     
    self.ST         = min(self.ST,0)                                            # cannot be positive                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
    # set the internal energy flux variables to 0 as long as the SWE is 0
    self.Ql         = ifthenelse(self.SWE==0,0,self.Ql)
    self.Qs         = ifthenelse(self.SWE==0,0,self.Qs)
    self.Qc         = ifthenelse(self.SWE==0,0,self.Qc)
    self.Qm         = ifthenelse(self.SWE==0,0,self.Qm)
    
           		
    ########################################## SNOWPACK ISOTOPE CONCENTRATION ##
    ############################################################################

    # snow concentration is cells with snow changes according to the new snow fall and is isotope concentration, and fractionation via evaporation
    # snowpack isotope concentration is thus a weighted average from snow precip inputs assuming INSTANTANEOUS MIXING
    self.concsn_o   = self.concsn
                                                                      
    self.concsn     = ifthenelse(pcror(self.Pice>0,self.Iunlm>0.001),           # if there is snow accumulation                     
                                ifthenelse(self.Wice>0.005,                     # For numerical stability for snowpacks less that 5 mm evaporation excluded
                                            (self.concsn * (self.Wice + self.Wliq) + self.PIthru * self.Pconc + self.Pliq * self.Pconc + self.Iunlm*self.concIst_o + (self.concsn - self.Efrac) * self.Eice)     
                                            / (self.Wice + self.Wliq + self.PIthru + self.Pliq + self.Iunlm + self.Eice),
                                            (self.concsn * (self.Wice + self.Wliq) + self.PIthru * self.Pconc + self.Pliq * self.Pconc + self.Iunlm*self.concIst_o) 
                                            / (self.Wice + self.Wliq + self.PIthru + self.Pliq + self.Iunlm)),
                                
                                ifthenelse(self.Wice==0,                        # no snow accumulation and no snowpack
                                            0,                                  # concentration of 0
                                            ifthenelse(self.Wice>0.005,         # For numerical stability for snowpacks less that 5 mm evaporation excluded
                                                        ((self.concsn * (self.Wice + self.Wliq) + (self.concsn - self.Efrac) * self.Eice) # if prior snowpack but no snowfall, sublimation changes concentration
                                                        / (self.Wice + self.Wliq + self.Eice)),
                                                        self.concsn)))                                      
    
    ################################################################ SNOW AGE ##
    ############################################################################
    if self.doAges:
        self.agesnow = ifthenelse(pcror(self.Pice>0,self.Iunlm>0.001),          # if there is snow accumulation
                                  ((self.agesnow + 1) * (self.Wice + self.Wliq) + self.PIthru + self.Pliq + self.Iunlm*self.ageIst_o)            
                                   / (self.Wice + self.Wliq + self.PIthru + self.Pliq + self.Iunlm),
                                  ifthenelse(self.Wice==0,                      # no snow accumulation and no snowpack
                                             1,                                 # age of 1
                                             self.agesnow+1))                   # if prior snowpack, snow ages one day
                          
    # calculate a "depletion coefficient" to drain the most depleted water out of snowpack first
    self.dSnow      = ifthenelse(self.Wice>0.010 ,self.dSnow+1,0)               # calculate how many days a given cell has had a snow cover of over 1 cm
    self.dMelt      = ifthenelse(self.dSnow>1,                                  # Days when snowpack has been releasing water. First: is there snow?
                                 ifthenelse(self.Wliqout>0.002,                 # yes, is there marked outflow (>2 mm) from snowpack?
                                             self.dMelt+1,                      # yes, one more day of water yield during a given snow accumulation period
                                             self.dMelt),                       # no, number of melt days remains the same
                                 0)                                             # if no prior snow cover, no melt days either                    
    
    # coefficient of how much the water outflow is depleted in relation to snowpack
    self.deplCoef   = ifthenelse(self.dMelt>0, -self.deplOffset/self.dMelt, 0)   # function of number of snowmelt days during snow season, first melts drain the most depleted water    
    
    self.concmelt   = ifthenelse((self.SWE_o - max(0,-self.deltaSWE))>0.01,     # is there a snowpack of over 1 cm, assume no preferential fractination for lower snowpacks
                                 self.concsn+self.deplCoef,                     # yes, determine melt rate to be more negative than snowpack
                                 min(self.concsn,self.concsn_o))                # no, any melt water concentration equals concentration of current or previous day snowpack
    
    self.concmelt   = ifthenelse(self.Wliqout>0,                                # is outflow from snowpack?
                                 self.concmelt,
                                 0)                                             # otherwise zero
            
                            
    self.concsn     = ifthenelse(pcrand((self.SWE_o - max(0,-self.deltaSWE))>0.01,self.Wliqout>0),     # is there an exisiting snowpack of over 1 cm and a snowmelt event??
                                 (self.SWE_o*self.concsn - max(0,-self.deltaSWE)*self.concmelt)/(self.SWE_o - max(0,-self.deltaSWE)), # yes, the concentration updated according to meltwater, if any
                                 self.concsn)                                   # no, stick to previous values
                                                                                                                    
               	                                                                                                                                      	                                                                                                                                                                                                                                                                	               	               	                                                                                                                                      	                                                                                                                                                                                                                                                                                                                                                                                          	               	               	                                                                                                                                      	                                                                                                                                                                                                                                                                                                                                                                                           	               	               	                                                                                                                                      	                                                                                                                                                                                                                                                                                                                                                                                           	               	               	                                                                                                                                      	                                                                                                                                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
    # WATER BALANCE CHECK
    ############################################################################
    # fluxes in
    waterBalanceCheck.waterBalanceCheck([self.Ptot*1000],\
    # fluxes out
    [self.Wliqout*1000, -self.Eice*1000,self.sIntE*1000],\
    # storage for previous timestep
    [self.Wice_o*1000, self.Wliq_o*1000,(self.Isnow_o/self.Isroo)*1000],\
    # storage for new timestep
    [self.Wice*1000, self.Wliq*1000,(self.Isnow/self.Isroo)*1000],\
    'Snow',True,\
    self.currentTimeStep(),threshold=1e-3,outputWB=False)
    ############################################################################
    
  






