

import subprocess
import datetime
import random
import os
import gc
import re
import math
import sys
import types

#import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import pcraster as pcr

import logging
logger = logging.getLogger(__name__)

def getMapTotal(mapFile):
    ''' outputs the sum of all values in a map file '''

    total, valid= pcr.cellvalue(pcr.maptotal(mapFile),1)
    return total

    
def waterBalanceCheck(fluxesIn,fluxesOut,preStorages,endStorages,processName,PrintOnlyErrors,dateStr,threshold=1e-5,outputWB=False,landmask=None):
    """ Returns the water balance for a list of input, output, and storage map files  """
    # modified by Edwin (22 Apr 2013)

    inMap   = pcr.spatial(pcr.scalar(0.0))
    outMap  = pcr.spatial(pcr.scalar(0.0))
    dsMap   = pcr.spatial(pcr.scalar(0.0))
    
    for fluxIn in fluxesIn:
        inMap   += fluxIn
    for fluxOut in fluxesOut:
        outMap  += fluxOut
    for preStorage in preStorages:
        dsMap   += preStorage
    for endStorage in endStorages:
        dsMap   -= endStorage
    
    errormap = pcr.spatial(inMap + dsMap- outMap)
    a,b,c = getMinMaxMean(errormap)
    
      
    if outputWB==True:
                     
        dirs = os.walk('.').next()[1]                                               # find all directories
        numdirs = [x for x in dirs if x.isdigit()]                             # get the number of directory where presently working
        
        # report MAPS of water balance error
        if dateStr<1000:
            pcr.report(errormap, str(numdirs[-1])+"/WBerr000."+str(dateStr).zfill(3))
        elif dateStr>=1000 and dateStr<2000:
            pcr.report(errormap, str(numdirs[-1])+"/WBerr001."+str(dateStr-1000).zfill(3))
        elif dateStr>=2000 and dateStr<3000:
            pcr.report(errormap, str(numdirs[-1])+"/WBerr002."+str(dateStr-2000).zfill(3))
        else:
            pcr.report(errormap, str(numdirs[-1])+"/WBerr003."+str(dateStr-2000).zfill(3))
        
        # report TIMESERIES of water balance error  
        parameter_file = open(str(numdirs[-1])+"/WBE_for_%s.txt"%processName, "a")  # append the water balance error for each timestep        
        parameter_values = float(a),float(b),float(c) 
        parameter_file.write(str(parameter_values)[1:-1]+'\n')                                 # write the actual output
        parameter_file.close()                                                      # close file


    if abs(a) > threshold or abs(b) > threshold:
        if PrintOnlyErrors: 
            
            msg  = "\n"
            msg += "\n"
            msg  = "\n"
            msg += "\n"
            msg += "##############################################################################################################################################\n"
            msg += "WARNING !!!!!!!! Water Balance Error %s Min %f Max %f Mean %f" %(processName,a,b,c)
            msg += "\n"
            msg += "##############################################################################################################################################\n"
            msg += "\n"
            msg += "\n"
            msg += "\n"
            
            print(msg)

    else:

        print "water balance OK"



def getMinMaxMean(mapFile,ignoreEmptyMap=False):
    mn = pcr.cellvalue(pcr.mapminimum(mapFile),1)[0]
    mx = pcr.cellvalue(pcr.mapmaximum(mapFile),1)[0]
    nrValues = pcr.cellvalue(pcr.maptotal(pcr.scalar(pcr.defined(mapFile))), 1 ) [0] #/ getNumNonMissingValues(mapFile)
    if nrValues == 0.0 and ignoreEmptyMap: 
        return 0.0,0.0,0.0
    else:
        return mn,mx,(getMapTotal(mapFile) / nrValues)


#~ example: 
            #~ vos.waterBalanceCheck([netLqWaterToSoil,\
                                   #~ self.irrGrossDemand,\
                                   #~ self.satExcess],\
                                  #~ [self.directRunoff,
                                   #~ self.openWaterEvap,
                                   #~ self.infiltration],\
                                  #~ [  preTopWaterLayer],\
                                  #~ [self.topWaterLayer],\
                                       #~ 'topWaterLayer',True,\
                                   #~ currTimeStep.fulldate,threshold=1e-4)
            
