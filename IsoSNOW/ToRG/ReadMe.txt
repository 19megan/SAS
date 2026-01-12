IsoSnow_v1.0: Code for simulating spatially distributed snow accumulation and melt and related isotopic ratios/fractionations 
Pertti Ala-aho, University of Aberdeen 5 July 2017

This distribution (IsoSnow_v1.0.zip) comes with input files that reproduce a synthetic simulation example in: Ala-aho, P., D. Tetzlaff, J. P. McNamara, H. Laudon, P. Kormos, and C. Soulsby (2017),
 Modelling the isotopic evolution of snowpack and snowmelt: testing a spatially distributed parsimonious approach, Water Resour. Res, doi: 10.1002/2017WR020650. 
For the problem description and model equations please see the publication (CC-BY), and the related supplementary material.

The IsoSnow_v1.0 is built with Python 2.7.9 64-bit AND PCraster 4.0.2_x86-64, which need to be installed to run the model. 
Available from: http://pcraster.geo.uu.nl/; https://www.python.org/downloads/ 
Spatial model output can be visualised with Aguila: http://pcraster.geo.uu.nl/projects/developments/aguila/

The distribution includes: 
•IsoSnow_setup_2017_07_04.py      -> script to setup the model run (using PCraster) 
•IsoSnow_initialPara_2017_07_04.py -> script to define parameter values and initial conditions 
•IsoSnow_reporting_2017_07_04.py  -> script to specify what model output is written to disk 
•IsoSnow_model_2016_12_22_PA.py  -> the actual model algoritms 
•waterBalanceCheck.py    -> script to check simulated water balance
•/Input     -> folder containing the input data 
•/1      -> folder where model output is written

For questions and comment please contact: pertti.ala-aho (at) abdc.ac.uk / pertti.ala-aho (at) oulu.fi
