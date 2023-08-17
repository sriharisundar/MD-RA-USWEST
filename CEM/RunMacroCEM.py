#Shell script for running macro CEM
#Order of inputs: interconn, co2cap

import sys,os
from MacroCEM import macroCEM

#Set working directory to location of this script
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

#Process inputs and call master function
inputData = sys.argv[1:] #exclude 1st item (script name)
interconn = inputData[0] #ERCOT, WECC, EI
co2Cap = int(inputData[1]) #integer giving % of CO2 emissions in final year relative to first year
wsGenFracOfDemand = int(inputData[2]) #integer giving % of demand that must be supplied by wind and solar generation
prm = int(inputData[3]) #integer giving % planning reserve margin
yr = int(inputData[4])
ret = True if int(inputData[5])==1 else False
cesmMember = [inputData[6]] if len(inputData)>6 else None #CESM Large Ensemble member (as list)

#Other inputs
climateChange = (cesmMember != None)

macroCEM(interconn,co2Cap,wsGenFracOfDemand,prm,yr,ret,cesmMember,climateChange)
