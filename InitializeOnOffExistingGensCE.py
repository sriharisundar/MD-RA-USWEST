#Michael Craig
#Dec 2, 2016

#Initialize on/off commitment state of existing units for CE model for first
#hour in each time block (seasons + demand day + ramp day).

import copy, pandas as pd
from operator import *
from GAMSAuxFuncs import *

#For first hour of each time block, gets gen stack, sorts gens by op cost, then
#meets demand and assigns dispatched gens to "on".
#Inputs: gen fleet, hours + hourly wind & solar gen + demand for CE model (1d list,
#slimmed to CE hours), dict mapping time block to hours in that time block
#Outputs: dict mapping block to dict mapping gen symbol to on/off in first hour of dict.
def initializeOnOffExistingGens(genFleetForCE,hoursForCE,netDemandOrig):
    #If have multi-index of ensemble members, just use a single member for net demand
    netDemand = netDemandOrig.copy()
    if type(netDemand.columns)==pd.core.indexes.multi.MultiIndex:
        firstEnsembleMember = netDemand.columns.get_level_values('ensembleMember')[0]
        netDemand = netDemand[firstEnsembleMember]

    #Create gen stack from cheapeast to most expensive (minus RE & storage)
    sortedSymbols,sortedCapacs,reGenSymbols = getGenStack(genFleetForCE)
    
    #For each block, find intersection b/wn block's initial net demand & gen stack to initialize on/off 
    blockToOnOff = dict()
    for bl in hoursForCE.unique():
        firstHr = hoursForCE.loc[hoursForCE == bl].index[0]
        hrNetDemand = netDemand.loc[firstHr].values[0]
        genToOnOff = getGenOnOffForHour(sortedSymbols,sortedCapacs,reGenSymbols,hrNetDemand)
        blockToOnOff[bl] = genToOnOff
    return blockToOnOff

#Get gen stack for fleet, returning sorted gen symbols, capacities, and RE gen symbols.
#Sorted in order of least -> most costly.
def getGenStack(genFleetForCE):
    noRE = genFleetForCE.loc[~genFleetForCE['FuelType'].isin(['Wind','Solar'])]
    re = genFleetForCE.drop(noRE.index)
    noRE = noRE.loc[~noRE['FuelType'].isin(['Energy Storage'])]
    noRE.sort_values(by="OpCost($/MWh)",inplace=True)
    noRESymbs = (noRE['ORIS Plant Code'].astype(str) + "+" + noRE['Unit ID'].astype(str))
    return noRE['GAMS Symbol'],noRE['Capacity (MW)'],re['GAMS Symbol']

#For given net demand value in hour, determine which gens should be on by 
#dispatching them using gen stack.
def getGenOnOffForHour(sortedSymbols,sortedCapacs,reGenSymbols,hrNetDemand):
    cumulativeCap = sortedCapacs.cumsum()
    online = cumulativeCap<hrNetDemand
    if False in online.values: online[online[online==False].index[0]] = True #set marginal generator as online
    onGens = sortedSymbols[online[online==True].index]
    genToOnOff = {**pd.Series(1,index=onGens).to_dict(),**pd.Series(1,index=reGenSymbols).to_dict()}
    if False in online.values: 
        offGens = sortedSymbols[online[online==False].index]
        genToOnOff = {**genToOnOff,**pd.Series(0,index=offGens).to_dict()}
    return genToOnOff
    