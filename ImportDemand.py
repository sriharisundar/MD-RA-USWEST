import operator, os, sys, pandas as pd

def importDemand(weatherYears,demandScen,currYear,transRegions,cesmMember):
    if weatherYears == [2012]: demand = importEFSDemand(currYear,transRegions,demandScen) #2012
    elif cesmMember != None: demand = importCCDemand(weatherYears,cesmMember) #future years so use CC demand
    else: demand = importERA5Demand(weatherYears) #historic non 2012 data so use reanalysis
    return demand

def importEFSDemand(currYear,transRegions,demandScen):
    #Initialize df
    if currYear > 2050: currYear = 2050
    dates = pd.date_range('1/1/'+str(currYear)+' 0:00','12/31/' + str(currYear) + ' 23:00',freq='H')
    dates = dates[~((dates.month == 2) & (dates.day == 29))] #ditch leap day
    demand = pd.DataFrame(index=dates)
    
    #Read EFS data
    filename = 'EP' + demandScen + '_FlexNONEload_hourly.csv'
    rawDemand = pd.read_csv(os.path.join('Data','REEDS', filename), delimiter=',',header=0)
    rawDemand = rawDemand.loc[rawDemand['year']==currYear]

    #Iterate through dict of zone:p regions (from REEDS) & aggregate demand for p-regions
    for zone,pRegions in transRegions.items():
        for p in pRegions:
            pDemand = rawDemand[p]
            if zone in demand.columns:
                demand[zone] += pDemand.values
            else:
                demand[zone] = pDemand.values
    return demand

def importCCDemand(weatherYears,cesmMembers):
    demands = list()
    #Import each ensemble member
    for c in cesmMembers:
        filename = 'demand_subregions_' + c + '.csv'
        demand = pd.read_csv(os.path.join('Data','CESM', filename), delimiter=',',header=0,index_col=0,parse_dates=True)

        #Slim down to current years
        demand = demand.loc[(demand.index.year >= weatherYears[0]) & (demand.index.year <= weatherYears[-1])]

        #If multiple ensemble members, create MultiIndex with 2 levels (for ensemble member & region)        
        if len(cesmMembers)>1: 
            idx = pd.MultiIndex.from_product([[c], demand.columns], names=['ensembleMember', 'region'])
            demand.columns = idx

        demands.append(demand)

    #Combine into 1 array
    demand = pd.concat(demands,axis=1)

    return demand

def importERA5Demand(weatherYears):
    demand = pd.read_csv(os.path.join('Data','ERA5','hourlydemand_subregions_historic.csv'), delimiter=',',header=0,index_col=0,parse_dates=True)
    
    #Slim down to current years
    demand = demand.loc[(demand.index.year >= weatherYears[0]) & (demand.index.year <= weatherYears[-1])]

    return demand