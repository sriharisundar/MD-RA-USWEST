#Michael Craig
#October 12, 2016
#Remove hydro (normal + pumped storage) units from fleet and subtract their monthly average generation
#from demand profile.

import copy, operator, os, numpy as np, pandas as pd, calendar
from CombinePlants import combinePlantsByRegion

#Inputs: gen fleet, demand and net demand dataframes
#Outputs: demand & net demand minus hydropower generation
def processHydro(genFleet, demand, netDemand, weatherYears, removeHydro, climateChange):
    weatherYears = weatherYears[0] if not climateChange else 2018
    if climateChange: print('In ProcessHydro, using 2018 weather year!')
    if 'Hydro' in genFleet['FuelType'].unique():
        #Get EIA data
        if weatherYears == 2012: netGenCols,rftCol = ['Netgen_'+calendar.month_name[i][:3] for i in range(1,13)],'Reported Fuel Type Code'
        elif weatherYears == 2019: netGenCols,rftCol = ['Netgen\n'+calendar.month_name[i] for i in range(1,13)],'Reported\nFuel Type Code'
        else: netGenCols,rftCol = ['Netgen\r\n'+calendar.month_name[i] for i in range(1,13)],'Reported\r\nFuel Type Code'

        yrGen = import923(weatherYears,netGenCols,rftCol)

        #Convert EIA data from monthly to hourly or daily
        # print('SKIPPIN HYDROPOWER') #for CESM runs, comment out next line 
        genPerStepAll = convert923ToSubMonthlyGen(yrGen,genFleet,netDemand,netGenCols)

        #If removing hydro gen from demand, remove hydro from fleet; otherwise aggregate hydro by zone
        if removeHydro:
            # print('SKIPPIN HYDROPOWER') #for CESM runs, uncomment next line
            # genPerStepAll = pd.DataFrame(0,index=netDemand.index,columns=netDemand.columns)
            demand -= genPerStepAll
            genFleet = genFleet.loc[genFleet['FuelType'] != 'Hydro']
        else: 
            for r in genFleet['region'].unique(): genFleet = combinePlantsByRegion(genFleet,'FuelType','Hydro',r)
    else: #create dummy copy of hourGenAll for later functions
        genPerStepAll = demand.copy()
        genPerStepAll *= 0

    return genFleet, genPerStepAll, demand

def convert923ToSubMonthlyGen(yrGen,genFleet,netDemand,netGenCols):
    #Determine time step (hourly or daily)
    hrsPerStep = 1 #if run total daily demand, use following code to scale hourly to daily capacity: # hrsPerStep = 24 if (netDemand.index[1]-netDemand.index[0]).days==1 else 1. Note that CC data is average, so don't scale now.
    #Pull out hydro units
    hydroUnits = genFleet.loc[genFleet['FuelType'] == 'Hydro']
    genPerStepAll = pd.DataFrame(index=netDemand.index,columns=netDemand.columns)
    for region in hydroUnits['region'].unique():
        #Aggregate hydro generators to plants, then get regional total capacity
        hydroRegion = hydroUnits.loc[hydroUnits['region']==region]
        initCap = hydroRegion['Capacity (MW)'].sum()
        capac = hydroRegion.groupby('ORIS Plant Code')['Capacity (MW)'].apply(lambda x: np.sum(x.astype(float))).reset_index()
        capac.index = capac['ORIS Plant Code']
        totalCapac = capac['Capacity (MW)'].sum()
        assert((initCap-totalCapac)<.01*initCap)
        
        #Match EIA data to hydro units
        genRegion = capac.merge(yrGen,how='left',left_index=True,right_index=True)

        #Get hourly generation based on regional net demand
        for mnth in range(1,13):
            #Get total month generation
            monthGen = genRegion[netGenCols[mnth-1]].sum()

            #Calculate hourly weights to determine hourly generation
            monthNetDemand = netDemand.loc[netDemand.index.month==mnth,region]
            monthNetDemand.loc[monthNetDemand<0] = 0 #sometimes negative, which messes w/ weights; so set these to zero

            #Calculate normalized weights so sum = 1, avoiding special case of all hours having negative net demand
            if int(monthNetDemand.max()) != 0: 
                wt = monthNetDemand/(monthNetDemand.max())
                wt = wt/wt.sum() 
            else: 
                wt = pd.Series(1/monthNetDemand.shape[0],index=monthNetDemand.index,name=region)

            #Estimate hourly generation using weights
            genPerStep = monthGen * wt * hrsPerStep
            assert((genPerStep.sum()-monthGen)<.01*monthGen)

            #If generation exceeds capacity in hours, reallocate that generation surplus to other hours
            hoursAboveCap,hoursBelowCap = genPerStep.loc[genPerStep>=(totalCapac*hrsPerStep)],genPerStep.loc[genPerStep<(totalCapac*hrsPerStep)]
            surplus = (hoursAboveCap - totalCapac).sum()
            while surplus > 0 and hoursBelowCap.shape[0] > 0:
                #Evenly spread surplus generation across hours below capacity
                genPerStep[hoursBelowCap.index] += surplus/hoursBelowCap.shape[0]
                #Cap generation at capacity
                genPerStep[hoursAboveCap.index] = totalCapac
                #Reassign hours and recalculate surplus
                hoursAboveCap,hoursBelowCap = genPerStep.loc[genPerStep>=(totalCapac*hrsPerStep)],genPerStep.loc[genPerStep<(totalCapac*hrsPerStep)]
                surplus = (hoursAboveCap - totalCapac).sum()

            #Might exit while loop w/ remaining surplus - if so, all hours are full, so truncate
            if surplus > 0: genPerStep.loc[genPerStep>=(totalCapac*hrsPerStep)] = totalCapac*hrsPerStep

            #Place gen for month & region into df w/ all months & regions
            genPerStepAll.loc[genPerStep.index,region] = genPerStep
        
            #Make sure all generation was allocated (or dropped as surplus)
            assert((monthGen - (genPerStep.sum()+surplus))<.0001*monthGen)
    return genPerStepAll

def import923(metYear,netGenCols,rftCol):
    #Import, skipping empty top rows
    yrGen = pd.read_csv(os.path.join('Data','EIA923','gen' + str(metYear) + '.csv'),skiprows=5,header=0,thousands=',')
    yrGen = yrGen[['Plant Id',rftCol]+netGenCols]

    #Data very rarely has a missing value - replace with a zero for simplicity
    yrGen = yrGen.replace('.',0)

    #Slim down to hydro facilities
    yrGen = yrGen.loc[yrGen[rftCol] == 'WAT']
    yrGen.drop([rftCol],axis=1,inplace=True)

    #Get rid of , text and convert to float (thousands=',' doesn't work above)
    for lbl in netGenCols:
        yrGen[lbl] = yrGen[lbl].astype(str).str.replace(',','')
        yrGen[lbl] = yrGen[lbl].astype(float)

    #Aggregate unit to plant level
    yrGen = yrGen.groupby('Plant Id').apply(lambda x: np.sum(x.astype(float))).reset_index(drop=True)

    #Reindex
    yrGen['Plant Id'] = yrGen['Plant Id'].astype(int)
    yrGen.index = yrGen['Plant Id']
    yrGen.drop(['Plant Id'],axis=1,inplace=True)

    return yrGen


