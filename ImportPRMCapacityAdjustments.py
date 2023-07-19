import os, copy, datetime, pandas as pd, geopandas as gpd, datetime as dt, numpy as np, xarray as xr
from os import path

def importPRMCapacityAdjustments(genFleetForCE, newTechsCE, demand, prmBasis, interconn, nonCCReanalysis, weatherYears, compressedGens, cesmMembers, defaultFOR=0.05):
    if cesmMembers != None and len(cesmMembers)>1: sys.exit('PRM capacity adjustments not built for multiple CESM members!')

    #Wind & solar not eligible for PRM if based on net demand because already accounted for in calculating planning reserve!
    prmEligWindSolar = 1 if prmBasis == 'demand' else 0

    #Wind and solar FORs equal default FOR
    windFOR,solarFOR = defaultFOR,defaultFOR

    #Initialize FORs
    fors = pd.DataFrame(defaultFOR,index=demand.index,columns=genFleetForCE['GAMS Symbol'])
    forsTechs = pd.DataFrame(defaultFOR,index=demand.index,columns=newTechsCE['GAMS Symbol'])

    if (nonCCReanalysis or cesmMembers != None) and interconn == 'WECC':
        #Get temperatures
        if cesmMembers != None: temps,tVar = importCESMTemperatures(weatherYears,cesmMembers[0]),'TREFHT'
        elif nonCCReanalysis: temps,tVar = importERA5Temperatures(weatherYears),'tas'

        #Calculate FORs for existing units
        fors = calculateFORs(temps,genFleetForCE,fors,compressedGens,tVar)

        #Calculate FORs for new units
        forsTechs = calculateNewTechFORs(forsTechs,fors,genFleetForCE,newTechsCE)
        print('Adding FORs for new techs only works with 1 CE run. In 2nd run, new techs do not have a lat/long, so cannot get location-specific TDFORs. To fix, add location or use original FORs.')

    return fors,windFOR,solarFOR,forsTechs,prmEligWindSolar

def importERA5Temperatures(weatherYears):
    #Load FORs file that contains temperatures (in C) and slice down to weather years
    temps = xr.open_dataset(os.path.join('Data','ERA5','wecc_FOR_ERA5_hourly_PST.nc'))
    temps = temps.sel(time=slice(str(weatherYears[0])+"-01-01", str(weatherYears[-1])+"-12-31"))
    return temps

def importCESMTemperatures(weatherYears,cesmMember):
    temps = xr.open_dataset(os.path.join('Data','CESM','wecc_derate_fields_' + cesmMember + '.nc'))
    temps = temps.sel(time=slice(str(weatherYears[0])+"-01-01", str(weatherYears[-1])+"-12-31")) #T in K    
    return temps

def calculateFORs(temps,fleet,fors,compressedGens,tVar):    
    #Import regression-based relationships for TDFORs. Taken from Murphy et al. (https://www.sciencedirect.com/science/article/pii/S0306261919311870)
    forsRegression = pd.read_excel(os.path.join('Data','TDFORRelationships.xlsx'),index_col=0,header=0)
    forsRegression /= 100 #want FORs as fractions, but given as percents in spreadsheet

    #Create dict mapping from fleet plant types to plant type codes in FOR sheet
    forPTMatching = {'Coal Steam':'ST','Biomass':'ST','O/G Steam':'ST','Geothermal':'ST',
        'Nuclear':'NU',
        'Landfill Gas':'CT','Combustion Turbine':'CT',
        'Combined Cycle':'CC',
        'Fossil Waste':'Other','Non-Fossil Waste':'Other','Fuel Cell':'Other','Municipal Solid Waste':'Other','Solar PV':'Other','Onshore Wind':'Other',
        'Hydro':'HD'}

    #For each generator w/ a TDFOR, replace temperatures w/ TDFOR
    for c in fors.columns:      
        genRow = fleet.loc[fleet['GAMS Symbol']==c].squeeze()

        if 'COMBINED' not in genRow['GAMS Symbol']: #If not a combined unit, get FOR that unit
            fors[c] = (calculateFORsForGen(temps,genRow,forsRegression,forPTMatching,tVar)).values #use values due to datetime format mismatch error
        else: #If a combined unit, get FORs for constituent units and average
            constituentUnits = compressedGens.loc[compressedGens['UnitCompressedInto']==c]
            constituentFORs = [calculateFORsForGen(temps,constituentUnits.iloc[i],forsRegression,forPTMatching,tVar) for i in range(constituentUnits.shape[0])]
            fors[c] = (pd.concat(constituentFORs,axis=1).mean(axis=1)).values #use values due to datetime format mismatch error
    return fors 

#Calculate FORs for given generator using temps
def calculateFORsForGen(temps,genRow,forsRegression,forPTMatching,tVar):
    pt,lat,lon = genRow['PlantType'],genRow['Latitude'],genRow['Longitude']

    #Get closest coords w/ non-NAN temperatures. If closest temperature has NANs, step towards center of region and try again
    latCenter,lonCenter = temps['lat'][int(temps.dims['lat']/2)].values,temps['lon'][int(temps.dims['lon']/2)].values
    latStep,lonStep = abs(temps['lat'][0].values-temps['lat'][1].values),abs(temps['lon'][0].values-temps['lon'][1].values)

    while np.isnan(temps.isel(time=-1).sel({'lat':lat,'lon':lon},method='nearest').compute()[tVar]):
        lat,lon = lat + ((latCenter-lat)/abs(latCenter-lat))*latStep,lon + ((lonCenter-lon)/abs(lonCenter-lon))*lonStep

    #Get temperatures & dt index from coordinate, then store in Series
    genTemps = temps.sel({'lat':lat,'lon':lon},method='nearest')
    genTemps,times = np.array(genTemps.variables[tVar][:]),genTemps['time'] #in time,lat,lon; swap so lat,lon,time
    genTemps = pd.Series(genTemps,index=times)

    #If CESM data, T given in K; convert to C and round to 5 degree intervals
    if np.nanmax(genTemps.values)>150: #screens for K values
        genTemps -= 273.15
        genTemps = 5*((genTemps/5).round())
        genTemps[genTemps>35],genTemps[genTemps<-15] = 35,-15

    #Get FOR type
    forType = forPTMatching[pt] if pt in forPTMatching else 'Other'

    #Replace temperatures w/ FORs        
    tToFOR = forsRegression[forType].to_dict()
    genFORs = genTemps.replace(tToFOR)
    return genFORs 

#Calculate FORs for new techs by taking average FORs of all existing generators
#of same plant type and region.
def calculateNewTechFORs(forsTechs,fors,genFleetForCE,newTechsCE):
    for tech in forsTechs.columns:
        #Get tech's plant type and region
        pt = newTechsCE.loc[newTechsCE['GAMS Symbol'] == tech,'PlantType'].values[0]
        region = newTechsCE.loc[newTechsCE['GAMS Symbol'] == tech,'region'].values[0]

        #Skip wind & solar techs, since don't have TDFORs
        if 'wind' not in tech and 'solar' not in tech:

            #Get existing generators of same pt and region
            gensOfPT = genFleetForCE.loc[genFleetForCE['PlantType']==pt]
            gensOfPTAndRegion = gensOfPT.loc[gensOfPT['region']==region]

            #If existing generators of same pt and region, get FORs and average them
            if gensOfPTAndRegion.shape[0]>0:
                gensFORs = fors[gensOfPTAndRegion['GAMS Symbol']]
                forsTechs[tech] = gensFORs.mean(axis=1)
    return forsTechs

