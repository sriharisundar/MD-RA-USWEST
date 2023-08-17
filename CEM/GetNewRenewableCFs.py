import os, copy, datetime, pandas as pd, datetime as dt, numpy as np
from os import path
from netCDF4 import Dataset
from GetRenewableCFs import *

#Output: dfs of wind and solar generation (8760 dt rows, arbitrary cols)
def getNewRenewableCFs(genFleet, tgtTz, weatherYears, currYear, pRegionShapes, 
            nonCCReanalysis, climateChange, cesmMembers, interconn, maxCapPerTech):
    #Importing single RE timeseries
    if not climateChange or len(cesmMembers)==1: 
        newCfs,maxCapPerTech = getSingleNewRenewableCFsTimeseries(genFleet, tgtTz, weatherYears, currYear, 
            pRegionShapes, nonCCReanalysis, climateChange, cesmMembers, interconn, maxCapPerTech)
    #Importing multiple RE timeseries for different ensemble members
    else: 
        newCfsAll = list()
        for cesmMember in cesmMembers:
            newCfs,maxCapPerTech = getSingleNewRenewableCFsTimeseries(genFleet, tgtTz, weatherYears, currYear, 
                pRegionShapes, nonCCReanalysis, climateChange, cesmMember, interconn, maxCapPerTech)
            newCfs.columns = pd.MultiIndex.from_product([[cesmMember], newCfs.columns], names=['ensembleMember', 'locs'])
            newCfsAll.append(newCfs)
        #Combine into 1 array
        newCfs = pd.concat(newCfsAll,axis=1)
    return newCfs,maxCapPerTech

#Output: dfs of wind and solar generation (8760 dt rows, arbitrary cols)
def getSingleNewRenewableCFsTimeseries(genFleet, tgtTz, weatherYears, currYear, 
        pRegionShapes, nonCCReanalysis, climateChange, cesmMember, interconn, maxCapPerTech):
    if currYear > 2050 and climateChange == False: currYear = 2050
    
    #Isolate wind & solar units
    windUnits, solarUnits = getREInFleet('Wind', genFleet), getREInFleet('Solar PV', genFleet)
    
    #Get list of wind / solar sites in region.
    lats,lons,cf,latlonRegion = loadData(weatherYears,pRegionShapes,cesmMember,interconn,climateChange,nonCCReanalysis)

    #Match existing gens to CFs (use this for figuring our spare capacity available at each coordinate given existing generators)
    get_cf_index(windUnits, lats, lons),get_cf_index(solarUnits, lats, lons)

    #Calculate new CFs for given met year, but (in some cases) setting dt index to currYear
    yrForCFs = weatherYears if (climateChange or (nonCCReanalysis and interconn == 'WECC')) else [currYear] #if not CC, relabel fixed met year to future year; if CC, have many years
    stateBounds = latlonRegion.reset_index(drop=True)
    stateBounds.columns = range(stateBounds.columns.size)
    cf = enforceStateBounds(cf, stateBounds)
    windCfs = calcNewCfs(windUnits, lats, lons, cf, 'wind', yrForCFs, maxCapPerTech['Wind'])
    solarCfs = calcNewCfs(solarUnits, lats, lons, cf, 'solar', yrForCFs, maxCapPerTech['Solar'])
    
    #Modify max capacity of wind and solar per grid cell based on size of grid cell (accounts for grid cells < 1 degree x 1 degree)
    for re in ['Wind','Solar']: maxCapPerTech[re] = maxCapPerTech[re] * abs(lats[1]-lats[0]) * abs(lons[1]-lons[0])

    #Downscale if desired - REPLACED 1/11/23 BY INSTEAD REMOVING LOWEST CF SITES PER REGION (SEE ADDWSSITESTONEWTECHS.PY)
    # windCfs, solarCfs = windCfs[windCfs.columns[::reDownFactor]],solarCfs[solarCfs.columns[::reDownFactor]]
    
    #Shift to target timezone
    print('Commented out shift timezone for Hari runs (GetNewRenewableCFs)! Ask Michael.')
    # if not climateChange: windCfs, solarCfs = shiftTz(windCfs, tgtTz, currYear, 'wind'), shiftTz(solarCfs, tgtTz, currYear, 'solar')
    return pd.concat([windCfs, solarCfs], axis=1),maxCapPerTech

#maxBuildPerLatLong: max MW per lat x long coordinate. 
def calcNewCfs(existingGens, lats, lons, cf, re, yrForCFs, maxBuildPerLatLong): 
    #Pull number of time steps from CF array (indexed by lat/lon/time, so time is idx 2)
    tSteps = cf[re].shape[2] 
    f = 'H' if tSteps >= 8760 else 'D'
    #For each lat/lon, check existing capacity and, if spare room for more renewables, add CFs
    cfs = dict()
    latDiff,lonDiff = abs(lats[1]-lats[0]),abs(lons[1]-lons[0])
    for latIdx in range(len(lats)):
        for lonIdx in range(len(lons)):
            lat, lon = lats[latIdx], lons[lonIdx]
            gensAtLoc = existingGens.loc[(existingGens['lat idx'] == latIdx) & (existingGens['lon idx'] == lonIdx)]
            existingCap = gensAtLoc['Capacity (MW)'].astype(float).sum()
            coordCfs = cf[re][latIdx, lonIdx, :]
            if existingCap > 0 or coordCfs.sum() > 0: #filter out coords w/ no gen
                if existingCap < (maxBuildPerLatLong*latDiff*lonDiff): 
                    cfs[re + 'lat' + str(round(lat, 3)) + 'lon' + str(round(lon, 3))] = coordCfs
                else:
                    cfs[re + 'lat' + str(round(lat, 3)) + 'lon' + str(round(lon, 3))] = 0
                    print('RE capacity has maxed out @ lat/long:', lat, lon, ', so CFs for that site are set to 0 to deter further investment.')
    #Add dt and set to Dataframe
    idx = pd.date_range('1/1/' + str(yrForCFs[0]) + ' 0:00','12/31/' + str(yrForCFs[-1]) + ' 23:00', freq=f)
#    idx = idx.drop(idx[(idx.month == 2) & (idx.day == 29)])
    return pd.DataFrame(cfs, index=idx)

def enforceStateBounds(cf, stateBounds):
   for re in cf:
        for row in stateBounds.index:
            for col in stateBounds.columns:
                cf[re][row,col] *= stateBounds.loc[row,col]
    # plotCFs(cf)
   return cf



# import matplotlib.pyplot as plt
# def plotCFs(cf):
#     avgCfs,lats,lons = np.zeros((23,23)),np.zeros(23),np.zeros(23)
#     for re in cf:
#         cfs = cf[re]
#         for lat in range(cfs.shape[0]):
#             for lon in range(cfs.shape[1]):
#                 avgCfs[lat,lon] = cfs[lat,lon].mean()
#                 # lats[lat] = cfs['lat'][lat]
#                 # lons[lon] = cfs['lon'][lon]

#         plt.figure()
#         ax = plt.subplot(111)
#         im = ax.contourf(avgCfs,cmap='plasma')#,extent = [np.min(lons),np.max(lons),np.min(lats),np.max(lats)])
#         cbar = ax.figure.colorbar(im, ax=ax)#, ticks=np.arange(vmin,vmax,int((vmax-vmin)/5)))
#         plt.title(re)
#     plt.show()
