import os, copy, datetime, pandas as pd, geopandas as gpd, datetime as dt, numpy as np, xarray as xr
from os import path
from statistics import mode
from netCDF4 import Dataset

#Output: dfs of wind and solar generation (8760 dt rows, arbitrary cols)
def getREGen(genFleet, tgtTz, weatherYears, currYear, pRegionShapes, nonCCReanalysis, climateChange, cesmMembers, interconn):
    #Importing single RE timeseries
    if not climateChange or len(cesmMembers)==1: 
        windGen, solarGen, windGenByRegion, solarGenByRegion, latlonRegion = getSingleREGenTimeseries(genFleet, 
            tgtTz, weatherYears, currYear, pRegionShapes, nonCCReanalysis, climateChange, cesmMembers, interconn)
    #Importing multiple RE timeseries for different ensemble members
    else: 
        windGens,solarGens,windGenByRegions,solarGenByRegions = list(),list(),list(),list()
        for cesmMember in cesmMembers:
            windGen, solarGen, windGenByRegion, solarGenByRegion, latlonRegion = getSingleREGenTimeseries(genFleet, 
                tgtTz, weatherYears, currYear, pRegionShapes, nonCCReanalysis, climateChange, cesmMember, interconn)

            windGen.columns = pd.MultiIndex.from_product([[cesmMember], windGen.columns.to_list()], names=['ensembleMember', 'locs'])
            solarGen.columns = pd.MultiIndex.from_product([[cesmMember], solarGen.columns.to_list()], names=['ensembleMember', 'locs'])
            windGenByRegion.columns = pd.MultiIndex.from_product([[cesmMember], windGenByRegion.columns], names=['ensembleMember', 'region'])
            solarGenByRegion.columns = pd.MultiIndex.from_product([[cesmMember], solarGenByRegion.columns], names=['ensembleMember', 'region'])

            windGens.append(windGen),solarGens.append(solarGen),windGenByRegions.append(windGenByRegion),solarGenByRegions.append(solarGenByRegion)

        #Combine into 1 array
        windGen,solarGen,windGenByRegion,solarGenByRegion = pd.concat(windGens,axis=1),pd.concat(solarGens,axis=1),pd.concat(windGenByRegions,axis=1),pd.concat(solarGenByRegions,axis=1)
    return windGen, solarGen, windGenByRegion, solarGenByRegion

def getSingleREGenTimeseries(genFleet, tgtTz, weatherYears, currYear, pRegionShapes, nonCCReanalysis, climateChange, cesmMember, interconn):
    if currYear > 2050 and climateChange == False: currYear = 2050
    #Isolate wind & solar units
    windUnits, solarUnits = getREInFleet('Wind', genFleet),getREInFleet('Solar', genFleet)
    #Get list of wind / solar sites in region and their CFs
    lats,lons,cf,latlonRegion = loadData(weatherYears,pRegionShapes,cesmMember,interconn,climateChange,nonCCReanalysis)
    #Match to CFs
    get_cf_index(windUnits, lats, lons),get_cf_index(solarUnits, lats, lons)
    #Get hourly generation (8760 x n df, n = num generators). Use given met year data but set dt index to currYear.
    yrForCFs = weatherYears if (climateChange or (nonCCReanalysis and interconn == 'WECC')) else [currYear] #if not CC, relabel fixed met year to future year; if CC, have many years
    windGen = convertCFToGeneration(windUnits,cf['wind'],yrForCFs)
    solarGen = convertCFToGeneration(solarUnits,cf['solar'],yrForCFs)
    #Shift into right tz
    print('Commented out shift timezone for Hari runs (GetRenewableCFs)! Ask Michael.')
    #windGen, solarGen = shiftTz(windGen, tgtTz, currYear, 'wind'),shiftTz(solarGen, tgtTz, currYear, 'solar')
    #Combine by region and fill in missing regions with zeros
    windGenByRegion = windGen.groupby(level='region', axis=1).sum()
    solarGenByRegion = solarGen.groupby(level='region', axis=1).sum()
    for genByRegion in [windGenByRegion, solarGenByRegion]:
        regionsNoGen = [r for r in genFleet['region'].unique() if r not in genByRegion.columns]
        for r in regionsNoGen: genByRegion[r] = 0
    return windGen, solarGen, windGenByRegion, solarGenByRegion, latlonRegion

def getREInFleet(reType, genFleet):
    reUnits = genFleet.loc[genFleet['FuelType'] == reType]
    return reUnits

#Get list of wind / solar sites in region and their CFs
def loadData(weatherYears,pRegionShapes,cesmMember,interconn,climateChange,nonCCReanalysis):
    if climateChange: lats, lons, cf, latlonRegion = loadClimateData(weatherYears, pRegionShapes, cesmMember, interconn)
    else: 
        if not nonCCReanalysis: lats, lons, cf, latlonRegion = loadNRELData(weatherYears, pRegionShapes)
        elif interconn == 'WECC': lats, lons, cf, latlonRegion = loadERA5Data(weatherYears, pRegionShapes,interconn) 
        else: lats, lons, cf, latlonRegion = loadMerraData(weatherYears, pRegionShapes) 
    return lats, lons, cf, latlonRegion

# Outputs: numpy arrays of lats and lons, and then a dictionary w/ wind and solar cfs
# as an np array of axbxc, where a/b/c = # lats/# longs/# hours in year
def loadClimateData(weatherYears, pRegionShapes, cesmMember, interconn, dataDir=os.path.join('Data','CESM')):
    print('Using CESM met data!')
    cesmMember = cesmMember[0] if type(cesmMember) == list else cesmMember
    #File and dir
    solarFile, windFile = path.join(dataDir, interconn.lower() + '_solarCF_' + cesmMember + '.nc'), path.join(dataDir, interconn.lower() + '_windCF_' + cesmMember + '.nc')
    # Error Handling
    if not (path.exists(solarFile) and path.exists(windFile)):
        error_message = 'Renewable Generation files not available:\n\t'+solarFile+'\n\t'+windFile
        raise RuntimeError(error_message)
    #Load data
    solarPowGen,windPowGen = xr.open_dataset(solarFile),xr.open_dataset(windFile) #assume solar and wind cover same geographic region

    #Convert long from 0-360 to -180-180
    if max(abs(np.array(solarPowGen.lon[:])))>180: 
        solarPowGen = solarPowGen.assign_coords(lon=(((solarPowGen.lon + 180) % 360) - 180))
        windPowGen = windPowGen.assign_coords(lon=(((windPowGen.lon + 180) % 360) - 180))

    #Slice weather years
    solarPowGen = solarPowGen.sel(time=slice(str(weatherYears[0]) + "-01-01", str(weatherYears[-1]) + "-12-31"))
    windPowGen = windPowGen.sel(time=slice(str(weatherYears[0]) + "-01-01", str(weatherYears[-1]) + "-12-31"))
    solarPowGen,windPowGen = solarPowGen.sel(member_id=cesmMember),windPowGen.sel(member_id=cesmMember)

    #Get lat and lons for both datasets
    latsPd,lonsPd = pd.DataFrame(solarPowGen['lat'], columns = ['lat']),pd.DataFrame(solarPowGen['lon'], columns=['lon'])
    latlonList = [(i, j) for i in latsPd.lat for j in lonsPd.lon]
    latlonPd = pd.DataFrame(data=latlonList, columns=['lat', 'lon'])
    latlonGpd = gpd.GeoDataFrame(latlonPd, geometry=gpd.points_from_xy(latlonPd.lon, latlonPd.lat))
    latlonPshapeJoin = gpd.sjoin(latlonGpd, pRegionShapes, how="inner", op='intersects')
    latlonPshapeJoin = latlonPshapeJoin.sort_values(by=['lat', 'lon'])
    latlonPshapeJoin = latlonPshapeJoin.reset_index()
    latlonRegion = pd.crosstab(latlonPshapeJoin['lat'],latlonPshapeJoin['lon']).reindex(index = solarPowGen['lat'],columns = solarPowGen['lon'],fill_value=0)

    #Store data
    cf = dict()
    cf["solar"] = np.array(solarPowGen.variables['Solar_CF'][:])
    cf["wind"] = np.array(windPowGen.variables['Wind_CF'][:])
    solarPowGen.close(), windPowGen.close()

    #Reshape arrays so index is lat/lon/time instead of time/lat/lon
    if cf['solar'].shape[0] > 360: #0 idx is time dimension
        cf['solar'],cf['wind'] = np.swapaxes(cf['solar'],0,1),np.swapaxes(cf['wind'],0,1)
        cf['solar'],cf['wind'] = np.swapaxes(cf['solar'],1,2),np.swapaxes(cf['wind'],1,2)

    return np.array(solarPowGen['lat'][:]),np.array(solarPowGen['lon'][:]),cf,latlonRegion

# Outputs: numpy arrays of lats and lons, and then a dictionary w/ wind and solar cfs
# as an np array of axbxc, where a/b/c = # lats/# longs/# hours in year
def loadERA5Data(weatherYears, pRegionShapes, interconn, dataDir=os.path.join('Data','ERA5')):
    print('Using ERA5 met data!')
    #File and dir
    # solarFile, windFile = path.join(dataDir, interconn.lower() + '_solarCF_ERA5.nc'), path.join(dataDir, interconn.lower() + '_windCF_ERA5.nc')
    solarFile, windFile = path.join(dataDir, interconn.lower() + '_solarCF_ERA5_hourly_PST.nc'), path.join(dataDir, interconn.lower() + '_windCF_ERA5_hourly_PST.nc')
    # Error Handling
    if not (path.exists(solarFile) and path.exists(windFile)):
        error_message = 'Renewable Generation files not available:\n\t'+solarFile+'\n\t'+windFile
        raise RuntimeError(error_message)
    #Load data
    solarPowGen,windPowGen = xr.open_dataset(solarFile),xr.open_dataset(windFile) #assume solar and wind cover same geographic region

    #Convert long from 0-360 to -180-180
    if max(abs(np.array(solarPowGen.lon[:])))>180: 
        solarPowGen = solarPowGen.assign_coords(lon=(((solarPowGen.lon + 180) % 360) - 180))
        windPowGen = windPowGen.assign_coords(lon=(((windPowGen.lon + 180) % 360) - 180))

    #Slice weather years
    solarPowGen = solarPowGen.sel(time=slice(str(weatherYears[0]) + "-01-01", str(weatherYears[-1]) + "-12-31"))
    windPowGen = windPowGen.sel(time=slice(str(weatherYears[0]) + "-01-01", str(weatherYears[-1]) + "-12-31"))

    #Round half-hour timestamps to hourly so align w/ demand
    solarPowGen['time'] = solarPowGen['time'].dt.floor('H')
    windPowGen['time'] = windPowGen['time'].dt.floor('H')

    #Get lat and lons for both datasets
    latsPd,lonsPd = pd.DataFrame(solarPowGen['lat'], columns = ['lat']),pd.DataFrame(solarPowGen['lon'], columns=['lon'])
    latlonList = [(i, j) for i in latsPd.lat for j in lonsPd.lon]
    latlonPd = pd.DataFrame(data=latlonList, columns=['lat', 'lon'])
    latlonGpd = gpd.GeoDataFrame(latlonPd, geometry=gpd.points_from_xy(latlonPd.lon, latlonPd.lat))
    latlonPshapeJoin = gpd.sjoin(latlonGpd, pRegionShapes, how="inner", op='intersects')
    latlonPshapeJoin = latlonPshapeJoin.sort_values(by=['lat', 'lon'])
    latlonPshapeJoin = latlonPshapeJoin.reset_index()
    latlonRegion = pd.crosstab(latlonPshapeJoin['lat'],latlonPshapeJoin['lon']).reindex(index = solarPowGen['lat'],columns = solarPowGen['lon'],fill_value=0)

    #Eliminate any negative CFs (solar has very small negative #s instead of zero generation)
    solarCFs = np.array(solarPowGen.variables['Solar_CF'][:])
    solarCFs[solarCFs<0] = 0
    windCFs = np.array(windPowGen.variables['Wind_CF'][:])
    windCFs[windCFs<0] = 0

    #Close xarray files
    solarPowGen.close(), windPowGen.close()

    #Store data in dictionary
    cf = dict()
    cf['solar'],cf['wind'] = solarCFs,windCFs

    #Reshape arrays so index is lat/lon/time instead of time/lat/lon
    if cf['solar'].shape[0] > 360: #0 idx is time dimension
        cf['solar'],cf['wind'] = np.swapaxes(cf['solar'],0,1),np.swapaxes(cf['wind'],0,1)
        cf['solar'],cf['wind'] = np.swapaxes(cf['solar'],1,2),np.swapaxes(cf['wind'],1,2)

    return np.array(solarPowGen['lat'][:]),np.array(solarPowGen['lon'][:]),cf,latlonRegion

# Get all necessary information from powGen netCDF files, VRE capacity factors and lat/lons
# Outputs: numpy arrays of lats and lons, and then a dictionary w/ wind and solar cfs
# as an np array of axbxc, where a/b/c = # lats/# longs/# hours in year
def loadMerraData(weatherYears, pRegionShapes, dataDir=os.path.join('Data','MERRA')):
    print('Using MERRA2 met data!')
    #File and dir
    solarFile, windFile = path.join(dataDir, str(weatherYears) + '_solar_generation_cf_US.nc'), path.join(dataDir, str(weatherYears) + '_wind_generation_cf_US.nc')
    # Error Handling
    if not (path.exists(solarFile) and path.exists(windFile)):
        error_message = 'Renewable Generation files not available:\n\t'+solarFile+'\n\t'+windFile
        raise RuntimeError(error_message)
    #Load data
    solarPowGen = Dataset(solarFile)
    windPowGen = Dataset(windFile) #assume solar and wind cover same geographic region

    #Get lat and lons for both datasets
    lats,lons = np.array(solarPowGen.variables['lat'][:]), np.array(solarPowGen.variables['lon'][:])
    latsPd = pd.DataFrame(lats, columns = ['lat'])
    lonsPd = pd.DataFrame(lons, columns=['lon'])

    latlonList = [(i, j) for i in latsPd.lat for j in lonsPd.lon]

    latlonPd = pd.DataFrame(data=latlonList, columns=['lat', 'lon'])
    latlonGpd = gpd.GeoDataFrame(latlonPd, geometry=gpd.points_from_xy(latlonPd.lon, latlonPd.lat))
    latlonPshapeJoin = gpd.sjoin(latlonGpd, pRegionShapes, how="inner", op='intersects')
    latlonPshapeJoin = latlonPshapeJoin.sort_values(by=['lat', 'lon'])
    latlonPshapeJoin = latlonPshapeJoin.reset_index()
    latlonRegion = (pd.crosstab(latlonPshapeJoin['lat'],latlonPshapeJoin['lon']).reindex(index = lats,columns = lons,fill_value=0))

    #Store data
    cf = dict()
    cf["solar"] = np.array(solarPowGen.variables['cf'][:])
    cf["wind"] = np.array(windPowGen.variables['cf'][:])
    solarPowGen.close(), windPowGen.close()

    # Error Handling
    if cf['solar'].shape != (lats.size, lons.size, 8760):
        print("powGen Error. Expected array of shape",lats.size,lons.size,8760,"Found:",cf['solar'].shape)
        return -1
    return lats,lons,cf,latlonRegion

# Get all necessary information from powGen netCDF files, VRE capacity factors and lat/lons
# Outputs: numpy arrays of lats and lons, and then a dictionary w/ wind and solar cfs
# as an np array of axbxc, where a/b/c = # lats/# longs/# hours in year
def loadNRELData(weatherYears, pRegionShapes):
    print('Using NREL met data!')
    # File and dir
    dataDir = 'Data\\RE'
    solarFile, windFile = path.join(dataDir, 'solar_cf_' + str(weatherYears) + '.csv'), path.join(dataDir, 'wind_cf_' + str(weatherYears) + '.csv')
    # Error Handling
    if not (path.exists(solarFile) and path.exists(windFile)):
        error_message = 'Renewable Generation files not available:\n\t'+solarFile+'\n\t'+windFile
        raise RuntimeError(error_message)
    # Load data
    solarPowGen = pd.read_csv(solarFile)
    windPowGen = pd.read_csv(windFile) #assume solar and wind cover same geographic region

    # Get lat and lons for both datasets
    lats_temp,lons_temp = solarPowGen['lat'],solarPowGen['lon']
    latsPd = pd.DataFrame(lats_temp, columns = ['lat'])
    latsPd = latsPd.drop_duplicates()
    latsPd = latsPd.sort_values(by=['lat'])
    latsPd = latsPd.reset_index()
    lats = latsPd['lat'].to_numpy()
    lats_temp = lats_temp.to_numpy()

    lonsPd = pd.DataFrame(lons_temp, columns=['lon'])
    lonsPd = lonsPd.drop_duplicates()
    lonsPd = lonsPd.sort_values(by=['lon'])
    lonsPd = lonsPd.reset_index()
    lons = lonsPd['lon'].to_numpy()
    lons_temp = lons_temp.to_numpy()

    latsAll = pd.Series(latsPd['lat'])
    lonsAll = pd.Series(lonsPd['lon'])

    latlonList = [(i, j)
               for i in latsPd.lat
               for j in lonsPd.lon]

    latlonPd = pd.DataFrame(data=latlonList, columns=['lat', 'lon'])
    latlonGpd = gpd.GeoDataFrame(latlonPd, geometry=gpd.points_from_xy(latlonPd.lon, latlonPd.lat))
    latlonPshapeJoin = gpd.sjoin(latlonGpd, pRegionShapes, how="inner", op='intersects')
    latlonPshapeJoin = latlonPshapeJoin.sort_values(by=['lat', 'lon'])
    latlonPshapeJoin = latlonPshapeJoin.reset_index()

    latlonRegion = (pd.crosstab(latlonPshapeJoin['lat'],
                                 latlonPshapeJoin['lon']).reindex(index = latsAll,
                                                                    columns = lonsAll,fill_value=0))
    # Store data
    cf = dict()
    cf["solar"] = np.zeros(shape=(len(latsPd), len(lonsPd), 8760))
    cf["wind"] = np.zeros(shape=(len(latsPd), len(lonsPd), 8760))

    for i in list(range(len(latsPd))):
        for j in list(range(len(lonsPd))):
            if (lats[i] in lats_temp and lons[j] in lons_temp):
                k_temp = np.where((lats_temp == lats[i]) & (lons_temp == lons[j]))
                try:
                    k = k_temp[0][0]
                    cf["solar"][i][j][:] = solarPowGen.iloc[k, 2:]
                    cf["wind"][i][j][:] = windPowGen.iloc[k, 2:]
                except:
                    cf["solar"][i][j][:] = solarPowGen.iloc[0, 2:]*0
                    cf["wind"][i][j][:] = windPowGen.iloc[0, 2:]*0

    # Error Handling
    if cf['solar'].shape != (lats.size, lons.size, 8760):
        print("powGen Error. Expected array of shape",lats.size,lons.size,8760,"Found:",cf['solar'].shape)
        return -1
    return lats,lons,cf,latlonRegion

# Convert the latitude and longitude of the vg into indices for capacity factor matrix,
#then save that index into df
# More detail: The simulated capacity factor maps are of limited resolution. This function
#               identifies the nearest simulated location for renewable energy generators
#               and replaces those generators' latitudes and longitudes with indices for 
#               for the nearest simulated location in the capacity factor maps
def get_cf_index(RE_generators, powGen_lats, powGen_lons):
    RE_generators.loc[:,"lat idx"] = find_nearest_impl(RE_generators["Latitude"].astype(float), powGen_lats).astype(int)
    RE_generators.loc[:,"lon idx"] = find_nearest_impl(RE_generators["Longitude"].astype(float), powGen_lons).astype(int)

# Find index of nearest coordinate. Implementation of get_RE_index
def find_nearest_impl(actual_coordinates, discrete_coordinates):
    indices = []
    for coord in actual_coordinates:
        indices.append((np.abs(coord-discrete_coordinates)).argmin())
    return np.array(indices)

# Find expected hourly capacity for RE generators. Of shape (8760 hrs, num generators)
# When running CC, using average demand, so don't scale up generation by 24
def convertCFToGeneration(RE_generators,cf,yrForCFs):
    #Pull number of time steps from CF array (indexed by lat/lon/time, so time is idx 2)
    tSteps = cf.shape[2] 
    f = 'H' if tSteps >= 8760 else 'D'
    # repeat nameplate capacity in array
    RE_nameplate = np.tile(RE_generators["Capacity (MW)"].astype(float),(tSteps,1))
    # multiply by variable hourly capacity factor
    times = np.tile(np.arange(tSteps),(RE_generators['Capacity (MW)'].size,1)).T 
    RE_capacity = np.multiply(RE_nameplate, cf[RE_generators["lat idx"], RE_generators["lon idx"], times])
    # convert to datetime index
    idx = pd.date_range('1/1/'+ str(yrForCFs[0]) + ' 0:00','12/31/' + str(yrForCFs[-1]) + ' 23:00',freq=f)
#    if (len(idx) > 365 and f == 'D') or (len(idx) > 8760 and f == 'H'): idx = idx.drop(idx[(idx.month==2) & (idx.day ==29)])
    reGen = pd.DataFrame(RE_capacity,index=idx,columns=[RE_generators['GAMS Symbol'].values,RE_generators['region'].values])
    reGen.columns.names = ['GAMS Symbol','region']
    return reGen

#shift tz (MERRA in UTC)
def shiftTz(reGen,tz,yr,reType):
    origIdx = reGen.index
    tzOffsetDict = {'PST':-7,'CST': -6,'EST': -5}
    reGen.index = reGen.index.shift(tzOffsetDict[tz],freq='H')
    reGen = reGen[reGen.index.year==yr]
    reGen = reGen.append([reGen.iloc[-1]]*abs(tzOffsetDict[tz]),ignore_index=True)
    if reType=='solar': reGen.iloc[-5:] = 0 #set nighttime hours to 0
    reGen.index=origIdx
    return reGen
