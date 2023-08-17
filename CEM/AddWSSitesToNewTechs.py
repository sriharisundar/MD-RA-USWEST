import pandas as pd, geopandas as gpd
from SetupTransmissionAndZones import assignGensToPRegions

def addWSSitesToNewTechs(newCfsOrig,newTechsCE,pRegionShapes,reDownFactor):
    #Create copy for isolating sites
    newCfs = newCfsOrig.copy()

    #If have multi-index of ensemble members, just use a single member for locations
    if type(newCfs.columns)==pd.core.indexes.multi.MultiIndex:
        firstEnsembleMember = newCfs.columns.get_level_values('ensembleMember')[0]
        newCfs = newCfs[firstEnsembleMember]

    sitesDfList = list()
    #For wind & solar, repeat tech row for each potential site, then remove original row
    for l,f in zip(['wind','solar'],['Wind','Solar']):
        re = newTechsCE.loc[newTechsCE['FuelType']==f]
        sites = [c for c in newCfs if l in c]
        sitesDf = pd.concat([re]*len(sites),ignore_index=True)
        sitesDf['PlantType'] = sites
        #Get lat/lon
        txt = sitesDf['PlantType'].str.split('lat',expand=True)[1]
        sitesDf[['Latitude','Longitude']] = txt.str.split('lon',expand=True).astype(float)
        newTechsCE.drop(re.index,inplace=True)
        sitesDfList.append(sitesDf)
    
    #Combine wind & solar rows into df, then map to regions
    sitesDf = pd.concat(sitesDfList)
    sitesDf = sitesDf.drop('region',axis=1)
    sitesDf = assignGensToPRegions(sitesDf,pRegionShapes)

    #Drop bottom percentile of wind and solar sites by annual average CF in each region
    allSitesToDrop = list()
    for region in sitesDf['region'].unique():
        regionGens = sitesDf.loc[sitesDf['region']==region]
        for re in ['wind','solar']:
            regionFTGens = regionGens.loc[regionGens['PlantType'].str.contains(re)]['PlantType'].values
            meanCfs = newCfs[regionFTGens].mean()
            meanCfs.sort_values(inplace=True) #ascending
            sitesToDrop = meanCfs.iloc[:int(meanCfs.shape[0]*reDownFactor)].index.values
            allSitesToDrop.extend(sitesToDrop)
    newCfsOrig.drop(allSitesToDrop,axis=1,inplace=True)
    sitesDf = sitesDf.loc[~sitesDf['PlantType'].isin(allSitesToDrop)]

    #Add remaining WS sites onto newTechsCE
    newTechsCE = pd.concat([newTechsCE,sitesDf],ignore_index=True)
    newTechsCE.reset_index(inplace=True,drop=True)

    #Create GAMS symbol as plant type + region
    newTechsCE['GAMS Symbol'] = newTechsCE['PlantType'] + newTechsCE['region']
    reRows = newTechsCE.loc[newTechsCE['ThermalOrRenewableOrStorage']=='renewable']
    reRows.index = reRows['PlantType']
    gamsDict = reRows['GAMS Symbol'].to_dict()

    #Relabel new CF columns to match GAMS symbols (needed for GAMS model later)
    if type(newCfsOrig.columns)==pd.core.indexes.multi.MultiIndex:
        newCfsOrig = newCfsOrig.T.loc[:,[k for k in gamsDict],:].rename(gamsDict,level='locs')
        newCfsOrig = newCfsOrig.T
    else:
        newCfsOrig = newCfsOrig[[k for k in gamsDict]].rename(gamsDict,axis=1)

    return newTechsCE,newCfsOrig
