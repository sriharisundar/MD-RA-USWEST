import sys, os, csv, operator, copy, time, random, warnings, numpy as np, datetime as dt, pandas as pd
from os import path; from netCDF4 import Dataset; from gams import *
from SetupGeneratorFleet import setupGeneratorFleet,compressAndAddSizeDependentParams
from AddCoolingTypes import addCoolingTypes
from ProcessHydro import processHydro
from UpdateFuelPriceFuncs import updateFuelPricesAndCosts
from ImportDemand import importDemand
from DemandFuncsCE import getHoursForCE
from IsolateDataForCE import isolateDataInCEHours,isolateDataInCEBlocks
from ImportNewTechs import getNewTechs
from RetireUnitsCFPriorCE import retireUnitsCFPriorCE
from CreateFleetForCELoop import createFleetForCurrentCELoop
from GetRenewableCFs import getREGen
from GetNewRenewableCFs import getNewRenewableCFs
from AddWSSitesToNewTechs import addWSSitesToNewTechs
from ProcessCEResults import saveCEBuilds,addNewGensToFleet,addNewLineCapToLimits
from ScaleRegResForAddedWind import scaleRegResForAddedWind
from CombinePlants import combineWindSolarStoPlants
from GAMSAddSetToDatabaseFuncs import *
from GAMSAddParamToDatabaseFuncs import *
from InitializeOnOffExistingGensCE import initializeOnOffExistingGens
from ReservesWWSIS import calcWWSISReserves
from GetIncResForAddedRE import getIncResForAddedRE
from SaveCEOperationalResults import saveCapacExpOperationalData
from WriteTimeDependentConstraints import writeTimeDependentConstraints
from WriteBuildVariable import writeBuildVariable
from CreateEmptyReserveDfs import createEmptyReserveDfs
from SetupTransmissionAndZones import setupTransmissionAndZones, defineTransmissionRegions
from DefineReserveParameters import defineReserveParameters
from CalculateDerates import calculateLineThermalDerates,calculatePlantCapacityDerates
from ImportPRMCapacityAdjustments import importPRMCapacityAdjustments

# SET OPTIONS
warnings.filterwarnings("ignore")
pd.set_option('display.max_rows', 10)
pd.set_option('display.max_columns', 10)

# SCALARS
mwToGW = 1000
lbToShortTon = 2000

# ###################################################################4############
# ##### UNIVERSAL PARAMETERS ####################################################
# ###############################################################################
def setKeyParameters(climateChange):
    # ### RUNNING ON SC OR LOCAL
    runOnSC = True                                     # whether running on supercomputer

    # ### START YEAR, END YEAR, AND STEPS FOR CE
    startYear, endYear, yearStepCE = 2022, 2031, 8

    # ### BUILD LIMITS AND RE UPSAMPLING
    reDownFactor = .4                           # fraction of wind & solar sites per region dropped; sites dropped with worst CFs
    buildLimitsCase = 1                               # 1 = reference case, 2 = limited nuclear, 3 = limited CCS and nuclear, 4 = limited hydrogen storage, 5 = limited transmission
    yearIncDACS = 2050                                  #year to include DACS - set beyond end period if don't want DACS

    # ### CE OPTIONS
    if climateChange: numBlocks, daysPerBlock, daysPerPeak = 1, 365*yearStepCE, 1    # num rep time blocks, days per rep block, and days per peak block in CE
    else: numBlocks, daysPerBlock, daysPerPeak = 4,7,1
    removeHydro = False                                  #whether to remove hydropower from fleet & subtract generation from demand, or to include hydro as dispatchable in CE w/ gen limit
    greenField = False                                  # whether to run greenField (set to True) or brownfield (False)
    stoInCE,seasStoInCE = True,False                    # whether to allow new storage,new seasonal storage in CE model

    # ### GENERIC DEMAND FLEXIBILITY PARAMETERS
    demandShifter = 0                                   # Percentage of hourly demand that can be shifted
    demandShiftingBlock = 4                             # moving shifting demand window (hours)
         
    return (runOnSC,startYear,endYear,yearStepCE,reDownFactor,yearIncDACS,
        numBlocks,daysPerBlock,daysPerPeak,removeHydro,greenField,stoInCE,seasStoInCE,
        demandShifter,demandShiftingBlock,buildLimitsCase)

def setNonCCWeatherData(yr):
    demandScen = 'REFERENCE'                        # NREL EFS demand scenario: 'REFERENCE','HIGH','MEDIUM' (ref is lower than med)
    nonCCReanalysis = True                                # == True: use reanalysis as renewable data source, == False: use NSDRB and WTK
    nonCCWeatherYear = [2012] if not nonCCReanalysis else [yr]   #EFS, WTK, NSRDB: 2012 weather year. Should use 2012 if using EFS.
    return nonCCWeatherYear,demandScen,nonCCReanalysis

def stateAssumptions(interconn,buildLimitsCase,yearStepCE,ret,reBuildRateMultiplier=10):
    # ### MAX BUILDS
    #Max builds (MW) per region or, for wind & solar, grid point
    areaPerLatLongBox = 9745 #km^2 per degree lat x long (https://www.usgs.gov/faqs/how-much-distance-does-a-degree-minute-and-second-cover-your-maps?qt-news_science_products=0#qt-news_science_products)
    windDens,solarDens = .9,5.7 #W/m^2 (equiv to MW/km^2); https://www.seas.harvard.edu/news/2018/10/large-scale-wind-power-would-require-more-land-and-cause-more-environmental-impact
    if ret:
        maxCapPerTech = {'Wind': areaPerLatLongBox * windDens, 'Solar': areaPerLatLongBox * solarDens, 
                         'Thermal': 999999, 'Combined Cycle': 50000,'Storage': 100000, 'Dac': -9999999, 'CCS': 999999, 
                         'Nuclear': 999999, 'Battery Storage': 10000, 'Hydrogen': 10000, 'Transmission': 10000} 
    else:
        maxCapPerTech = {'Wind': areaPerLatLongBox * windDens, 'Solar': areaPerLatLongBox * solarDens, 
                         'Thermal': 0, 'Combined Cycle': 0,'Storage': 99999, 'Dac': 0, 'CCS': 0, 
                         'Nuclear': 0, 'Battery Storage': 99999, 'Hydrogen': 0, 'Transmission': 99999} 
    
    if buildLimitsCase == 2: maxCapPerTech['Nuclear'] = 9000 
    elif buildLimitsCase == 3: maxCapPerTech['CCS'],maxCapPerTech['Nuclear'] = 1500,9000
    elif buildLimitsCase == 4: maxCapPerTech['Hydrogen'] = 0
    elif buildLimitsCase == 5: maxCapPerTech['Transmission'] = 10
    #Max wind & solar builds (MW) per interconnection & region
    if interconn == 'WECC': #see SetupTransmissionAndZones for list of zones
        maxREPerZone = {'Wind': {'NWPP_NE':99999,'CAMX':99999,'Desert_Southwest':99999,'NWPP_Central':99999,'NWPP_NW':99999},
                        'Solar':{'NWPP_NE':99999,'CAMX':99999,'Desert_Southwest':99999,'NWPP_Central':99999,'NWPP_NW':99999}} 
        maxREInInterconn = {'Wind':6700/2*yearStepCE*reBuildRateMultiplier,'Solar':7500/1*yearStepCE*reBuildRateMultiplier} #max build per interconnection per CE run. WECC 2022 WARA: 7.5 GW solar expected to be added in 2023; EIA 860: 6.7 GW & 5.8 GW wind & solar added in 2020 & 2021 across WECC
    else:
        sys.exit('Need to set up max RE per zone! MacroCEM, stateAssumptions function')
    
    # ### CO2 EMISSION CAPS [https://www.eia.gov/environment/emissions/state/, table 3]
    if interconn == 'ERCOT': co2EmsInitial =  130594820     #METRIC TONS. Initial emission for ERCOT: 130594820.
    elif interconn == 'EI': co2EmsInitial =  1274060000
    elif interconn == 'WECC': co2EmsInitial =  248800000    #2019. METRIC TONS. wa,or,ca,nm,az,nv,ut,co,wy,id,mt

    # ### CE AND UCED/ED OPTIONS
    balAuths = 'full'                                   # full: run for all BAs in interconn. TODO: add selection of a subset of BAs. [10/31 MC note: Papa Yaw has this code]
    compressFleet = True                                                # whether to compress fleet
    tzAnalysis = {'ERCOT':'CST','EI':'EST','WECC':'PST'}[interconn]     # timezone for analysis
    fuelPrices = importFuelPrices('Reference case')                     # import fuel price time series
    transmissionEff = 0.95                                              # efficiency of transmission between zones (https://ars.els-cdn.com/content/image/1-s2.0-S2542435120305572-mmc1.pdf)
    stoEff = 0.81

    # ### CE OPTIONS
    runCE,ceOps = True,'ED'                           # ops are 'ED' or 'UC' (econ disp or unit comm constraints)
    includeRes = False                                  # whether to include reserves in CE & dispatch models (if False, multiplies reserve timeseries by 0)
    retireByAge = ret                                  # whether to retire by age or not
    prmBasis = 'demand'                      # whether basis for planning reserve margin is peak demand ('demand') or net demand ('netdemand')

    retirementCFCutoff = .3                             # retire units w/ CF lower than given value
    ptEligRetCF = ['Coal Steam']                        # which plant types retire based on capacity factor (economics)
    discountRate = 0.07 

    # ### WARNINGS OR ERRORS
    if ceOps == 'UC': sys.exit('CEwithUC.gms needs to be updated for DACS operations - add DACS constraints and include gentechs set')

    return (balAuths,co2EmsInitial,compressFleet,tzAnalysis,fuelPrices,transmissionEff,stoEff,
        runCE,ceOps,includeRes,retireByAge,prmBasis,retirementCFCutoff,ptEligRetCF,discountRate,maxCapPerTech,
        maxREPerZone,maxREInInterconn)

def storageAssumptions():
    stoMkts = 'energy'                            # energy,res,energyAndRes - whether storage participates in energy, reserve, or energy and reserve markets
    stoFTLabels = ['Energy Storage','Pumped Storage']
    stoDuration = {'Energy Storage':'st','Hydrogen':'lt','Battery Storage':'st','Flywheels':'st','Batteries':'st','Pumped Storage':'st'} # mapping plant types to short-term (st) or long-term (lt) storage
    stoPTLabels = [pt for pt in stoDuration ]
    initSOCFraction = {pt:{'st':.1,'lt':.05}[dur] for pt,dur in stoDuration.items()} # get initial SOC fraction per st or lt storage
    stoMinSOC = 0     # min SOC
    return stoMkts,stoFTLabels,stoDuration,stoPTLabels,initSOCFraction,stoMinSOC

def importFuelPrices(fuelPriceScenario):
    fuelPrices = pd.read_csv(os.path.join('Data', 'Energy_Prices_Electric_Power.csv'), skiprows=4, index_col=0)
    fuelPrices = fuelPrices[[col for col in fuelPrices if fuelPriceScenario in col]]
    fuelPrices.columns = [col.split(':')[0] for col in fuelPrices.columns]
    return fuelPrices    
# ###############################################################################
# ###############################################################################
# ###############################################################################

# ###############################################################################
# ##### MASTER FUNCTION #########################################################
# ###############################################################################
#Main function to call. 
#Inputs: interconnection (EI, WECC, ERCOT); CO2 emissions in final year as fraction
#of initial CO2 emissions.
def macroCEM(interconn,co2EndPercent,wsGenFracOfDemand,prm,yr,ret,cesmMembers,climateChange):
    # Set key parameters
    (runOnSC,startYear,endYear,yearStepCE,reDownFactor,yearIncDACS,
        numBlocks,daysPerBlock,daysPerPeak,removeHydro,greenField,stoInCE,seasStoInCE,
        demandShifter,demandShiftingBlock,buildLimitsCase) = setKeyParameters(climateChange)

    # Set assumptions
    nonCCWeatherYear,demandScen,nonCCReanalysis = setNonCCWeatherData(yr)
    (balAuths,co2EmsInitial,compressFleet,tzAnalysis,fuelPrices,transmissionEff,stoEff,
        runCE,ceOps,includeRes,retireByAge,prmBasis,retirementCFCutoff,ptEligRetCF,discountRate,
        maxCapPerTech,maxREPerZone,maxREInInterconn) = stateAssumptions(interconn,buildLimitsCase,yearStepCE,ret)
    stoMkts,stoFTLabels,stoDuration,stoPTLabels,initSOCFraction,stoMinSOC = storageAssumptions()
    (regLoadFrac, contLoadFrac, regErrorPercentile, flexErrorPercentile, regElig, contFlexInelig, regCostFrac,
        rrToRegTime, rrToFlexTime, rrToContTime) = defineReserveParameters(stoMkts, stoFTLabels)

    # Create results directory
    resultsDirAll = 'Results'+interconn+'C'+str(co2EndPercent)+'RE'+str(wsGenFracOfDemand)+'PRM'+str(prm)+('EM'+cesmMembers[0] if climateChange else '')+'Yr'+str(yr)+'Ret'+('1' if ret else '0')
    if not os.path.exists(resultsDirAll): os.makedirs(resultsDirAll)
    pd.Series(co2EmsInitial).to_csv(os.path.join(resultsDirAll,'initialCO2Ems.csv'))

    # Setup initial fleet and demand
    (genFleet, compressedGens, transRegions, pRegionShapes, lineLimits, lineDists, lineCosts) = getInitialFleetAndTransmission(startYear, 
        fuelPrices, compressFleet, resultsDirAll, regElig, regCostFrac, stoMinSOC, greenField, interconn, balAuths, 
        contFlexInelig, stoFTLabels, stoPTLabels, stoEff, stoInCE)

    # Run CE and/or ED/UCED
    for currYear in range(startYear, endYear, yearStepCE):
        # Set CO2 cap
        currCo2Cap = co2EmsInitial + (co2EndPercent/100*co2EmsInitial - co2EmsInitial)/((endYear-1) - startYear) * (currYear - startYear)
        print('Entering year ', currYear, ' with CO2 cap (million tons):', round(currCo2Cap/1e6))

        # Create results directory
        resultsDir = os.path.join(resultsDirAll,str(currYear))
        if not os.path.exists(resultsDir): os.makedirs(resultsDir)
        
        # Get weather years
        weatherYears = list(range(currYear-yearStepCE+1,currYear+1)) if climateChange else nonCCWeatherYear #+1 so, e.g. 2021-2030

        # Get electricity demand profile
        demand = importDemand(weatherYears,demandScen,currYear,transRegions,cesmMembers)
        demand.to_csv(os.path.join(resultsDir,'demandInitial'+str(currYear)+'.csv'))

        # Run CE
        if currYear > startYear and runCE:
            print('Starting CE')
            #Initialize results & inputs
            if currYear == startYear + yearStepCE: priorCEModel, priorHoursCE, genFleetPriorCE = None, None, None
            (genFleet, genFleetPriorCE, lineLimits,priorCEModel, priorHoursCE) = runCapacityExpansion(genFleet, demand, currYear, weatherYears, prm, prmBasis,
                                                                discountRate, fuelPrices, currCo2Cap, numBlocks, daysPerBlock, daysPerPeak,
                                                                retirementCFCutoff, retireByAge, tzAnalysis, resultsDir,
                                                                maxCapPerTech, regLoadFrac, contLoadFrac, regErrorPercentile, flexErrorPercentile,
                                                                rrToRegTime, rrToFlexTime, rrToContTime, regElig, regCostFrac, ptEligRetCF,
                                                                genFleetPriorCE, priorCEModel, priorHoursCE, stoInCE, seasStoInCE,
                                                                ceOps, stoMkts, initSOCFraction, includeRes, reDownFactor, demandShifter,
                                                                demandShiftingBlock, runOnSC, interconn, yearIncDACS, transRegions,pRegionShapes,
                                                                lineLimits, lineDists, lineCosts, contFlexInelig, 
                                                                nonCCReanalysis, stoFTLabels, transmissionEff, removeHydro, climateChange, cesmMembers,
                                                                maxREPerZone, maxREInInterconn, compressedGens, wsGenFracOfDemand)
# ###############################################################################
# ###############################################################################
# ###############################################################################

# ###############################################################################
# ###### SET UP INITIAL FLEET AND DEMAND ########################################
# ###############################################################################
def getInitialFleetAndTransmission(startYear, fuelPrices, compressFleet, resultsDir, regElig, regCostFrac, 
        stoMinSOC, greenField, interconn, balAuths, contFlexInelig, stoFTLabels, stoPTLabels, stoEff, stoInCE):
    # GENERATORS
    genFleet = setupGeneratorFleet(interconn, startYear, fuelPrices, stoEff, stoMinSOC, stoFTLabels, stoInCE)

    # ADD COOLING TYPES TO GENERATORS TO CAPTURE DERATINGS WHEN RUNNING CLIMATE SCENARIOS
    genFleet = addCoolingTypes(genFleet, interconn)

    # DEFINE TRANSMISSION REGIONS
    transRegions = defineTransmissionRegions(interconn, balAuths)

    # TRANSMISSION
    genFleet, transRegions, limits, dists, costs, pRegionShapes = setupTransmissionAndZones(genFleet, transRegions, interconn, balAuths)
    for df, l in zip([limits, dists, costs],['Limits', 'Dists', 'Costs']): df.to_csv(os.path.join(resultsDir, 'transmission' + l + 'Initial.csv'))
    genFleet.to_csv(os.path.join(resultsDir, 'genFleetInitialPreCompression.csv'))

    # COMBINE GENERATORS FOR SMALLER GEN FLEET AND ADD SIZE DEPENDENT PARAMS (COST, REG OFFERS, UC PARAMS)
    genFleet,compressedGens = compressAndAddSizeDependentParams(genFleet, compressFleet, regElig, contFlexInelig, regCostFrac, stoFTLabels, stoPTLabels)

    # IF GREENFIELD, ELIMINATE EXISTING GENERATORS EXCEPT TINY NG, WIND, & SOLAR PLANT (TO AVOID CRASH IN LATER FUNCTIONS)
    if greenField: genFleet = stripDownGenFleet(genFleet, greenField)

    # SAVE FILES
    genFleet.to_csv(os.path.join(resultsDir,'genFleetInitial.csv')),compressedGens.to_csv(os.path.join(resultsDir,'compressedUnitsFromGenFleet.csv'))

    return genFleet, compressedGens, transRegions, pRegionShapes, limits, dists, costs
# ###############################################################################
# ###############################################################################
# ###############################################################################

# ###############################################################################
# ###### RUN CAPACITY EXPANSION #################################################
# ###############################################################################
def runCapacityExpansion(genFleet, demand, currYear, weatherYears, prm, prmBasis, discountRate, fuelPrices, currCo2Cap, numBlocks,
                         daysPerBlock, daysPerPeak, retirementCFCutoff, retireByAge, tzAnalysis, resultsDirOrig, maxCapPerTech,
                         regLoadFrac,contLoadFrac, regErrorPercentile, flexErrorPercentile, rrToRegTime, rrToFlexTime,  rrToContTime,
                         regElig, regCostFrac, ptEligRetCF, genFleetPriorCE, priorCEModel, priorHoursCE, stoInCE, seasStoInCE,
                         ceOps, stoMkts, initSOCFraction, includeRes, reDownFactor, demandShifter, demandShiftingBlock, runOnSC,
                         interconn, yearIncDACS, transRegions, pRegionShapes, lineLimits, lineDists, lineCosts, contFlexInelig,
                         nonCCReanalysis, stoFTLabels, transmissionEff, removeHydro, climateChange, cesmMembers,
                         maxREPerZone, maxREInInterconn, compressedGens, wsGenFracOfDemand):
    # ###############CREATE RESULTS DIRECTORY FOR CE RUN AND SAVE INITIAL INPUTS
    resultsDir = os.path.join(resultsDirOrig, 'CE')
    if not os.path.exists(resultsDir): os.makedirs(resultsDir)
    print('Entering CE loop for year ' + str(currYear))
    lineLimits.to_csv(os.path.join(resultsDir,'lineLimitsForCE' + str(currYear) + '.csv'))
    pd.Series(maxCapPerTech).to_csv(os.path.join(resultsDir,'buildLimitsForCE' + str(currYear) + '.csv'))
    pd.Series(currCo2Cap).to_csv(os.path.join(resultsDir,'co2CapCE' + str(currYear) + '.csv'))

    # ###############PREPARE INPUTS FOR CEM
    # Update new technology and fuel price data    
    newTechsCE = getNewTechs(regElig, regCostFrac, currYear, stoInCE, seasStoInCE,
                             fuelPrices, yearIncDACS, transRegions, contFlexInelig, weatherYears)
    genFleet = updateFuelPricesAndCosts(genFleet, currYear, fuelPrices, regCostFrac)

    # Retire units and create fleet for current CE loop
    if priorCEModel != None:                    # if not in first CE loop
        genFleet = retireUnitsCFPriorCE(genFleet, genFleetPriorCE, retirementCFCutoff,
            priorCEModel, priorHoursCE, ptEligRetCF, currYear)
    genFleet, genFleetForCE = createFleetForCurrentCELoop(genFleet, currYear, retireByAge)
    genFleetForCE.to_csv(os.path.join(resultsDir, 'genFleetForCEPreRECombine' + str(currYear) + '.csv'))
    
    # Combine wind, solar, and storage plants by region
    genFleetForCE = combineWindSolarStoPlants(genFleetForCE)
    
    # Get renewable CFs by plant and region and calculate net demand by region
    print('Loading RE data')
    windGen, solarGen, windGenRegion, solarGenRegion = getREGen(genFleet, 
        tzAnalysis, weatherYears, currYear, pRegionShapes, nonCCReanalysis, climateChange, cesmMembers, interconn)
    netDemand = demand - windGenRegion - solarGenRegion

    # Remove hydropower generation from demand using net-demand-based heuristic
    genFleetForCE,hydroGen,demand = processHydro(genFleetForCE, demand, netDemand, weatherYears, removeHydro, climateChange) 
    genFleetForCE.to_csv(os.path.join(resultsDir, 'genFleetForCE' + str(currYear) + '.csv'))

    # Get hours included in CE model (representative + special blocks)
    (hoursForCE, planningReserve, blockWeights, socScalars, planningReserveHour, blockNamesChronoList, 
        lastRepBlockNames, specialBlocksPrior) = getHoursForCE(demand, netDemand, windGenRegion, solarGenRegion,
        daysPerBlock, daysPerPeak, currYear, resultsDir, numBlocks, prm, prmBasis, climateChange)

    # Get CFs for new wind and solar sites and add wind & solar sites to newTechs
    newCfs,maxCapPerTech = getNewRenewableCFs(genFleet, tzAnalysis, weatherYears, currYear, 
        pRegionShapes, nonCCReanalysis, climateChange, cesmMembers, interconn, maxCapPerTech)
    newTechsCE,newCfs = addWSSitesToNewTechs(newCfs, newTechsCE, pRegionShapes, reDownFactor)

    # Calculating thermal power plant & thermal line capacity deratings, FORs, and capacity eligibilities towards PRM
    print('Calculating deratings and capacity adjustments')
    capDerates,capDeratesTechs,newTechsCE = calculatePlantCapacityDerates(genFleetForCE, newTechsCE, demand, weatherYears, cesmMembers, compressedGens)
    lineDerates = calculateLineThermalDerates(lineLimits, demand) 
    fors,windFOR,solarFOR,forsTechs,prmEligWindSolar = importPRMCapacityAdjustments(genFleetForCE, 
                            newTechsCE, demand, prmBasis, interconn, nonCCReanalysis, weatherYears, compressedGens, cesmMembers)

    # Initialize which generators are on or off at start of each block of hours (useful if CE has UC constraints)
    onOffInitialEachPeriod = initializeOnOffExistingGens(genFleetForCE, hoursForCE, netDemand)

    # Set reserves for existing and incremental reserves for new generators
    print('Calculating reserves')
    if includeRes:
        cont, regUp, flex, regDemand, regUpSolar, regUpWind, flexSolar, flexWind = calcWWSISReserves(windGenRegion, solarGenRegion, demand, regLoadFrac,
                                                                                                     contLoadFrac, regErrorPercentile, flexErrorPercentile)
        regUpInc, flexInc = getIncResForAddedRE(newCfs, regErrorPercentile, flexErrorPercentile)
    else:
        cont, regUp, flex, regDemand, regUpSolar, regUpWind, flexSolar, flexWind, regUpInc, flexInc = createEmptyReserveDfs(windGenRegion, newCfs)

    # Get timeseries hours for CE (demand, wind, solar, new wind, new solar, reserves) & save dfs
    (demandCE, windGenCE, solarGenCE, newCfsCE, contCE, regUpCE, flexCE, regUpIncCE, 
        flexIncCE, forsCE, forsTechsCE, capDeratesCE, capDeratesTechsCE, lineDeratesCE) = isolateDataInCEHours(hoursForCE, 
        demand, windGenRegion, solarGenRegion, newCfs, cont, regUp, flex, regUpInc, flexInc, fors, forsTechs, capDerates, capDeratesTechs, lineDerates)
    # Get total hydropower generation potential by block for CE
    [hydroGenCE] = isolateDataInCEBlocks(hoursForCE,hydroGen)

    # Save CE inputs
    for df, n in zip([windGen, solarGen, windGenRegion, solarGenRegion, newCfs, demand, netDemand, cont, regUp, flex, regUpInc, flexInc, regDemand, regUpSolar, regUpWind, flexSolar, flexWind, hydroGen, fors, forsTechs, capDerates, capDeratesTechs, lineDerates],
                     ['windGen','solarGen','windGenRegion','solarGenRegion','windSolarNewCFs','demand','netDemand','contRes','regUpRes','flexRes','regUpInc','flexInc','regUpDemComp','regUpSolComp','regUpWinComp','flexSolComp','flexWinComp','hydroGen','fors','forTechs','capDerates','capDeratesTechs','lineDerates']):
        df.to_csv(os.path.join(resultsDir, n + 'FullYr' + str(currYear) + '.csv'))
    for df, n in zip([demandCE, windGenCE, solarGenCE, newCfsCE, newTechsCE, contCE, regUpCE, flexCE, regUpIncCE, flexIncCE, hydroGenCE, forsCE, forsTechsCE, capDeratesCE, capDeratesTechsCE, lineDeratesCE, hoursForCE],
                     ['demand', 'windGen', 'solarGen','windAndSolarNewCFs','newTechs','contRes','regUpRes','flexRes','regUpInc','flexInc','hydroGen','fors','forTechs','capDerates','capDeratesTechs','lineDerates','hoursByBlock']):
        df.to_csv(os.path.join(resultsDir, n + 'CE' + str(currYear) + '.csv'))
    for scalar, n in zip([windFOR,solarFOR,prmEligWindSolar,planningReserve],['windFOR','solarFOR','windSolarPRMElig','planningReserveCE']): 
        pd.Series(scalar).to_csv(os.path.join(resultsDir,n + str(currYear) + '.csv'))   
    pd.DataFrame([[k, v] for k, v in socScalars.items()],columns=['block','scalar']).to_csv(os.path.join(resultsDir,'socScalarsCE' + str(currYear) + '.csv'))

    # ###############SET UP CAPACITY EXPANSION
    # Create GAMS workspace and database. Parameters and sets are put into database below.
    ws, db, gamsFileDir = createGAMSWorkspaceAndDatabase(runOnSC)

    # Write .gms files for inclusion in CE. These files vary based on CE parameters, e.g. treatment of time.
    writeTimeDependentConstraints(blockNamesChronoList, stoInCE, seasStoInCE, gamsFileDir, ceOps, lastRepBlockNames, specialBlocksPrior, removeHydro)
    writeBuildVariable(ceOps, gamsFileDir)

    # Enter sets and parameters into database
    genSet, hourSet, hourSymbols, zoneOrder, lineSet, zoneSet = edAndUCSharedFeatures(db, genFleetForCE, hoursForCE, demandCE, contCE,regUpCE,flexCE,
                                                                             demandShifter, demandShiftingBlock, rrToRegTime, rrToFlexTime, rrToContTime,
                                                                             solarGenCE, windGenCE, transRegions, lineLimits, transmissionEff, capDeratesCE,
                                                                             lineDeratesCE)  
    stoGenSet, stoGenSymbols = storageSetsParamsVariables(db, genFleetForCE, stoMkts, stoFTLabels)
    stoTechSet, stoTechSymbols = ceSharedFeatures(db, planningReserveHour, genFleetForCE, newTechsCE, planningReserve, discountRate, currCo2Cap,
                                      genSet, hourSet, hourSymbols, newCfsCE, maxCapPerTech, maxREPerZone, maxREInInterconn, regUpIncCE, flexIncCE, 
                                      stoMkts, lineDists, lineCosts, lineSet, zoneOrder, ceOps, interconn, stoFTLabels, zoneSet, 
                                      forsCE, forsTechsCE, capDeratesTechsCE, windFOR, solarFOR, prmEligWindSolar, wsGenFracOfDemand)
    if ceOps == 'UC': ucFeatures(db, genFleetForCE, genSet),
    ceTimeDependentConstraints(db, hoursForCE, blockWeights, socScalars, ceOps, onOffInitialEachPeriod, 
                genSet, genFleetForCE, stoGenSet,stoGenSymbols, newTechsCE, stoTechSet, stoTechSymbols, initSOCFraction,
                hydroGenCE, zoneSet)

    # Run CE model
    print('Running CE for ' + str(currYear))
    capacExpModel, ms, ss = runGAMS('CEWith{o}.gms'.format(o=ceOps), ws, db)

    # ########## SAVE AND PROCESS CE RESULTS
    pd.Series([ms,ss],index=['ms','ss']).to_csv(os.path.join(resultsDir, 'msAndSsCE' + str(currYear) + '.csv'))
    saveCapacExpOperationalData(capacExpModel, genFleetForCE, newTechsCE, hoursForCE, transRegions, lineLimits, resultsDir, 'CE', currYear)
    newGens,newStoECap,newStoPCap,newLines = saveCEBuilds(capacExpModel, resultsDir, currYear)
    genFleet = addNewGensToFleet(genFleet, newGens, newStoECap, newStoPCap, newTechsCE, currYear)
    lineLimits = addNewLineCapToLimits(lineLimits, newLines)
    genFleet.to_csv(os.path.join(resultsDir, 'genFleetAfterCE' + str(currYear) + '.csv'))
    lineLimits.to_csv(os.path.join(resultsDir, 'lineLimitsAfterCE' + str(currYear) + '.csv'))

    return (genFleet, genFleetForCE, lineLimits, capacExpModel, hoursForCE)

# ###############################################################################
# ###############################################################################
# ###############################################################################

# ###############################################################################
# ################## GAMS FUNCTIONS #############################################
# ###############################################################################
def createGAMSWorkspaceAndDatabase(runOnSC):
    # currDir = os.getcwd()
    if runOnSC:
        gamsFileDir = 'GAMS'
        gamsSysDir = '/home/mtcraig/gams40_3'
    else:
        gamsFileDir = 'C:\\Users\\mtcraig\\Desktop\\Research\\Models\\MacroCEM\\GAMS'
        gamsSysDir = 'C:\\GAMS\\40'
        # gamsFileDir = r"C:\Users\atpha\Documents\Postdocs\Projects\NETs\Model\EI-CE\GAMS"
        # gamsSysDir = r"C:\GAMS\win64\30.2"
    ws = GamsWorkspace(working_directory=gamsFileDir, system_directory=gamsSysDir)
    db = ws.add_database()
    return ws, db, gamsFileDir

def runGAMS(gamsFilename, ws, db):
    t0 = time.time()
    model = ws.add_job_from_file(gamsFilename)
    opts = GamsOptions(ws)
    opts.defines['gdxincname'] = db.name
    model.run(opts, databases=db)
    ms, ss = model.out_db['pModelstat'].find_record().value, model.out_db['pSolvestat'].find_record().value
    if (int(ms) != 8 and int(ms) != 1) or int(ss) != 1: print('*********Modelstat & solvestat:', ms, ' & ', ss, ' (ms1 global opt, ms8 int soln, ss1 normal)')
    print('Time (mins) for GAMS run: ' + str(round((time.time()-t0)/60)))
    return model, ms, ss

def edAndUCSharedFeatures(db, genFleet, hours, demand, contRes, regUpRes, flexRes, demandShifter, demandShiftingBlock, rrToRegTime, rrToFlexTime,
                          rrToContTime, hourlySolarGen, hourlyWindGen, transRegions, lineLimits, transmissionEff, capDeratesCE, lineDeratesCE, cnse=10000, co2Price=0):
    # SETS
    genSet = addGeneratorSets(db, genFleet)
    hourSet, hourSymbols = addHourSet(db, hours)
    zoneSet,zoneSymbols,zoneOrder = addZoneSet(db, transRegions)
    lineSet,lineSymbols = addLineSet(db, lineLimits)

    # PARAMETERS
    # Demand and reserves
    addDemandParam(db, demand, hourSet, zoneSet, demandShifter, demandShiftingBlock, mwToGW)
    addReserveParameters(db, contRes, regUpRes, flexRes, rrToRegTime, rrToFlexTime, rrToContTime, hourSet, zoneSet, mwToGW)

    # CO2 cap or price
    addCo2Price(db, co2Price)

    # Generators
    addGenParams(db, genFleet, genSet, mwToGW, lbToShortTon, zoneOrder)
    addCapacityDerates(db, genSet, hourSet, capDeratesCE)
    addExistingRenewableMaxGenParams(db, hourSet, zoneSet, hourlySolarGen, hourlyWindGen, mwToGW)
    addSpinReserveEligibility(db, genFleet, genSet)
    addCostNonservedEnergy(db, cnse)

    # Transmission lines
    addLineParams(db,lineLimits, transmissionEff, lineSet, zoneOrder, mwToGW)
    addLineDerates(db, lineSet, hourSet, lineDeratesCE)
    return genSet, hourSet, hourSymbols, zoneOrder, lineSet, zoneSet

def storageSetsParamsVariables(db, genFleet, stoMkts, stoFTLabels):
    (stoGenSet, stoGenSymbols) = addStoGenSets(db, genFleet, stoFTLabels)
    addStorageParams(db, genFleet, stoGenSet, stoGenSymbols, mwToGW, stoMkts)
    return stoGenSet, stoGenSymbols

def ed(db, socInitial, stoGenSet):
    addStorageInitSOC(db, socInitial, stoGenSet, mwToGW)

def ucFeatures(db, genFleet, genSet):
    addGenUCParams(db, genFleet, genSet, mwToGW)
    
def uc(db, stoGenSet, genSet, socInitial, onOffInitial, genAboveMinInitial, mdtCarriedInitial):
    addStorageInitSOC(db, socInitial, stoGenSet, mwToGW)
    addEguInitialConditions(db, genSet, onOffInitial, genAboveMinInitial, mdtCarriedInitial, mwToGW)

def ceSharedFeatures(db, planningReserveHour, genFleet, newTechs, planningReserve, discountRate, co2Cap, 
        genSet, hourSet, hourSymbols, newCfs, maxCapPerTech, maxREPerZone, maxREInInterconn, regUpInc, 
        flexInc, stoMkts, lineDists, lineCosts, lineSet, zoneOrder, ceOps, interconn, stoFTLabels, zoneSet,
        forsCE, forsTechsCE, capDeratesTechsCE, windFOR, solarFOR, prmEligWindSolar, wsGenFracOfDemand):
    # Sets
    addPlanningReserveHourSubset(db, planningReserveHour)
    addStorageSubsets(db, genFleet, stoFTLabels)
    (techSet, renewTechSet, stoTechSet, stoTechSymbols, thermalSet, dacsSet, CCSSet) = addNewTechsSets(db, newTechs)

    # Long-term planning parameters
    addPlanningReserveParam(db, planningReserve, mwToGW)
    addDiscountRateParam(db, discountRate)
    addCO2Cap(db, co2Cap)

    # New tech parameters
    addGenParams(db, newTechs, techSet, mwToGW, lbToShortTon, zoneOrder, True)
    addCapacityDerates(db, techSet, hourSet, capDeratesTechsCE, True)
    addHourlyFORs(db, forsCE, genSet, hourSet)
    addHourlyFORs(db, forsTechsCE, techSet, hourSet, True)
    addFORScalars(db, windFOR, solarFOR, prmEligWindSolar)
    addTechCostParams(db, newTechs, techSet, stoTechSet, mwToGW)
    addRenewTechCFParams(db, renewTechSet, hourSet, newCfs)
    addMaxNewBuilds(db, newTechs, thermalSet, stoTechSet, dacsSet, CCSSet, zoneSet, maxCapPerTech, maxREPerZone, maxREInInterconn, mwToGW)
    addWSMinGen(db, wsGenFracOfDemand)
    if ceOps == 'UC': addGenUCParams(db, newTechs, techSet, mwToGW, True)
    addResIncParams(db, regUpInc, flexInc, renewTechSet, hourSet)
    addSpinReserveEligibility(db, newTechs, techSet, True)
    addStorageParams(db, newTechs, stoTechSet, stoTechSymbols, mwToGW, stoMkts, True)
    addNewLineParams(db, lineDists, lineCosts, lineSet, maxCapPerTech, zoneOrder, interconn, mwToGW)
    return stoTechSet, stoTechSymbols

def ceTimeDependentConstraints(db, hoursForCE, blockWeights, socScalars, ceOps, onOffInitialEachPeriod,
        genSet, genFleet, stoGenSet, stoGenSymbols, newTechs, stoTechSet, stoTechSymbols, 
        initSOCFraction, hydroGenCE, zoneSet):
    addHourSubsets(db, hoursForCE)
    addSeasonDemandWeights(db, blockWeights)
    addBlockSOCScalars(db, socScalars)
    if ceOps == 'UC': addInitialOnOffForEachBlock(db, onOffInitialEachPeriod, genSet)
    addStoInitSOCCE(db, genFleet, stoGenSet, stoGenSymbols, mwToGW, initSOCFraction)
    addStoInitSOCCE(db, newTechs, stoTechSet, stoTechSymbols, mwToGW, initSOCFraction, True)
    addHydroGenLimits(db, hydroGenCE, zoneSet, mwToGW)

# ###############################################################################
# ###############################################################################
# ###############################################################################
