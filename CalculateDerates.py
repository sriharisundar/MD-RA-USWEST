import os, pandas as pd, numpy as np, xarray as xr, csv

#Calculate deratings for THERMAL POWER PLANTS ONLY. Set derates to 0 for other plants.
#Derates are expressed as fraction of capacity; available capacity in GAMS = capacity * (1 - derate)
#Outputs: dfs w/ datetime rows, gen columns
def calculatePlantCapacityDerates(genFleet, newTechsCE, demand, weatherYears, cesmMembers, compressedGens):
    #Initialize derates to zero
    capDerates = pd.DataFrame(0,index=demand.index,columns=genFleet['GAMS Symbol'])
    capDeratesTechs = pd.DataFrame(0,index=capDerates.index,columns=newTechsCE['GAMS Symbol'])

    #Derate units if running climate change analysis
    if cesmMembers != None: 
        #Derate existing units
        if len(cesmMembers) == 1:
            capDerates = derateExistingUnits(genFleet,capDerates,weatherYears,cesmMembers[0],compressedGens)
        else:
            sys.exit('HAVE NOT TESTED THE MULTIPLE ENSEMBLE MEMBERS IN PLANT DERATES')
            capDeratesAll = list()
            for cesmMember in cesmMembers:
                capDerates = pd.DataFrame(0,index=demand.index,columns=genFleet['GAMS Symbol'])
                capDerates = derateExistingUnits(genFleet,capDerates,weatherYears,cesmMember,compressedGens)
                capDerates.columns = pd.MultiIndex.from_product([[cesmMember], capDerates.columns], names=['ensembleMember', 'gens'])
                capDeratesAll.append(capDerates)
            #Combine into 1 array
            capDerates = pd.concat(capDeratesAll,axis=1)

        #Create new technologies & use regional average derate
        capDeratesTechs,newTechsCE = calculateNewTechDerates(capDerates,genFleet,newTechsCE)
        
    return capDerates,capDeratesTechs,newTechsCE

def calculateLineThermalDerates(lines, demand):
    print('set line deratings to zero')
    lineDerates = pd.DataFrame(0,index=demand.index,columns=lines['GAMS Symbol'])
    return lineDerates

################################################################################
#CREATE DERATES FOR NEW TECHNOLOGIES
#Calculate derates for new techs by taking average derates of all existing generators
#of same plant type & cooling type in region. Need to first diversify derated
#thermal techs by cooling type.
def calculateNewTechDerates(capDerates,genFleet,newTechsCE):
    #Split derated techs by cooling type
    newTechsCE = splitTechsByCT(newTechsCE)
    #Reinitialize derates to zero with new techs
    capDeratesTechs = pd.DataFrame(0,index=capDerates.index,columns=newTechsCE['GAMS Symbol'])
    
    for tech in capDeratesTechs.columns:
        #Get tech row
        techRow = newTechsCE.loc[newTechsCE['GAMS Symbol'] == tech].squeeze()

        #Get tech's plant type, region, and cooling type
        pt,region,cool = techRow['PlantType'],techRow['region'],techRow['CoolingType']

        #Drop CCS from plant type when looking up existing since no CCS exist
        pt = pt.replace(' CCS','') if 'CCS' in pt else pt

        #Skip wind & solar techs, since don't have TDFORs
        if techRowEligForDerates(techRow):
            #Get existing generators of same pt, region, cool type
            gens = genFleet.loc[genFleet['PlantType']==pt]
            gens = gens.loc[gens['region']==region]
            gens = gens.loc[gens['CoolingType']==cool]

            #Get derates of matching existing gens and average them
            if gens.shape[0]>0:
                genDerates = capDerates[gens['GAMS Symbol']]
                capDeratesTechs[tech] = genDerates.mean(axis=1)
            else: #if no other gens of that type in region, drop that gen as a new tech option
                capDeratesTechs.drop(columns=tech,inplace=True)
                newTechsCE.drop(techRow.name,inplace=True)
    newTechsCE.reset_index(drop=True,inplace=True)
    return capDeratesTechs,newTechsCE

def splitTechsByCT(newTechsCE,deratePTs=['Combined Cycle','Coal Steam CCS','Combined Cycle CCS']):
    newTechsCE['CoolingType'] = ''

    #Copy rows that are eligible for derates to make DC cooling
    techDerateRows = newTechsCE.loc[newTechsCE['PlantType'].isin(deratePTs)].copy()
    techDerateRows['CoolingType'],techDerateRows['GAMS Symbol'] = 'DC',techDerateRows['GAMS Symbol']+'DC'

    #Make original rows RC cooling
    newTechsCE.loc[newTechsCE['PlantType'].isin(deratePTs),'CoolingType'] = 'RC'
    newTechsCE.loc[newTechsCE['PlantType'].isin(deratePTs),'GAMS Symbol'] += 'RC'

    #Add in tech derate rows and reset index
    newTechsCE = pd.concat([newTechsCE,techDerateRows],ignore_index=True)
    newTechsCE.sort_values('PlantType',inplace=True)
    newTechsCE.reset_index(drop=True,inplace=True)
    return newTechsCE

def techRowEligForDerates(genRow):
    return ((genRow['PlantType'] in ['Combined Cycle','Coal Steam CCS','Combined Cycle CCS'] and genRow['CoolingType'] in ['RC','DC']))
################################################################################

############### DERATE THERMAL POWER PLANTS ####################################
def derateExistingUnits(genFleet,capDerates,weatherYears,cesmMember,compressedGens):
    #Set recirculating and dry cooling design assumptions
    rcDesign = '95-75' #['90-70','95-75','100-80'] 90-70 more vulnerable than 100-80
    dcDesign = 'ITD 45' #['ITD 25','ITD 35','ITD 45','ITD 55'] #55 more vulnerable than 25

    #Load derate coefficients for coal & NG RC & DC units (using Aviva Loew paper)
    coalRC,coalDC,ngRC,ngDC = loadCoeffs()

    #Load weather
    metVars = xr.open_dataset(os.path.join('Data','CESM','wecc_derate_fields_' + cesmMember + '.nc'))
    metVars = metVars.sel(time=slice(str(weatherYears[0])+"-01-01", str(weatherYears[-1])+"-12-31"))

    #Calculate derates for units that are not a combination of other units
    derates = [calculateGenDerates(metVars,genFleet.iloc[i],rcDesign,dcDesign,coalRC,coalDC,ngRC,ngDC) 
                for i in range(genFleet.shape[0]) if ('COMBINED' not in genFleet.iloc[i]['GAMS Symbol'] and rowEligForDerates(genFleet.iloc[i]))]

    #For units that are a combination of other units, calculate derates for constituent units, then take capacity-weighted average
    combinedRows = genFleet.loc[genFleet['GAMS Symbol'].str.contains('COMBINED')]
    for row in range(combinedRows.shape[0]):
        combinedRow = combinedRows.iloc[row]
        constituentUnits = compressedGens.loc[compressedGens['UnitCompressedInto']==combinedRow['GAMS Symbol']]

        if rowEligForDerates(combinedRow):
            #Calculate derates of each constituent unit
            constituentDerates = [calculateGenDerates(metVars,constituentUnits.iloc[i],rcDesign,dcDesign,coalRC,coalDC,ngRC,ngDC) for i in range(constituentUnits.shape[0])
                                    if rowEligForDerates(constituentUnits.iloc[i])]
            constituentDerates = pd.concat(constituentDerates,axis=1)

            #Take capacity-weighted average 
            wts = constituentUnits['Capacity (MW)']/constituentUnits['Capacity (MW)'].sum()
            constituentDerates = (constituentDerates*wts).sum(axis=1)

            #Save derates
            constituentDerates.name = combinedRow['GAMS Symbol']
            derates.append(constituentDerates)

    #Combine all derates and then replace placeholder derates
    derates = pd.concat(derates,axis=1)
    for c in derates.columns: capDerates[c] = derates[c].values #can't use pd replacement b/c mismatch in date format
    return capDerates

def rowEligForDerates(genRow):
    return ((genRow['PlantType'] in ['Coal Steam','Combined Cycle'] and genRow['CoolingType'] in ['RC','DC']) or (genRow['PlantType'] == 'Combustion Turbine'))

def calculateGenDerates(metVars,genRow,rcDesign,dcDesign,coalRC,coalDC,ngRC,ngDC):
    #Get generator characteristics
    cool,pt,lat,lon = genRow['CoolingType'],genRow['PlantType'],genRow['Latitude'],genRow['Longitude']

    #Get closest coords w/ non-NAN met vars. If closest met var has NANs, step towards center of region and try again
    latCenter,lonCenter = metVars['lat'][int(metVars.dims['lat']/2)].values,metVars['lon'][int(metVars.dims['lon']/2)].values
    latStep,lonStep = abs(metVars['lat'][0].values-metVars['lat'][1].values),abs(metVars['lon'][0].values-metVars['lon'][1].values)
    while np.isnan(metVars.isel(time=-1).sel({'lat':lat,'lon':lon},method='nearest').compute()['TREFHT']):
        lat,lon = lat + ((latCenter-lat)/abs(latCenter-lat))*latStep,lon + ((lonCenter-lon)/abs(lonCenter-lon))*lonStep

    #Get temperatures & dt index from coordinate, then store in Series
    genVars = metVars.sel({'lat':lat,'lon':lon},method='nearest')
    rh = np.array(genVars.variables['RHREFHT'][:]) #relative humidity in %
    t = np.array(genVars.variables['TREFHT'][:]) #temperature in K
    ps = np.array(genVars.variables['PS'][:]) #air pressure in Pa

    #Units needed for regressions: temperature (F), rh (%), air pressure, so convert T from K to C and pressure from Pa to psia.
    t -= 273.15
    ps *= .000145

    #Combine into df
    unitMet = pd.DataFrame({'t':t,'rh':rh,'ps':ps},index=genVars['time']) 

    #Do derating
    if pt == 'Coal Steam':
        if cool == 'RC': availCapac = derateRC(coalRC,unitMet,rcDesign)
        elif cool == 'DC': availCapac = derateDC(coalDC,unitMet,dcDesign)
    elif pt == 'Combined Cycle':
        if cool == 'RC': availCapac = derateRC(ngRC,unitMet,rcDesign)
        elif cool == 'DC': availCapac = derateDC(ngDC,unitMet,dcDesign)
    elif pt == 'Combustion Turbine': availCapac = derateCTs(unitMet)
    
    #Rename series and return fraction of capacity derated (not available capacity)
    if 'GAMS Symbol' in genRow.index: availCapac.name = genRow['GAMS Symbol'] #constituent units do not have GAMS Symbol
    else: availCapac.name = genRow.name
    return 1 - availCapac 

def loadCoeffs():
    baseDir = os.path.join('Data','ThermalDerateCoefficients')
    coalRC = readCSVto2dList(os.path.join(baseDir,'pcrecap.csv'))
    coalDC = readCSVto2dList(os.path.join(baseDir,'pcdccap.csv'))
    ngRC = readCSVto2dList(os.path.join(baseDir,'ngrecap.csv'))
    ngDC = readCSVto2dList(os.path.join(baseDir,'ngdccap.csv'))
    return coalRC,coalDC,ngRC,ngDC 

#Derate GTs (i.e., CTs not at a combined cycle plant) using equation from 
#Bartos and Chester 2015, SI, Equation 48.
def derateCTs(unitMet):
    cct = .0083 #units of 1/degC
    availCapac = -cct*unitMet['t'] + 1.15
    availCapac[availCapac<0],availCapac[availCapac>1] = 0,1
    return availCapac

#Table S10 = coal RC, S15 = NGCC RC. Designs: 90-70, 95-75, 100-80.
def derateRC(coeffs,unitMet,rcDesign):
    varNames = ['Air Temperature (F)','Relative Humidity (%)','Interaction Term',
                'Intercept']
    coeffVals = getCoeffs(coeffs,varNames,rcDesign)
    availCapac = (cToF(unitMet['t'])*coeffVals[0] + unitMet['rh']*coeffVals[1] + 
            cToF(unitMet['t'])*unitMet['rh']*coeffVals[2] + coeffVals[3])
    availCapac[availCapac<0],availCapac[availCapac>1] = 0,1
    return availCapac

#Table S12 = coal DC, S17 = NGCC DC. Designs: ITD 25, ITD 35, ITD 45, ITD 55.
def derateDC(coeffs,unitMet,dcDesign):
    varNames = ['Air Temperature (F)','Air Pressure (psia)','Intercept']
    coeffVals = getCoeffs(coeffs,varNames,dcDesign)
    availCapac = (cToF(unitMet['t'])*coeffVals[0] + unitMet['ps']*coeffVals[1] + coeffVals[2])
    availCapac[availCapac<0],availCapac[availCapac>1] = 0,1
    return availCapac

def cToF(c):
    return c*9/5 + 32

#Extract coeff vals from CSVs
def getCoeffs(coeffs,varNames,design):
    col = coeffs[0].index(design)
    rowLabels = [row[0] for row in coeffs]
    vals = list()
    for v in varNames: vals.append(float(coeffs[rowLabels.index(v)][col].split(' ')[0]))
    return vals    

def readCSVto2dList(fileNameWithDir):
    with open(fileNameWithDir,'r') as f:
        f = csv.reader(f)
        f = list(f)
    return f
################################################################################







################################################################################
#BELOW CODE IS FROM E3/IEC PROJECT
def calcTransmission(temps,fname,outputDir):
    #Interface ratings
    ratings = pd.read_csv('E3TransmissionInterfaceRatings.csv',header=0)
    ratings.index = ratings['zone interface'].str.upper()

    #Convert F to K
    tempsK = ((temps-32)/1.8) + 273

    #Get ratings for each interface
    l1,l2 = list(),list()
    for c in tempsK.columns:
        #Get hourly ampacities at ambient air temperatures
        kV = ratings.loc[c,'kv rating'] #look up kV of line
        amps = calcRatedAmpacity(tempsK[c],kV) #calculate ampacity
        amps.name = c

        #Get ampacities @ standard test conditions
        stdAmps = getStandardTransmissionAmpacities(c,tempsK.index,kV)
        
        #Calculate percent reduction from standard ampacities
        ampRedux = calculateAmpReduction(amps,stdAmps)

        #Save percent reduction from standard ampacities
        l1.append(ampRedux),l2.append(amps)
    ampRedux,amp = pd.concat(l1,axis=1),pd.concat(l2,axis=1)
    ampRedux.to_csv(os.path.join(outputDir,'ampacities','ampacityRedux-' + str(kV) + '-' + fname+'.csv'))
    amp.to_csv(os.path.join(outputDir,'ampacities','ampacity-' + str(kV) + '-' + fname+'.csv'))

    return lineDerates

#Data from Bartos et al., Impacts of rising air temperatures on electric transmission ampacity and peak electricity load in the United States
#https://iopscience.iop.org/article/10.1088/1748-9326/11/11/114008/meta
#SI: https://cfn-live-content-bucket-iop-org.s3.amazonaws.com/journals/1748-9326/11/11/114008/revision1/erl114008_suppdata.pdf?AWSAccessKeyId=AKIAYDKQL6LTV7YY2HIK&Expires=1664123799&Signature=cEIyLjlR1KuqsWZNjAbkVBGGzPQ%3D
#Tamb: K. 
def calcRatedAmpacity(Tamb,kV): 
    a_s = .9 #absorptivity of conductor; Table S3
    delta = 1000 #W/m^2; solar radiation
    Tcond = 75 + 273 #K; max allowable conductor T is 75 C
    eps = .7 #emissivity of conductor; Table S3

    D,resistance = conductorTable()
    resistance = resistance.loc[75,kV] #assume conductor T = 75 C
    D = D.loc[D['V(kV)']==kV,'Diam(m)'].values[0]

    #Get average heat transfer coefficient 
    h = calcH(Tamb,kV,Tcond,D)

    term1 = np.pi * h * D * (Tcond - Tamb) #W/m
    term2 = np.pi * eps * constants.Stefan_Boltzmann * D * (Tcond**4 - Tamb**4) #W/m
    term3 = delta*D*a_s #W/m

    return ((term1+term2-term3)/resistance)**.5 #A

#Calculate heat transfer coefficient using Bartos SI equations. Although SI is ambiguous,
#use film T (not ambient T) for calculations of Pr, v, & k. 
#Output units of h are W/m^2/K.
def calcH(Tamb,kV,Tcond,D):
    Tamb = calcFilmT(Tamb,Tcond)

    #Get dynamic viscosity (v), thermal conductivity (k), and Prandtls number (Pr) as func of T
    df = transmissionTable() #df of T (as index), v (dyn visc), k (therm cond), and Pr (Prandtls)
    #Downscale and linearly interpolate
    df2 = pd.DataFrame(index=list(range(200,601)))
    df2['PrandtlNumber'] = df['PrandtlNumber']
    df2['DynamicViscosity(m^2/s)'] = df['DynamicViscosity(m^2/s)']
    df2['ThermalConductivity(W/m-K)'] = df['ThermalConductivity(W/m-K)']
    df = df2.interpolate()

    #Get k,v,Pr for each T value
    l = [Tamb.copy()]
    for c in df.columns:
        s = Tamb.copy()
        s = s.replace(df[c].to_dict())
        s.name = c
        l.append(s)
    df = pd.concat(l,axis=1)

    Pr,v,k = df['PrandtlNumber'],df['DynamicViscosity(m^2/s)'],df['ThermalConductivity(W/m-K)']
    
    #Define other constants
    V = 0.61 #m/s; wind speed; Table S3
    Re = V*D/v #Reynold's number

    #Return heat transfer coefficient
    term1 = .62 * ((Re)**.5) * (Pr**(1/3))
    term2 = (1+((.4/Pr)**(2/3)))**(1/4)
    term3 = (1+(((Re)/282000)**(5/8)))**(4/5)
    term4 = k/D
    return .3 + term1/term2*term3*term4

#Table S1
def transmissionTable():
    df = pd.DataFrame({'Temperature':list(range(200,601,50)),'DynamicViscosity(m^2/s)':[7.59,11.44,15.89,20.92,26.41,32.39,38.79,45.57,52.69],
        'ThermalConductivity(W/m-K)':[18.1,22.3,26.3,30,33.8,37.3,40.7,43.9,46.9],'PrandtlNumber':[.737,.72,.707,.7,.69,.686,.684,.683,.685]})
    df['DynamicViscosity(m^2/s)'] *= (10**-6)
    df['ThermalConductivity(W/m-K)'] *= (10**-3)
    df.index = df['Temperature']
    return df

#Table S2 - function of ambient air temperature
def conductorTable():
    D = pd.DataFrame({'V(kV)':[500,345,230,115,69],'Diam(m)':[.0304,.0304,.0362,.0277,.0183]}) #convert diam from cm to m so Reynold's # is dimensionless (Eq. 12 of SI)
    #cols = nominal voltage (kV), rows = temperature (C), values = AC resistance (Ohms/km)
    resistance = pd.DataFrame({500:[.061,.067,.073],345:[.061,.067,.073],230:[.044,.048,.052],115:[.073,.08,.087],69:[.17,.186,.203]},index={25,50,75})
    resistance /= 1000 #convert Ohms/km to Ohms/m
    return D,resistance

#Eqn 8 in SI
def calcFilmT(Tcond,Tamb):
    Tfilm = .5*(Tcond+Tamb) 
    Tfilm = Tfilm.astype(int)
    return Tfilm

#Check if running ampacity model @ standard conditions (75 C conductor T, 25 C ambient T)
#yields correct ampacity values (per https://www.prioritywire.com/specs/ACSR.pdf).
#ACSR specs: 69: Linnet, 529A; 115: Condor, 889A; 230: Martin, 1232A; 345 & 500: Cardinal, 995A
#Results:    69:         513A; 115:         871A; 230:         1207A; 345 & 500:           974A
def getStandardTransmissionAmpacities(c,idx,kV,testTK=25+273): #use K temperature
    tempsK = pd.Series(testTK,index=idx,name=c)
    return calcRatedAmpacity(tempsK,kV)

def calculateAmpReduction(amps,stdAmps):
    ampRedux = amps.copy()
    #Get reduction for ampacity decreases
    ampRedux[ampRedux<stdAmps] = (stdAmps - ampRedux)/stdAmps
    #Ignore ampacity increases
    ampRedux[ampRedux>=stdAmps] = 0
    return ampRedux