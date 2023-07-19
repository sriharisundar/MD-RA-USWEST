#WECC LINES: 500: ALL CAMX LINES, NW:NE. 345: CENTRAL:DSW. 220-287: CENTRAL-NE, CENTRAL-NW
#PER https://atlas.eia.gov/apps/all-energy-infrastructure-and-resources/explore

#Calculate temperature-dependent ampacity ratings of lines, then convert to percent
#reduction in ampacity (and, per P=VI, power rating), which ignores ampacity increases. 
#temps: F.
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