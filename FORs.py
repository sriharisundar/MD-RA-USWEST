#Calculate temperature-dependent FORs at thermal generators

fors = pd.read_excel(forFile,index_col=0,header=0)    

def calcTDFORs(fors,fleet,temps,fname,outputDir,aggregateGens):  
    #Define whether working w/ counties or zones
    region = 'Zone' if aggregateGens else 'county'

    #Convert to C
    temps = (temps-32)/1.8

    #Round temperatures to nearest 5th degree to match TDFORs
    temps = 5*((temps/5).round())
    temps[temps>35] = 35
    temps[temps<-15] = -15
    
    #For each generator w/ a TDFOR, replace temperatures w/ TDFOR
    tdfors = pd.DataFrame(index=temps.index,columns=fleet['line ref no'])
    for c in tdfors.columns:      
        forType,location = fleet.loc[fleet['line ref no']==c,'forType'].values[0],fleet.loc[fleet['line ref no']==c,region].values[0]
        locTemps = temps[location]
        if forType in fors.columns:
            tToFOR = fors[forType].to_dict()
            tdfors[c] = locTemps.replace(tToFOR)
    
    #Save file
    tdfors.to_csv(os.path.join(outputDir,'fors','fors-'+('agg-' if aggregateGens else '')+fname+'.csv'))