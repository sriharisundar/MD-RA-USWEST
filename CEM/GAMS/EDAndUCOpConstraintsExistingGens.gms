*Michael Craig 14 May 2020

Positive Variables
*Nonserved energy (GWh)
                 vNse(h)
;

Equations
*Meet demand and reserves
         meetdemand(h,z)                   meet supply with demand
         meetflexreserve(h)              meet spin reserve requirement
         meetcontreserve(h)              meet contingency reserve requirement
         meetregupreserve(h)             meet reg up reserve requirement
*Storage operational constraints
                 genplusuprestosoc(storageegu,h)                limit generation plus up reserves based on SOC
         defstateofcharge(storageegu,h)          define relationship between state of charge and charging and discharging for storage units
                ;

******************SYSTEM-WIDE GENERATION AND RESERVE CONSTRAINTS*******************
*Demand = generation by new and existing plants
meetdemand(h,z)..          sum(egu$[pEguzone(egu)=ORD(z)],vGen(egu,h)) + vNse(h) + sum(l$[pLinesink(l)=ORD(z)],vLineflow(l,h)) =g= (pDemand(h,z) + vShiftedDemand(h,z)) + sum(storageegu$[pEguzone(storageegu)=ORD(z)],vCharge(storageegu,h)) + sum(l$[pLinesource(l)=ORD(z)],vLineflow(l,h));

*Meet spinning and regulation reserve requirements
meetflexreserve(h)..     sum(egu,vFlex(egu,h)) =g= pFlexreserves(h);
meetcontreserve(h)..     sum(egu,vCont(egu,h)) =g= pContreserves(h);
meetregupreserve(h)..    sum(egu,vRegup(egu,h)) =g= pRegupreserves(h);
***********************************************************************************

******************STORAGE SOC TRACKER****************
*Limit gen plus up reserves to SOC
genplusuprestosoc(storageegu,h) .. vGen(storageegu,h) + vRegup(storageegu,h) + vFlex(storageegu,h) + vCont(storageegu,h) =l= (vStateofcharge(storageegu,h-1)$[ORD(h)>1]+ pInitsoc(storageegu)$[ORD(h)=1]);

*Link state of charge to charging and discharging
defstateofcharge(storageegu,h) .. vStateofcharge(storageegu,h) =e= (vStateofcharge(storageegu,h-1)$[ORD(h)>1] + pInitsoc(storageegu)$[ORD(h)=1] - 1/sqrt(pEfficiency(storageegu)) * vGen(storageegu,h) + sqrt(pEfficiency(storageegu)) * vCharge(storageegu,h));
****************************************************
