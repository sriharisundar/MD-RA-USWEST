#!/bin/bash

for interconn in "WECC"; do
    for co2cap in "100"; do
	for wsgen in "45"; do
	    for prm in "13"; do
		for yr in "2016" "2017" "2018" "2019"; do
		    for ret in "0"; do
       			sbatch RunMacroCEMJob.sbat $interconn $co2cap $wsgen $prm $yr $ret
			sleep 15s
		    done
		done
	    done
	done
    done
done
