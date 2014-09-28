#!/bin/bash 

#create the file of v1 on MIC
echo "NP, M, Time, Ninit, Nlast" > ../CSV-Scale/BSMPI-scale.csv
for np in 2 3 4 5 6 7 
do
    cat ../../exp-mpi/par-M-10000-np-$np.out | while read oneline 
do
	echo $oneline >> ../CSV-Scale/BSMPI-scale.csv
done 
done



