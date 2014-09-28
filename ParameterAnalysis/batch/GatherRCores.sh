#!/bin/bash 

#create the file of v1 on MIC
echo "Time, Error, Count, M, Prob" > ../CSV-core/bv1-mic-cores.csv
for core in 60c 120c 180c 240c 
do
    cat ../../exp-cores/par-ompv1-mic-$core.out | while read oneline 
do
	echo $oneline >> ../CSV-core/bv1-mic-cores.csv
done 
done

#create the file of v2 on MIC
echo "Time, Error, Count, M, Prob" > ../CSV-core/bv2-mic-cores.csv
for core in 60c 120c 180c 240c 
do
    cat ../../exp-cores/par-ompv2-mic-$core.out | while read oneline 
do
	echo $oneline >> ../CSV-core/bv2-mic-cores.csv
done 
done


