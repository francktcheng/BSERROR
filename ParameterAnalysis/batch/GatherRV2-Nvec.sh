#!/bin/bash 

#create the file of v1 on MIC
printf "Nvec, Time,Error,Count, M, Prob\n" > ../CSV-V2-Vec/bsV2-Nvec.csv

for((Nvec=8; Nvec<=480; Nvec+=8)) 
do
	printf ""$Nvec", " >> ../CSV-V2-Vec/bsV2-Nvec.csv
	cat ../../exp-v2-nvec/par-ompv2-mic-Nvec-$Nvec.out >> ../CSV-V2-Vec/bsV2-Nvec.csv   

done


