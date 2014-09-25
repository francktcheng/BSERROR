#!/bin/bash 

Ncache=0
Chunck=0
Nvec=0

#create the file of v2 on MIC
echo "Nvec,Chunck,NN,Time,Error,Prob" > ./CSV/bsV2-Mic-Chunck-0-Data.csv
for((Nvec=8; Nvec<=80; Nvec+=8))
do
	cat ../bsResult-Mis2/par-ompv2-mic-$Nvec-0.out | while read oneline
do 
	echo $oneline >> ./CSV/bsV2-Mic-Chunck-0-Data.csv
done
done

#create the file of v2 on CPU
echo "Nvec,Chunck,NN,Time,Error,Prob" > ./CSV/bsV2-CPU-Chunck-0-Data.csv
for((Nvec=8; Nvec<=80; Nvec+=8))
do
#	cat ../bsResult-Misc/par-ompv2-cpu-$Nvec-1.out >> ./CSV/bsV2-CPU-Chunck-$Chunck-Data.csv
#	printf "\n" >>./CSV/bsV2-CPU-Chunck-1-Data.csv
	cat ../bsResult-Mis2/par-ompv2-cpu-$Nvec-0.out | while read oneline
do 
	echo $oneline >> ./CSV/bsV2-CPU-Chunck-0-Data.csv
done
done


