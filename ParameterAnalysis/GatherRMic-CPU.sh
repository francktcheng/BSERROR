#!/bin/bash 

Ncache=0
Chunck=0
Nvec=0

#create the file of v1 on MIC
for((Chunck=10; Chunck<=100; Chunck+=10)) 
do
	echo "Ncache,Chunck,NN,Time,Error,Prob" > ./CSV/bsV1-Mic-Chunck-$Chunck-Data.csv
	for((Ncache=100; Ncache<=60000; Ncache+=100))
	do
		cat ../bsResult-mic-cpu/par-ompv1-mic-$Ncache-$Chunck.out | while read oneline
	do 
		echo $oneline >> ./CSV/bsV1-Mic-Chunck-$Chunck-Data.csv
	done
	done
done

#create the file of v1 on CPU
for((Chunck=10; Chunck<=100; Chunck+=10)) 
do
	echo "Ncache,Chunck,NN,Time,Error,Prob" > ./CSV/bsV1-CPU-Chunck-$Chunck-Data.csv
	for((Ncache=100; Ncache<=60000; Ncache+=100))
	do
		cat ../bsResult-mic-cpu/par-ompv1-cpu-$Ncache-$Chunck.out >> ./CSV/bsV1-CPU-Chunck-$Chunck-Data.csv
		printf "\n" >>./CSV/bsV1-CPU-Chunck-$Chunck-Data.csv
	done
done

#create the file of v2 on MIC
for((Chunck=10; Chunck<=100; Chunck+=10)) 
do
	echo "Nvec,Chunck,NN,Time,Error,Prob" > ./CSV/bsV2-Mic-Chunck-$Chunck-Data.csv
	for((Nvec=8; Nvec<=4800; Nvec+=8))
	do
		cat ../bsResult-mic-cpu/par-ompv2-mic-$Nvec-$Chunck.out | while read oneline
	do 
		echo $oneline >> ./CSV/bsV2-Mic-Chunck-$Chunck-Data.csv
	done
	done
done

#create the file of v2 on CPU
for((Chunck=10; Chunck<=100; Chunck+=10)) 
do
	echo "Nvec,Chunck,NN,Time,Error,Prob" > ./CSV/bsV2-CPU-Chunck-$Chunck-Data.csv
	for((Nvec=8; Nvec<=4800; Nvec+=8))
	do
		cat ../bsResult-mic-cpu/par-ompv2-cpu-$Nvec-$Chunck.out >> ./CSV/bsV2-CPU-Chunck-$Chunck-Data.csv
		printf "\n" >>./CSV/bsV2-CPU-Chunck-$Chunck-Data.csv
	done
done


