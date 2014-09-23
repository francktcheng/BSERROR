#!/bin/bash 

Ncache=0
Chunck=0


for((Chunck=10; Chunck<=100; Chunck+=10)) 
do
	echo "Ncache,Chunck,Time,Error,Prob" > bserrorV1-120core-Chunck-$Chunck-Data.csv
	for((Ncache=100; Ncache<=60000; Ncache+=100))
	do
		cat ../bserrorV1-120core-Result/par-omp-$Ncache-$Chunck.out | while read oneline
	do 
		echo $oneline >> bserrorV1-120core-Chunck-$Chunck-Data.csv
	done
	done
done




