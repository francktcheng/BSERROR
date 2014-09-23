#!/bin/bash 

Ncache=0
Chunck=0


for((Chunck=10; Chunck<=100; Chunck+=10)) 
do
	echo "Ncache,Chunck,Time,Error,Prob" > bserror-Chunck-$Chunck-Data.csv 
	for((Ncache=10000; Ncache<=20000; Ncache+=100))
	do
		cat ../bserrorTestResult/par-omp-$Ncache-$Chunck.out | while read oneline
	do 
		echo $oneline >> bserror-Chunck-$Chunck-Data.csv
	done
	done
done




