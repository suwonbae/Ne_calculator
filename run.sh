#!/bin/bash

# change num_chains and num_procs depending on model and resource
num_chains=27
num_procs=544

num=1
while [ $num -le $num_chains ]
do
	echo $num > chain.txt
	ibrun -np $num_procs a.out > stdout.$num &&

	num=$(($num+1))
	rm chain.txt
done
