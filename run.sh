#!/bin/bash

# num_chains and num_procs according to model and resource
num_chains=27
num_procs=544

# do not change the following lines
num=1
while [ $num -le $num_chains ]
do
	echo $num > chain.txt
	ibrun -np $num_procs a.out > stdout.$num &&

	num=$(($num+1))
	rm chain.txt
done
