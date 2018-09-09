num=1
num_chains=27
num_procs=544

while [ $num -le $num_chains ]
do
	echo $num > chain.txt
	ibrun -np $num_procs a.out > stdout.$num &&

	num=$(($num+1))
done
