num=1
num_chain=27
while [ $num -le $num_chain ]
do
	echo $num > chain.txt
	ibrun -np 544 a.out > stdout.$num &&

	num=$(($num+1))
done
