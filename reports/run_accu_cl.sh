#!/bin/bash

#$ -N morpheus
#$ -wd /workdir2/auder/morpheus/reports
#$ -m abes
#$ -M benjamin@auder.net
#$ -pe make 50
rm -f .output
#$ -o .output
#$ -j y

module load R/3.6.1

N=100
nc=50

for n in "5000" "10000" "100000" "500000" "1000000"; do
	for d in 2 5 10; do
		for link in "logit" "probit"; do
			R --slave --args N=$N n=$n nc=$nc d=$d link=$link <accuracy.R >out_${n}_${link}_${d} 2>&1
		done
	done
done
