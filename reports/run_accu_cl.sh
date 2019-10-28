#!/bin/bash

# Lancement: qsub -o .output -j y run_accu_cl.sh

#$ -N morpheus
#$ -wd /workdir2/auder/morpheus/reports
#$ -m abes
#$ -M benjamin@auder.net
#$ -pe make 10
#$ -l h_vmem=1G
rm -f .output

module load R/3.6.1

N=100
n=1e5
nc=10

for d in 2 5 10; do
	for link in "logit" "probit"; do
		for weights in "1,1,0"; do
			R --slave --args N=$N n=$n nc=$nc d=$d link=$link weights=$weights <accuracy.R >out_${n}_${link}_${d}_${weights} 2>&1
		done
	done
done

#for d in 2 5; do
#	for n in 5000 10000 100000 500000 1000000; do
#		for link in "logit" "probit"; do
#			R --slave --args N=$N n=$n nc=$nc d=$d link=$link <accuracy.R >out_$n$link$d 2>&1
#		done
#	done
#done
