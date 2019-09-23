#!/bin/bash

#$ -N morpheus
#$ -m abes
#$ -M benjamin@auder.net
#$ -pe make 5
#$ -l h_vmem=1G
#$ -j y
#$ -o .output
rm -f .output

WORKDIR=/workdir2/auder/morpheus/reports
cd $WORKDIR

module load R/3.6.0

N=1000
n=1e5
nc=50

link=logit
# and disable d=20 to run faster

# arg --vanilla maybe possible on cluster
for d in 2 5 10; do
	#for link in "logit" "probit"; do
		R --slave --args N=$N n=$n nc=$nc d=$d link=$link <accuracy.R >out_$n$link$d 2>&1
	#done
done

#for d in 2 5; do
#	for n in 5000 10000 100000 500000 1000000; do
#		for link in "logit" "probit"; do
#			R --slave --args N=$N n=$n nc=$nc d=$d link=$link <accuracy.R >out_$n$link$d 2>&1
#		done
#	done
#done
