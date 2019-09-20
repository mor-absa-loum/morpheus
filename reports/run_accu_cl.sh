#!/bin/bash

#PBS -l nodes=1:ppn=15,mem=8gb,pmem=512mb
#PBS -j oe

#PBS -o .output
rm -f .output

WORKDIR=/workdir2/auder/morpheus/reports
cd $WORKDIR

module load R

# arg --vanilla maybe possible on cluster
for d in 2 5 10 20; do
	for link in "logit" "probit"; do
		R --slave --args N=1000 n=1e5 nc=15 d=$d link=$link <accuracy.R >out$d$link 2>&1
	done
done

#for d in 2 5; do
#	for n in 5000 10000 100000 500000 1000000; do
#		for link in "logit" "probit"; do
#			R --slave --args N=1000 n=$n nc=64 d=$d link=$link <accuracy.R >out_$n$link$d 2>&1
#		done
#	done
#done
