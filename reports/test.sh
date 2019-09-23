#!/bin/bash

# arg --vanilla maybe possible on cluster
for d in 2 5; do
	for link in "logit" "probit"; do
		R --slave --args N=10 n=1e3 nc=3 d=$d link=$link <accuracy.R >out$d$link 2>&1
	done
done

#for d in 2 5; do
#	for n in 5000 10000 100000 500000 1000000; do
#		for link in "logit" "probit"; do
#			R --slave --args N=1000 n=$n nc=64 d=$d link=$link <accuracy.R >out_$n$link$d 2>&1
#		done
#	done
#done
