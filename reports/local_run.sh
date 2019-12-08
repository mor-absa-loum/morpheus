#!/bin/bash

N=100
n=1e5
nc=3
nstart=5

for d in 2 5 10; do
	for link in "logit" "probit"; do
		R --slave --args N=$N n=$n nc=$nc d=$d link=$link nstart=$nstart <multistart.R >out_${n}_${link}_${d}_${nstart} 2>&1
	done
done
