#!/bin/bash

N=100
n=5e3
nc=3
#nstart=5

for n in "5000" "10000" "100000"; do
  for d in 2 5 10; do
    for link in "logit" "probit"; do
      #R --slave --args N=$N n=$n nc=$nc d=$d link=$link nstart=$nstart <multistart.R >out_${n}_${link}_${d}_${nstart} 2>&1
      R --slave --args N=$N n=$n nc=$nc d=$d link=$link <accuracy.R >out_${n}_${link}_${d} 2>&1
    done
  done
done
