#!/bin/bash

#PBS -l nodes=1:ppn=16,mem=8gb,pmem=512mb
#PBS -j oe

#PBS -o .output
rm -f .output

WORKDIR=/workdir2/auder/morpheus/reports
cd $WORKDIR

module load R

# arg --vanilla maybe possible on cluster
R --slave --args N=1000 nc=16 link=logit <timings.R >out 2>&1
