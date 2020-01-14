#!/bin/bash

#$ -N morpheus
#$ -wd /workdir2/auder/morpheus/reports
#$ -m abes
#$ -M benjamin@auder.net
#$ -pe make 50
rm -f .output
#$ -o .output
#$ -j y

module load R

N=1000
nc=50
link="logit"

# arg --vanilla maybe possible on cluster
R --slave --args N=$N nc=$nc link=$link <timings.R >out 2>&1
