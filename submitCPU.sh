#!/bin/bash
# $1 = number of CPU cores desired
# $2 = number of nodes desired
# wayness = number of cores per node = (1 if $1 == 1, else 2) = min(2, $1)
wayness=`eval 'python -c "print min($1,2)"'`
#cores=`eval 'python -c "from math import ceil; print int(ceil(4*$1/8.0)*8*3)"'`
cores=`eval 'python -c "from math import ceil; print int($2 * 8)"'`
#! todo create data dir for each run
set -x
name=`date +%Y_%m_%d_%H_%M_%S`
oldDir=`pwd`
mkdir $name
cd $name
env > env
hostname >> info
echo $name >> info
uname -a >> info
git log -n 1 --oneline >> info

qsub -V -cwd -pe $(echo -n $wayness)way $cores -q normal -P hpc -l h_rt=0:12:00 -A TG-ASC100041 ../run.sh -c -o log -e errlog

