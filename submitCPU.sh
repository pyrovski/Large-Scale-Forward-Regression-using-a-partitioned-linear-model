#!/bin/bash
# $1 = number of CPU cores desired
# $2 = number of nodes desired
# $1 must be a >=1 multiple of $2
# wayness = number of cores per node
wayness=`eval 'python -c "print int($1/$2)"'`
#cores=`eval 'python -c "from math import ceil; print int(ceil(4*$1/8.0)*8*3)"'`
cores=`eval 'python -c "from math import ceil; print int($2 * 8)"'`
#! todo create data dir for each run
set -x

name=`date +%Y_%m_%d_%H_%M_%S_%N`
oldDir=`pwd`
mkdir $name
cd $name
env > env
hostname >> info
echo $name >> info
uname -a >> info
git log -n 1 --oneline >> info
git branch | grep '*' >> info
echo 'cores: '$1 >> info
echo 'nodes: '$2 >> info

#echo $(echo -n $wayness)way $cores
qsub -V -cwd -pe $(echo -n $wayness)way $cores -q normal -P hpc -l h_rt=0:12:00 -A TG-ASC100041 -o log -e errlog ../run.sh -c > sublog 2>suberrlog

