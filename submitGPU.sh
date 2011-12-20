#!/bin/bash
# $1 = number of GPUs desired
# $2 = number of nodes desired
# number of cores to request = 8 * ceil(4 * $1 / 8)
# wayness = number of cores per node
# for 1 core per socket
cores=`eval 'python -c "print int($2 * 8)"'`
#wayness=`eval 'python -c "print min($1,2)"'`
wayness=`eval 'python -c "print int($1/$2)"'`

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
echo 'cores: $1' >> info
echo 'nodes: $2' >> info

qsub -V -cwd -pe $(echo -n $wayness)way $cores -q normal -P gpgpu -l h_rt=0:12:00 -A TG-ASC100041 -o log -e errlog ../run.sh
