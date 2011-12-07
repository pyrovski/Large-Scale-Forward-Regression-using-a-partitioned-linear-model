#!/bin/bash
# $1 = number of CPU cores desired
# $2 = number of nodes desired
# wayness = number of cores per node = (1 if $1 == 1, else 2) = min(2, $1)
wayness=`eval 'python -c "print min($1,2)"'`
#cores=`eval 'python -c "from math import ceil; print int(ceil(4*$1/8.0)*8*3)"'`
cores=`eval 'python -c "from math import ceil; print int($2 * 8)"'`
#! todo create data dir for each run
set -x
qsub -V -cwd -pe $(echo -n $wayness)way $cores -q normal -P hpc -l h_rt=0:12:00 -A TG-ASC100041 ./run.sh -c
# development instead of normal for up to 8 nodes, 1hr

