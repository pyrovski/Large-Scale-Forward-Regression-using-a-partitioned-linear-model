#!/bin/bash
# $1 = number of GPUs desired
# number of cores to request = 8 * ceil(4 * $1 / 8)
# wayness = number of cores per node = (1 if $1 == 1, else 2) = min(2, $1)
# for 1 core per socket
cores=`eval 'python -c "from math import ceil; print int(ceil(4.0*$1/8)*8)"'`
wayness=`eval 'python -c "print min($1,2)"'`
#! todo create data dir for each run
set -x
qsub -V -cwd -pe $(echo -n $wayness)way $cores -q normal -P gpgpu -l h_rt=0:12:00 ./run.sh
# development instead of normal for up to 8 nodes, 1hr
