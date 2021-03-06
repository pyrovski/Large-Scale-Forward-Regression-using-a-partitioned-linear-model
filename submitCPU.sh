# 
# 
# Copyright (c) 2011, The Arizona Board of Regents on behalf of 
# The University of Arizona
# 
# All rights reserved.
# 
# Developed by Peter Bailey, Tapasya Patki, and Greg Striemer with
# support from the iPlant Collaborative as a collaboration between
# participants at BIO5 at The University of Arizona (the primary hosting
# institution), Cold Spring Harbor Laboratory, The University of Texas
# at Austin, and individual contributors. Find out more at
# http://www.iplantcollaborative.org/.
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are
# met:
# 
#  * Redistributions of source code must retain the above copyright 
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright 
#    notice, this list of conditions and the following disclaimer in the 
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of the iPlant Collaborative, BIO5, The University 
#    of Arizona, Cold Spring Harbor Laboratory, The University of Texas at 
#    Austin, nor the names of other contributors may be used to endorse or 
#    promote products derived from this software without specific prior 
#    written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
# TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
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
#TG-MCB110022
#TG-ASC100041
qsub -V -cwd -pe $(echo -n $wayness)way $cores -q normal -P hpc -l h_rt=0:12:00 -A TG-MCB110022 -o log -e errlog ../run.sh -c > sublog 2>suberrlog

