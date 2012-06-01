#!/bin/bash
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
rm -f .tmp

#id
echo $1 >> .tmp

#branch
cat $1/info | grep '*' | cut -d' ' -f2 >> .tmp

#nodes
cat $1/info | grep 'nodes:'|cut -d':' -f2 >> .tmp

#cores
cat $1/info | grep 'cores:'|cut -d':' -f2 >> .tmp

#rank
cat $1/id_0.log|./exRank.sh >> .tmp

#cpu_gpu
echo $2 >> .tmp

#iterations
cat $1/id_0.log|./exIter.sh >> .tmp

#time
cat $1/id_0.log|./exTotal.sh >> .tmp

#comp time
cat $1/id_0.log|./exComp.sh >> .tmp

#SNPs on rank 0
snp=`cat $1/id_0.log|grep 'has SNPs'| cut -d' ' -f5|cut -d'-' -f2`
python -c "print $snp + 1" >> .tmp

#SNP length
egrep -o 'num_r[[:space:]]+[[:digit:]]+' $1/id_0.log|tr -s ' ' | cut -d' ' -f2 >> .tmp

#regurgitate
cat .tmp|tr -d ' '|tr '\n' ' '
echo
