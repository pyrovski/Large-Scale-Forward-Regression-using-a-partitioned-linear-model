#!/bin/bash
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
cat $1/id_0.log|exRank.sh >> .tmp

#cpu_gpu
echo $2 >> .tmp

#iterations
cat $1/id_0.log|exIter.sh >> .tmp

#time
cat $1/id_0.log|exTotal.sh >> .tmp

#comp time
cat $1/id_0.log|exComp.sh >> .tmp

#SNPs on rank 0
snp=`cat $1/id_0.log|grep 'has SNPs'| cut -d' ' -f5|cut -d'-' -f2`
python -c "print $snp + 1" >> .tmp

#regurgitate
cat .tmp|tr -d ' '|tr '\n' ' '
echo
