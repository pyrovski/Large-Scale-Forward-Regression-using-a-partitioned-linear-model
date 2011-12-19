#!/bin/bash
rm -f .tmp
echo $1 >> .tmp
cat $1/info | grep '*' | cut -d' ' -f2 >> .tmp
cat $1/id_0.log|exRank.sh >> .tmp
echo $2 >> .tmp
cat $1/id_0.log|exIter.sh >> .tmp
cat $1/id_0.log|exTotal.sh >> .tmp
cat $1/id_0.log|exComp.sh >> .tmp
snp=`cat $1/id_0.log|grep 'has SNPs'| cut -d' ' -f5|cut -d'-' -f2`
python -c "print $snp + 1" >> .tmp

cat .tmp|tr '\n' ' '
echo
