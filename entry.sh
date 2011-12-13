#!/bin/bash
rm -f .tmp
echo $1 >> .tmp
cat $1|exRank.sh >> .tmp
echo $2 >> .tmp
cat $1|exIter.sh >> .tmp
cat $1|exTotal.sh >> .tmp
cat $1|exComp.sh >> .tmp
cat .tmp|tr '\n' ' '
echo
