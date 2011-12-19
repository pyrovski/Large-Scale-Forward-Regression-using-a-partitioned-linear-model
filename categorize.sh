#!/bin/bash
find . -name 'id_0.log*' -print|sort -n|grep -v '~'>.list
#cpu=`eval "echo $list|grep cpu"`
cat .list|xargs -I{} grep -Hci gpu {}>.gpucount
cat .gpucount|egrep -v ':0'|cut -d':' -f1|cut -d'/' -f2>.gpu
cat .gpucount|grep ':0'|cut -d':' -f1|cut -d'/' -f2>.cpu

echo id branch rank cpu_gpu iterations time comp SNPs_on_rank_0 >table

cat .cpu|xargs -I{} ./entry.sh {} cpu >> table
cat .gpu|xargs -I{} ./entry.sh {} gpu >> table
