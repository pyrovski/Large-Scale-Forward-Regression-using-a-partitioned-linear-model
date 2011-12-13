#!/bin/bash
find . -name 'run*.sh.o*' -print|grep -v '~'>.list
#cpu=`eval "echo $list|grep cpu"`
cat .list|xargs -I{} grep -Hci gpu {}>.gpucount
cat .gpucount|egrep -v ':0'|cut -d':' -f1>.gpu
cat .gpucount|grep ':0'|cut -d':' -f1>.cpu

echo id rank cpu_gpu iterations time comp>table

cat .cpu|xargs -I{} ./entry.sh {} cpu >> table
cat .gpu|xargs -I{} ./entry.sh {} gpu >> table