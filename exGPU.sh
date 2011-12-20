#!/bin/bash
cat|grep 'GPU comp'|grep -v per|cut -d':' -f2|tr -d 's'|Rscript -e 'data=read.table(pipe("cat /dev/stdin"), header=F);cat(sum(data),"\n")'
