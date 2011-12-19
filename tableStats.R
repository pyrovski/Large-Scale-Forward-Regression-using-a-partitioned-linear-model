#!/usr/bin/env Rscript
a = read.table("table",header=T)

#some examples:

master_gpu_100000 = intersect(which(a$branch == "master"), which(a$cpu_gpu == "gpu"))
master_gpu_100000 = intersect(master_gpu_100000, which(a$SNPs_on_rank_0 == 100000))
summary(a$comp[master_gpu_100000])

master_cpu_100000 = intersect(which(a$branch == "master"), which(a$cpu_gpu == "cpu"))
master_cpu_100000 = intersect(master_cpu_100000, which(a$SNPs_on_rank_0 == 100000))
summary(a$comp[master_cpu_100000])
