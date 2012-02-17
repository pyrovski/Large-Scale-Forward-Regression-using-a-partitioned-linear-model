#!/usr/bin/env Rscript
a = read.table('filteredTable', header=T)
gpu = which(a$cpu_gpu == 'gpu')
cpu = which(a$cpu_gpu == 'cpu')
pdf('lenSNP.pdf')
plot(a$SNP_len[cpu], a$comp[cpu]/a$SNPs_on_rank_0[cpu], ylim=c(0, max(a$comp[cpu]/a$SNPs_on_rank_0[cpu])))
points(a$SNP_len[gpu], a$comp[gpu]/a$SNPs_on_rank_0[gpu])
dev.off()