#!/usr/bin/env Rscript
a = read.table('filteredTable', header=T)
p = order(a$comp)
a = a[p,]
gpu = which(a$cpu_gpu == 'gpu')
cpu = which(a$cpu_gpu == 'cpu')
pdf('lenSNP.pdf')
labcex=1.2
maincex=1.2
plot(a$SNP_len[cpu], a$comp[cpu]/a$SNPs_on_rank_0[cpu], 
    ylim=c(0, max(a$comp[cpu]/a$SNPs_on_rank_0[cpu])), 
    main='Time per SNP vs SNP length', 
    xlab='SNP length', ylab='Computation time per SNP', type='l', col='red', lty=1, lwd=2,
    sub='1M double-precision SNPs, 1 CPU/1 GPU', cex.lab=labcex, cex.main=maincex)
lines(a$SNP_len[gpu], a$comp[gpu]/a$SNPs_on_rank_0[gpu], col='black', lty=2, lwd=2)
legend(x='topleft', legend=c('cpu', 'gpu'), col=c('red', 'black'), lwd=2, lty=1:2)
dev.off()
