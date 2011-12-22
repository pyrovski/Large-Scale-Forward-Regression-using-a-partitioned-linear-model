#!/usr/bin/env Rscript

mergeSels = function(selList){
  result = unlist(selList[1])
  for(i in 2:(length(selList))){
    result = intersect(result, unlist(selList[i]))
  }
  return(result)
}

a = read.table("filteredTable",header=T)
#cat(names(a), '\n')

#some examples:

#master_gpu_100000 = intersect(which(a$branch == "master"), which(a$cpu_gpu == "gpu"))
#master_gpu_100000 = intersect(master_gpu_100000, which(a$SNPs_on_rank_0 == 100000))
#summary(a$comp[master_gpu_100000])

#master_cpu_100000 = intersect(which(a$branch == "master"), which(a$cpu_gpu == "cpu"))
#master_cpu_100000 = intersect(master_cpu_100000, which(a$SNPs_on_rank_0 == 100000))
#summary(a$comp[master_cpu_100000])

master_gpu_1000000_29i = mergeSels(list(which(a$branch == "master"), which(a$cpu_gpu == "gpu"), 
		       which(a$SNPs_on_rank_0 == 1000000), which(a$iterations == 29)))
#summary(a$comp[master_gpu_1000000_29i])


#################################################################
# generate plot for small GPU vs large GPU vs cpu implementations
#################################################################
branchSelGPU_small = which(a$branch == 'test_GPU_all')
branchSelGPU_large = which(a$branch == 'master')
branchSelCPU = which(a$branch == 'master')

gpuSel = which(a$cpu_gpu == 'gpu')
cpuSel = which(a$cpu_gpu == 'cpu')

coreSel = which(a$cores.gpus == 1)
selGPU_small = mergeSels(list(branchSelGPU_small, gpuSel, coreSel))
selGPU_large = mergeSels(list(branchSelGPU_large, gpuSel, coreSel))
selCPU = mergeSels(list(branchSelCPU, cpuSel, coreSel))

conf = a$SNPs_on_rank_0[selGPU_small]
uconf = sort(unique(conf))

# merge entries matching selGPU_small and uconf
cat('branch cores/gpus cpu_gpu count SNPs_on_rank_0 mean_time std_time\n')
mGPU_small = c()
mGPU_large = c()
mCPU = c()

stdGPU_small = c()
stdGPU_large = c()
stdCPU = c()
for(i in 1:length(uconf)){
  valsGPU_small = a$comp[intersect(which(a$SNPs_on_rank_0 == uconf[i]), selGPU_small)]
  mGPU_small[i] = mean(valsGPU_small)
  stdGPU_small[i] = sd(valsGPU_small)

  valsGPU_large = a$comp[intersect(which(a$SNPs_on_rank_0 == uconf[i]), selGPU_large)]
  mGPU_large[i] = mean(valsGPU_large)
  stdGPU_large[i] = sd(valsGPU_large)

  valsCPU = a$comp[intersect(which(a$SNPs_on_rank_0 == uconf[i]), selCPU)]
  mCPU[i] = mean(valsCPU)
  stdCPU[i] = sd(valsCPU)
  
  cat('test_GPU_all 1 gpu', uconf[i], ' ', mGPU_small[i], ' ', stdGPU_small[i], '\n')
  cat('master 1 gpu', uconf[i], ' ', mGPU_large[i], ' ', stdGPU_large[i], '\n')
  cat('master 1 cpu', uconf[i], ' ', mCPU[i], ' ', stdCPU[i], '\n')
}

yl = range(c(mGPU_small, mCPU))
pdf('smallDataGPU.pdf')
p = plot(uconf, mGPU_small, log='xy', xlab='log(SNP count)', ylab='log (computation time) (s)', main='computation time vs. SNP count', sub='small datasets', type='l', col='black', ylim=yl, lwd='2', lty=1)
lines(uconf, mGPU_large, col='blue', lwd='2', lty=2)
lines(uconf, mCPU, col='red', lwd='2', lty=3)
legend(x='topleft', legend=c('gpu small', 'gpu large', 'cpu'), col=c('black', 'blue', 'red'), lwd=2, lty=1:3)
