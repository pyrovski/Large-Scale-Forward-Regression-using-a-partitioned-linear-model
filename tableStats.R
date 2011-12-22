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
plot(uconf, mGPU_small, log='xy', xlab='log(SNP count)', ylab='log (computation time) (s)', main='computation time vs. SNP count', sub='small datasets', type='l', col='black', ylim=yl, lwd='2', lty=1)
lines(uconf, mGPU_large, col='blue', lwd='2', lty=2)
lines(uconf, mCPU, col='red', lwd='2', lty=3)
legend(x='topleft', legend=c('gpu small', 'gpu large', 'cpu'), col=c('black', 'blue', 'red'), lwd=2, lty=1:3)


iterSel = which(a$iterations == 50)

################################################################################
#
# compare small problem size instances across node counts, implementations
#
################################################################################

SNPSel100k = which(a$SNPs_on_rank_0 == 100000)

sel = mergeSels(list(gpuSel, branchSelGPU_small, iterSel, SNPSel100k))

selL = mergeSels(list(gpuSel, branchSelGPU_large, iterSel, SNPSel100k))

conf = a[sel,c('nodes','cores.gpus')]
p = order(conf$nodes,conf$cores.gpus)
uconf = unique(conf[p,])

m = c()
std = c()
count = c()

mL = c()
stdL = c()
countL = c()

totalSNPs = c()
gpuCount = c()
for(i in 1:length(uconf[,1])){
  if(uconf[i,'nodes'] == uconf[i,'cores.gpus'] && uconf[i,'cores.gpus'] > 1){
    next()
  }

  gpuCount[i] = paste(uconf[i, 'cores.gpus'], 'GPUs')
  
  # summarize matching configurations;
  # they should have equal numbers of total SNPs
  iSel = mergeSels(list(sel, which(a$nodes == uconf[i,'nodes']), which(a$cores.gpus == uconf[i,'cores.gpus'])))
  vals = a$comp[iSel];
  m[i] = mean(vals);
  std[i] = sd(vals);
  count[i] = length(vals)
  
  iSelL = mergeSels(list(selL, which(a$nodes == uconf[i,'nodes']), which(a$cores.gpus == uconf[i,'cores.gpus'])))
  valsL = a$comp[iSelL];
  mL[i] = mean(valsL);
  stdL[i] = sd(valsL);
  countL[i] = length(valsL)

  # this is not entirely precise; the last rank could have fewer than SNPs_on_rank_0 SNPs
#  totalSNPs[i] = paste(format(uconf[i,'cores.gpus'] * 100000, scientific=T), 'SNPs')
}
pdf('smallDataGPUMPI.pdf');
#plot(uconf[,'cores.gpus'], m)
barplot(t(matrix(c(m,mL),ncol=2)), names.arg=gpuCount, beside=T,main=paste('weak scaling across GPUs via MPI'), lwd=2, ylab='computation time (s)', legend.text=c('gpu small','gpu large'), sub='100k SNPs/GPU', ylim=c(0,1.2*max(m,mL)))
print(count)

################################################################################
#
# compare 1 SNP/thread block and 2 SNP/thread block, double and single precision
#
################################################################################

gpuSel = which(a$cpu_gpu == 'gpu')
D2Sel = which(a$branch == 'twoSNP')
D1Sel = which(a$branch == 'master')
S2Sel = which(a$branch == 'twoSNP_single')
S1Sel = which(a$branch == 'single')

# 1M, 1 gpu
SNPSel1 = which(a$SNPs_on_rank_0 == 1000000)

# 1M, 2 gpu2
SNPSel2 = which(a$SNPs_on_rank_0 == 500000)

nodeSel = which(a$nodes == 1)
core1Sel = which(a$cores == 1)
core2Sel = which(a$cores == 2)
D2core1 = mergeSels(list(gpuSel,D2Sel,core1Sel))
D2core2 = mergeSels(list(gpuSel,D2Sel,core2Sel))
D1core1 = mergeSels(list(gpuSel,D1Sel,core1Sel))
D1core2 = mergeSels(list(gpuSel,D1Sel,core2Sel))
S2core1 = mergeSels(list(gpuSel,S2Sel,core1Sel))
S2core2 = mergeSels(list(gpuSel,S2Sel,core2Sel))
S1core1 = mergeSels(list(gpuSel,S1Sel,core1Sel))
S1core2 = mergeSels(list(gpuSel,S1Sel,core2Sel))
