# 
# 
# Copyright (c) 2011, The Arizona Board of Regents on behalf of 
# The University of Arizona
# 
# All rights reserved.
# 
# Developed by Peter Bailey, Tapasya Patki, and Greg Striemer with
# support from the iPlant Collaborative as a collaboration between
# participants at BIO5 at The University of Arizona (the primary hosting
# institution), Cold Spring Harbor Laboratory, The University of Texas
# at Austin, and individual contributors. Find out more at
# http://www.iplantcollaborative.org/.
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are
# met:
# 
#  * Redistributions of source code must retain the above copyright 
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright 
#    notice, this list of conditions and the following disclaimer in the 
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of the iPlant Collaborative, BIO5, The University 
#    of Arizona, Cold Spring Harbor Laboratory, The University of Texas at 
#    Austin, nor the names of other contributors may be used to endorse or 
#    promote products derived from this software without specific prior 
#    written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
# TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
#!/usr/bin/env Rscript

mergeSels = function(selList){
  result = unlist(selList[1])
  for(i in 2:(length(selList))){
    result = intersect(result, unlist(selList[i]))
  }
  return(result)
}


labcex=1.2
maincex=1.2

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
plot(uconf, mGPU_small, log='xy', xlab='log(SNP count)', ylab='log (computation time) (s)', main='Computation time vs. SNP count', sub='Small datasets', type='l', col='black', ylim=yl, lwd='2', lty=1, cex.lab=labcex, cex.main=maincex)
lines(uconf, mGPU_large, col='blue', lwd='2', lty=2)
lines(uconf, mCPU, col='red', lwd='2', lty=3)
legend(x='topleft', legend=c('gpu small', 'gpu large', 'cpu'), col=c('black', 'blue', 'red'), lwd=2, lty=1:3)


iterSel = which(a$iterations == 50)

################################################################################
#
# Compare small problem size instances across node counts, implementations.
# Weak scaling.
#
################################################################################

SNPSel100k = which(a$SNPs_on_rank_0 == 100000)

sel = mergeSels(list(gpuSel, branchSelGPU_small, iterSel, SNPSel100k))
selL = mergeSels(list(gpuSel, branchSelGPU_large, iterSel, SNPSel100k))
selC = mergeSels(list(cpuSel, branchSelGPU_large, iterSel, SNPSel100k))

conf = a[sel,c('nodes','cores.gpus')]
p = order(conf$nodes,conf$cores.gpus)
uconf = unique(conf[p,])

m = c()
std = c()
count = c()

mL = c()
stdL = c()
countL = c()

mC = c()
stdC = c()
countC = c()

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

  iSelC = mergeSels(list(selC, which(a$nodes == uconf[i,'nodes']), which(a$cores.gpus == uconf[i,'cores.gpus'])))
  valsC = a$comp[iSelC];
  mC[i] = mean(valsC);
  stdC[i] = sd(valsC);
  countC[i] = length(valsC)

  # this is not entirely precise; the last rank could have fewer than SNPs_on_rank_0 SNPs
#  totalSNPs[i] = paste(format(uconf[i,'cores.gpus'] * 100000, scientific=T), 'SNPs')
}
pdf('smallDataGPUMPI.pdf');
#plot(uconf[,'cores.gpus'], m)
data = c(mC/m,mC/mL)*100.0 - 100.0
# gpu only vs selective
barplot(t(matrix(data,ncol=2)),
	names.arg=gpuCount, beside=T,
	main=paste('Weak scaling across GPUs via MPI'), 
	lwd=2, 
	ylab='Improvement over CPU MPI (%)', 
	legend.text=c('GPU only','Selective GPU'), 
	sub='100k SNPs/GPU', 
	ylim=c(0,1.2*max(data)), 
	cex.lab=labcex, 
	cex.main=maincex)
cat('weak scaling counts:\n')
print(count)
print(countL)
print(countC)

################################################################################
#
# Compare small problem size instances across node counts, implementations.
# Strong scaling.
#
################################################################################

sel = mergeSels(list(gpuSel, branchSelGPU_small, iterSel))
selL = mergeSels(list(gpuSel, branchSelGPU_large, iterSel))
selC = mergeSels(list(cpuSel, branchSelGPU_large, iterSel))

SNPSel = which(a$SNPs_on_rank_0 * a$cores.gpus == 100000)

conf = a[sel,c('nodes','cores.gpus')]
p = order(conf$nodes,conf$cores.gpus)
uconf = unique(conf[p,])

m = c()
std = c()
count = c()
cih = c()
cil = c()

mL = c()
stdL = c()
countL = c()
cihL = c()
cilL = c()

mC = c()
stdC = c()
countC = c()
cihC = c()
cilC = c()

totalSNPs = c()
gpuCount = c()
for(i in 1:length(uconf[,1])){
  if(uconf[i,'nodes'] == uconf[i,'cores.gpus'] && uconf[i,'cores.gpus'] > 1){
    next()
  }

  gpuCount[i] = paste(uconf[i, 'cores.gpus'], 'GPUs')
  
  # summarize matching configurations;
  # they should have equal numbers of total SNPs
  uconfSel = mergeSels(list(which(a$nodes == uconf[i,'nodes']), which(a$cores.gpus == uconf[i,'cores.gpus']), SNPSel))

  iSel = intersect(sel, uconfSel)
  vals = a$comp[iSel];
  m[i] = mean(vals);
  std[i] = sd(vals);
  count[i] = length(vals)

  iSelL = intersect(selL, uconfSel)
  valsL = a$comp[iSelL];
  mL[i] = mean(valsL);
  stdL[i] = sd(valsL);
  countL[i] = length(valsL)

  iSelC = intersect(selC, uconfSel)
  valsC = a$comp[iSelC];
  mC[i] = mean(valsC);
  stdC[i] = sd(valsC);
  countC[i] = length(valsC)

  tmp = t.test(valsC)
  cilC[i] = tmp$conf.int[1]
  cihC[i] = tmp$conf.int[2]

  tmp = t.test(mC[i]/vals*100-100)
  cil[i] = tmp$conf.int[1]
  cih[i] = tmp$conf.int[2]

  tmp = t.test(mC[i]/valsL*100-100)
  cilL[i] = tmp$conf.int[1]
  cihL[i] = tmp$conf.int[2]
  cat('test_GPU_all ', uconf[i,'cores.gpus'], ' gpus ', m[i], ' ', std[i], '\n')
  cat('master ', uconf[i,'cores.gpus'], ' gpus ', mL[i], ' ', stdL[i], '\n')
  cat('master ', uconf[i,'cores.gpus'], ' cpus ', mC[i], ' ', stdC[i], '\n')
}
pdf('smallDataGPUMPIStrong.pdf');
#plot(uconf[,'cores.gpus'], m)
data = c(mC/m,mC/mL)*100.0 - 100.0
barplot(t(matrix(data,ncol=2)),
	 names.arg=gpuCount, 
	 beside=T,
	 main=paste('Strong scaling across GPUs via MPI (small datasets)'), 
	 lwd=2, 
	 ylab='Improvement over CPU MPI (%)', 
	 legend.text=c('GPU only','selective GPU'), 
	 sub='100k SNPs total; 2 CPUs/GPUs per node', 
	 ylim=c(1.2*min(cil,cilL),1.2*max(data)), 
	 cex.lab=labcex, 
	 cex.main=maincex)

# add error bars
#for(i in 1:length(uconf[,1])){
#  lines(3*(c(i,i)-.5), c(cil[i],cih[i]), lwd=3, col='red')
#  lines(3*(i-.5)+c(-.3,.3), c(cil[i],cil[i]),lwd=3, col='red')
#  lines(3*(i-.5)+c(-.3,.3), c(cih[i],cih[i]),lwd=3, col='red')
#
#  lines(3*(c(i,i)-.5)+1, c(cilL[i],cihL[i]), lwd=3, col='red')
#  lines(3*(i-.5)+1+c(-.3,.3), c(cilL[i],cilL[i]),lwd=3, col='red')
#  lines(3*(i-.5)+1+c(-.3,.3), c(cihL[i],cihL[i]),lwd=3, col='red')
#}
cat('strong scaling counts:\n')
print(count)
print(countL)
print(countC)

################################################################################
#
# Summarize small problem size instances across node counts for GPU only.
# Single precision.
# Strong scaling.
#
################################################################################

SNPSel = which(a$SNPs_on_rank_0 * a$cores.gpus == 100000)
confSel = union(which(a$cores.gpus == 1), which(a$cores.gpus == 2*a$nodes))
branchSelGPUSmallSingle = which(a$branch == 'GPU_all_single')
sel = mergeSels(list(gpuSel, branchSelGPUSmallSingle, iterSel, SNPSel, confSel))
conf = a[sel,c('nodes','cores.gpus')]
p = order(conf$nodes,conf$cores.gpus)
uconf = unique(conf[p,])

m = c()
sockets = c()
count = c()
cih = c()
cil = c()
for(i in 1:length(uconf[,1])){
  uconfSel = mergeSels(list(which(a$nodes == uconf[i,'nodes']), 
  	     which(a$cores.gpus == uconf[i,'cores.gpus']), sel))

  vals = a$comp[uconfSel]
  m[i] = mean(vals)
  count[i] = length(vals)
  sockets[i] = uconf[i,'cores.gpus']
}
#pdf('cpuStrong.pdf')
#plot(sockets, m, type='l', lwd=2, main='CPU strong scaling via MPI',   sub='1M SNPs', xlab='Cores (2 cores/node)', ylab='Computation time (s)', log='xy', axes=F, cex.lab=labcex, cex.main=maincex)
#axis(at=sockets, side=1)
#axis(at=m, labels=format(m,digits=3), side=2)
#box()

#points(sockets, m, pch=19)

#pdf('cpuStrongSNPsPerSecond.pdf')
#plot(sockets, 1000000/m, type='l', lwd=2, main='CPU strong scaling via MPI',  sub='1M SNPs', xlab='Cores (2 cores/node)', ylab='SNPs/s', log='xy', #cex.lab=labcex, cex.main=maincex, axes=F)
#axis(at=sockets, side=1)
#axis(at=1000000/m, labels=format(m,digits=3), side=2)
#box()
#points(sockets, 1000000/m, pch=19)

cat('GPU small strong scaling counts:\n')
print(count)
print(sockets)
print(m)

################################################################################
#
# Compare large problem size instances across node counts for CPU only.
# Strong scaling.
#
################################################################################

SNPSel = which(a$SNPs_on_rank_0 * a$cores.gpus == 1000000)
confSel = union(which(a$cores.gpus == 1), which(a$cores.gpus == 2*a$nodes))
sel = mergeSels(list(cpuSel, branchSelCPU, iterSel, SNPSel, confSel))
conf = a[sel,c('nodes','cores.gpus')]
p = order(conf$nodes,conf$cores.gpus)
uconf = unique(conf[p,])

m = c()
sockets = c()
count = c()
cih = c()
cil = c()
for(i in 1:length(uconf[,1])){
  uconfSel = mergeSels(list(which(a$nodes == uconf[i,'nodes']), 
  	     which(a$cores.gpus == uconf[i,'cores.gpus']), sel))

  vals = a$comp[uconfSel]
  m[i] = mean(vals)
  count[i] = length(vals)
  sockets[i] = uconf[i,'cores.gpus']
}
pdf('cpuStrong.pdf')
plot(sockets, m, type='l', lwd=2, main='CPU strong scaling via MPI', 
     sub='1M SNPs', xlab='Cores (2 cores/node)', ylab='Computation time (s)', log='xy', axes=F, cex.lab=labcex, cex.main=maincex)
axis(at=sockets, side=1)
axis(at=m, labels=format(m,digits=3), side=2)
box()

points(sockets, m, pch=19)

pdf('cpuStrongSNPsPerSecond.pdf')
plot(sockets, 1000000/m, type='l', lwd=2, main='CPU strong scaling via MPI', 
     sub='1M SNPs', xlab='Cores (2 cores/node)', ylab='SNPs/s', log='xy', cex.lab=labcex, cex.main=maincex, axes=F)
axis(at=sockets, side=1)
axis(at=1000000/m, labels=format(m,digits=3), side=2)
box()
points(sockets, 1000000/m, pch=19)

cat('CPU strong scaling counts:\n')
print(count)
print(sockets)
print(m)


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
