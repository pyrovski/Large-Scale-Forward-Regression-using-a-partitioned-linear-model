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

branchSel = which(a$branch == 'test_GPU_all')
coreSel = which(a$cores.gpus == 1)
sel = mergeSels(list(branchSel, coreSel))

conf = a$SNPs_on_rank_0[sel]
uconf = sort(unique(conf))

# merge entries matching sel and uconf
cat('branch cores/gpus count SNPs_on_rank_0 mean_time std_time\n')
m = c()
std = c()
for(i in 1:length(uconf)){
  vals = a$comp[intersect(which(a$SNPs_on_rank_0 == uconf[i]), sel)]
  m[i] = mean(vals)
  std[i] = sd(vals)
  
  cat('test_GPU_all 1 ', uconf[i], ' ', m[i], ' ', std[i], '\n')
}
