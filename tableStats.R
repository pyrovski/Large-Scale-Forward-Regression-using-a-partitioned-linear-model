#!/usr/bin/env Rscript
cat('nodes cores adagio log count mean_time std_time\n')
a = read.table('filteredTable.dat',header=T)
conf = matrix(c(a$nodes,a$cores),ncol=2)
uniqueConfs = unique(conf)
rows = nrow(uniqueConfs)
for(i in 1:rows){
  nodes = uniqueConfs[i,1]
  cores = uniqueConfs[i,2]
  match1 = (conf[,1] == uniqueConfs[i,1])
  match2 = (conf[,2] == uniqueConfs[i,2])
  matchRows = which(match1 & match2)
  alg = matrix(c(a[matchRows,]$adagio, a[matchRows,]$log),ncol=2)
  uniqueAlgMatchRow = unique(alg)
  for(j in 1:nrow(uniqueAlgMatchRow)){
    adagio = uniqueAlgMatchRow[j,1]
    log = uniqueAlgMatchRow[j,2]
    # for all rows that match this configuration of nodes and cores, adagio and log,
    # generate the mean time and stdev time
    #matchAlg = (alg == uniqueAlgMatchRow[j,])
    matchAlg1 = (alg[,1] == uniqueAlgMatchRow[j,1])
    matchAlg2 = (alg[,2] == uniqueAlgMatchRow[j,2])
    # rows of matchRows
    matchAlgRows = which(matchAlg1 & matchAlg2)
    sel = matchRows[matchAlgRows]
    count = length(sel)
    write(paste(nodes, cores, adagio, log, count, mean(a[sel,]$runtime), sd(a[sel,]$runtime)), file = "", sep = '\t')
  }
}
