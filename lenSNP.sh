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
