#!/usr/bin/env Rscript
require(igraph)
require(WGCNA) # for fast cor function

ll <- load('data.RData')
X <- datExpr_FC

#############################################
## Computing the full co-expression matrix ##
#############################################
n <- ncol(X) 
C <- abs(cor(X))
C1 <- C
    

################################################
## Connected components measured using igraph ##
################################################

fout <- 'res_clusters_init.RData'
if (!file.exists(fout)){

    thrs <- seq(0.3,0.9,0.02)
    ccs <- list()
    counter <- 1
    for (thr in thrs){
        print(thr)
        C1[C1 < thr] <- 0
        G <- graph.adjacency(C1, mode='upper', weighted=T, diag=F)
        cc <- clusters(G)
        print(cc$no)
        ccs[[counter]] <- cc
        counter <- counter + 1
    }

    save(thrs, ccs, file=fout)
} else {
    load(fout)
}


nb_clus <- unlist(lapply(ccs, function(cc) cc$no))

pdf('plotNumberOfClusters.pdf',width=8)
par(mfrow=c(1,2))
plot(thrs,nb_clus,
     ylim = c(0,nrow(C)),
     pch=19,
     type='o',
     xlab='Correlation cutoff',
     ylab='Number of clusters'
     )
abline(h = nrow(C),
       lty=2,
       col='red')
plot(thrs[-1],diff(nb_clus),
     pch=19,
     type='o',
     xlab='Correlation cutoff',
     ylab='Difference in number of clusters'
     )
dev.off()

LCC_size <- unlist(lapply(ccs, function(cc) max(cc$csize)))

pdf('plotLCCSize.pdf',width=10)
par(mfrow=c(1,2), mar=c(6,6,5,5), cex.lab=1.5,cex.axis=1.5)
plot(thrs,LCC_size,
     ylim = c(0,nrow(C)),
     pch=19,
     type='o',
     xlab='Correlation cutoff',
     ylab='LCC size'
     )
abline(h = nrow(C),
       lty=2,
       col='red')
plot(thrs[-1],diff(LCC_size),
     pch=19,
     type='o',
     xlab='Correlation cutoff',
     ylab='Difference in LCC size'
     )
dev.off()


##############################
## Getting the core cluster ##
##############################
thr <- 0.5
C1 <- C
C1[C1<thr] <- 0
G = graph.adjacency(C1, mode='upper', weighted=T, diag=F)
cc <- clusters(G)

pdf('plotClusterSizes.pdf')
plot.0(cc$csize[order(-cc$csize)],
       type='o',
       log='xy',
       lwd=3,
       xlab='Connected component (ranked by size)',
       ylab='Size'
       )
dev.off()

inds_filt = which(cc$membership==which.max(cc$csize))

save(inds_filt,file='res_inds_filter.RData')

X <- datExpr_FC
Y <- datTraits_FC[,1]
cc1 <- cor.0(X[,inds_filt],Y)
cc2 <- cor.0(X[,-inds_filt],Y)
cc_r <- sapply(1:1000,function(i) cor.0(X,sample(Y)))

pdf(paste0('plotHist_filt.pdf'))
hist2(
      abs(cc1),abs(cc_r), 
      legend1='Kept', 
      legend2='random',
      legend.pos='topright',
      nbreaks=15,
      #main=paste0('-log10(p) = ',signif(-log10(p),3)),
      xlab='Correlation with phenotype (abs)',
      ylab='Frequency',
      lwd=1
      )
dev.off()

pdf(paste0('plotHist_filtered_out.pdf'))
hist2(
      abs(cc2),abs(cc_r), 
      legend1='Filtered out', 
      legend2='random',
      legend.pos='topright',
      nbreaks=15,
      #main=paste0('-log10(p) = ',signif(-log10(p),3)),
      xlab='Correlation with phenotype (abs)',
      ylab='Frequency',
      lwd=1
      )
dev.off()

