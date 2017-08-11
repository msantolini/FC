#!/usr/bin/env Rscript
require(WGCNA) # for faster "cor" function
set.seed(1) # for reproducibility


##########
## DATA ##
##########
# X is the expression personal fold-change matrix, rows are N strains, columns are M genes
# Y is the phenotype (hypertrophy), a vector of length M
# inds_filt is the indices of the genes kept after a pre-filtering

load('data.RData')
load('res_inds_filter.RData')
X <- datExpr_FC[,inds_filt]
Y <- datTraits_FC


###############
## FUNCTIONS ##
###############
# Displays 2 histograms on a same plot

hist2 <- function(dat1, dat2, nbreaks=50, log.breaks=F,main="", xlab="Value", ylab="", legend1="Data 1", legend2="Data 2",legend.pos='topleft',xlim=NULL,ylim=NULL,cex.axis=1.3,cex.lab=1.3,cex.leg=1.3,cex.main=1.3,col1=rgb(0,0,1,1/4),col2 = rgb(1,0,0,1/4), border=F, freq=FALSE,bty='n',...){
    "Produces two superimposed histograms for two sets of data,
    with alpha color blending."

    par(mar=c(6,6,5,5))

    if (log.breaks==TRUE){
        breaks <- 10^seq(from=log10(min(dat1,dat2)),
                         to=log10(max(dat1,dat2)),
                         length.out=nbreaks)
    } else {
        breaks <- seq(from=min(dat1,dat2),
                      to=max(dat1,dat2),
                      length.out=nbreaks
                      )
    }
    h1=hist(dat1,breaks=breaks, plot = FALSE)
    h2=hist(dat2,breaks=breaks, plot = FALSE)
    if (!length(xlim)) xlim <- c(min(h1$mids,h2$mids), max(h1$mids, h2$mids))
    if (!length(ylim)){
        if (freq==FALSE){
            ylim <- c(min(h1$density,h2$density), max(h1$density, h2$density))
        } else {
            ylim <- c(min(h1$counts,h2$counts), max(h1$counts, h2$counts))
        }
    }

    hist(dat1,
         breaks=breaks,
         freq = freq,
         col=col1,
         xlim = xlim,  
         ylim = ylim,  
         main=main,
         ylab=ylab,
         xlab=xlab,
         cex.main = cex.main,
         cex.lab = cex.lab,
         cex.axis = cex.axis, 
         border=border,
         ...
         )
    hist(dat2, 
         breaks=breaks,
         freq = freq,         
         col=col2, 
         add=TRUE,
         border=border,
         ...
         )
    legend(legend.pos, 
           legend = c(legend1, legend2), 
           border=border,
           bty=bty,
           fill = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)),
           cex = cex.leg) 
}

##########
## MAIN ##
##########

#==========================================================#
# 1. Compute gene-trait correlation, and a randomized case #
#==========================================================#
# Note that in the randomization, I keep the gene-gene correlation structure
# constant, so I really test how that given gene set would correlate with
# a random trait by chance. It is also better for comparing with WGCNA: you
# can see how WGCNA eigengenes would correlate with random traits and compute a 
# similar p-value.

cc <- abs(cor(X,Y,use="pairwise.complete.obs"))
cc_r <- abs(sapply(1:1000,function(i) cor(X,sample(Y),use="pairwise.complete.obs")))

pdf('plotHistograms.pdf')
hist2(cc, cc_r,
      legend1='Observed',
      legend2='Randomized',
      nbreaks=20,
      legend.pos='topright',
      xlab='Correlation with trait (abs.)',
      ylab='Frequency'
      )
dev.off()


#=============================#
# 2. Then, compute enrichment #
#=============================#

breaks <- sort(cc)
props_obs <- 1-ecdf(cc)(breaks)
props_r <- 1-ecdf(cc_r)(breaks)
s2n <- props_obs / props_r

png(paste0('plotEnrichment.png'),600,600,bg='transparent')
plot.0(rev(s2n),
       log='x',
       ylab='Enrichment',
       xlab='Gene rank',
       type='o',
       lwd=3,
       cex.axis=2.5,
       cex.lab=2.5,
       cex.main=3,
       col=rgb.0('blue'),
       pch=19
       )
abline(h=1,
       lty=2,
       col='red',
       lwd=5
       )
dev.off()


#===========================#
# 3. Output order gene list #
#===========================#

ind_cutoff <- which.max(rev(s2n))
write(symbols[inds_filt[order(-cc)]][1:ind_cutoff], 
      file='genes_ordered.txt')

