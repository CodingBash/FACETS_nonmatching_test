#!/usr/bin/env Rscript
argv = commandArgs(trailingOnly=TRUE)

library(facets)

readSnpMat <- function(filename, skip = 0L, err.thresh = Inf, del.thresh = Inf,
                       perl.pileup = FALSE){
    pileup <- read.csv(filename, stringsAsFactors = FALSE,
                       colClasses = rep(c("character", "numeric", "character",
                                          "numeric"), c(1, 1, 2, 8)))
    colnames(pileup) <- c(colnames(pileup)[1:4], "File1R", "File1A",
                          "File1E", "File1D", "File2R", "File2A", "File2E", "File2D")
    if (pileup$Chromosome[1] == "chr1") {
        pileup$Chromosome <- gsub("chr", "", pileup$Chromosome)
    }
    ii <- which(pileup$File1E <= err.thresh & pileup$File1D <=
                    del.thresh & pileup$File2E <= err.thresh & pileup$File2D <=
                    del.thresh)
    rcmat <- pileup[ii, 1:2]
    rcmat$NOR.DP <- pileup$File1R[ii] + pileup$File1A[ii]
    rcmat$NOR.RD <- pileup$File1R[ii]
    rcmat$TUM.DP <- pileup$File2R[ii] + pileup$File2A[ii]
    rcmat$TUM.RD <- pileup$File2R[ii]
    return(rcmat)
}

rcmat <- readSnpMat(paste0("snpTmpG5_", argv[1], "_", argv[2], ".csv"))

set.seed(3224)
xx <- preProcSample(rcmat, ndepth=20, ndepthmax=1000)
oo <- procSample(xx, cval=150)
fit<-emcncf(oo)

saveRDS(xx, paste0("facetsG5XX_", argv[1], "_", argv[2], ".rds"))
saveRDS(oo, paste0("facetsG5OO_", argv[1], "_", argv[2], ".rds"))
saveRDS(fit, paste0("facetsG5Fit_", argv[1], "_", argv[2], ".rds"))
q()


