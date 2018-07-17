#!/usr/bin/env Rscript
argv = commandArgs(trailingOnly=TRUE)

total_samps <- 7 # We have 7 total samples (and 14 BAMs)
TUMOR <- as.numeric(argv[1]) # ID of the target tumor BAM
NORM <- seq(from = 1, to = total_samps, by = 1) 

library(facets)

createSnpMat <- function(snpMatrix, normalId, tumorId, skip = 0L, err.thresh = Inf, del.thresh = Inf,
                       perl.pileup = FALSE){
	tumorCol <- 4 + (1 + (8 * (tumorId - 1)))
	normalCol <- 4 + (1 + (8 * (normalId - 1))) + 4
	
	pileup <- snpMatrix[,c(1,2,3,4,
				normalCol, normalCol + 1, normalCol + 2, normalCol + 3,
				tumorCol, tumorCol + 1, tumorCol + 2, tumorCol + 3)]
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

snpMatrix <- as.data.frame(read.csv("filteredSnpTCGATest20G5.csv", stringsAsFactors = FALSE, 
	colClasses = rep(c("character", "numeric", "character", "numeric"), c(1, 1, 2, 14 * 4))))
runFacets <- function(rcmat, tumorId, normalId){
	print(paste("Running FACETS on tumorId:", tumorId, "and normalID:", normalId))
	set.seed(3224)
	xx <- preProcSample(rcmat, ndepth=20, ndepthmax=1000)
	oo <- procSample(xx, cval=150)
	fit<-emcncf(oo)

	saveRDS(xx, paste0("facetsG5XX_", tumorId, "_", normalId, ".rds"))
	saveRDS(oo, paste0("facetsG5OO_", tumorId, "_", normalId, ".rds"))
	saveRDS(fit, paste0("facetsG5Fit_", tumorId, "_", normalId, ".rds"))
}

for(NORM.i in NORM){ # Run for non-matching normal
	runFacets(createSnpMat(snpMatrix, NORM.i, TUMOR), TUMOR, NORM.i)
}

q()


