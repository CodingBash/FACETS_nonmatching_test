#!/usr/bin/env Rscript

coverageFilteredSnpMatrix <- as.data.frame(read.csv("filteredSnpTCGATest20G5.csv", stringsAsFactors = FALSE))
col_names <- colnames(coverageFilteredSnpMatrix)

print(nrow(coverageFilteredSnpMatrix))

delete_rows <- c()
delete_rows.i <- 1
for(coverageFilteredSnpMatrix.row in seq(1, nrow(coverageFilteredSnpMatrix))){ # For each row
  for(normal.col in seq(from=4+1, to = ncol(coverageFilteredSnpMatrix), by=8)){ # For each normal column (starting at FileXR)
	print(paste("row=", coverageFilteredSnpMatrix.row, "col=", normal.col))
	logOR <- log2(coverageFilteredSnpMatrix[coverageFilteredSnpMatrix.row, normal.col + 1] / coverageFilteredSnpMatrix[coverageFilteredSnpMatrix.row, normal.col + 4 + 1])
	print(paste("logOR=", logOR, "tv=", coverageFilteredSnpMatrix[coverageFilteredSnpMatrix.row, normal.col + 4 + 1], "nv=", coverageFilteredSnpMatrix[coverageFilteredSnpMatrix.row, normal.col + 1]))
    if(is.nan(logOR)) next # TODO: Should nan be removed? let FACETS decide
    if(logOR > 4 || logOR < -4){ # If the FileXR + FileXA is less than coverage requirement
      delete_rows[delete_rows.i] <- coverageFilteredSnpMatrix.row
      delete_rows.i <- delete_rows.i + 1
      break
    }
  }
}

logORFilteredSnpMatrix <- coverageFilteredSnpMatrix[-delete_rows,]

write.table(logORFilteredSnpMatrix, "logORFilteredSnpTCGATest20G5.csv", row.names = F, sep = ",", quote = FALSE)

