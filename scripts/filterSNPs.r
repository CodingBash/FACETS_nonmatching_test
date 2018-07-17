#!/usr/bin/env Rscript

snpMatrix <- as.data.frame(read.csv("snpTCGATest20G5.csv", stringsAsFactors = FALSE))
col_names <- colnames(snpMatrix)

print(nrow(snpMatrix))
cov <- 20

delete_rows <- c()
delete_rows.i <- 1
for(snpMatrix.row in seq(1, nrow(snpMatrix))){ # For each row
  for(normal.col in seq(from=4+1, to = ncol(snpMatrix), by=8)){ # For each normal column (starting at FileXR)
	print(paste("row=", snpMatrix.row, "col=", normal.col))
    if(snpMatrix[snpMatrix.row, normal.col] + snpMatrix[snpMatrix.row, normal.col + 1] < cov){ # If the FileXR + FileXA is less than coverage requirement
      delete_rows[delete_rows.i] <- snpMatrix.row
      delete_rows.i <- delete_rows.i + 1
      break
    }
  }
}

filteredSnpMatrix <- snpMatrix[-delete_rows,]

write.table(filteredSnpMatrix, "filteredSnpTCGATest20G5.csv", row.names = F, sep = ",", quote = FALSE)

