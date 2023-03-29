cat("** Start script\n")

library(mixOmics)

args <- commandArgs(trailingOnly = TRUE)
qq_path <- args[1]
mm_path <- args[2]
qm_path <- args[3]


#qq_path <- "S:\\U_Proteomica\\UNIDAD\\software\\MacrosRafa\\data\\Metabolomics\\PESA_Integromics\\Data\\Analysis\\02-Correlations\\ALDH4\\myRData\\qq_a.tsv"
#mm_path <- "S:\\U_Proteomica\\UNIDAD\\software\\MacrosRafa\\data\\Metabolomics\\PESA_Integromics\\Data\\Analysis\\02-Correlations\\ALDH4\\myRData\\mm_a.tsv"
#qm_path <- "S:\\U_Proteomica\\UNIDAD\\software\\MacrosRafa\\data\\Metabolomics\\PESA_Integromics\\Data\\Analysis\\02-Correlations\\ALDH4\\myRData\\qm_a.tsv"

cat("** Reading qq table\n")
qq <- read.csv(qq_path, sep='\t')

cat("** Reading mm table\n")
mm <- read.csv(mm_path, sep='\t')

cat("** Applying rCCA\n")
qqmm <- rcc(qq, mm, method='shrinkage')

# Get equiangular vectors
cat("** Calculating equiangular vectors\n")
Z <- qqmm$variates$X + qqmm$variates$Y

# Protein projection
cat("** qq projection\n")
x <- rbind(
  cor(
    Z[,1],
    qq  
  ),
  cor(
    Z[,2],
    qq  
  )
)

# Metabolite projection
cat("** mm projection\n")
y <- rbind(
  cor(
    Z[,1],
    mm  
  ),
  cor(
    Z[,2],
    mm  
  )
)

# Inner product between protein and metabolite
cat("** Inner product\n")
M <- t(x) %*% y


# Write output
cat("** Writing output table\n")
write.table(M, qm_path, sep="\t")

cat("** End Script\n")
