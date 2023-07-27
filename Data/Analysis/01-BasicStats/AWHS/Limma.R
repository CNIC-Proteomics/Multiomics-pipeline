library(limma)

setwd("S:\\U_Proteomica\\UNIDAD\\software\\MacrosRafa\\data\\Metabolomics\\PESA_Integromics\\Data\\Analysis\\01-BasicStats\\AWHS")

xq_path <- "S:\\U_Proteomica\\UNIDAD\\software\\MacrosRafa\\data\\Metabolomics\\PESA_Integromics\\Data\\Proteomics\\AWHS\\WorkingFiles\\Xq_minus_X_norm.tsv"
xm_path <- "S:\\U_Proteomica\\UNIDAD\\software\\MacrosRafa\\data\\Metabolomics\\PESA_Integromics\\Data\\Metabolomics\\AWHS\\WorkingFiles\\Xm_norm_MS2.tsv"
mdata_path <- "S:\\U_Proteomica\\UNIDAD\\software\\MacrosRafa\\data\\Metabolomics\\PESA_Integromics\\Data\\Metadata\\AWHS\\WorkingFiles\\main_metadata.tsv"

xq <- read.csv(xq_path, sep="\t", row.names = 'Seqn')
xq <- t(xq)

xm <- read.csv(xm_path, sep="\t", row.names = 'fid')
xm <- t(xm)

mdata <- read.csv(mdata_path, sep="\t", row.names = 'Seqn')


#
# Proteomics
#

xq_design <- mdata[colnames(xq), 'Group']
xq_design <- cbind(C=1, DvsC=as.integer(xq_design=='D'))

fit <- lmFit(xq, xq_design)
fit <- eBayes(fit)
qLimma <- topTable(fit, coef="DvsC", adjust="BH", number=nrow(xq))

write.table(qLimma, file='qLimma.tsv', sep='\t', row.names = T, col.names=T)

#
# Metabolomics
#

xm_design <- mdata[colnames(xm), 'Group']
xm_design <- cbind(C=1, DvsC=as.integer(xm_design=='D'))

fit <- lmFit(xm, xm_design)
fit <- eBayes(fit)
mLimma <- topTable(fit, coef="DvsC", adjust="BH", number=nrow(xm))

write.table(mLimma, file='mLimma.tsv', sep='\t', row.names = T, col.names=T)
