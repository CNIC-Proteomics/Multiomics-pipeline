library(sva)

args <- commandArgs(trailingOnly = TRUE)
Rpath <- args[1]

#Rpath <- "S:\\U_Proteomica\\UNIDAD\\software\\MacrosRafa\\data\\Metabolomics\\PESA_Integromics\\Data\\Proteomics\\ALDH4\\WorkingFiles\\myRData"

dat <- read.csv(paste0(Rpath, "\\dat.tsv"), sep = "\t", row.names = 1)
batch <- read.csv(paste0(Rpath, "\\batch.tsv"), sep = "\t", row.names = 1)
mod <- read.csv(paste0(Rpath, "\\mod.tsv"), sep = "\t", row.names = 1)

dfa <- ComBat(
    dat = dat,
    batch = batch$batch,
    mod = mod,
    par.prior = FALSE,
    prior.plots = FALSE
)

write.table(dfa, file = paste0(Rpath, "\\dfa.tsv"), sep = "\t")