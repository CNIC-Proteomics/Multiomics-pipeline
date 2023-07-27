#
# Import libraries
#

library(statTarget)
library(readxl)

#
# Constants
#
RLOESS <- T

basePath <- "S:\\U_Proteomica\\UNIDAD\\software\\MacrosRafa\\data\\Metabolomics\\PESA_Integromics\\Data\\Metabolomics\\PESA_V2\\OriginalFiles"
mode <- "HILP" # C18N, C18P, HILN, HILP

for (mode in c('C18N', 'C18P', 'HILP')) {
setwd(file.path(basePath, 'LOESS', mode))

#
# Prepare samFile
#

samFile <- read_excel(
  path=file.path(basePath, 'RBR_PESA_V2.xlsx'),
  sheet=mode
)


#
# Prepare samPeno
#

samPeno <- read_excel(
  path=file.path(basePath, 'RBR_SamPheno.xlsx'),
  sheet=mode
)[, c('sample', 'batch', 'class', 'order')]

samPeno$class[samPeno$class == 'QC'] <- NA


#
# Write csv files
#

samFile_path <- file.path(basePath, 'LOESS', mode, 'samFile.csv')
write.csv(samFile, samFile_path, quote=FALSE, row.names = F)

samPeno_path <- file.path(basePath, 'LOESS', mode, 'samPeno.csv')
write.csv(samPeno, samPeno_path, quote=FALSE, row.names = F)

#
# Execute shiftCor
#

if (RLOESS){
shiftCor(
  samPeno = samPeno_path,
  samFile = samFile_path,
  Frule = 0.8,
  MLmethod = "QCRLSC",
  degree=2,
  QCspan = 0.2,
  imputeM = "KNN",
  coCV = 30
)
}

}

