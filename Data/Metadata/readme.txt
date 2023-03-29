############
# METADATA #
############

#      #
# PESA #
#      #

#
# Original Files
#

- "RBR_main_metadata.xlsx": Metadatos extraídos de la tabla de metabolómica "V1_C18_NEG_LOESS.xlsx"


#
# Working Files
#

Generado usando WorkingFilesGenerator.ipynb

- "main_metadata.tsv": Igual que "RBR_main_metadata.xlsx" pero poniendo variables en columnas



#      #
# AWHS #
#      #

#
# Original Files
#

- "Cohortes_AWHS.xlsx": Metadatos copiados de "S:\LAB_JVC\RESULTADOS\AWHS\Metadatos\Cohortes_AWHS.xlsx"
- "STAT-LOG-NEG_xEstefi(AutoRecovered).xlsx": Matriz con datos de metabolómica donde Alessia
  ha añadido información sobre cada batch y orden de pinchazo. De aquí obtenemos la tabla
  "RBR_MetabolomicsBatch.xlsx" que asocia a cada Seqn el batch y el orden de pinchazo.

#
# Working Files
#

- "main_metadata.tsv": Hemos añadido una columna Group que indica si el individuo es control (C) o enfermo (D)