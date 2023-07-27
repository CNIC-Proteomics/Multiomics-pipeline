##############
# PROTEOMICS #
##############


#      #
# PESA #
#      #

#
# Original Files
# 

- "All_Xq_Wq_V1.xlsx":

Fichero copiado de:
"S:\U_Proteomica\PROYECTOS\PESA_omicas\Results_tables_V1\Con_fraccionamiento\Xq\All_Xq_Wq_V1.xlsx"

De la sheet "Xq_CV" extraemos los valores Xq-Xgrand_mean para cada proteína y en cada muestra.
En concreto, de las columnas QI a AHN las hemos copiado en el fichero "RBR_Xq_minus_X.xlsx".


- "Seqn_TMT_Cohorte.xlsx":

Fichero copiado de:
"S:\U_Proteomica\PROYECTOS\PESA_omicas\Results_tables_V1\Con_fraccionamiento\Seqn_TMT_Cohorte.xlsx"


#
# Working Files
#

Creados usando WorkingFilesGenerator.ipynb

- "TMT2seqn.tsv": Mapear cada TMT con el identificador del paciente (Not necessary)

- "q2info.tsv": Mapear cada proteína con Np (número de péptidos ¿?)

- "Xq_minus_X.tsv": Valores de X'q sin filtrar

- "Xq_minus_X_norm.tsv": Valores de Xq aplicando:
	- Mantener proteínas con presencia superior al 80%
	- Imputar missing values mediante KNN n=35
	- Centrado y escalado de los valores


#      #
# AWHS #
#      #

#
# Original Files
#

- "All_Xq_Wq_absoluta.xlsx":

Fichero copiado de:
"S:\LAB_JVC\RESULTADOS\AWHS\All_Zq\All_Xq_Wq_absoluta.xlsx"

De este fichero copiamos la Sheet='Xq' en el fichero "RBR_Xq_minus_X.xlsx", donde restamos la grand mean
en la Sheet2 (Xq-GranMean(X)).


- "SPSS_all.xlsx":

Fichero copiado de:
"S:\LAB_JVC\RESULTADOS\AWHS\All_Zq\Sin_Fraccionamiento\All_Zq_AWHS.xlsx"

Del header de este fichero extraemos el mapeo de Cohorte-TMT a Seqn en "RBR_Seqn_TMT_Cohorte"


#
# Working Files
#

Mismos ficheros que los de PESA, pero creados con otro WorkingFilesGenerator.ipynb



#       #
# ALDH4 #
#       #

#
# Original Files
#

- Los ficheros de proteómica corresponientes a este experimento se encuentran en:
	"S:\U_Proteomica\LABS\LAB_ARR\ClonesAb-atherosclerosis\Higados\proteomica"

- Para generar la tabla Xq_minus_X(_norm).tsv, con la estructura NxQ necesitamos extraer
  los X'inf de las carpetas TMT[12]/msf/[tag]/data/Q2A_lowerNormW.xls

- Para conocer la identidad de cada tag debemos emplear el siguiente fichero:
	. "S:\U_Proteomica\LABS\LAB_ARR\ClonesAb-atherosclerosis\Higados\proteomica\tags.xlsx"
	. "S:\U_Proteomica\LABS\LAB_ARR\ClonesAb-atherosclerosis\Higados\proteomica\protocolo.docx"

- Para obtener el número de scanes/péptidos de cada proteína empleamos los idq de cada TMT:
	. "S:\U_Proteomica\LABS\LAB_ARR\ClonesAb-atherosclerosis\Higados\proteomica\TMT1\msf\ID_Q_XV.txt"
	. "S:\U_Proteomica\LABS\LAB_ARR\ClonesAb-atherosclerosis\Higados\proteomica\TMT2\msf\ID_Q_XV.txt"

Cada observación proviene de un ratón con un identificador asociado. Necesitamos la relación (ratón, TMT, tratamiento)

#
# Working Files
#

Empleando WorkingFilesGenerator generamos los siguientes ficheros:

- Xq_minus_X.tsv: Tabla NxQ con los valores de Xq'inf
- Xq_minus_X_norm.tsv: Tabla NxQ con los valores Xq'inf 
	- Estandarizados por Q
	- 
- q2info.tsv: Información de proteínas (id, desc, ScanFreq)
- n2info.tsv: 



#         #
# PESA V2 #
#         #

#
# Original Files
# 

- Los datos cuantitativos se obtuvieron directamente de los ficheros Q2A_outStats.xls del SanXoT.
  La estructura del path es la siguiente:
  "S:\U_Proteomica\PROYECTOS\PESA_omicas\2a_Cohorte_120_V2\Proteomics\TMT_Fraccionamiento\TMT1\SanXoT\127_C\data\Q2A_outStats.xls"

- Para calcular el LowerNorm (X') hemos restado en estos ficheros la columna Xinf-Xsup.