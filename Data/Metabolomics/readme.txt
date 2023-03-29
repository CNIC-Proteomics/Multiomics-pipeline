################
# METABOLOMICS #
################


#      #
# PESA #
#      #

#
# Original Files
#

Los siguientes ficheros:

- "PESA_V1_C18_NEG_Filtered.xlsx"
- "PESA_V1_C18_POS_Filtered.xlsx": 	. Trabajar con datos que van desde la celda 36 a 733.
						. A la derecha en azul el batch que hubo que excluir
- "PESA_V1_HILIC_NEG_Filtered.xlsx"
- "PESA_V1_HILIC_POS_Filtered.xlsx"
- "V1_Identifications.xlsx"

han sido copiados del siguiente path:
S:\U_Proteomica\PROYECTOS\PESA_omicas\METABOLOMICS\V1_Multiomics

Los cuatro primeros ficheros contienen las intensidades obtenidas para cada feature (fila)
a lo largo de las diferentes muestras (columnas), incluyendo el QC.
El quinto fichero contiene las identificaciones de las features más relevantes.

RBR_f2i.xlsx y RBR_PESA_V1.xlsx contiene en cada sheet la info de cada plataforma.

Utilizando los metadatos generamos el "RBR_samPeno.xlsx" con el
formato del Pheno File de StatTarget, que necesitaremos para aplicar el QC-RLSC:
	. Name = sample
	. Group = class
	. Global Order = order
	. Global Batch = batch
Para evitar posibles errores, cada matriz debe tener su sheet samPeno usando como sample
la primera fila (no da lugar a errores). Observaciones:
	. En C18N intercambiamos 527 y 528 para que el QC quede al final (de otra forma da error).
	. La muestra C18_C1809_2n está asociada al Global batch 9, pero su Global order pertenece al
	  10. Le reasignamos el Global batch 10 para la corrección de señal.
	. En HILIC algunas de las muestras con Global batch 2, pertenecen al 3. Lo corregimos.

* En C18P faltan algunas observaciones debido a caída de señal en el equipo para estos batches.
  Construimos matriz global eliminando estas observaciones

* LOESS lo hacemos fijando el smoothing a 0.2 (Cross Validation genera efecto batch)

#
# Working Files
#

Creados usando WorkingFilesGenerator.ipynb

- "Xm.tsv": ## No la hemos sacado, pues no la usaremos
- "Xm_norm.tsv": Aplicar logaritmo, centrado y escalado a la salida del LOESS
  Los missing values fueron imputados empleando valores aleatorios de una Normal(0,1)
  tras el centrado y escalado.
- "f2info.tsv": Mapeo feature con información (columna-carga, m/z, rt)



#      #
# AWHS #
#      #


#
# Original Files
#

- "STAT-LOG-POS_xEstefi.xlsx":

Fichero copiado de: "S:\U_Proteomica\PROYECTOS\PESA_omicas\LIPIDOMICS\AWHS_Lipidomics\AWHS_Statistics\STAT-LOG-POS_xEstefi.xlsx"
De este fichero obtenemos "RBR_LOG_POS_INFO.xlsx" con las siguientes Sheets:
  - Sheet1: Datos originales, aunque hemos añadido un identificador a cada feature
  - fid2log: Matriz Px2N con Log10(Concentracion) (P=Numero de features, N=Número de muestras) 
    La concentración fue medida con precursor MS1 y fragmento MS2 (Alessia aconseja empezar con MS2)
    Las muestras están identificadas con Seqn
  - fid2LipidInfo: Mapear fid con información de los lípidos
  - sample2Seqn: Mapear sample_id con Seqn (no creo que sea necesario utilizarlo)

*** Cuidado con los missing values. Distinguimos 3 casos:
	- Casos en los que no se detectó el compuesto con el equipo (pero sí el estándar interno)
	  Lo que se hizo fue imputar la concentración asignando 0.00001. Al aplicar log es -5.
	- Casos en los que no se produjo la fragmentación en ningún caso (None). Los eliminamos de la tabla MS2.
	- Casos en los que no se ha detectado el estándar interno (NaN). Imputamos por KNN. 

- "STAT-LOG-NEG_xEstefi.xlsx":

Fichero copiado de: "S:\U_Proteomica\PROYECTOS\PESA_omicas\LIPIDOMICS\AWHS_Lipidomics\AWHS_Statistics\STAT-LOG-NEG_xEstefi.xlsx"
De este fichero obtenermos "RBR_LOG_NEG_INFO.xlsx" con las mismas Sheets descritas arriba


#
# Working Files
#

En las tablas siguientes unificamos positivo y negativo.

- "Xm_MS1.tsv": Valores de la tabla original. Mediciones con el precursor
- "Xm_MS2.tsv": Valores de la tabla original. Mediciones con la fragmentación
- "Xm_norm_MS1.tsv": Valores originales estandarizados (precursor)
- "Xm_norm_MS2.tsv": Valores originales estandarizados (fragmentación)
	. Para obtener la tabla normalizada:
		. Eliminado features con más de un 20% de missing values
		. Imputar missing values con KNN (n=3)
		. Centrado y escalado
- "m2info.tsv": Mapear fid a información lipídica