#
# Factor Analysis
#

- Utilizamos el modelo Multi-Omics Factor Analysis (MOFA), una versión de los modelos estadísticos
Factor Analysis adaptado para el estudio de datos ómicos. El modelo presupone la existencia de variables
latentes, de las que dependen las features analizadas. Features que dependen del mismo factor
tienen un comportamiento similar. El modelo permite identificar conjuntos de features con un comportamiento
similar u opuesto entre las muestras.

- Paper del modelo:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6010767/pdf/MSB-14-e8124.pdf

- GitHub:
https://github.com/bioFAM/MOFA2
https://biofam.github.io/MOFA2/installation.html

- Web desde la que explorar los resultados del modelo
http://www.ebi.ac.uk/shiny/mofa/


- ¿Qué podemos ver?

	1) Identificar conjuntos de features asignados de manera significativa a cada uno de los factores.
	Nos centramos únicamente en factores con una fracción de varianza explicada superior al 1.5% en 
	algunas de las ómicas (en el paper ponen 2%, lo reducimos para retener más factores).

	Para determinar la significatividad podríamos emplear p-valores empíricos y ajustados, pero al
	generar datos aleatorios el modelo no es capaz de encontrar patrones y genera error.

	2) Cuáles son estos conjuntos de features en cada factor y qué procesos biológicos tienen asociados.
	Esto en principio no da información sobre procesos/features asociados a la aterosclerosis subclínica,
	tan solo una descripción de las features que aparecen asociados conjuntamente a un proceso biológico.

	3) Para emplear el MOFA en el estudio de la aterosclerosis subclínica, o cualquiera de los metadatos
	asociados, debemos usar la proyección de las observaciones sobre cada uno de los factores y estudiar
	mediante tests estadísticos la asociación de cada uno de estos factores con las diferentes condiciones
	de las observaciones (Grupo, glucosa, colesterol...). Considerar tests no paramétricos.

	Otra opción es fijarnos en features asociadas a ese factor y que son significativas
	en aterosclerosis.

	Comparar estos resultados con los obtenidos en el análisis de grafo.


#
# AWHS
#

- El primer factor gobierna y explica el 30.60% de la variabilidad observada en los lípidos.
  En concreto, el primer factor factor "absorbe" o explica la variación conjunta de los TAG.
  Aquí sería interesante ver, no tanto qué metabolitos están controlados por este factor, pues
  vemos que son la inmensa mayoría (fundamentalmente triglicéridos), 
  sino comprobar qué proteínas cambian conjuntamente con
  el conjunto de los TAG. Hacemos para ello un enriquecimiento con las proteínas del primer
  factor.

- El segundo factor resulta más interesante al mostrar de manera más clara la correlación
  entre proteínas y diferentes lípidos, indicando agrupamiento y clusterización entre pacientes.
  
	-> Separación entre control y tratamiento por el segundo factor
	-> Identificar clusters de paciente y buscar en metadatos en qué se diferencian del resto.

- El tercer factor es relevante únicamente a nivel metabolómico, pues apenas explica variabilidad
  de los datos de proteómica. En concreto, este factor está asociado a diferentes tipos de TAG
  donde se identifican dos grupos que se comportan de manera contraria. Un grupo son TAG [4x:1o2]
  y el otro grupo son TAG [5x: >5]. 

	-> Separación entre control y tratamiento por el tercer factor
	-> Identificar clusters de pacientes y sus metadato ¿?

- A partir del tercer factor ya dejamos de ver patrones interesantes.