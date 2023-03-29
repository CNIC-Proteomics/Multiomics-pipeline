#
# Correlations
#

- En esta sección calculamos correlaciones intra/inter-ómica utilizando los siguientes métodos:
	- Media de Pearson, Spearman y Kendall* (PSK): QxQ & QxM & MxM
	- Partial Correlation (PC): QxQ & MxM
	- Graphical Lasso (intra) (RPC): QxQ & MxM
	- CCA (inter): QxM
	- rCCA (inter): QxM

- rCCA y Graphical Lasso son las que usaremos para la construcción de los grafos. Pero cuidado,
  pretendemos construir un grafo de correlación diferencial. No pretendemos ver features
  cuya media varía entre control y caso (para eso bastaría con ttest...). Pretendemos identificar
  parejas/grupos de biomoléculas cuya correlación/asociación se ve significativamente alterada
  como consecuencia del tratamiento/enfermedad. 

- Aspectos a considerar:
	- Cada uno de los métodos se aplicará a 3 subconjuntos de observaciones:
		. Usando todas las observaciones
		. Usando observaciones Control
		. Usando observaciones Enfermos

	- Además, obtendremos una cuarta matriz "delta-Correlation" con la diferencia
	  de correlaciones Enfermos-Casos. El objetivo es poder obtener una red 
	  de coexpresión/correlación diferencial (Differential Correlation Network).

	- Calcularemos pvalores empíricos para las correlaciones, lo que implica
	  aplicar el mismo esquema a un número determinado de matrices de datos
	  "barajadas"/aleatorizadas. Por tanto, cada matriz de correlaciones obtenida tendrá asociada 
	  su matriz de pvalores

	- A partir de las matrices de correlaciones generaremos Clustered Heatmaps:
		. QxQ & MxM & QxM
	  El heatmap de QxQ y MxM no muestra patrones relevantes o son difíciles de ver
	  El heatmap de QxM sí muestra una clusterización visual.

	- Además, empleando los valores Xq y Xm obtendremos un Clustered Heatmap esta vez
	  basado en la distancia euclídea de las observaciones y las features:
		. QxN & MxN & (Q U M)xN

- Conclusiones
	AWHS: 
	. La pipeline permite obtener para cada (método, ómica, muestra) tres tablas que permiten
	  construir el grafo de la siguiente sección:
		. Correlaciones
		. P-values empíricos (Monte Carlo Simulation/Permutation test)
		. Adjusted p-values

	. Para la construcción del grafo aplicaremos un pvalor ajustado (BH) de 0.01:	
		. QxQ : Graphical LASSO
		. MxM : Graphical LASSO
		. QxM : rCCA

#
# ALDH4
#

- En este caso trabajamos con 3 comparativas:
	. A12 - PBS
	. A12 - B1-8
	. B1-8 - PBS

- Dado que son múltiples comparaciones, es conveniente construir un único grafo que permita
  extraer conclusiones biológicas sin necesidad de hacer comaparaciones entre redes. Dos
  formas de contruir un grafo que sintetice la información de interés:

	. Para la construcción del grafo consideraremos las asociaciones entre features significativas
  en la comparación A12vsPBS siempre que no sean significativas (sustrayendo) en B1-8vsPBS.
  Esta es una manera de asegurarnos que la significatividad es derivada del efecto específico
  de A12 y no de la inoculación del anticuerpo.

	. Otra posibilidad sería considerar la intersección de las asociaciones significativas A12vsPBS y A12vsB1-8. Tal vez esta opción genere un grafo más pobre. Aunque es la mejor opción,
de hecho, podría ser una excusa para no aplicar la FDR, que no nos sale significativa.

- Es interesante obtener un heatmap clusterizado con las correlaciones rCCA entre q y m usando
  datos de A12vsPBs. Manteniendo el clusterizado, podemos obtener el mismo heatmap pero esta
  vez con los datos de A12vsB1-8. Las zonas que conservan el patrón revelan reproducibilidad.
  Es decir, correlaciones entre biomoléculas que aumentan o disminuyen como consecuencia del
  efecto específico de A12.

  Este par de heatmaps introduce la construcción del grafo fijándonos de manera específica
  en las asociaciones que se comportan igual. 
