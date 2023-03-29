#
# Identifications
#

- El cálculo de las correlaciones emplea algoritmos costosos computacionalmente. Especialmente,
  Graphical Lasso y rCCA. Por ello es conveniente aplicar un paso previos de feature selection.

- Como primer método para la selección de features emplearemos un F-test derivado de la construcción
  de un modelo lineal (ANOVA). Tomamos aquellos elementos significativos a un pvalor de 0.05.

- Si solo usáramos el F-test podríamos
  perder proteínas-metabolitos asociados cuya correlación cambia entre control y tratamiento.
  Piensa que pueden alternar o variar su correlación sin ser buenos clasificadores.

- Por ello, aplicamos un método de selección de features no supervisado que preserve
  la variabilidad de la muestra. La PCA calcula las componentes t.q. la proyección de las
  observaciones sobre ellas maximiza la varianza. Usamos una versión del PCA conocida
  como sparse PCA, que penaliza loadings distintos de 0.

- La idea es emplear el mayor alpha que retiene un número de features inferior a 1000 (en nuestro caso)
  Por debajo de 1000 features las correlaciones se calculan en un tiempo razonable. 

Refs:
	. https://scikit-learn.org/stable/modules/decomposition.html#sparsepca
	. https://hastie.su.domains/Papers/spc_jcgs.pdf
	. https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.SparsePCA.html

- Cuidado! Las componentes no son ortogonales y no se puede calcular la varianza explicada por componente.
  Existe una aproximación para el cálculo de la variabilidad retenida en el espacio latente: 
  https://arxiv.org/pdf/1907.03989.pdf

- Esta fracción de variabilidad explicada puede ser útil para asegurarnos de que la sPCA retiene
  variabilidad en el espacio latente. Obviamente, la variabilidad retenida preservando las features
  que no se anulan en el sPCA es mayor y ligeramente superior a su fracción asociada.