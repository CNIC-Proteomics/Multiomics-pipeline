#
# Basic Stats
#

- El objetivo de esta sección es realizar cálculos estadísticos sobre cada ómica de manera independiente.
- De esta manera, aplicaremos diferentes tipos de tests, contrastes y modelos, que asignen a cada feature
  un p-valor, coeficiente o indicador asociado a su significatividad.
- Todos estos modelos nos indicarán, para cada feature, la capacidad predictiva sobre una métrica
  asociada a la aterosclerosis subclínica. Esta métrica será la variable dependiente en los modelos.
- Por tanto, aplicaremos diferentes modelos y diferentes variables dependientes, que procedemos a enumerar.

- Variable dependiente:
	. Grupo Control/Enfermo
	. Tamaño de placa
	. HDL/LDL
	. Diabetes
	. Tension/Hipertensión
	. Calcio
	. Glucosa
	. Edad

- Los modelos aplicados serán:
	. Logistic Regression: Simple / Multiple / Multiple Regularised
	. Linear Regression: Simple / Multiple / Multiple Regularised*
	. (O)PLS (DA): Regularised*
	. Correlation: Pearson / Spearman / Kendall
	. t-test & U-test

# Resultados
- El resultado será un (multi-index) dataframe que asocie a cada feature el valor correspondiente.
- Creamos una librería en Python que paralelice.

# Utilidad
- Esta tabla permitirá señalar en el grafo los nodos/features significativos según diferentes criterios,
  algo que puede ser útil al comparar grafo Control, Disease y Aleatorio.
- Además de señalar los nodos, la idea es centrar la comparación y el análisis de métrica en ellos.

# ALDH4
- En el caso de ALDH4 aplicaremos ANOVA y Kruskal-Wallis con post-hocs y corregidos por FDR. 
- Esto nos permitirá identificar features que difieren entre los 3 grupos y específicamente
  aquellos que distinguen A12 de B1-8 y PBS. 
