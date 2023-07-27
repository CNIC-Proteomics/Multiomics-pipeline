- Anotación putativa de las m/z empleando CMM y simplificado con TurboPutative.
- Exploratory Analysis and Quality Control Check
	o PCA independiente y combinada y análisis de varianza entre la proyección de las componentes y los metadatos
	o UMAP
	o Histograma & Boxplot
	o Scatter con distribución de logFC por clase lipídica
- MOFA
	o Porcentaje de explicación de cada factor
	o Análisis de varianza entre la proyección y los metadatos
	o Features con coeficientes más significativos (y su LogFC) / Distribución de coeficientes
	o Heatmap con distribución de valores considerando los N coeficientes
- rCCA
	o Cálculo de las correlaciones considerando todas las muestras
		. Tabla para explorar todas las correlaciones
		. Representar correlaciones entre proteína y clase lipídica (diapos)
	o Cálculo de las correlaciones diferenciales
		. Tabla para explorar todas las correlaciones
		. Representar correlaciones diferenciales para proteína y clase lipídica
- Graphical Lasso
	o Cálculo de correlaciones diferenciales en cada ómica & Detección de comunidades & Significatividad de las comunidades mediante proyección y análisis de varianza & Enriquecimiento de las comunidades & Conexión entre comunidades de proteómica y metabolómica (esta parte no está clara)
- Elastic Net
	o Construcción de Elastic Net con proteómica, metabolómica y ambas
		. AUC / -log(Spearman_pvalue)
		. Detección de features con elementos significativos y boxplot
		. Comparar modelo conjunto con modelos individuales
		. Evolución mediante eliminación iterada
