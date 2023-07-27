# TAREAS
Tareas:
1) Datos originales de PESA y ALDH4. Ambos son experimentos
Untargeted, por lo que hay que aplicar el preprocesado correspondiente. 
	. Utilizar la función shiftCor de statTarget en R:
		- https://rdrr.io/bioc/statTarget/man/shiftCor.html
		- https://stattarget.github.io/docs/quick-start/
		  Formato de tablas
		- https://www.sciencedirect.com/science/article/abs/pii/S0003267018309395?via%3Dihub
		  De este paper sacar que el valor corregido se obtiene mediante un ratio
		- Version de escritorio: https://github.com/statTarget/statTarget2

2) Análisis de grafos
	. Generar grafos siguiendo las instrucciones del readme
	. Jugar con NetworkX y Pyvis para la representación interactiva
	. Generar grafos aleatorios
	. Analizar grafos y comparar con aleatorios.
	  (orden, tamaño, componentes conexas, densidad, distribución de grado, cluster idx, camino característico)
	  Nodos más importantes a partir de su centralidad (betweenness, cercanía, lejanía, grado)
	. Representación cuqui del grafo

3) Enriquecimiento con GO terms asociados a proteínas

----
4) Complementar lo observado en el grafo con el correlation circle plot 
   (los datos para esto igual hay que sacarlo, si no, lo hacemos en R),
   cluster image map (X vs Y --> la matriz de correlación coloreada y agrupada) 
   y heatmap (X/Y vs sample --> los valores de las features coloreados y agrupados)

5) Los clusters/proteínas-metabolitos observados deben interpretarse por qué aparecen conectados y 
   si hay diferencias entre control/enfermo, correlación con placa...
   En el grafo podemos señalar vertices potencialmente significativos siguiendo algún criterio:
	- ttest con corrección sobre C/D
	- Regresión logística (regularizada)
	- Regresión lineal múltiple (regularizada) sobre tamaño de placa (statsmodel para sacar pvalores)
	- rPLS-reg

  Por tanto, tenemos un grafo de biomoléculas correlacionadas. 
  Sobre el grafo podemos aplicar diferentes criterios de significatividad a los vértices y buscar clusters.
  

#
# RECURSOS DE OTRAS BASES DE DATOS
#


----

ConsensusPathDB data access: http://cpdb.molgen.mpg.de/

Esta base de datos incluye tres tipos de datos que podemos descargar en tablas:	
	- Interacciones entre proteínas
	- Conjuntos de proteínas asociadas a biological pathways (UniProt)
	- Conjuntos de metabolitos asociados a biological pathways (Kegg / ChEBI / PubChem)
	  Muchas de estas asociaciones son extraídas de REACTOME, que vemos en la siguiente sección

* Para poder utilizar las asociaciones metabolitos-pathways vamos a necesitar mapear identificador
de kegg a chebi. Como vemos más adelante, a partir del identificador de chebi podemos sacar el nombre del compuesto.

El mapeo kegg-->chebi lo podemos hacer de dos formas:
	- https://rest.kegg.jp/conv/chebi/compound
	  Lo hemos sacado de las KEGG API (https://www.kegg.jp/kegg/rest/keggapi.html); Sección CONV
	- Paquete BioServices de Python: https://bioservices.readthedocs.io/en/latest/compound_tutorial.html

----

REACTOME: https://reactome.org/download-data

Recursos para mapear proteínas (UniprotID) y metabolitos (CHEBI) a rutas bioquímicas.

----

CHEBI: https://www.ebi.ac.uk/chebi/downloadsForward.do

REACTOME mapea los metabolitos utilizando el id de CHEBI, por lo que es preciso una
relación entre dicho ID y el nombre de los compuestos.

----

STITCH:

STITCH almacena relaciones/interacciones entre proteínas y metabolitos. 
Las proteínas se indican con identificador de ENS
Los metabolitos se indican con identificador de PubChem (precedido de un prefijo CIDm/CIDs/CID1)

http://stitch.embl.de/cgi/download.pl?UserId=rIxM13o39d3A&sessionId=i3tLi5pt8P2o

----

PubChem: https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest-tutorial#_Toc458584416

PubChem nos ofrece servicios REST para obtener todos los sinónimos de un compuesto químico
https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/glucose/synonyms/TXT

REST API: cualquiera de los sinónimos a PubChem id
https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/CHEBI:4167/cids/TXT

REST API: Del PubChem id a todos los sinónimos
https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/10000/synonyms/TXT

----

* El enfoque sería el siguiente: En la generación de los grafos de interacciones
los compuestos químicos deben estar representados por el PubChem ID:
	- Todos los nombres/identificadores (de kegg, de chebi...) se pueden expresar en pubchem con su REST
	- Los features se identifican putativamente, por lo que no tienen un identificador concreto de pubchem
	  sino que aglutinan varios. La identificación deberían ser múltiples PubChem plausibles (en un nodo)
	  y unidos a todos los pathways.

	- Esto apunta a la construcción de un grafo para la interpretación y confirmación biológica de posibles biomoléculas relevantes siguiendo este esquema:
		- Features significativas (por un proceso previo de selección) e identificadas putativamente
		- Las features (con sus pubchem) conectados a sus rutas según las bases de datos vistas
		- Las features conectadas entre sí por correlación (u otro criterio)
		- Proteínas significativas (por un proceso previo de seleción)
		- Proteínas conectadas a pathways y/o GO usando bases de datos
		- Proteínas conectadas entre sí por correlación (u otro criterio) y considerando STITCH
		- Proteínas y Metabolitos conectados por correlación (u otro criterio) y teniendo en cuenta la información y puntuación de STITCH

		--> Las componentes conexas 'completas' de este grafo podrían ofrecer una interpretación biológica de los elementos significativos
		--> También podría ayudar en la priorización de metabolitos candidatos

	- Previo a lo anterior, habría que aplicar diferentes modelos para la identificación de posibles biomarcadores o 
	  (elementos predictores: proteínas o metabolitos). Para ello podríamos aplicar los modelos de machine learning 
 	  que vemos en el siguiente apartado.




#
# MODELOS/WF PARA EL ANÁLISIS
#

---
Normalización de intensidades Proteómica:

Se debe realizar la integración s2p2q2all. Hay dos posibilidades:
1) Trabajar con Xq' extraídas del LowerNorm. Estos valores son X'q = Xq - Xall
por lo que podemos hacer estandarizar las X'q entre muestras.
X'q ya incluye la transformación logarítmica, y el denominador del ratio se anula
en la estandarización. Estos datos se asemejan a lo que se obtiene en metabolómica

2) Otra posibilidad es trabajar con Zq (preferencia de Jesús) por aprovechar toda la
potencia del modelo. Implica transformar los datos de metabolómica para normalizar
dentro de cada muestra entre proteínas.


El valor asignado a una proteína q en un individuo I, sería el sumatorio de las intensidades
de los scanes s_i que integran a q. Hay tres fuentes de variación entre individuos que habría que corregir:
	- Diferencias entre scanes que integran a la misma q (posición en el pico del scan, scan de péptidos que ionizan más...)
	- Diferencia de sistemática de intensidades entre canales de un mismo experimento/batch
	- Diferencia sistemática entre intensidades de diferentes experimentos

---
Missing Data Inputation
- Mantener features presentes en más de 80% control y 80% disease
- La estrategias principales para la imputación son:
	- KNN
	- Tomar datos aleatorios de una distribución normal (0,1), tras el estandarizado

---
ANÁLISIS CONJUNTO E INDEPENDIENTE DE PREDICTORES

Seguir esquema de mix-omics (aunque en Python¿?): 
https://mixomicsteam.github.io/Bookdown/intro.html

---

Selección de variables

Se debe distinguir Feature Selection de Feature Extraction.
Feature Selection selecciona las variables con mayor poder predictivo (depende del método)
	- Simplifica proceso de integración
	- Reduce complejidad, ruido y riesgo de overfitting en modelos de ML
	- Más interpretable

Feature Extraction calcula variables ficticias como una combinación lineal de
las variables originales que maximizan la información (depende del método)
	- Reduce complejidad, redundancia y ruido, y mejora rendimiento en modelos ML
	- Clustering
	- Permite visualizar los datos (EDA)
	- Permite mostrar features potencialmente significativas (medir significatividad estadística mediante pvalores empíricos)

En relación con este último punto, y para facilitar la interpretación de la 
reducción de dimensionalidad, son útiles los siguientes plots:
	- Correlation Circle Plots: Representar cada variable sobre un eje de coordenadas (cada eje representa un componente).
 	  El valor en cada eje es la correlación entre la variable y el componente.
	  En caso de tener variables centradas y escaladas, la correlación es la proyección de la variable sobre el componente.
	  El producto escalar está asociado a la correlación de variables. Como u.v = |u||v|cos(uv)
	  el ángulo ofrece información de la correlación entre variables (Esto permite comparar ómicas).
	  Ángulo agudo --> Correlación positiva
	  Ángulo obtuso --> Correlación negativa
	- Relevant Network: Grafo donde las variables son nodos y las aristas conexiones a partir de cierta correlación
	- Clustered Image Map: Heatmap con features de ambas ómicas. Se colorea en función de correlación y se agrupa/clusteriza por distancia.

Estas figuras/análisis precisan el cálculo de una matriz de correlaciones (aproximadas) entre las features
de una y otra ómica. El paper https://biodatamining.biomedcentral.com/articles/10.1186/1756-0381-5-19 explica
varias formas de obtener esta matriz mediante (r)CCA y (r)PLS(reg/can). Para ello, lo más inmediato es emplear
el paquete MixOmics, de R. Tal vez debas hacer un script para generar esta matriz de asociaciones en R usando este paquete.
De ahí, te los llevas a Python para representar y dibujar.
También sería útil calcular pvalores empíricos para determinar significatividad de las correlaciones.
- rCCA: http://mixomics.org/case-studies/rcca-nutrimouse-case-study/
- sPLS: http://mixomics.org/methods/spls/
	- sPLS-reg: http://mixomics.org/case-studies/spls-liver-toxicity-case-study/ 
	  --> Permite incluir múltiples datos clínicos de los pacientes

¿Para qué sirve esto?
	- Obtener un ranking de proteínas-metabolitos altamente correlacionados. Esto es, identificar subconjuntos de proteínas
	  y metabolitos altamente correlacionados. Esto es algo que solo se puede hacer con análisis multiómico, por lo que no es
	  comparable con un análisis independiente.
	- Construir un grafo bipartito que conecta Proteómica y Metabolómica (Hay que pensar otra forma de hacer conexiones intra)


Algunos muy usados son
- (Sparse) PCA
- (Sparse) CCA
- (Sparse) PLS: Esto se puede hacer sobre otras variables de los pacientes (tamaño de placa, fumador...)
- (Reg.) GLM: Ídem.
- Random Forest & Boosting (Ada-Boost): Ídem.


Los Feature Selection pueden ser:
- Filter-based: Filtro mediante análisis univariante (chi2, F-test)
- Wrapper: Construir un modelo ML con distintas features hasta encontrar la mejor combinación 
(Recursive feature elimination)
- Embedded: Obtener modelo ML que en su propia construcción selecciona features. Los LASSO, Ridge entran aquí


----

ANÁLISIS INTEGRAL DE AMBAS MATRICES ÓMICAS

Debemos tener cuidado con el desequilibrio generado por un número distinto de variables en cada ómica.

- (r)CCA: Para la comparación entre el bloque de features en Proteómica y Metabolómica.
Puedes seguir esquema y representación de mixomics:
http://mixomics.org/case-studies/rcca-nutrimouse-case-study/
Existe un PLS-Canonical (en lugar de regression), pero es muy parecido a CCA.

- Aplicar Multi-Block PLS (DIABLO en mix-omics)
Este esquema de librería en Python puede ser muy útil:
https://github.com/DTUComputeStatisticsAndDataAnalysis/MBPLS/blob/master/examples/real_world_applications/Carbohydrate_Microarray_PLS.ipynb

Puede ser interesante la sparsed version de DIABLO (seleccionar de antemano número de variables es problemático. Better en python tal vez¿?).

---
ANÁLISIS DEL GRAFO CONJUNTO Y MÉTRICAS ASOCIADAS


- Cómo construir el grafo:
	- Empezar la construcción del grafo empleando una combinación de correlaciones (Pearson, Spearman, Kendall)
	- Usar matrices de correlación generadas por modelos de mixomics. rCCA & rPLS-can
	- Usar también paper de MOTA (tras haber usado mixomics previamente)

- Con qué comparar el grafo:
	- Grafos aleatorios con orden y tamaño igual
	- Grafos generados con el mismo criterio pero realizando un shuffling de los datos

Construcción del grafo siguiendo esquema:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7241240/

Apoyarte en el paquete Mixomics para cálculo de rCCA
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3630015/

Paquete en Python para cálculo de rgCCA
https://cca-zoo.readthedocs.io/en/latest/index.html


- NetworkX para cálculo de métricas
- Pyvis para la representación interactiva


#
# PATHS A LOS DATOS
#

- Metabolómica (Incluye metadatos)
S:\U_Proteomica\PROYECTOS\PESA_omicas\METABOLOMICS\V1_afterLOESS


#
# OBJETIVOS DEL ESTUDIO / VALIDACIÓN DEL MODELO
#

Podemos definir dos objetivos en el estudio:

- Descripción del metabolismo global:
	- Descripción en individuos sanos de asociaciones metabolito-proteína-procesos biológicos
	- Alteración de estas asociaciones en individuos enfermos. (Influencia de otro factores hdl, placa... en metabolismo)
	. Este análisis se puede realizar construyendo grafo para control (referencia) y otro para enfermos.
	. MOTA genera un único grafo, pero reteniendo relaciones que se ven invertidas en ambos grupos.
	. Otro punto interesante podría ser generar grafos aleatorios que permita realizar contrastes en hipótesis nula.
	. Es interesante tener conocimiento de features asociadas a grupo o factores de riesgo, de cara a la comparación de grafos.

- Identificación de biomarcadores:
	- En este caso el objetivo sería identificar proteínas-metabolitos que permiten predecir la aterosclerosis
	  mejorando el resto de biomarcadores ya existentes. 
	- Aplicar modelos de regresión corrigiendo por los factores de riesgo establecidos e identificar significativos.
	- En este punto sería interesante identificar proteínas-metabolitos que interaccionen (efecto multiplicativo?).

- Debes concretar qué se entiende por estadística tradicional aplicada al análisis
  de Proteómica y Metabolómica para la identificación de biomoléculas con cambios
  estadísticamente significativos (no necesariamente biomarcadores).
  Puede ser ttest, ANOVA, modelo predictivo (regresión, clasificación), OPLS-DA

  . Estefanía: Filtra proteínas mediante análisis de correlación con tamaño de placa.
  Para cada una de las proteínas filtradas aplica regresión logística (predecir caso/control), de modo
  que el pvalor derivado de la regresión logística local permite identificar biomoléculas
  significativas.
  Este mismo criterio se puede aplicar a features de Metabolómica.
 

- Por otro lado, fijar uno o varios pipelines de integración (ya descritos,
  o combinación de los ya descritos) que conduzca a la identificación de biomoléculas
  con cambios estadísticamente significativos.
  No debemos restringirnos al estudio de relaciones entre la condición Sano/Enfermo 
  y el resto de features. Pueden emplear otras variables como
	. Glucosa
	. Diabetes
	. Fumador
	. HDL
	. Tamaño de placa
	. Tension (hipertension)
  De hecho, podría ser interesante analizar independencia de estas variables con 
  test chi-cuadrado + bonferroni

-->	Un aspecto clave del estudio será aplicar una comparativa que permita
	demostrar la bondad/alcance/virtud de la integración frente al análisis
	independiente. Algunas ideas para fundamentar/sostener esta demostración
	aparecen en el paper de MOTA:

	. Identificar subconjuntos de proteínas-metabolitos correlacionados entre sí. Tratar 
	  de interpretar biológicamente estas asociaciones. Se puede realizar este análisis sin
	  considerar la condición control/enfermo (aunque a posteriori se puede analizar su comportamiento
	  en diferentes factores de las observaciones). MOTA permitiría identificar correlaciones 
	  alteradas por la enfermedad directamente. Otra posibilidad es construir grafo con features
	  filtradas por significatividad u otro criterio.
	  No hay comparación posible con el análisis independiente.

	. Biomoléculas significativas exclusivas de la integración (y del análisis independiente) --> 
	  Por qué sale significativa y relacionar/buscar en la Literatura

	. Estabilidad de la integración entre diferentes cohortes -->
	  Comprobar si las biomóleculas significativas en la integración se conservan
	  en diferentes cohortes (en nuestro caso PESA y ARAGON).

	. Capacidad predictiva de las biomóleculas significativas -->
	  Comparar la capacidad predictiva, medida en AUC (o con otra métrica)
	  usando biomoléculas significativas en integración y en análisis independiente.
	  Estabilidad de la capacidad predictiva entre cohortes.

	. Análisis de enriquecimiento (GO en proteínas y Pathways en proteínas y metabolitos) 
	  empleando las biomoléculas de los diferentes modelos --> Se prueba la bondad 
	  del modelo de integración si aparecen GO/Pathways asociados a la enfermedad
	  de estudio (AT). Mirar topGO.


- Estrategias de validación del algoritmo de integración multiómica MOTA
	- Paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6986235/pdf/nihms-1026053.pdf
	- El objetivo es comparar el algoritmo de integración con métodos tradicionales
	  empleados en estadística, en este caso ttest (la comparación pretende 
	  resaltar virtudes del nuevo algoritmo, en este caso su mayor estabilidad).
	  
	  1) En este caso, el algoritmo de integración permite obtener un 
	  una puntuación y un ranking para los metabolitos, así como la construcción 
	  de un grafo inter-ómico, relacionado con la puntuación. Una forma que tienen de demostrar
	  la bondad del modelo es identificar metabolitos que son significativos en el modelo
	  pero no en los análisis tradicionales. La significatividad en el modelo
	  se justifica y razona mediante las conexiones del grafo, que a su vez se validan
	  con la bibliografía.
	  
	  2) Se aplican los diferentes modelos a comparar (ttest, idingo, MOTA) en diferentes cohortes.
	  Se comparan cuántas features son significativas en ambas cohortes para cada modelo. Se comprueba
	  que MOTA tiene 4, frente a 2 de las otras. Se intenta explicar el porqué de dicha estabilidad.
	  Se consideran solo features que se detectaron en ambas cohortes. 

	  3) Otra forma de realizar la comparación es tomar los 10 metabolitos 
	  más significativos en cada método.
	  A continuación se mide la capacidad predictiva de estos metabolitos 
	  mediante el cálculo del AUC asociado a un clasificador generado mediante Random Forest en cada caso.
	  Esto se realiza con 2 conjuntos de datos (cohortes) distintas (validación).

	  Llevado a nuestro estudio, un objetivo podría ser comparar los
	  biomarcadores obtenidos al analizar cada ómica de manera independiente empleando
	  los métodos estadísticos comúnmente utilizados (qué hacen en el paper de PESA)
	  Acto seguido, comparamos su capacidad predictiva mediante el AUC de algún clasificador.