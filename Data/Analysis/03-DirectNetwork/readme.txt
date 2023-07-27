#
# Network Analysis
#

. El objetivo de este apartado es construir los grafos Proteina-Metabolitos utilizando
  los datos de correlación y significatividad (pvalores) obtenidos en el apartado anterior.

. Todo grafo es un par de conjuntos (V, E), donde V es el conjunto de nodos y E el conjuntos
  de aristas. En este caso el conjunto V serán las proteínas y los metabolitos, mientras que
  una arista conecta dos nodos si las biomoléculas correspondientes están
  correlacionadas. En concreto, existe una arista entre dos biomoléculas si presentan una
  correlación diferencial significativa entre el grupo Control y Enfermo.

  El grafo contiene información sobre las correlaciones entre biomoléculas que se han visto
  alteradas por la condición Enfermo.

. Construimos primero un grafo correspondiente al proteoma y otro correspondiente al metaboloma.
  Las correlaciones empleadas para ella son las obtenidas a partir del Graphical Lasso.
  A continuación unimos estos dos grafos e incluimos las aristas entre proteínas y metabolitos
  a partir de las correlaciones obtenidas mediante el modelo rCCA.


Procedimiento para el análisis de grafo:

1) Métricas generales en cada uno de los grafos por separado: qq / mm / qm
    . Orden / Tamaño / Densidad / Índice de clusterización
    . Distribución de grado / Número de componentes conexas / Número de cliqués maximales
    --> Cuando se pueda, se comparan estas métricas con las obtenidas para grafos aleatorios.

2) Identificación de comunidades/clústeres estables en grafos qq y mm.
    . La identificación de comunidades la realizamos empleando el algoritmo de Leiden
      que se trata de una versión mejorada del algoritmo de Louvain.

    . Para garantizar la estabilidad de la comunidad aplicamos el método 
      conocido como "consenso de particiones" (Consensus Clustering). Consiste en aplicar
      el algoritmo para la detección de comunidades empleando diferentes semillas, de tal
      manera que retendremos aquellas comunidades que aparezcan en todas las ejecuciones.
      Para determinar el número de iteraciones, aplicamos el algoritmo a grafos aleatorios
      hasta que el >95% de los grafos no retengan ninguna comunidad con un tamaño superior
      al establecido.
      --> Este método tiene como principal inconveniente la pérdida de numerosas comunidades 
      potencialmente relevantes. Aunque tendremos la seguridad de que los clústeres obtenidos
      están interconectados de manera estable.
      https://www.nature.com/articles/srep00336
    
    * Debemos comprobar si las proteínas del paper de PESA están contenidas en
    alguno de los clústeres

3) Para cada uno de los clústeres calcularemos lo que en WGCNA llaman eigengene. Aquí lo 
denominamos autofeature. A partir de aquí podemos calcular la conexión entre clústeres
mediante la correlación (contraste de hipótesis incluido) de sus autofeatures.
A su vez, para cada clúster calculamos la asociación con los traits de los pacientes,
prestando especial atención a la condición Caso/Control. Esto genera tabla relevante.
Esta información permitirá construir una "meta-network", es decir, un grafo donde los
nodos sean clústeres (size del nodo proporcional a su asociación a la condición Caso/Control)
y las aristas representan correlaciones entre clústeres (width proprocional a la significatividad
de la correlación)

  * Es conveniente anotar para cada clúster el tipo de metabolitos/proteínas
  que aparecen y su asociación.
  * Sacar tablas con datos de clusters y pvalores
    
    Metabolomica 
      . Cluster 3: Asociado a smoke y hay una mayor presencia de compuestos halogenados
    
    Proteomica: Asociación directa a aterosclerosis subclínica los clusters 1,3,4
      . Cluster 1: Apolipoproteínas A,C e Ig de distinto tipo. Asociado a metabolismo
      de lípidos en STRING y complemento.
      . Cluster 3: Aparecen la haptoglobina, la PIGR y la C4, entre otras y está
      muy asociado a la condición Caso/Control, y a otros muchos metadatos. En STRING
      da complemento y respuesta inmune
      . Cluster 4: Es menos relevante. No da nada en STRING, y contiene algunas
      Ig.

4) Del apartado anterior nos interesan las conexiones entre clústeres qq y mm:
    . Seleccionar conexiones significativas (u otro criterio como número aristas...)
      Test hipergeométrico donde N es total de aristas posibles, D son aristas posibles
      entre ambos clústeres, n es número de aristas existentes y x es número de aristas
      existentes entre ambos clústeres.
    . Para cada "bi-clúster" calcular autofeature
      . Asociación con metadatos
      . Grado de pertenencia de cada nodo al bi-clúster (es por correlación, sacar pvalor)
        Analizar las listas obtenidas (aterrizaje)
      . De entre los nodos anteriores, indicar aquellos que están asociados.
    
5) Density plots que separen caso y control ¿?