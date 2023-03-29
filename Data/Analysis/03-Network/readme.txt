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