# Introduction

## Abstract
Spatial transcriptomics revolutionizes transcriptomics by incorporating positional information. However, an emergency problem is to find out the gene expression pattern which can reveal the special region in tissue and find out the genes only expression in those regions. 

 Here we propose “STMiner” based on the Gaussian mixture model to solve this problem. STMiner is a bottom-up methodology algorithm. It is initiated by fitting a parametric model of gene spatial distributions and constructing a distance array between them utilizing the Hellinger distance. Genes are clustered, thereby recognizing spatial co-expression patterns across distinct gene classes.

STMiner is implemented as an open-source Python package and is available for use at [STMiner](https://github.com/PSSUN/STMiner).


## Algorithm

### Fit the Guassian Mixture Model

### Build th Distance Array

### Multidimensional scaling


<div><iframe allowtransparency="yes" frameborder="0" width="420" height="400" src="../_static/mds.html"/></iframe></div>


**Referance**:
 - [“Modern Multidimensional Scaling - Theory and Applications”](https://link.springer.com/book/10.1007/0-387-28981-X)
 - [“Nonmetric multidimensional scaling: a numerical method”](http://cda.psych.uiuc.edu/psychometrika_highly_cited_articles/kruskal_1964b.pdf)
 - [“Multidimensional scaling by optimizing goodness of fit to a nonmetric hypothesis”](http://cda.psych.uiuc.edu/psychometrika_highly_cited_articles/kruskal_1964a.pdf)
 - [wikipedia](https://en.wikipedia.org/wiki/Multidimensional_scaling)

### Cluster the genes


