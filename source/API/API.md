# API

## Attribute of STMiner

| Attribute           | Type         | Description                        |
|---------------------|--------------|------------------------------------|
| adata               | Anndata      | Anndata for loaded spatial data    |
| genes_patterns      | dict         | GMM model for each gene            |
| genes_distance_aray | pd. DataFrame | Distance between each GMM          |
| genes_labels        | pd. DataFrame | Gene name and their pattern labels |

## Main Methods of STMiner 

### load data

* read_h5ad()
* read_stereo()
* read_bmk()

### preprocess

* get_genes_csr_array()

### fit model & clustering

* fit_pattern()
* build_distance_array()
* cluster()

### visualization

* plot_pattern()
* plot_intersection()
* plot_heatmap()
* plot_genes()
* plot_gene()
