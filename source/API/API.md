# API

## Main Attributes of STMiner

| Attributes           | Type          | Description                            |
| -------------------- | ------------- | -------------------------------------- |
| adata                | Anndata       | Anndata for loaded spatial data        |
| patterns             | dict          | Spatial distributions pattern of genes |
| genes_patterns       | dict          | GMM model for each gene                |
| global_distance      | pd. DataFrame | Distances between genes and background |
| mds_features         | array         | embedding features of genes            |
| genes_distance_array | pd. DataFrame | Distance between each GMM              |
| genes_labels         | pd. DataFrame | Gene name and their pattern labels     |
| plot                 | Object        | Call plot to visualization             |

---

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

* plot.plot_pattern()
* plot.plot_intersection()
* plot.plot_heatmap()
* plot.plot_genes()
* plot.plot_gene()
