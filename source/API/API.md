# API

## Attribute of STMiner

| Attribute           | Type         | Description                        |
|---------------------|--------------|------------------------------------|
| adata               | Anndata      | Anndata for loaded spatial data    |
| genes_patterns      | dict         | GMM model for each gene            |
| genes_distance_aray | pd.DataFrame | Distance between each GMM          |
| genes_labels        | pd.DataFrame | Gene name and their pattern labels |

## Methods of STMiner

### load data
- read_h5ad()
- read_gem()
- merge_bin()
### preprocess
- fit_pattern()
- normalize()
- log1p()
### fit model
- fit_pattern()
### build distance array & clustering
- cluster()
### visualization
- plot_pattern()
- plot_heatmap()