# Tutorial

## Zebrafish

### import package

```python
from STMiner.SPFinder import SPFinder
```

### Load data

```python
file_path = 'I://10X_Visium_hunter2021spatially_sample_C_data.h5ad'
sp = SPFinder()
sp.read_h5ad(file=file_path, bin_size=1)
```

 - The parameter **bin_size** specifies the size of merged cells (spots). If not specified, no merging is performed. If set to 50, 50x50 cells/spots will be merged into a single cell/spot. Due to low sequencing depth in some datasets, cells/spots are often merged during analysis (e.g., stereo-seq). However, 10x data typically does not require merging.
 - The ST datasets was storaged in **.adata** object of **sp**, you can use **sp.adata** to check them:

```python
sp.adata
```

Besides, **sp** Obj has many useful attributes which can be used for visualization or integrated into other pipelines (such as scanpy). See [API](https://stminerdoc.readthedocs.io/en/latest/API/API.html) for more details.

### Find SVG

```python
sp.get_genes_csr_array(min_cells=50, log1p=False)
sp.spatial_high_variable_genes()
```
 - The parameter **min_cells** was used to filter genes that are too sparse to generate a reliable spatial distribution.
 - The parameter **log1p** was used to avoid extreme values affecting the results. For most open-source h5ad files, log1p has already been executed, so the default value here is False.
 - You can perform STMiner in your interested gene sets. Use parameter **gene_list** to input the gene list to STMiner. Then, STMiner will only calculate the given gene set of the dataset.
</br>
You can see output while computing as follows:

```bash
Parsing distance array...: 100%|██████████| 10762/10762 [01:12<00:00, 149.11it/s]
Computing ot distances...:  10%|▉         | 1069/10762 [03:04<31:11,  6.12it/s]  
```


### Fit GMM

```python
sp.fit_pattern(n_comp=20, gene_list=list(sp.global_distance[:2000]['Gene']))
```

You can see output while computing as follows:

```bash
Fitting GMM...:  10%|▉         | 190/2000 [00:42<04:36,  6.54it/s] 
```

 - **n_comp**： Number of components for each GMM model
 - **gene_list**: Gene list to fit GMM model
</br>

### Build distance array

```python
sp.build_distance_array()
```
This step calculates the distance between genes' spatial distributions.
You can visualize the distance array by:

```python
import seaborn as sns
sns.clustermap(sp.genes_distance_array)
```

<div align=center><img src="../_static/heatmap.png" width = "400"/></div>

### build distance matrix & clustering

```python
sp.cluster(n_clusters=6)
```

**n_clusters**: Number of cluster

### Result & Visualization

The result are stored in **genes_labels**:

```python
spf.genes_labels
```

The output looks like the following:

|     | gene_id        | labels |
|-----|----------------|--------|
| 0   | Cldn5          | 2      |
| 1   | Fyco1          | 2      |
| 2   | Pmepa1         | 2      |
| 3   | Arhgap5        | 0      |
| 4   | Apc            | 5      |
| ..  | ...            | ...    |
| 95  | Cyp2a5         | 0      |
| 96  | X5730403I07Rik | 0      |
| 97  | Ltbp2          | 2      |
| 98  | Rbp4           | 4      |
| 99  | Hist1h1e       | 4      |

To visualize the patterns by heatmap:

```python
sp.get_pattern_array()
sp.plot.plot_pattern(heatmap=False,
                     s=5,
                     rotate=True,
                     reverse_y=True,
                     reverse_x=True,
                     vmax=95,
                     cmap='Spectral_r',
                     output_path='./')
```

 - **heatmap**: If True, plot a heatmap. If False, plot a scatterplot. False is the default.
 - **s**: Spot size 
 - **rotate\reverse_y\reverse_x**: Adjust the axis of plot.
 - **cmap**: cmap of plot
 - **vmax**: The percentage of the highest value of plots. Avoid the effect of large values for visualization.
 - **output_path**: If set, save the figure to path
To visualize the genes by labels:

```python
sp.plot.plot_genes(label=0, n_gene=8, s=5, reverse_y=True, reverse_x=True)
```

**n_gene**: Number of genes to visualize


To visualize the specific gene (such as ***BRAFhuman***):
```python
hcc2l.plot.plot_gene('BRAFhuman', 
                     spot_size=10,
                     global_matrix_spot_size=10,
                     rotate=True, 
                     reverse_y=True, 
                     reverse_x=True, 
                     vmax=95, 
                     cmap='Spectral_r',
                     figsize=(5,5),
                     save_path='./',
                     format='png')
```

**reverse_y, reverse_x, rotate** is optional, they are used to adjust coordinate here.

---

## Hepatocellular carcinoma (HCC)

The HCC ST data can be download [here](http://lifeome.net/supp/livercancer-st/data.htm)

### import package

```python
from STMiner.SPFinder import SPFinder
```
### Load data

```python
file_path = 'I://HCC-1L.h5ad'
sp = SPFinder()
sp.read_h5ad(file=file_path)
```

### Find SVG

```python
hcc1l.fit_pattern(n_comp=20, gene_list=imm_a)
hcc1l.build_distance_array()
```

### Custom analysis (load interested gene set)
STMiner allows to input the genes or gene sets of interest and calculated the distance between all genes and the given gene/genes.

```python
import numpy as np
from STMiner.Algorithm.distance import compare_gmm_distance
from sklearn import mixture

# Input interested gene sets
imm_genes = ['CCL2','CCL3','CCL4','CCL5','CCL8','CCL18','CCL19','CCL21','CXCL9','CXCL10','CXCL11','CXCL13']
imm_genes_in_hcc1l = []
```

### Custom analysis (get patterns of interested gene set)
```python
for i in imm_genes:
    if i in list(hcc1l.adata.var.index):
        imm_genes_in_hcc1l.append(i)

# Get interested pattern
hcc1l.fit_pattern(n_comp=20, gene_list=imm_genes_in_hcc1l)
hcc1l.build_distance_array()
hcc1l.cluster_gene(n_clusters=1, mds_components=2)
hcc1l.get_pattern_array(vote_rate=0.2)
```

### Custom analysis (get patterns of all genes)
```python
def array_to_list(matrix) -> np.array:
    coords = np.column_stack(np.where(matrix > 0))
    counts = matrix[matrix > 0].flatten()
    result = np.repeat(coords, counts, axis=0)
    return result
gmm = mixture.GaussianMixture(n_components=20)
gmm.fit(array_to_list(np.round(hcc1l.patterns_matrix_dict[0]).astype(np.int32)))
```

### Custom analysis (compare all genes with interested gene set)
```python
df = compare_gmm_distance(gmm, hcc1l.patterns)
```
