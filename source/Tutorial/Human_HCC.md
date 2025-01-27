## Human hepatocellular carcinoma

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

### Load interested gene set
STMiner allows to input the genes or gene sets of interest and calculated the distance between all genes and the given gene/genes.

```python
import numpy as np
from STMiner.Algorithm.distance import compare_gmm_distance
from sklearn import mixture

# Input interested gene sets
imm_genes = ['CCL2','CCL3','CCL4','CCL5','CCL8','CCL18','CCL19','CCL21','CXCL9','CXCL10','CXCL11','CXCL13']
imm_genes_in_hcc1l = []
for i in imm_genes:
    if i in list(hcc1l.adata.var.index):
        imm_genes_in_hcc1l.append(i)
```

### Custom analysis (get patterns of interested gene set)
```python
# Get interested pattern
hcc1l.fit_pattern(n_comp=20, gene_list=imm_genes_in_hcc1l)
hcc1l.build_distance_array()
hcc1l.cluster_gene(n_clusters=1, mds_components=2)
hcc1l.get_pattern_array(vote_rate=0.2)
```

### Get patterns of all genes
```python
def array_to_list(matrix) -> np.array:
    coords = np.column_stack(np.where(matrix > 0))
    counts = matrix[matrix > 0].flatten()
    result = np.repeat(coords, counts, axis=0)
    return result
gmm = mixture.GaussianMixture(n_components=20)
gmm.fit(array_to_list(np.round(hcc1l.patterns_matrix_dict[0]).astype(np.int32)))
```

### Cmpare all genes with interested gene set
```python
df = compare_gmm_distance(gmm, hcc1l.patterns)
```

{bdg-link-primary}`https://example.com`