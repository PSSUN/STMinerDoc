# Human hepatocellular carcinoma

The HCC ST data can be download {bdg-link-primary}`here <http://lifeome.net/supp/livercancer-st/data.htm>`.


## Import package

```python
import numpy as np
import pandas as pd
import scanpy as sc

from sklearn import mixture
from STMiner.SPFinder import SPFinder
from STMiner.Algorithm.distance import compare_gmm_distance
```

## Load data

```python
data = sc.read_10x_h5("I://HCC-5A/filtered_feature_bc_matrix.h5")  # Replace with your h5 file path
position=pd.read_csv("I://HCC-5A/spatial/tissue_positions_list.csv", header=None, index_col=0) # Replace with your tissue_positions_list.csv file path
position.columns = ['in_tissue','x','y','px','py']
data.obs = pd.merge(data.obs, position, left_index=True, right_index=True)
sc.pp.filter_genes(data, min_cells=50)
hcc = SPFinder(data) # Load anndata to STMiner
```

## Get patterns of interested gene set

STMiner allows to input the genes or gene sets of interest and calculated the distance between all genes and the given gene/genes.

```python
imm_genes =  ['CCL2','CCL3','CCL4','CCL5','CCL8','CCL18','CCL19','CCL21','CXCL9','CXCL10','CXCL11','CXCL13']
hcc.get_pattern_of_given_genes(gene_list=imm_genes)
```

## Cmpare all genes with interested gene set

```python
hcc.fit_pattern(n_comp=20) # Fit patterns of all genes
df = compare_gmm_distance(hcc.custom_pattern, hcc.patterns) # Compare the distance between all genes and the given gene set
```

## Optional pattern enrichment and spatial pathways

The custom-gene-set comparison above does not itself create gene clusters.
After fitting patterns, build the distance matrix and cluster genes before using
pattern-specific enrichment:

```python
hcc.build_distance_array()
hcc.cluster_gene(n_clusters=6)
enrichment = hcc.enrich_patterns(organism="hsapiens", sources=["GO:BP", "KEGG"])
fig, axes, spot_scores, terms = hcc.plot.plot_pathways(
    enrichment, bandwidth=1.5, contours=3,
    output_path="hcc_spatial_pathways.svg", show=False,
)
```

Use normalized nonnegative expression, or select an existing normalized layer
with `layer=`. The continuous view is an expression-weighted spatial distribution;
the returned spot scores retain unsmoothed mean expression. For the statistical
background, ID mapping, and full parameter descriptions, see
[Pattern enrichment and spatial pathway maps](Pathway_spatial_maps.md).
