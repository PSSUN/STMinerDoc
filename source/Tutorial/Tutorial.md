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
sp.read_h5ad(file=file_path)
```

### Find SVG

```python
sp.fit_pattern(n_comp=20, n_top_genes=200, min_cells=200)
```

### Fit GMM

```python
sp.fit_pattern(n_comp=10, gene_list=list(sp.global_distance[:2000]['Gene']))

```
**n_comp**： Number of components for each GMM model
</br>
**gene_list**: Gene list to fit GMM model
</br>


### Build distance array

```python
sp.build_distance_array()
```

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
sp.plot.plot_pattern(vmax=99, reverse_y=True, reverse_x=True, s=5, output_path='./Adult.eps')
```

**s**: Spot size 
**output_path**: If set, save the figure to path

To visualize the genes by labels:

```python
sp.plot.plot_genes(label=0, n_gene=8, s=5, reverse_y=True, reverse_x=True)
```
**n_gene**: Number of genes to visualize


### Additional function - Marked region

```python
# Open the GUI of STMiner
sp.app.run()
```
Load the image in UI:

<div style="text-align: center"><img src="../_static/t1.png" width="300" height="300" title="STMiner UI"><p align="center">Load the image</p></div>

Cut the image:

<div style="text-align: center"><img src="../_static/t2.png" width="300" height="300" title="STMiner UI"><p align="center">Cut the image</p></div>

Mark the image:
