# Benchmark for SVG detections

## Evaluate SVG
**NOTE:** To install other packages such as SpatialDE/hotspot/SpaGFT, please refer to their official site. 

### SpaGFT
```python
import pandas as pd
import numpy as np
import scanpy as sc
import SpaGFT as spg

def process_h5ad_file(file):
    spagft_svg(file)


def spagft_svg(file):
    adata = sc.read_h5ad(file)    
    sc.pp.filter_genes(adata, min_cells=50)
    name = os.path.basename(file)
    pre = name.split('.')[0]
    adata.obs.loc[:, ['array_row', 'array_col']] = adata.obsm['spatial']
    (ratio_low, ratio_high) = spg.gft.determine_frequency_ratio(adata, ratio_neighbors=1)    
    df_res = spg.detect_svg(adata,
                        spatial_info=['array_row', 'array_col'],
                        ratio_low_freq=ratio_low,
                        ratio_high_freq=ratio_high,
                        ratio_neighbors=1,
                        filter_peaks=True,
                        S=6)

    df_res = df_res.loc[adata.var_names]
    df = df_res.sort_values(by='svg_rank', ascending=True)
    df.to_csv('./benchmark/spagft_' + pre + '_50.csv')


if __name__ == '__main__':
    current_dir = 'I://mm10'
    for root, _, files in os.walk(current_dir):
        for file in files:
            if file.endswith('h5ad'):
                h5ad_file_path = os.path.join(root, file)
                # Process the H5AD file here
                process_h5ad_file(h5ad_file_path)
```

### scGCO
```python
import pandas as pd
import numpy as np
import scanpy as sc
from scGCO import *

def process_h5ad_file(file):
    adata = sc.read_h5ad(file)
    name = os.path.basename(file)
    pre = name.split('.')[0]
    sc.pp.filter_genes(adata, min_cells=50)
    data = pd.DataFrame(
        adata.X.todense(), columns=adata.var_names, index=adata.obs_names
    )
    data_norm = normalize_count_cellranger(data)
    exp = data.iloc[:, 0]
    locs = adata.obsm['spatial'].copy()
    cellGraph = create_graph_with_weight(locs, exp)
    gmmDict= gmm_model(data_norm)

    df_res = identify_spatial_genes(locs, data_norm, cellGraph,gmmDict)
    df_res = df_res.loc[adata.var_names]
    df = pd.concat([df_res, adata.var], axis=1)
    df.to_csv('/data/home/pssun/scgco_' + pre + '_50.csv')


if __name__ == '__main__':
    current_dir = '/data/home/pssun/mm10/'
    for root, _, files in os.walk(current_dir):
        for file in files:
            if file.endswith('h5ad'):
                h5ad_file_path = os.path.join(root, file)
                # Process the H5AD file here
                process_h5ad_file(h5ad_file_path)
```

### nnSVG
```R
library(nnSVG)
library(anndata)
library(SpatialExperiment)
library(scran)

reticulate::use_condaenv("base", required = TRUE)


csv_files <- list.files(path = '/data/home/st/', pattern = "\\.h5ad$", full.names = TRUE)
out_dir <- "/data/home/result/"

for (file in csv_files) {
    file_name <- tools::file_path_sans_ext(basename(file))
    print(file_name)
    adata <- anndata::read_h5ad(file)
    counts <- t(as.matrix(adata$X))
    colnames(counts) <- adata$obs_names
    rownames(counts) <- adata$var_names
    loc <- as.data.frame(adata$obsm[['spatial']])
    row_data = adata$var
    row_data$gene_id = rownames(row_data)
    row_data$feature_type = "Gene Expression"
    colnames(loc) <- c("x", "y")
    rownames(loc) <- colnames(counts)
    spe <- SpatialExperiment(
      assays = list(counts = counts),
      rowData = row_data,
      colData = loc,
      spatialCoordsNames = c("x", "y"))
    spe <- computeLibraryFactors(spe)
    spe <- logNormCounts(spe)

    # filter any new zeros created after filtering low-expressed genes
    # remove genes with zero expression
    # remove spots with zero expression
    ix_zero_spots <- colSums(counts(spe)) == 0
    table(ix_zero_spots)
    if (sum(ix_zero_spots) > 0) {
      spe <- spe[, !ix_zero_spots]
    }
    dim(spe)
    ix_zero_genes <- rowSums(counts(spe)) <= 50
    table(ix_zero_genes)

    if (sum(ix_zero_genes) > 0) {
      spe <- spe[!ix_zero_genes, ]
    }
    dim(spe)

    print('Start')
    spe <- nnSVG(spe, n_threads=48)
    df <- rowData(spe)
    output_file <- file.path(out_dir, paste0("nnsvg_", file_name, ".csv"))
    write.csv(df, file=output_file, quote=FALSE)
}
```

### SPARK
```R
library(SPARK)
library(Matrix)


dir_path <- "/data/home/st/"
csv_files <- list.files(path = dir_path, pattern = "\\.csv$", full.names = TRUE)
out_dir <- "/data/home/result/"

for (file in csv_files) {
  count <- read.csv(file, row.names = 1, check.names = FALSE)
  file_name <- tools::file_path_sans_ext(basename(file))
  print(file_name)
  tmp <- as.matrix(count)
  info <- cbind.data.frame(
    x = as.numeric(sapply(strsplit(colnames(tmp), split = "x"), "[", 1)),
    y = as.numeric(sapply(strsplit(colnames(tmp), split = "x"), "[", 2)),
    total_counts = apply(tmp, 2, sum)
  )
  rownames(info) <- colnames(tmp)
  location <- as.matrix(info)
  counts_sparse <- as.matrix(tmp)
  spark <- CreateSPARKObject(counts=count, percentage = 0, min_total_counts = 0, location=info[, 1:2])
  spark@lib_size <- apply(spark@counts, 2, sum)
  spark <- spark.vc(spark,
                    covariates = NULL,
                    lib_size = spark@lib_size,
                    num_core = 24,
                    verbose = F)
  spark <- spark.test(spark,
                      check_positive = T,
                      verbose = F)
  df <- as.data.frame(spark@res_mtest)
  output_file <- file.path(out_dir, paste0("spark_", file_name, ".csv"))
  write.csv(df, file=output_file, quote=FALSE)
  df[is.na(df)] <- 1
  print('-----------------------------------------')
}

```

### trendsceek
```R
library(trendsceek)
library(Matrix)


dir_path <- "G://trendsceek/"
csv_files <- list.files(path = dir_path, pattern = "\\.csv$", full.names = TRUE)
out_dir <- "E://benchmark"

for (file in csv_files) {
  count <- read.csv(file, row.names = 1, check.names = FALSE)
  file_name <- tools::file_path_sans_ext(basename(file))
  print(file_name)
  counts_filt <- genefilter_exprmat(count,
  min.expr = 3,
    min.ncells.expr = 50)
  vargenes_stats <- calc_varstats(counts_filt,
    counts_filt,
    quant.cutoff = 0.5,
   method = "glm"
  )
  topvar.genes <- rownames(vargenes_stats[["real.stats"]])[1:2000]
  output_file <- file.path(out_dir, paste0("trendsceek_", file_name, ".csv"))
  write.csv(topvar.genes, output_file)
}
```

### SPARKX
```R
library(SPARK)
library(Matrix)


dir_path <- "G://sparkx/"
csv_files <- list.files(path = dir_path, pattern = "\\.csv$", full.names = TRUE)
out_dir <- "E://benchmark"

for (file in csv_files) {
  count <- read.csv(file, row.names = 1, check.names = FALSE)
  file_name <- tools::file_path_sans_ext(basename(file))
  print(file_name)

  tmp <- as.matrix(count)
  info <- cbind.data.frame(
    x = as.numeric(sapply(strsplit(colnames(tmp), split = "x"), "[", 1)),
    y = as.numeric(sapply(strsplit(colnames(tmp), split = "x"), "[", 2)),
    total_counts = apply(tmp, 2, sum)
  )
  rownames(info) <- colnames(tmp)
  location <- as.matrix(info)
  counts_sparse <- as.matrix(tmp)
  sparkX <- sparkx(counts_sparse, location, numCores = 1, option = "mixture")
  df <- as.data.frame(sparkX$res_mtest)
  df <- df[order(df$adjustedPval), ]
  output_file <- file.path(out_dir, paste0("sparkx_", file_name, ".csv"))
  write.csv(df, output_file)
}
```



### SpatialDE/SpatialDE2/Hotspot/STMiner
```python
import os
import NaiveDE
import SpatialDE
import hotspot
import scanpy as sc
import pandas as pd
import numpy as np
from STMiner import SPFinder

def process_h5ad_file(file):
    hotspot_svg(file)
    spatial_de_svg(file)
    stminer_svg(file)

def stminer_svg(file):
    name = os.path.basename(file)
    print(name)
    pre = name.split('.')[0]
    sp = SPFinder()
    sp.read_h5ad(file, bin_size=1)
    sp.get_genes_csr_array(min_cells=50, log1p=True)
    sp.spatial_high_variable_genes()
    df = sp.global_distance
    df.to_csv('E://benchmark/stminer_' + pre + '_50.csv')

def hotspot_svg(file):
    adata = sc.read_h5ad(file)
    name = os.path.basename(file)
    pre = name.split('.')[0]
    sc.pp.filter_genes(adata, min_cells=50)
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata.var['Gene'] = list(adata.var.index)
    adata.var.query('highly_variable==True').loc[:, ['Gene', 'dispersions_norm']].to_csv(
        'E://benchmark/benchmark/seurat_' + pre + '_50.csv')

    hs = hotspot.Hotspot(
        adata,
        model='danb',
        latent_obsm_key="X_pca",
        umi_counts_obs_key="total_counts"
    )

    hs.create_knn_graph(weighted_graph=False, n_neighbors=50)
    hs_results = hs.compute_autocorrelations()
    hs_results.to_csv('E://benchmark/hotspots_' + pre + '_50.csv')

def spatial_de_svg(file):
    name = os.path.basename(file)
    pre = name.split('.')[0]
    adata = sc.read_h5ad(file)
    sc.pp.filter_genes(adata, min_cells=50)
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.obs['x'] = adata.obsm['spatial'][:, 0]
    adata.obs['y'] = adata.obsm['spatial'][:, 1]

    exp_df = pd.DataFrame(adata.X.todense(), index=adata.obs.index, columns=adata.var.index)
    norm_expr = NaiveDE.stabilize(exp_df.T).T
    resid_expr = NaiveDE.regress_out(adata.obs, norm_expr.T, 'np.log(total_counts)').T
    sample_resid_expr = resid_expr.sample(n=len(adata.var), axis=1, random_state=1)
    X = adata.obs[['x', 'y']]
    results = SpatialDE.run(np.array(X), sample_resid_expr)
    sign_results = results.query('qval < 0.05')
    sign_results = sign_results.sort_values(by=['FSV'], ascending=False)

    sign_results.to_csv('E://benchmark/spatialde_' + pre + '_50.csv')

def spatial2_de_svg(file):
    name = os.path.basename(file)
    pre = name.split('.')[0]
    adata = sc.read_h5ad(file)
    sc.pp.filter_genes(adata, min_cells=500)
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.obs['x'] = adata.obsm['spatial'][:, 0]
    adata.obs['y'] = adata.obsm['spatial'][:, 1]
    total_counts = sc.get.obs_df(adata, keys=["total_counts"])
    exp_df = pd.DataFrame(adata.X.todense(), index=adata.obs.index, columns=adata.var.index)
    norm_expr = NaiveDE.stabilize(exp_df.T).T
    adata.X = NaiveDE.regress_out(total_counts, norm_expr.T, 'np.log(total_counts)').T
    X = adata.obs[['x', 'y']]
    df_res = SpatialDE.fit(adata, normalized=True, control=None)
    df_res.set_index("gene", inplace=True)
    df_res = df_res.loc[adata.var_names]

    df = pd.concat([df_res, adata.var], axis=1)
    df.to_csv('E://benchmark/spatialde2_' + pre + '_50.csv')


if __name__ == '__main__':
    # current_dir is the h5ad file directory.
    current_dir = 'E://benchmark_data'
    for root, _, files in os.walk(current_dir):
        for file in files:
            if file.endswith('h5ad'):
                h5ad_file_path = os.path.join(root, file)
                # Process the H5AD file here
                process_h5ad_file(h5ad_file_path)
                
```

## Compare the distance between SVGs and non-SVGs

```python
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from tqdm import tqdm
from statannot import add_stat_annotation
from sklearn.metrics.pairwise import (
    cosine_similarity, euclidean_distances, manhattan_distances, additive_chi2_kernel
)
top = 2000
p = 3
for tag in [
    # File names in last step.
]:
    print(tag)
    hotspot = list(pd.read_csv(
        f'E://benchmark\hotspots_{tag}_50.csv')[:top]['Gene'])
    seurat = list(pd.read_csv(
        f'E://benchmark\seurat_{tag}_50.csv').sort_values(
        by='dispersions_norm', ascending=False)[:top]['Gene'])
    stminer = list(
        pd.read_csv(f'E://benchmark//stminer_{tag}_50.csv')[:top]['Gene'])
    spatialde = list(pd.read_csv(
        f'E://benchmark\spatialde_{tag}_50.csv')[:top]['g'])
    tot = list(set(hotspot + seurat + stminer + spatialde))
    sp = SPFinder()
    sp.read_h5ad(f"E://benchmark_data//{tag}.h5ad", bin_size=10)
    sp.get_genes_csr_array(min_cells=50)

    total_dict = {}
    for i in tqdm(tot, desc='Get csr...'):
        total_dict[i] = csr_matrix(np.array(scale_matrix(get_exp_array(sp.adata, i)).flatten().reshape(1, -1)))


    global_arr = scale_matrix(sp.plot.get_global_matrix(False, False, False).todense())
    X = np.array(global_arr.flatten().reshape(1, -1))

    ###############################################################################################################
    hotspot_result_df = pd.DataFrame(index=list(set(hotspot + seurat + stminer + spatialde)),
                                     columns=['1/Cosine Similarity', 'Euclidean Distances', 'Manhattan Distances',
                                              '-Additive Chi-Square Kernel', 'Minkowski Distance'])
    for gene in tqdm(hotspot, desc='hotspot'):
        Y = total_dict[gene].toarray()
        
        cosine_sim = cosine_similarity(X, Y)[0][0]
        hotspot_result_df['1/Cosine Similarity'][gene] = 1 / cosine_sim
        
        euclidean_dist = euclidean_distances(X, Y)[0][0]
        hotspot_result_df['Euclidean Distances'][gene] = euclidean_dist
        
        manhattan_dist = manhattan_distances(X, Y)[0][0]
        hotspot_result_df['Manhattan Distances'][gene] = manhattan_dist
        
        additive_chi2_ker = additive_chi2_kernel(X, Y)[0][0]
        hotspot_result_df['-Additive Chi-Square Kernel'][gene] = -additive_chi2_ker
        
        minkowski_dist = np.sum(np.abs(X - Y) ** p) ** (1 / p)
        hotspot_result_df['Minkowski Distance'][gene] = minkowski_dist

    ###############################################################################################################
    print('stminer')
    stminer_result_df = pd.DataFrame(index=list(set(hotspot + seurat + stminer + spatialde)),
                                     columns=['1/Cosine Similarity', 'Euclidean Distances', 'Manhattan Distances',
                                              '-Additive Chi-Square Kernel', 'Minkowski Distance'])
    for gene in stminer:
        Y = total_dict[gene].toarray()
        
        cosine_sim = cosine_similarity(X, Y)[0][0]
        stminer_result_df['1/Cosine Similarity'][gene] = 1 / cosine_sim
        
        euclidean_dist = euclidean_distances(X, Y)[0][0]
        stminer_result_df['Euclidean Distances'][gene] = euclidean_dist
        
        manhattan_dist = manhattan_distances(X, Y)[0][0]
        stminer_result_df['Manhattan Distances'][gene] = manhattan_dist
        
        additive_chi2_ker = additive_chi2_kernel(X, Y)[0][0]
        stminer_result_df['-Additive Chi-Square Kernel'][gene] = -additive_chi2_ker
        
        minkowski_dist = np.sum(np.abs(X - Y) ** p) ** (1 / p)
        stminer_result_df['Minkowski Distance'][gene] = minkowski_dist

    ###############################################################################################################
    print('seurat')
    seurat_result_df = pd.DataFrame(index=list(set(hotspot + seurat + stminer + spatialde)),
                                    columns=['1/Cosine Similarity', 'Euclidean Distances', 'Manhattan Distances',
                                             '-Additive Chi-Square Kernel', 'Minkowski Distance'])
    for gene in seurat:
        Y = total_dict[gene].toarray()
        
        cosine_sim = cosine_similarity(X, Y)[0][0]
        seurat_result_df['1/Cosine Similarity'][gene] = 1 / cosine_sim
        
        euclidean_dist = euclidean_distances(X, Y)[0][0]
        seurat_result_df['Euclidean Distances'][gene] = euclidean_dist
        
        manhattan_dist = manhattan_distances(X, Y)[0][0]
        seurat_result_df['Manhattan Distances'][gene] = manhattan_dist
        
        additive_chi2_ker = additive_chi2_kernel(X, Y)[0][0]
        seurat_result_df['-Additive Chi-Square Kernel'][gene] = -additive_chi2_ker
        
        minkowski_dist = np.sum(np.abs(X - Y) ** p) ** (1 / p)
        seurat_result_df['Minkowski Distance'][gene] = minkowski_dist

    ###############################################################################################################
    spatialde_result_df = pd.DataFrame(index=list(set(hotspot + seurat + stminer + spatialde)),
                                       columns=['1/Cosine Similarity', 'Euclidean Distances', 'Manhattan Distances',
                                                '-Additive Chi-Square Kernel', 'Minkowski Distance'])
    for gene in spatialde:
        Y = total_dict[gene].toarray()
        
        cosine_sim = cosine_similarity(X, Y)[0][0]
        spatialde_result_df['1/Cosine Similarity'][gene] = 1 / cosine_sim
        
        euclidean_dist = euclidean_distances(X, Y)[0][0]
        spatialde_result_df['Euclidean Distances'][gene] = euclidean_dist
        
        manhattan_dist = manhattan_distances(X, Y)[0][0]
        spatialde_result_df['Manhattan Distances'][gene] = manhattan_dist
        
        additive_chi2_ker = additive_chi2_kernel(X, Y)[0][0]
        spatialde_result_df['-Additive Chi-Square Kernel'][gene] = -additive_chi2_ker
        
        minkowski_dist = np.sum(np.abs(X - Y) ** p) ** (1 / p)
        spatialde_result_df['Minkowski Distance'][gene] = minkowski_dist

    hotspot_result_df = hotspot_result_df.dropna()
    stminer_result_df = stminer_result_df.dropna()
    seurat_result_df = seurat_result_df.dropna()
    spatialde_result_df = spatialde_result_df.dropna()


    fig, axes = plt.subplots(1, 5, figsize=(12, 3.5))
    colors = ['#a07f64', '#a9bc76', '#d99b3f', '#5862ab']

    sns.boxplot([list(stminer_result_df['Euclidean Distances']),
                 list(hotspot_result_df['Euclidean Distances']),
                 list(seurat_result_df['Euclidean Distances']),
                 list(spatialde_result_df['Euclidean Distances'])], showfliers=False, ax=axes[0], width=.5,
                palette=colors)
    axes[0].set_ylabel('Euclidean Distances', fontname="Arial", fontsize=10)
    # axes[0].set_title('Euclidean Distances')
    # axes[0].set_ylim(0, 250)
    data_eu = pd.DataFrame({
        'Method': ['stMINER'] * len(stminer_result_df) +
                  ['Hotspot'] * len(hotspot_result_df) +
                  ['Seurat'] * len(seurat_result_df) +
                  ['SpatialDE'] * len(spatialde_result_df),
        'Euclidean Distances': list(stminer_result_df['Euclidean Distances']) +
                               list(hotspot_result_df['Euclidean Distances']) +
                               list(seurat_result_df['Euclidean Distances']) +
                               list(spatialde_result_df['Euclidean Distances'])
    })
    add_stat_annotation(axes[0], data=data_eu, x='Method', y='Euclidean Distances',
                        box_pairs=[('stMINER', 'Hotspot'),
                                   ('stMINER', 'Seurat'),
                                   ('stMINER', 'SpatialDE')],
                        test='Mann-Whitney', text_format='star', loc='outside', verbose=2)
    #######################################################################################################

    sns.boxplot([list(stminer_result_df['1/Cosine Similarity']),
                 list(hotspot_result_df['1/Cosine Similarity']),
                 list(seurat_result_df['1/Cosine Similarity']),
                 list(spatialde_result_df['1/Cosine Similarity'])], showfliers=False, ax=axes[1], width=.5,
                palette=colors)
    axes[1].set_ylabel('1/Cosine Similarity', fontname="Arial", fontsize=10)
    data = pd.DataFrame({
        'Method': ['stMINER'] * len(stminer_result_df) +
                  ['Hotspot'] * len(hotspot_result_df) +
                  ['Seurat'] * len(seurat_result_df) +
                  ['SpatialDE'] * len(spatialde_result_df),
        '1/Cosine Similarity': list(stminer_result_df['1/Cosine Similarity']) +
                               list(hotspot_result_df['1/Cosine Similarity']) +
                               list(seurat_result_df['1/Cosine Similarity']) +
                               list(spatialde_result_df['1/Cosine Similarity'])
    })
    add_stat_annotation(axes[1], data=data, x='Method', y='1/Cosine Similarity',
                        box_pairs=[('stMINER', 'Hotspot'),
                                   ('stMINER', 'Seurat'),
                                   ('stMINER', 'SpatialDE')],
                        test='Mann-Whitney', text_format='star', loc='outside', verbose=2)
    #######################################################################################################
    sns.boxplot([list(stminer_result_df['Manhattan Distances']),
                 list(hotspot_result_df['Manhattan Distances']),
                 list(seurat_result_df['Manhattan Distances']),
                 list(spatialde_result_df['Manhattan Distances'])], showfliers=False, ax=axes[3], width=.5,
                palette=colors)
    axes[3].set_ylabel('Manhattan Distances', fontname="Arial", fontsize=10)
    data_Manhattan = pd.DataFrame({
        'Method': ['stMINER'] * len(stminer_result_df) +
                  ['Hotspot'] * len(hotspot_result_df) +
                  ['Seurat'] * len(seurat_result_df) +
                  ['SpatialDE'] * len(spatialde_result_df),
        'Manhattan Distances': list(stminer_result_df['Manhattan Distances']) +
                               list(hotspot_result_df['Manhattan Distances']) +
                               list(seurat_result_df['Manhattan Distances']) +
                               list(spatialde_result_df['Manhattan Distances'])
    })
    add_stat_annotation(axes[3], data=data_Manhattan, x='Method', y='Manhattan Distances',
                        box_pairs=[('stMINER', 'Hotspot'),
                                   ('stMINER', 'Seurat'),
                                   ('stMINER', 'SpatialDE')],
                        test='Mann-Whitney', text_format='star', loc='outside', verbose=2)
    ######################################################################################################################
    plt.grid(False)
    sns.boxplot([list(stminer_result_df['-Additive Chi-Square Kernel']),
                 list(hotspot_result_df['-Additive Chi-Square Kernel']),
                 list(seurat_result_df['-Additive Chi-Square Kernel']),
                 list(spatialde_result_df['-Additive Chi-Square Kernel'])], showfliers=False, ax=axes[4], width=.5,
                palette=colors)
    axes[4].set_ylabel('-Additive Chi-Square Kernel', fontname="Arial", fontsize=10)

    data_Manhattan = pd.DataFrame({
        'Method': ['stMINER'] * len(stminer_result_df) +
                  ['Hotspot'] * len(hotspot_result_df) +
                  ['Seurat'] * len(seurat_result_df) +
                  ['SpatialDE'] * len(spatialde_result_df),
        '-Additive Chi-Square Kernel': list(stminer_result_df['-Additive Chi-Square Kernel']) +
                                       list(hotspot_result_df['-Additive Chi-Square Kernel']) +
                                       list(seurat_result_df['-Additive Chi-Square Kernel']) +
                                       list(spatialde_result_df['-Additive Chi-Square Kernel'])
    })
    add_stat_annotation(axes[4], data=data_Manhattan, x='Method', y='-Additive Chi-Square Kernel',
                        box_pairs=[('stMINER', 'Hotspot'),
                                   ('stMINER', 'Seurat'),
                                   ('stMINER', 'SpatialDE')],
                        test='Mann-Whitney', text_format='star', loc='outside', verbose=2)
    ##################################################################################################
    sns.boxplot([list(stminer_result_df['Minkowski Distance']),
                 list(hotspot_result_df['Minkowski Distance']),
                 list(seurat_result_df['Minkowski Distance']),
                 list(spatialde_result_df['Minkowski Distance'])], showfliers=False, ax=axes[2], width=.5,
                palette=colors)
    axes[2].set_ylabel('Minkowski Distance', fontname="Arial", fontsize=10)
    data_Manhattan = pd.DataFrame({
        'Method': ['stMINER'] * len(stminer_result_df) +
                  ['Hotspot'] * len(hotspot_result_df) +
                  ['Seurat'] * len(seurat_result_df) +
                  ['SpatialDE'] * len(spatialde_result_df),
        'Minkowski Distance': list(stminer_result_df['Minkowski Distance']) +
                              list(hotspot_result_df['Minkowski Distance']) +
                              list(seurat_result_df['Minkowski Distance']) +
                              list(spatialde_result_df['Minkowski Distance'])
    })
    add_stat_annotation(axes[2], data=data_Manhattan, x='Method', y='Minkowski Distance',
                        box_pairs=[('stMINER', 'Hotspot'),
                                   ('stMINER', 'Seurat'),
                                   ('stMINER', 'SpatialDE')],
                        test='Mann-Whitney', text_format='star', loc='outside', verbose=2)
    sns.set_style("white")
    for ax in axes.flatten():
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_xticklabels(['STMiner', 'Hotspot', 'Seurat', 'SpatialDE'], fontname="Arial", fontsize=8)
    axes[0].text(-0.45, 0.5, tag.split('_spaceranger')[0], rotation=90, transform=axes[0].transAxes, va='center', fontsize=10)
    plt.tight_layout()
    plt.savefig(f"{tag}.eps", format='eps')
    plt.show()
```
