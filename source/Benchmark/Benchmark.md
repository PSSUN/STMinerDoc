# Benchmark for SVG detections

## Evaluate SVG

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
    df.to_csv('E://benchmark/benchmark/stminer_' + pre + '_50.csv')

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
    hs_results.to_csv('E://benchmark/benchmark/hotspots_' + pre + '_50.csv')

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

    print(pre)
    results = SpatialDE.run(np.array(X), sample_resid_expr)
    sign_results = results.query('qval < 0.05')
    sign_results = sign_results.sort_values(by=['FSV'], ascending=False)

    sign_results.to_csv('E://benchmark/spatialde_' + pre + '_50.csv')

if __name__ == '__main__':
    # current_dir is the h5ad file dir
    current_dir = 'E://benchmark_data'
    for root, _, files in os.walk(current_dir):
        for file in files:
            if file.endswith('h5ad'):
                h5ad_file_path = os.path.join(root, file)
                # Process the H5AD file here
                process_h5ad_file(h5ad_file_path)
                
```

## Compare SVG and non-SVG 

```python

top = 2000

for tag in [
    # file name in last step
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
    import numpy as np
    from tqdm import tqdm
    from sklearn.metrics.pairwise import (
        cosine_similarity, euclidean_distances, manhattan_distances, additive_chi2_kernel
    )

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
        
        p = 3
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
        
        p = 3
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
        
        p = 3
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
        
        p = 3
        minkowski_dist = np.sum(np.abs(X - Y) ** p) ** (1 / p)
        spatialde_result_df['Minkowski Distance'][gene] = minkowski_dist

    hotspot_result_df = hotspot_result_df.dropna()
    stminer_result_df = stminer_result_df.dropna()
    seurat_result_df = seurat_result_df.dropna()
    spatialde_result_df = spatialde_result_df.dropna()
    import seaborn as sns
    import matplotlib.pyplot as plt
    from statannot import add_stat_annotation

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
