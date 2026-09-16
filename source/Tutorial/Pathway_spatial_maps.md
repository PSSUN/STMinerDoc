# Pattern enrichment and spatial pathway maps

This workflow enriches the gene sets of spatial patterns, then maps their actual
GO/KEGG-hit genes back onto the tissue. Start with a loaded `SPFinder` object
named `sp`, after `sp.cluster_gene()`. `get_pattern_array()` is not required.
These interfaces describe the current source checkout; older PyPI packages may
not include them. See [Installation](../Install/Install.md).

## Enrich each pattern

```python
enrichment = sp.enrich_patterns(
    organism="drerio",
    sources=["GO:BP", "KEGG"],
)
```

Choose the organism matching your data: `drerio` for zebrafish, `hsapiens` for
human, or `mmusculus` for mouse. These are g:Profiler names, distinct from KEGG
entry-query codes such as `hsa` and `mmu`.

Enrichment uses the existing `scanpy.queries.enrich()` / `gprofiler-official`
interface. It requires internet access and submits gene identifiers and the
background to g:Profiler. The default background is `sp.adata.var_names`; pass
`background=` with the genes eligible for selection if your tested universe is
narrower. Default correction is `g_SCS`; `correction_method="fdr"` or
`"bonferroni"` are also supported. Default sources are GO:BP, GO:MF, GO:CC, KEGG.

The result includes `pattern`, `query`, `source`, `native` (term ID), `name`,
adjusted `p_value`, and `intersections` (gene lists). Only significant terms are
returned unless `all_results=True`. Preserve the lists in `intersections`;
CSV round trips may turn them into strings requiring explicit parsing.

## Plot continuous spatial densities

```python
fig, axes, spot_scores, terms = sp.plot.plot_pathways(
    enrichment,
    top_n=2,
    min_genes=2,
    bandwidth=1.5,
    contours=3,
    output_path="spatial_pathways.svg",
    show=False,
)
```

The default view uses continuous color gradients and thin contours. Titles show
the pattern, term name/ID, adjusted P, and genes used/returned hits. Colors show
an **expression-weighted spatial density**, not a calibrated probability of
pathway activation. P values belong to pattern-level enrichment, not a spatial
significance test.

At each observation, the starting score is the arithmetic mean expression of
genes belonging to the returned term intersection, the current pattern, and
`adata.var_names`. Repeated grid positions are averaged for display. A Gaussian
kernel smooths this grid; masking and unit-sum normalization follow. Density
normalization removes absolute expression magnitude: compare spatial
concentration rather than activity strength. The same pathway in two patterns
may use different intersection genes.

Use normalized nonnegative expression in `adata.X`, or select an existing layer
with `layer="log_normalized"`. No automatic normalization is performed.
The source matrix, its precision, pattern assignments, and existing pattern
matrices are not modified.

| Parameter | Default | Meaning |
|---|---|---|
| `labels` | All patterns | Selected labels, e.g. `[0, 2]`. |
| `sources` | All returned sources | A string or list, e.g. `"KEGG"`. |
| `term_ids` | No restriction | Exact returned IDs, e.g. `["KEGG:04115"]`. |
| `top_n` | `2` | Maximum qualifying terms per pattern across selected sources, ranked by adjusted P; still applies with `term_ids`. |
| `p_value_cutoff` | `0.05` | Maximum adjusted P. |
| `min_genes` | `2` | Minimum mapped pattern-hit genes. |
| `bandwidth` | `1.5` | Gaussian width in grid units; larger values give broader features. |
| `contours` | `3` | Relative-peak contours (25%, 50%, 75% for three), not confidence or probability-mass contours. Zero hides them. |
| `support_radius` | `1.5` | Maximum distance from measured positions for automatic support, in grid units. |
| `tissue_mask` | `None` | Registered boolean `[x, y]` array matching the grid; overrides automatic support. |
| `shared_scale` | `True` | Shared color scale; False emphasizes individual distributions. |
| `display` | `"density"` | Use `"expression"` for unsmoothed mean-expression spots or heatmaps. |

Automatic support is not histological segmentation. Supply `tissue_mask` for a
precise outline. `image_path`, rotations, `num_cols`, `figsize`, `cmap`, `colorbar`,
`dpi`, and export follow `plot_pattern`. Images must already be aligned to the
untransformed grid. PDF, SVG, PNG, TIFF, and EPS suffixes are supported.

```python
# Compare unsmoothed expression with a common numerical scale.
fig, axes, spot_scores, terms = sp.plot.plot_pathways(
    enrichment, sources="KEGG", display="expression",
    heatmap=True, shared_scale=True, show=False,
)
```

## Inspect scores, genes, and missing mappings

`spot_scores` always contains unsmoothed mean-expression scores with the original
observation index. Its columns match `terms.index`. `terms` records `genes`,
`unmapped_genes`, `n_genes`, `n_hits`, `layer`, pattern, source, term ID/name, and
adjusted P. Attributes record omitted terms/patterns and display settings.

```python
terms[["pattern", "native", "genes", "unmapped_genes"]]
terms.attrs["skipped_terms"]
terms.attrs["unrepresented_patterns"]

# Calculate scores without creating a figure.
spot_scores, terms = sp.score_pathways(enrichment, sources="KEGG")
```

Patterns without qualifying terms are omitted and reported. If no terms remain,
an explicit error is raised. Rerun enrichment after changing cluster membership;
the stored membership snapshot is checked when available.

IDs are matched exactly, without changing case or guessing orthologues. For a
different namespace, use `gene_map={"returned_id": "adata_var_name"}`. For an
existing g:Profiler result without a `pattern` column, provide an explicit
`pattern_map={"query_name": pattern_label}` when scoring or plotting.

## Enrichment overview and KEGG queries

```python
from STMiner import KEGGFinder

finder = KEGGFinder(timeout=30)
fig, axes = finder.plot_enrichment(
    enrichment, sources=["GO:BP", "KEGG"], top_n=5,
    save_path="enrichment_overview.svg", show=False,
)

# Human entry query, independent of the spatial dataset above.
finder.find("hsa00010")
genes = finder.get_gene_dataframe(strict=True)
finder.find("map00010")
empty_genes = finder.get_gene_dataframe()
compounds = finder.get_section_dataframe("COMPOUND")
```

The overview bubble plot uses adjusted P for color and hit count for size; its
x-axis shows gene ratio for one query or cluster names for multiple queries.
Generic KEGG reference entries may lack a species-specific GENE section: the
default gene table is empty, whereas `strict=True` raises an explanatory error.

`KEGGFinder` also exposes `enrich_go`, `enrich_kegg`, `cluster_gene_set`,
`get_go_annotations`, `get_pathways_by_gene`, `get_pathway_network`,
`search_uniprot`, and `map_identifiers`. See the [API](../API/API.rst).
Close returned figures with `plt.close(fig)` after batch export.
