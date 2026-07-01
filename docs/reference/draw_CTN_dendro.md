# Draw CTN dendrogram Dendrogram visualization for the hierarchical clustering results on the pairwise Jaccard distance between CTNs.

Draw CTN dendrogram Dendrogram visualization for the hierarchical
clustering results on the pairwise Jaccard distance between CTNs.

## Usage

``` r
draw_CTN_dendro(
  cl,
  h,
  k = NULL,
  CTN_annotation = NULL,
  celltype_palette = NULL,
  dot_fontsize = 12,
  CTN_group_palette = NULL,
  show_text = TRUE,
  x_expand = c(0.001, 0.001),
  y_lim = c(-2, 1),
  rename_celltype = NULL
)
```

## Arguments

- cl:

  An object of class "hclust" which describes the tree produced by the
  clustering process. Output of the
  [`jaccard_dist_hclust()`](jaccard_dist_hclust.md) function.

- h:

  Numeric scalar or vector with heights where the tree should be cut.
  [`stats::cutree()`](https://rdrr.io/r/stats/cutree.html) for details.
  At least one of `k` or `h` must be specified, `k` overrides `h` if
  both are given.

- k:

  An integer scalar or vector with the desired number of groups. See
  [`stats::cutree()`](https://rdrr.io/r/stats/cutree.html) for details.

- CTN_annotation:

  (Optional) A data frame with the original CTN names (`CTN`), the
  corresponding cluster label (`cluster`), and the new name for the
  consolidated CTN (`annotation_final`) after merging.

- celltype_palette:

  Input of the [`color_CTN_names()`](color_CTN_names.md) function. A
  named vector with colors as values and cell type labels as names. Must
  contain all cell types.

- dot_fontsize:

  Input of the [`color_CTN_names()`](color_CTN_names.md) function. Size
  of the dots. Default is 16pt.

- CTN_group_palette:

  (Optional) A named vector with colors as values and CTN categories as
  names. Only used when `CTN_annotation` is provided.

- show_text:

  Whether to display CTN names alongside the dendrogram (default:
  `TRUE`).

- x_expand:

  A numeric vector of length two of the multiplicative range expansion
  factors for the x axis. See
  [`ggplot2::expansion()`](https://ggplot2.tidyverse.org/reference/expansion.html)
  for details.

- y_lim:

  A numeric vector of length two providing limits of the scale. Use NA
  to refer to the existing minimum or maximum.

- rename_celltype:

  A function for renaming cell type labels.

## Value

A list with components:

- p:

  Dendrogram

- label_df:

  Label data. See
  [`ggdendro::dendro_data()`](https://andrie.github.io/ggdendro/reference/dendro_data.html)
  for details.
