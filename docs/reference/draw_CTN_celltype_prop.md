# CTN cell type proportion dotplot

`cowplot` and `ggtext` packages are required for the function.

## Usage

``` r
draw_CTN_celltype_prop(
  CTN_merged_celltype_prop,
  CTN_dend,
  celltype_palette,
  rename_celltype = waiver(),
  facet_var = "clust",
  x_expand = c(0.05, 0.05)
)
```

## Arguments

- CTN_merged_celltype_prop:

  A data frame representing the cell type proportion for each CTN with
  the following columns: `CTN` (the CTN names), `Celltype` (cell type
  labels), and `prop` (cell type proportion of the CTN).

- CTN_dend:

  Output of the [`draw_CTN_dendro()`](draw_CTN_dendro.md) function

- celltype_palette:

  A named vector with colors as values and cell type labels as names.
  Must contain all cell types.

- rename_celltype:

  A function for renaming cell type labels.

- facet_var:

  Faceting variables. Default to `"clust"` from `CTN_dend$label_df`

- x_expand:

  A numeric vector of length two of the multiplicative range expansion
  factors for the x axis. See
  [`ggplot2::expansion()`](https://ggplot2.tidyverse.org/reference/expansion.html)
  for details.

## Details

Dotplot visualization of cell type abundance within each cell-type triad
niche (CTN).

Dotplot visualization of the cell type abundance of cell-type triad
niches (CTN). Bubbles are colored by cell type labels, and their sizes
represent relative cell type proportions with each CTN.
