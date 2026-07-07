# CTN cell type proportion dotplot

Dotplot visualization of cell type abundance within each cell-type triad
niche (CTN).

[cowplot](https://wilkelab.org/cowplot/reference/cowplot-package.html)
and [ggtext](https://wilkelab.org/ggtext/reference/ggtext.html) packages
are required for the function.

## Usage

``` r
draw_CTN_celltype_prop(
  CTN_merged_celltype_prop,
  CTN_dend,
  facet_var = "clust",
  celltype_palette,
  radius_range = c(0.5, 8),
  rename_celltype = waiver(),
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

- facet_var:

  Faceting variables, usually CTN category. Default to `"clust"` from
  `CTN_dend$label_df`

- celltype_palette:

  A named vector with colors as values and cell type labels as names.
  Must contain all cell types.

- radius_range:

  A numeric vector of length two that specifies the minimum and maximum
  size of the bubbles. See
  [`ggplot2::scale_radius()`](https://ggplot2.tidyverse.org/reference/scale_size.html)
  for details.

- rename_celltype:

  A function for renaming cell type labels.

- x_expand:

  A numeric vector of length two of the multiplicative range expansion
  factors for the x axis. See
  [`ggplot2::expansion()`](https://ggplot2.tidyverse.org/reference/expansion.html)
  for details.

## Details

Dotplot visualization of the cell type abundance of cell-type triad
niches (CTN). Bubbles are colored by cell type labels, and their sizes
represent relative cell type proportions with each CTN.
