# Spatial Visualization of CTNs and Cell Types

Generates plots to visualize the spatial distribution of Celltype Triad
Niches (CTNs) and individual cell types across tissue sections.

[ggtext](https://wilkelab.org/ggtext/reference/ggtext.html) package is
required for the function.

## Usage

``` r
draw_CTN_celltype(
  Full_CTN_table,
  CTN = "Bcell_CD4T_CD8T",
  img_list,
  celltype_palette,
  CTN_palette = c(CTN = "#E69F00", Unassigned = "grey90"),
  rename_celltype = NULL,
  scales = "fixed",
  ncol = length(img_list),
  legend.ncol = NULL,
  aspect.ratio = NULL
)
```

## Arguments

- Full_CTN_table:

  A data frame containing spatial coordinates (`X`, `Y`), image IDs
  (`ImageID`), cell type labels (`Celltype`), and niche labels for all
  cell-type triads

- CTN:

  CTN label; e.g. B_CD4T_CD8T.

- img_list:

  A list of images to be visualized.

- celltype_palette:

  A named vector with colors as values and cell type labels as names.
  Must contain all cell types.

- CTN_palette:

  A named vector with colors as values and CTN labels (e.g. CTN vs.
  Unassigned) as names.

- rename_celltype:

  A function for renaming cell type labels.

- scales:

  `scales` parameter
  [`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html).

- ncol:

  Number of columns passed to
  [`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html).

- legend.ncol:

  Number of columns of legend keys, passed to
  [`ggplot2::guide_legend()`](https://ggplot2.tidyverse.org/reference/guide_legend.html).

- aspect.ratio:

  Aspect ratio of the panel, passed directly to
  [`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html).

## Value

A list of
[ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
objects:

- p_celltype:

  A spatial plot showing individual cell type distributions.

- p_CTN:

  A spatial plot showing the distribution of the specified CTN.
