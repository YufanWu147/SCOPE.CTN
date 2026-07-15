# Visualize CTN cell type composition similarity using UMAP

Conduct dimensionality reduction on cell-type triad niche (CTN) cell
type composition using the Uniform Manifold Approximation and Projection
(UMAP) method (McInnes et al., 2018).

The first two UMAP embeddings are visualized using a bubble plot, where
the sizes reflect the total cell count within each CTN, and the colors
denote their CTN category.

[ggrepel](https://ggrepel.slowkow.com/reference/ggrepel.html),
[cowplot](https://wilkelab.org/cowplot/reference/cowplot-package.html),
and [uwot](https://jlmelville.github.io/uwot/reference/umap.html)
packages are required for the function.

## Usage

``` r
draw_CTN_umap(
  Full_CTN_table,
  CTN_merged_celltype_prop,
  CTN_ls,
  CTN_annotation = NULL,
  CTN_group_var = "cluster",
  CTN_group_palette = NULL,
  radius_range = c(1, 8),
  umap_n_neighbors = 10,
  umap_seed = 42,
  ...
)
```

## Arguments

- Full_CTN_table:

  A data frame containing clustering labels for all cell-type triads

- CTN_merged_celltype_prop:

  A data frame representing the cell type proportion for each CTN with
  the following columns: `CTN` (the CTN names), `Celltype` (cell type
  labels), and `prop` (cell type proportion of the CTN).

- CTN_ls:

  A character vector containing the names of the CTNs.

- CTN_annotation:

  A data frame with the original CTN names (`CTN`), the corresponding
  cluster label (`cluster`), and the new name for the consolidated CTN
  (specified in `annot_var`) after merging.

- CTN_group_var:

  A character string specifying the column in `CTN_annotation` on CTN
  grouping. Will be used to color the dots.

- CTN_group_palette:

  A named vector with colors as values and CTN categories as names.

- radius_range:

  A numeric vector of length two that specifies the minimum and maximum
  size of the bubbles. See
  [`ggplot2::scale_radius()`](https://ggplot2.tidyverse.org/reference/scale_size.html)
  for details.

- umap_n_neighbors:

  `n_neighbors` parameter
  [`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html).
  The size of local neighborhood (in terms of number of neighboring
  sample points) used for manifold approximation.

- umap_seed:

  `seed` parameter
  [`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html).
  Integer seed to use to initialize the random number generator state.

- ...:

  Other parameters to pass on to
  [`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html)

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.
