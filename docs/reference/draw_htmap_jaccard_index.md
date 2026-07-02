# Jaccard similarity index heatmap

`ComplexHeatmap` and `grid` packages are required for the function.

## Usage

``` r
draw_htmap_jaccard_index(jaccard_mat, cl, num_clusters)
```

## Arguments

- jaccard_mat:

  A numeric matrix containing the Jaccard similarity indices between all
  pairs of CTNs. Output of the
  [`jaccard_dist_hclust()`](jaccard_dist_hclust.md) function.

- cl:

  An object of class "hclust" which describes the tree produced by the
  clustering process. Output of the
  [`jaccard_dist_hclust()`](jaccard_dist_hclust.md) function.

- num_clusters:

  An integer specifying the target number of clusters to slice or
  highlight on the heatmap.

## Value

A
[ComplexHeatmap::Heatmap](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap-class.html)
object.

## Details

Heatmap visualization of the Jaccard similarity indices between all
pairs of cell-type triad Niches (CTNs) using an ordered heatmap based on
hierarchical clustering.
