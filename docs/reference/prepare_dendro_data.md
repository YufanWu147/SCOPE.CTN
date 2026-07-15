# Prepare input data CTN dendrogram

Internal function for preparing input data to
[`draw_CTN_dendro()`](https://github.com/YufanWu147/SCOPE.CTN/reference/draw_CTN_dendro.md)
function

[ggdendro](https://andrie.github.io/ggdendro/reference/ggdendro-package.html)
package is required for the function.

## Usage

``` r
prepare_dendro_data(cl, h, k = NULL)
```

## Arguments

- cl:

  An object of class "hclust" which describes the tree produced by the
  clustering process. Output of the
  [`jaccard_dist_hclust()`](https://github.com/YufanWu147/SCOPE.CTN/reference/jaccard_dist_hclust.md)
  function.

- h:

  Numeric scalar or vector with heights where the tree should be cut.
  [`stats::cutree()`](https://rdrr.io/r/stats/cutree.html) for details.
  At least one of `k` or `h` must be specified, `k` overrides `h` if
  both are given.

- k:

  An integer scalar or vector with the desired number of groups. See
  [`stats::cutree()`](https://rdrr.io/r/stats/cutree.html) for details.

## Value

A list with components (See
[`ggdendro::dendro_data()`](https://andrie.github.io/ggdendro/reference/dendro_data.html)
for details):

- segments:

  Line segment data

- labels:

  Label data

## References

Atrebas. Dendrograms in R, a lightweight approach.
<https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/>
