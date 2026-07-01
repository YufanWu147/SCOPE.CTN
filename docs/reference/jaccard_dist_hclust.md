# Perform hierarchical clustering on Jaccard distance

Computes the Jaccard distance matrix between all pairs of cell-type
triad niches (CTNs) and performs hierarchical clustering analysis on the
result.

## Usage

``` r
jaccard_dist_hclust(Full_CTN_table, CTN_ls, num_cores = 1)
```

## Arguments

- Full_CTN_table:

  A data frame containing clustering labels for all cell-type triads

- CTN_ls:

  A vector of strings containing all CTNs

- num_cores:

  The number of cores to use in
  [`pbmcapply::pbmclapply()`](https://rdrr.io/pkg/pbmcapply/man/pbmclapply.html),
  i.e. at most how many child processes will be run simultaneously. Can
  only be set at 1 on Windows.

## Value

A list contains:

- jaccard_mat_mgcv:

  A matrix containing the Jaccard similarity indices (1-Jaccard
  distance) between every pair of CTNs

- cl:

  An object of class "hclust" which describes the tree produced by the
  clustering process. See
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html) for details.

- hclust_silhouette:

  A dataframe containing the average and median silhouette scores for
  evaluating cluster quality. See
  [`cluster::silhouette()`](https://rdrr.io/pkg/cluster/man/silhouette.html)
  for details.
