# Subclustering cells labeled as "CTN" in initial analysis

Subclustering cells labeled as "CTN" in initial analysis

## Usage

``` r
run_SCOPE_subcluster(
  Full_CTN_table,
  cell_neighbor_table,
  CTN_annotation,
  annot_var = "annotation",
  subclusternum = 2,
  output_dir = getwd(),
  verbose = TRUE
)
```

## Arguments

- Full_CTN_table:

  A data frame containing spatial coordinates (`X`, `Y`), image IDs
  (`ImageID`), cell type labels (`Celltype`), and niche labels for all
  cell-type triads

- cell_neighbor_table:

  Same cell_neighbor_table used as the input for the initial clustering

- CTN_annotation:

  Data frame with CTN annotations. Must contain columns `CTN` (initial
  CTN names) and `annot_var` (consolidated CTN names)

- annot_var:

  A character string specicying the column name representing the
  consolidated CTN names in `CTN_annotation`

- subclusternum:

  Number of subclusters

- output_dir:

  Output directory used in
  [`run_SCOPE()`](https://github.com/YufanWu147/SCOPE.CTN/reference/run_SCOPE.md).

- verbose:

  Whether progress is printed during mini-batch k-means clustering.
  Passed on to
  [ClusterR::MiniBatchKmeans](https://rdrr.io/pkg/ClusterR/man/MiniBatchKmeans.html).
  Defaults to `TRUE`.

## Value

A list with the following attributes:

- Full_CTN_table_reclustered:

  A data frame containing the subcluster labels with the same structure
  as `Full_CTN_table`

- WCSS:

  Within-cluster sum of squares
