# Multi-view cell type triad niche screening

Cell type triad niche (CTN) detection using SCOPE.

## Usage

``` r
run_SCOPE(
  combination_table,
  cell_neighbor_table,
  dat_in,
  clusternum = 2,
  min.prop = 0.3,
  linear = FALSE,
  mgcv_df = 15,
  iter.max = 100,
  epsilon = 1e-04,
  alpha = 0.5,
  num_cores = 1,
  rename_CTN = TRUE,
  save_results = TRUE,
  output_dir = getwd(),
  suffix = ""
)
```

## Arguments

- combination_table:

  Data frame containing all cell type combinations for CTN screening.
  Each row is a combination.

- cell_neighbor_table:

  \\N\times T\\ neighborhood cell type proportion table. Column names
  must contain all cell type labels in `combination_table`.

- dat_in:

  Input data frame with \\N\\ rows. Columns must contain `CellID` and
  `Celltype`.

- clusternum:

  The number of clusters \\K\\. The default value is 2.

- min.prop:

  Minimum neighborhood cell type proportion threshold. Proportions below
  this value will all be set to 0. The default value is 0.3.

- linear:

  Whether to run linear version of SCOPE (i.e. not applying thin plate
  regression spline smoothing). Default is `FALSE`.

- mgcv_df:

  The dimension of the basis used to represent the thin plate regression
  spline smooth term. The default is 15. See
  [mgcv::s](https://rdrr.io/pkg/mgcv/man/s.html) for details. Disabled
  when `linear=TRUE`.

- iter.max:

  The maximum number of iterations allowed. The default is 100.

- epsilon:

  Convergence threshold. The default is \\10^{-4}\\.

- alpha:

  Decay rate \\\alpha\\ used for computing the exponential moving
  average of weight vector \\\theta\\. In interation \\l\\,
  \\{\theta^{(l)}}^{EMA} = \alpha \theta^{(l)} + (1-\alpha)
  {\theta^{(l-1)}}^{EMA}, \alpha\in\[0, 1\]\\. The default value is 0.5.

- num_cores:

  The number of cores to use in
  [pbmcapply::pbmclapply](https://rdrr.io/pkg/pbmcapply/man/pbmclapply.html),
  i.e. at most how many child processes will be run simultaneously. Can
  only be set at 1 on Windows.

- rename_CTN:

  Whether to rename niches based on the proportion of core cell types.
  The cluster with the highest proportion will be named as `CTN` when
  `clusternum = 2` and `CTN1` otherwise. The cluster with the lowest
  core cell type proportion will be renamed as `Unassigned`. Other
  clusters will be named as `CTN2, ...` accordingly when
  `clusternum > 2`.

- save_results:

  Whether to save the niche screening results to the output directory.
  If `TRUE` (default) then the results for each CTN will be saved into a
  separate file. If `FALSE` then will all CTN screening results will be
  combined into a single data frame and returned as function output.

- output_dir:

  Output directory. Will automatically create a new directory when the
  target directory doesn't exist.

- suffix:

  File name suffix.

## Value

- res:

  CTN screening results. Only returns when `save_results` is `FALSE`.
