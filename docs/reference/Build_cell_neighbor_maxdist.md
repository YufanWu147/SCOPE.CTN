# Constructing neighborhood cell type count and proportion matrices

Define spatial neighborhood using
[`FNN::get.knn()`](https://rdrr.io/pkg/FNN/man/get.knn.html) and compute
neighborhood cell counts and proportions by all available cell types.

## Usage

``` r
Build_cell_neighbor_maxdist(
  dat_slide,
  knn_num,
  cell_types_available,
  min.radius = 100,
  max.radius = 500,
  min.n.neighbors = 10,
  return.nn.ids = FALSE,
  prob.list = seq(0, 1, 0.25)
)
```

## Arguments

- dat_slide:

  Input data frame for a single tissue section. Must contain columns
  including cell ids (`CellID`), spatial coordinates (`X`, `Y`), and
  cell type labels (`Celltype`).

- knn_num:

  The number of nearest neighbors. The default value is 20.

- cell_types_available:

  A vector of length \\T\\ containing all available cell type labels.

- min.radius:

  Lower threshold for neighborhood distance boundary. The default value
  is 100.

- max.radius:

  Upper threshold for neighborhood distance boundary. The default value
  is 500.

- min.n.neighbors:

  Minimum number of nearest neighbors within the neighborhood distance
  boundary. Cells with fewer than `min.n.neighbors` (default: 10)
  neighbors will be excluded from downstream analysis.

- return.nn.ids:

  Whether to return cell ids of the nearest neighbors of each cell.

- prob.list:

  Probabilities between \\\[0, 1\]\\ to produce quantiles for
  neighborhood distance distribution

## Value

A list contains:

- count.df:

  An \\N \times T\\ neighborhood cell type count matrix. Rows with all
  `NaN` values are cells with fewer than fewer than `min.n.neighbors`
  neighbors within the distance boundary.

- prop.df:

  An \\N \times T\\ neighborhood cell type proportion matrix. Rows with
  all `NaN` values are cells with fewer than fewer than
  `min.n.neighbors` neighbors within the distance boundary.

- dist.quantiles:

  Quantiles (including the 99th percentile) of the of the distances to
  all nearest neighbors corresponding to the probabilities in
  `prop.list`. Can be used to examine the distribution of the distances.

- neighbor.cell.ids:

  A data frame containing cell ids of the nearest neighbors of each
  focal cell.

## Details

The neighborhood distance boundary `max.dist` is initially set as the
99th percentile of all distances from focal cells to their the nearest
neighbors. If it falls outside the predefined range (default: \\\[100,
500\]\\), it will be reset to the minimum (`min.radius`) or maximum
threshold (`max.radius`).

For cells with fewer than `min.n.neighbors` neighbors within the
boundary, their corresponding neighborhood cell type counts and
proportions will all be set to `NaN` and they will be subsequently
excluded from niche clustering.
