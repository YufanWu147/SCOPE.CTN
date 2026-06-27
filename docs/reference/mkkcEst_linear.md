# Robust multiple kernel K-means clustering on a multiview data

An internal function adapted from MKKC::mkkcEst for performing robust
multiple kernel K-means clustering on a multiview data.

## Usage

``` r
mkkcEst_linear(
  X_m_list,
  centers,
  iter.max = 100,
  epsilon = 1e-04,
  theta = rep(1/length(X_m_list), length(X_m_list)),
  alpha = 0.5
)
```

## Arguments

- X_m_list:

  A list of length \\P\\ (default: 4) containing \\P\\ kernel matrices.

- centers:

  The number of clusters \\K\\.

- iter.max:

  The maximum number of iterations allowed. The default is 100.

- epsilon:

  Convergence threshold. The default is \\10^{-4}\\.

- theta:

  Initial values for kernel coefficients. The default is 1/P for all
  views.

- alpha:

  Decay rate used for computing the exponential moving average of
  \\\theta\\. In interation \\\scriptl\\, \\{\theta^{(\scriptl)}}^{EMA}
  = \alpha \theta^{(\scriptl)} + (1-\alpha)
  {\theta^{(\scriptl-1)}}^{EMA}, \alpha\in\[0, 1\]\\. The default value
  is 0.5.

## Value

A list contains:

- cluster:

  A vector of integers (from `1:K`) indicating the cluster to which each
  point is allocated.

- totss:

  The total sum of squares.

- withinss:

  Matrix of within-cluster sum of squares by cluster, one row per view.

- withinsscluster:

  Vector of within-cluster sum of squares, one component per cluster.

- withinssview:

  Vector of within-cluster sum of squares, one component per view.

- tot.withinss:

  Total within-cluster sum of squares, i.e. `sum(withinsscluster)`.

- betweenssview:

  Vector of between-cluster sum of squares, one component per view.

- tot.betweenss:

  The between-cluster sum of squares, i.e. `totss-tot.withinss`.

- clustercount:

  The number of clusters `K`.

- coefficients:

  The kernel coefficients

- H:

  The continuous clustering assignment

- size:

  The number of points, one component per cluster.

- iter:

  The number of iterations.

- MbatchKm_WCSS_per_cluster:

  Within-cluster sum of squares of each cluster when recovering cluster
  labels from the final H using minibatch K-means.

- MbatchKm_centroids:

  Final cluster centroids.
