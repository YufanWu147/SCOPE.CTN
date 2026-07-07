# Visualize silhouette score

Evaluating hierarchical clustering performance for selecting optimal
cluster number.

[cowplot::cowplot-package](https://wilkelab.org/cowplot/reference/cowplot-package.html)
package is required for the function.

## Usage

``` r
draw_silhoutte_score(hclust_silhouette, num_clusters)
```

## Arguments

- hclust_silhouette:

  A data frame containing the average (`avg_si`) and median Silhouette
  score (`median_si`) under cluster number `K`. Output of the
  [`jaccard_dist_hclust()`](jaccard_dist_hclust.md) function.

- num_clusters:

  Optimal number of clusters.
