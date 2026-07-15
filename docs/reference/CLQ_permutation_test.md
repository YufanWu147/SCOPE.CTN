# Compute colocation quotient (CLQ) on permutated data

Randomly shuffles the cell type labels within each tissue section and
computes permutated CLQs between every cell type pair, creating a null
distribution that is compared with the observed CLQ to compute a
permutation test p-value.

## Usage

``` r
CLQ_permutation_test(CLQ_obs, CLQ_perm)
```

## Arguments

- CLQ_obs:

  Observed CLQs computed using
  [`compute_CLQ_observed()`](https://github.com/YufanWu147/SCOPE.CTN/reference/compute_CLQ_observed.md).

- CLQ_perm:

  Permutated CLQs computed using
  [`compute_CLQ_permuted()`](https://github.com/YufanWu147/SCOPE.CTN/reference/compute_CLQ_permuted.md).

## Value

A data frame containing the following columns:

- `Celltype`:

  Focal cell type.

- `Celltype_neighbor`:

  Neighboring cell type.

- `ImageID`:

  Tissue section ID.

- `p_perm`:

  Nominal p-value of the permutation test.

- `p_perm.adj`:

  Benjamini-Hochberg-adjusted p-value of the permutation test.

## Details

The permutation test p-value of tissue section \\l\\ is defined as

\$\$p\_{b\rightarrow a}^{(l)} =
\sum\_{k=1}^{\mbox{iter.num}}I(CLQ\_{b\rightarrow a}^{(l), k} \>
CLQ\_{b\rightarrow a}^{(l), obs}),\$\$

\\CLQ\_{b\rightarrow a}^{(l), obs}\\ is the observed CLQ computed using
[`compute_CLQ_observed()`](https://github.com/YufanWu147/SCOPE.CTN/reference/compute_CLQ_observed.md)
on the true cell type labels.

\\CLQ\_{b\rightarrow a}^{(l), k}\\ is the permuted CLQ computed using
[`compute_CLQ_permuted()`](https://github.com/YufanWu147/SCOPE.CTN/reference/compute_CLQ_permuted.md)
on the randomly shuffled cell type labels, where `iter.num` is the
number of permutations (defaults: 500).

## References

Bouchard, G. et al. A quantitative spatial cell-cell colocalizations
framework enabling comparisons between in vitro assembloids and
pathological specimens. *Nat Commun* 16, 1392 (2025).
<https://doi.org/10.1038/s41467-024-55129-6>

## See also

[compute_CLQ_observed](https://github.com/YufanWu147/SCOPE.CTN/reference/compute_CLQ_observed.md),
[compute_CLQ_permuted](https://github.com/YufanWu147/SCOPE.CTN/reference/compute_CLQ_permuted.md)
