# Computing colocation quotient (CLQ) on permutated data

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
  [`compute_CLQ_observed()`](compute_CLQ_observed.md).

- CLQ_perm:

  Permutated CLQs computed using
  [`compute_CLQ_permutated()`](compute_CLQ_permutated.md).

## Value

A data frame containing the nominal and Benjamini-Hochberg-adjusted CLQ
permutation test p-values.

## Details

The permutation test p-value of tissue section \\l\\ is defined as

\$\$p\_{b\rightarrow a}^{(l)} =
\sum\_{k=1}^{\mbox{iter.num}}I(CLQ\_{b\rightarrow a}^{(l), k} \>
CLQ\_{b\rightarrow a}^{(l), obs}),\$\$

where \\CLQ\_{b\rightarrow a}^{(l), k}\\ is the permutated CLQ computed
using [`compute_CLQ_permutated()`](compute_CLQ_permutated.md) on the
randomly shuffled cell type labels, and \\CLQ\_{b\rightarrow a}^{(l),
obs}\\ is the observed CLQ computed using
[`compute_CLQ_observed()`](compute_CLQ_observed.md) on the true cell
type labels.

## References

Bouchard, G. et al. A quantitative spatial cell-cell colocalizations
framework enabling comparisons between in vitro assembloids and
pathological specimens. Nat Commun 16, 1392 (2025).
<https://doi.org/10.1038/s41467-024-55129-6>
