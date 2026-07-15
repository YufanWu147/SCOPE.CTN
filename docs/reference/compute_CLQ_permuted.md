# Compute colocation quotient (CLQ) on permuted data

Randomly shuffle cell type labels within each tissue section and compute
the permuted CLQs between every cell type pair.

This creates a null distribution for comparison with the observed CLQ
([`compute_CLQ_observed()`](https://github.com/YufanWu147/SCOPE.CTN/reference/compute_CLQ_observed.md))
to compute a permutation p-value.

## Usage

``` r
compute_CLQ_permuted(
  cell_neighbor_ids,
  dat_in,
  min.ncell = 100,
  num_cores = 5,
  iter.num = 500,
  seed = 42
)
```

## Arguments

- cell_neighbor_ids:

  Cell IDs of the nearest neighbors obtained using
  [`Build_cell_neighbor_maxdist()`](https://github.com/YufanWu147/SCOPE.CTN/reference/Build_cell_neighbor_maxdist.md).

- dat_in:

  An input data frame with \\N\\ rows. Columns must contain `CellID`,
  `ImageID` and `Celltype`.

- min.ncell:

  An integer specifying the minimum cell count threshold. Tissue
  sections with fewer than `min.ncell` cells will be excluded from the
  analysis.

- num_cores:

  Number of cores to use in
  [`pbmcapply::pbmclapply()`](https://rdrr.io/pkg/pbmcapply/man/pbmclapply.html),
  which determines how many permutations are processed simultaneously.
  Windows environments only support `1` core.

- iter.num:

  Number of permutations (default: 500).

- seed:

  Random seed passed to
  [`base::set.seed()`](https://rdrr.io/r/base/Random.html) to ensure
  reproducibility when shuffling the cell type labels.

## Value

A data frame of \\T\times T\times L\times\\ `iter.num` rows containing
the permutated CLQs between focal (`Celltype`) and neighboring cell
types (`Celltype_neighbor`) in each tissue section (`ImageID`), where
\\T\\ is the number of available cell types and \\L\\ is the number of
tissue sections.

## Details

In each tissue section, the cell type labels are randomly shuffled for
`iter.num` (default: 500) times to generate permutated CLQ values

\$\$CLQ\_{b\rightarrow a}^{(l), k} = \left\\\begin{array}{ll} 0, &
N_a^{(l)} \leq 5 \mbox{ or} \\ N_b^{(l)} \leq 5 \\
\frac{C\_{b\rightarrow a}^{(l), k} / N_a^{(l)}}{N_b^{(l)} / (N^{(l)} -
1)}, & \text{otherwise} \end{array}\right.\$\$

where \\C\_{b\rightarrow a}^{(l), k}\\ is the number of cells of cell
type \\b\\ within the neighborhood of cell type \\a\\ in iteration \\k,
k=\\1,2,...,`iter.num`, and \\N_a^{(l)}\\, \\N_b^{(l)}\\, \\N^{(l)}\\
are the number of cells for cell type \\a\\, \\b\\ and the total number
of cells in tissue section \\l\\, respectively.

## References

Bouchard, G. et al. A quantitative spatial cell-cell colocalizations
framework enabling comparisons between in vitro assembloids and
pathological specimens. *Nat Commun* 16, 1392 (2025).
<https://doi.org/10.1038/s41467-024-55129-6>

## See also

[Build_cell_neighbor_maxdist](https://github.com/YufanWu147/SCOPE.CTN/reference/Build_cell_neighbor_maxdist.md),
[compute_CLQ_observed](https://github.com/YufanWu147/SCOPE.CTN/reference/compute_CLQ_observed.md),
[CLQ_permutation_test](https://github.com/YufanWu147/SCOPE.CTN/reference/CLQ_permutation_test.md)
