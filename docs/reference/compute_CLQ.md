# Computing colocation quotient (CLQ)

Computes colocation quotient (CLQ) between available cell type pairs.

## Usage

``` r
compute_CLQ(interact_matrix, dat_in)
```

## Arguments

- interact_matrix:

  Input

- dat_in:

  Input data frame with \\N\\ rows. Columns must contain `CellID` and
  `Celltype`.

## Value

A \\T\times T\\ matrix containing CLQs between every cell type pair,
where \\T\\ is the number of available cell types. The rows of the
matrix are the focal cell types and the columns are the neighboring cell
types.

## Details

The colocalization of a pair of cell types \\a\\ and \\b\\ in tissue
section \\l\\ can be evaluated with the CLQ metric

\$\$CLQ\_{b\rightarrow a} = \left\\\begin{array}{ll}
\frac{C\_{b\rightarrow a} / N_a}{N_b / (N - 1)} & N_a \> 5 \mbox{and} \\
N_b \> 5 \\ 0, & \text{otherwise} \end{array}\right.\$\$

where \\C\_{b\rightarrow a}\\ is the number of cells of cell type \\b\\
within the neighborhood of cell type \\a\\, and \\N_a\\, \\N_b\\, \\N\\
are the number of cells for cell type \\a\\, \\b\\ and the total number
of cells in tissue section \\l\\, respectively.

## References

Bouchard, G. et al. A quantitative spatial cell-cell colocalizations
framework enabling comparisons between in vitro assembloids and
pathological specimens. Nat Commun 16, 1392 (2025).
<https://doi.org/10.1038/s41467-024-55129-6>
