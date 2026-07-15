# Compute colocation quotient (CLQ)

Compute colocation quotient (CLQ) between available cell type pairs.

## Usage

``` r
compute_CLQ(interact_matrix, dat_in)
```

## Arguments

- interact_matrix:

  Input matrix \\C^{(l)} \in T \times T\\, where \\T\\ is the number of
  available cell types. The element in the \\a\\-th row and \\b\\-th
  column \\C\_{b\rightarrow a}^{(l)}\\ represents the number of cells of
  cell type \\b\\ within the neighborhood of cell type \\a\\ in tissue
  section \\l\\.

- dat_in:

  An input data frame with \\N^{(l)}\\ rows, where \\N^{(l)}\\
  represents the number of cells in tissue section \\l\\. The columns
  must include `CellID` and `Celltype`.

## Value

A \\T\times T\\ matrix containing CLQs of tissue section \\l\\ between
every cell type pair, where \\T\\ is the number of available cell types.
The rows of the matrix are the focal cell types and the columns are the
neighboring cell types.

## Details

The colocalization of a pair of cell types \\a\\ and \\b\\ in tissue
section \\l\\ can be evaluated with the CLQ metric

\$\$CLQ\_{b\rightarrow a}^{(l)} = \left\\\begin{array}{ll}
\frac{C\_{b\rightarrow a}^{(l)} / N_a^{(l)}}{N_b^{(l)} / (N^{(l)} - 1)}
& N_a^{(l)} \> 5 \mbox{ or} \\ N_b^{(l)} \> 5 \\ 0, & \text{otherwise}
\end{array}\right.\$\$

where \\C\_{b\rightarrow a}^{(l)}\\ is the number of cells of cell type
\\b\\ within the neighborhood of cell type \\a\\, and \\N_a^{(l)}\\,
\\N_b^{(l)}\\, \\N^{(l)}\\ are the number of cells for cell type \\a\\,
\\b\\ and the total number of cells in tissue section \\l\\,
respectively.

## References

Bouchard, G. et al. A quantitative spatial cell-cell colocalizations
framework enabling comparisons between in vitro assembloids and
pathological specimens. Nat Commun 16, 1392 (2025).
<https://doi.org/10.1038/s41467-024-55129-6>
