# Computing observed colocation quotient (CLQ)

Computes the observed colocation quotient (CLQ) between available cell
type pairs in each tissue section in a spatial omics dataset.

## Usage

``` r
compute_CLQ_observed(
  cell_neighbor_counts,
  dat_in,
  min.ncell = 100,
  num_cores = 1
)
```

## Arguments

- cell_neighbor_counts:

  \\N \times T\\ neighborhood cell type count matrix obtained using
  [`Build_cell_neighbor_maxdist()`](Build_cell_neighbor_maxdist.md).

- dat_in:

  Input data frame with \\N\\ rows. Columns must contain `CellID`,
  `ImageID` and `Celltype`.

- min.ncell:

  Minimum cell number threshold. Tissue sections with fewer than
  `min.ncell` cells will be excluded from the analysis.

- num_cores:

  The number of cores to use in
  [`pbmcapply::pbmclapply()`](https://rdrr.io/pkg/pbmcapply/man/pbmclapply.html),
  i.e. at most how many child processes will be run simultaneously. Can
  only be set at 1 on Windows.

## Value

A data frame of \\T\times T\times L\\ rows containing observed CLQs
between focal (`Celltype`) and neighboring cell types
(`Celltype_neighbor`) in each tissue section (`ImageID`), where \\T\\ is
the number of available cell types and \\L\\ is the number of tissue
sections.

## Details

The colocalization of a pair of cell types \\a\\ and \\b\\ in tissue
section \\l\\ can be evaluated with the observed CLQ

\$\$CLQ\_{b\rightarrow a}^{(l), obs} = \left\\\begin{array}{ll} 0, &
N_a^{(l)} \leq 5 \mbox{ or} \\ N_b^{(l)} \leq 5 \\
\frac{C\_{b\rightarrow a}^{(l), obs} / N_a^{(l)}}{N_b^{(l)} / (N^{(l)} -
1)}, & \text{otherwise} \end{array}\right.\$\$

where \\C\_{b\rightarrow a}^{(l), obs}\\ is the number of cells of cell
type \\b\\ within the neighborhood of cell type \\a\\, and
\\N_a^{(l)}\\, \\N_b^{(l)}\\, \\N^{(l)}\\ are the number of cells for
cell type \\a\\, \\b\\ and the total number of cells in tissue section
\\l\\, respectively.

## References

Bouchard, G. et al. A quantitative spatial cell-cell colocalizations
framework enabling comparisons between in vitro assembloids and
pathological specimens. Nat Commun 16, 1392 (2025).
<https://doi.org/10.1038/s41467-024-55129-6>
