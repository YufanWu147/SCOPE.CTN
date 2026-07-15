# Rename consolidated CTNs based on dominant cell types

Renaming the consolidated Cell-type Triads Niches (CTNs) using the most
frequent core cell types within the original triads. If a tie occurs,
the function resolves it by identifying the cell type with the highest
overall proportion.

## Usage

``` r
CTN_annotate(CTN_cluster_label, CTN_merged_celltype_prop, celltype_levs)
```

## Arguments

- CTN_cluster_label:

  A data frame containing hierarchical clustering results for merging
  the highly overlapped CTNs. It must include the following two columns:
  `CTN` (the CTN names) and `cluster` (the corresponding cluster label).

- CTN_merged_celltype_prop:

  A data frame representing the cell type proportion for each CTN with
  the following columns: `CTN` (the CTN names), `Celltype` (cell type
  labels), and `prop` (cell type proportion of the CTN).

- celltype_levs:

  A character vector defining the levels or ordering of the cell types.

## Value

A list contains:

- annotation:

  A data frame with the original CTN names (`CTN`), the corresponding
  cluster label (`cluster`), and the new name for the consolidated CTN
  (`annotation`) after merging.

- core_celltypes_summary:

  A data frame summarizing clusters with at least 2 CTNs. It includes
  the number of appearances of each core cell type in the original CTN
  names (`n_appearance`) and their relative proportion within the
  consolidated CTNs (`prop`).

## Details

If the top two most frequent cell types only appear twice in the
original CTN names, the function flags the consolidated CTN with an
asterisk (\*). Duplicated new CTN names are flagged with an plus sign
("+"). We suggest double-checking the flagged CTNs and manually renaming
them if necessary.
