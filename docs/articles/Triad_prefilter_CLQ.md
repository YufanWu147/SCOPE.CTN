# Prefiltering cell type triads using CLQ permutation test

``` r

library(SCOPE.CTN)
library(tidyverse)
library(cowplot)
```

SCOPE identifies niches characterized by cell-type triads. However, as
the number of cell types grows, the combinations of three cell types
increase exponentially. To reduce the computational burden, we utilize a
colocation quotient (CLQ) permutation test to prioritize cell type
triads that exhibit significant spatial proximity.

``` r

data.frame(n_celltype = 3:50, n_combn = choose(3:50, 3)) %>%
  ggplot(aes(x = n_celltype, y = n_combn)) +
  geom_point(size = 1) +
  theme_cowplot() + 
  labs(x = "Number of Cell Types", y = "Number of triads")
```

![](Triad_prefilter_CLQ_files/figure-html/load%20packages-1.png)

Following [Bouchard et al.](https://doi.org/10.1038/s41467-024-55129-6),
the colocalization of a pair of cell types \\a\\ and \\b\\ in tissue
section \\l\\ can be evaluated with metric

\\CLQ\_{b\rightarrow a}^{(l)} = \begin{cases} \frac{C\_{b\rightarrow
a}^{(l)} / N_a^{(l)}}{N_b^{(l)} / (N^{(l)} - 1)} & N_a^{(l)} \> 5 \mbox{
and} \\ N_b^{(l)} \> 5 \\ 0, & \text{otherwise}\end{cases},\\

where \\CLQ\_{b\rightarrow a}^{(l)}\\ is the number of cells of cell
type \\b\\ within the neighborhood of cell type \\a\\; \\N_a^{(l)}\\,
\\N_b^{(l)}\\, and \\N^{(l)}\\ are the number of cells for cell type
\\a\\, \\b\\ and the total number of cells in tissue section \\l\\,
respectively.

Here we use the 45-sample pulmonary fibrosis (PF) Xenium dataset from
[Vannan et al. *Nat Genet*
2025](https://www.nature.com/articles/s41588-025-02080-x) to demonstrate
how to perform CLQ permutation test. This dataset contains 47 cell
types, which can form 16,215 potential triads.

For each tissue section, we first calculate an observed CLQ value
\\CLQ\_{b\rightarrow a}^{(l),obs}\\ using the true cell type labels.
Please refer to the [Detecting TLS-like CTNs in PF Xenium
data](https://github.com/YufanWu147/SCOPE.CTN/articles/PF_Xenium.md) to
see how `PF_Xenium_dat_in` and `PF_Xenium_cell_neighbor_results` were
obtained.

``` r

dataset_name <- "PF_Xenium"
data_path <- paste0("vignette_data/", dataset_name, "/")   # !!! Change the directory

# Load data
load(paste0(data_path, dataset_name, "_dat_in.RData"))
load(paste0(data_path, dataset_name, "_PF_Xenium_cell_neighbor_counts_nnids.RData"))

# Neighborhood cell counts
cell_neighbor_counts <- lapply(
  cell_neighbor_results_maxdist, "[[", "count.df") %>%
  do.call("rbind", .)
# Nearest neighbor cell ids
cell_neighbor_ids <- lapply(
  cell_neighbor_results_maxdist, "[[", "neighbor.cell.ids") %>%
  do.call("rbind", .)
```

``` r

CLQ_obs <- compute_CLQ_observed(
  cell_neighbor_counts = cell_neighbor_counts,
  dat_in = Xenium_PF_dat_in, min.ncell = 500)
```

We then randomly shuffle the cell type labels for 500 times to create a
null distribution of CLQs assuming that the cell types are distributed
randomly. \\CLQ\_{b\rightarrow a}^{(l),k},k = 1,2,…,500\\.

``` r

CLQ_perm <- compute_CLQ_permuted(
  cell_neighbor_ids = cell_neighbor_ids, 
  dat_in = Xenium_PF_dat_in, min.ncell = 500,
  iter.num = 500, seed = 42, n_cores = 5)
```

The permutation test p-value is defined as

\\p\_{b\rightarrow a}^{(l)} =
\sum\_{k=1}^{500}{I\left(CLQ\_{b\rightarrow a}^{(l), k} \>
CLQ\_{b\rightarrow a}^{(l), obs}\right)},\\

where \\I\\ is the indicator function.

``` r

CLQ_perm_pval <- CLQ_permutation_test(CLQ_obs = CLQ_obs, CLQ_perm = CLQ_perm)
```

The function returns both nominal and Benjamini-Hochberg (BH) adjusted
p-values. We prioritize cell type triads where all constituent pairs
exhibit significant CLQs (adjusted p-value \< 0.05) in at least 10
tissue sections.

``` r

## Exclude cell types with too few cells
celltype_to_exclude <- table(dat_in$Celltype)[table(dat_in$Celltype) < 500]

## All cell type combinations
celltype_combs <- combn(setdiff(levels(dat_in$Celltype), celltype_to_exclude), 3) %>%
  t() %>% as.data.frame()
celltype_combs_signif <- celltype_combs %>%
  left_join(CLQ_perm_pval, by = c("V1" = "Celltype", "V2" = "Celltype_neighbor")) %>%
  left_join(CLQ_perm_pval, by = c("V1" = "Celltype", "V3" = "Celltype_neighbor", "ImageID"),
            suffix = c(".1", "")) %>%
  left_join(CLQ_perm_pval, by = c("V2" = "Celltype", "V3" = "Celltype_neighbor", "ImageID"),
            suffix = c(".2", ".3")) 

combination_table_filtered = celltype_combs_signif %>%
  filter(p_perm.adj.1 < 0.05 & p_perm.adj.2 < 0.05 & p_perm.adj.3 < 0.05) %>%
  group_by(V1, V2, V3) %>%
  summarise(n = n()) %>%
  filter(n >= 10) %>%
  select(-n)

print("Number of cell type triads: ", nrow(combination_table_filtered))
```

``` r

sessionInfo()
#> R version 4.4.1 (2024-06-14)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS 26.5.2
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: America/Chicago
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] cowplot_1.2.0      lubridate_1.9.5    forcats_1.0.1      stringr_1.6.0     
#>  [5] dplyr_1.2.0        purrr_1.2.1        readr_2.2.0        tidyr_1.3.2       
#>  [9] tibble_3.3.1       ggplot2_4.0.2.9000 tidyverse_2.0.0    SCOPE.CTN_0.0.1   
#> 
#> loaded via a namespace (and not attached):
#>  [1] gmp_0.7-5.1         sass_0.4.10         generics_0.1.4     
#>  [4] stringi_1.8.7       lattice_0.22-6      hms_1.1.4          
#>  [7] digest_0.6.39       magrittr_2.0.4      timechange_0.4.0   
#> [10] evaluate_1.0.5      grid_4.4.1          RColorBrewer_1.1-3 
#> [13] fastmap_1.2.0       jsonlite_2.0.0      Matrix_1.7-0       
#> [16] survival_3.6-4      mgcv_1.9-3          scales_1.4.0       
#> [19] ClusterR_1.3.3      textshaping_1.0.5   jquerylib_0.1.4    
#> [22] cli_3.6.5           rlang_1.1.7         pbmcapply_1.5.1    
#> [25] splines_4.4.1       withr_3.0.2         cachem_1.1.0       
#> [28] yaml_2.3.12         FNN_1.1.4.1         otel_0.2.0         
#> [31] tools_4.4.1         parallel_4.4.1      tzdb_0.5.0         
#> [34] assertthat_0.2.1    vctrs_0.7.2         R6_2.6.1           
#> [37] lifecycle_1.0.5     fs_1.6.7            htmlwidgets_1.6.4  
#> [40] ragg_1.5.1          irlba_2.3.7         pkgconfig_2.0.3    
#> [43] desc_1.4.3          pkgdown_2.2.0       pillar_1.11.1      
#> [46] bslib_0.10.0        gtable_0.3.6        data.table_1.18.2.1
#> [49] glue_1.8.0          Rcpp_1.1.1          systemfonts_1.3.1  
#> [52] xfun_0.57           tidyselect_1.2.1    rstudioapi_0.18.0  
#> [55] knitr_1.51          farver_2.1.2        htmltools_0.5.9    
#> [58] nlme_3.1-164        labeling_0.4.3      rmarkdown_2.30     
#> [61] compiler_4.4.1      S7_0.2.1
```
