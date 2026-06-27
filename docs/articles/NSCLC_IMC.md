# Systematic CTN screening in NSCLC IMC data

This tutorial shows an example of systematic screening of cell-type
triad niches (CTN) in the non-small cell lung cancer (NSCLC) imaging
mass cytometry (IMC) data using SCOPE. The dataset comes from [Cords et
al. *Cancer Cell*
2024](https://www.cell.com/cancer-cell/fulltext/S1535-6108(23)00449-X)
and is composed of 2,072 tissue samples from a cohort of 1,070 patients.

First, we load the packages that will be used for the analysis and
visualization.

``` r
library(SCOPE.CTN)
library(data.table)
library(tidyverse)
library(ggtext)
library(cowplot)
dataset_name <- "NSCLC_IMC"
data_path <- paste0("vignette_data/", dataset_name, "/")   # !!! Change the directory
source(paste0(data_path, dataset_name, "_celltype_palette.R")) # cell type color palette
```

## 1. Preparing input data

The SingleCellExperiment object and clinical data for this dataset can
be downloaded from [Zenodo](https://zenodo.org/records/7961844).

- `sce_all_annotated.rds` can be found in
  `SingleCellExperiment Objects.zip`.
- `clinical_data_ROI.csv` can be found in `comp_csv_files.zip`.

``` r
library(SingleCellExperiment)
NSCLC_obj <- readRDS(paste0(data_path, "sce_all_annotated.rds"))
```

The input data frame for SCOPE should contain the following columns:

- `CellID`: cell id
- `ImageID`: sample id
- `X`: X coordinate
- `Y`: Y coordinate
- `Celltype`: cell type labels

``` r
# Extract cell metadata
NSCLC_IMC_dat_in = colData(NSCLC_obj)[,c("Tma_ac","CellID","cell_type","Center_X","Center_Y")]
NSCLC_IMC_dat_in <- as.data.frame(NSCLC_IMC_dat_in)

# Rename the columns
colnames(NSCLC_IMC_dat_in) <- c("ImageID","CellID","Celltype","X","Y")
NSCLC_IMC_dat_in$Celltype <- factor(NSCLC_IMC_dat_in$Celltype, NSCLC_IMC_celltype_levs)
```

We exclude the images with low cell counts (\< 500) from further
analysis.

``` r
ncell_by_slide <- table(NSCLC_IMC_dat_in$ImageID)
slide_list = unique(names(ncell_by_slide)[which(ncell_by_slide>500)]) # 1984 images left
NSCLC_IMC_dat_in <- NSCLC_IMC_dat_in %>%
  filter(ImageID %in% slide_list)
save(NSCLC_IMC_dat_in, file = paste0(data_path, dataset_name, "_dat_in.RData"))
```

## 2. Create neighborhood cell proportion matrix

Next, we define the spatial neighborhood using use the k-nearest
neighbors and compute the neighborhood cell type composition matrix.
Here the number of nearest neighbors $`k`$ is set at default value 20.

In addition, we enforce a distance boundary for the nearest neighbors
and exclude the cells falling outside this radius. The boundary is
initially defined as the 99th percentile of all distances from focal
cells to their nearest neighbors. If it falls outside a predefined
range, we reset it to either the minimum (`min.radius`) or maximum
threshold (`max.radius`).

``` r
cell_neighbor_results_maxdist = pbmcapply::pbmclapply(
  1:length(slide_list), function(i) {
    Build_cell_neighbor_maxdist(
      dat_slide=dat_in[dat_in$ImageID==slide_list[i], ], knn_num = 20, 
      min.radius = 100, max.radius = 120, min.n.neighbors = 10,
      cell_types_available = levels(dat_in$Celltype),
      return.nn.ids = TRUE)},
  mc.cores = 5)


## cell type proportion within the neighborhood
cell_neighbor_table_maxdist = lapply(
  cell_neighbor_results_maxdist, "[[", "prop.df") %>% 
  do.call(rbind, .)
## cell counts by cell type within the neighborhood
cell_neighbor_counts = lapply(
  cell_neighbor_results_maxdist, "[[", "count.df") %>% 
  do.call(rbind, .)
## cell ids of the nearest neighbors
cell_neighbor_ids = lapply(
  cell_neighbor_results_maxdist, "[[", "neighbor.cell.ids") %>% 
  do.call(rbind, .)
## distance to nearest neighbors
dist.quantiles = lapply(
  cell_neighbor_results_maxdist, "[[", "dist.quantiles") %>% 
  do.call(rbind, .)

save(cell_neighbor_table_maxdist, 
     file = paste0(data_path, dataset_name, "_cell_neighbor_table_maxdist.RData"))
save(cell_neighbor_counts, 
     file = paste0(data_path, dataset_name, "_cell_neighbor_counts.RData"))
save(cell_neighbor_ids, 
     file = paste0(data_path, dataset_name, "_cell_neighbor_ids.RData"))
save(dist.quantiles, 
     file = paste0(data_path, dataset_name, "_dist.quantiles.RData"))
rm(cell_neighbor_results_maxdist)
```

``` r
load(paste0(data_path, dataset_name, "_cell_neighbor_table_maxdist.RData"))
load(paste0(data_path, dataset_name, "_dist.quantiles.RData"))
```

We can examine the distribution of 99th percentile of all distance to
the nearest neighbors to select the threshold parameters min.radius and
max.radius accordingly.

``` r
dist.quantiles %>%
  as.data.frame() %>% 
  ggplot(aes(x = `99%`)) +
  geom_histogram(color = "white", fill = "skyblue", bins = 20) +
  labs(y = "99th percentile of Dist.\nto Nearest Neighbors", x = NULL) +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45))
```

![](NSCLC_IMC_files/figure-html/Examine%20distance%20distribution-1.png)

Cells that retain fewer than `min.n.neighbors` neighbors within the
distance boundary will have their neighborhood cell-type proportions set
to `NaN` and will be excluded from subsequent niche clustering. Here we
set `min.n.neighbors` to 10.

``` r
cells_to_exclude = apply(cell_neighbor_table_maxdist, 1, function(x) {any(is.na(x))}) %>%
  as.data.frame() %>%
  `colnames<-`("is_na") %>%
  rownames_to_column("CellID") %>%
  left_join(NSCLC_IMC_dat_in[, c("CellID", "X", "Y", "ImageID")], by = "CellID")

# Number of cells to exclude
table(cells_to_exclude$is_na)
#> 
#>   FALSE    TRUE 
#> 5960794     944

cells_to_exclude %>%
  group_by(ImageID) %>%
  mutate(n_NaN = sum(is_na)) %>%
  filter(n_NaN > 10) %>% ungroup %>%
  ggplot(aes(x = X, y = Y, color = is_na)) +
  geom_point(shape = 16, size= 0.1) +
  scale_color_manual(values = c("brown2", "grey90"), 
                     breaks = c(TRUE, FALSE), name = "Cells To\nExclude") +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  facet_wrap(~ImageID, scales = "free", ncol = 7) +
  theme_void() +
  coord_cartesian(expand = FALSE) +
  theme(legend.key.size = unit(2.5, "mm"),
        panel.background = element_blank(),
        aspect.ratio = 1)
```

![](NSCLC_IMC_files/figure-html/Examine%20loose%20cells-1.png)

## 3. Prefilter cell type triads using CLQ permutation test

SCOPE identifies niches characterized by cell-type triads. However, as
the number of cell types grows, the combinations of three cell types can
increase exponentially. Therefore, we utilize a colocation quotient
(CLQ) permutation test to prioritize cell-type triads (CTs) that exhibit
significant spatial proximity.

Please refer to the [Prefiltering cell type triads using CLQ permutation
test](Triad_prefilter_CLQ.md) vignette for details.

``` r
load(paste0(data_path, dataset_name, "_cell_neighbor_counts.RData"))
load(paste0(data_path, dataset_name, "_cell_neighbor_ids.RData"))

# Compute observed CLQ
CLQ_obs <- compute_CLQ_observed(
  cell_neighbor_counts = cell_neighbor_counts,
  dat_in = NSCLC_IMC_dat_in)

# Compute permutated CLQ by randomly shuffling cell type labels for 500 times
CLQ_perm <- compute_CLQ_permutated(
  cell_neighbor_ids = cell_neighbor_ids, 
  dat_in = NSCLC_IMC_dat_in,
  iter.num = 500, seed = 42, n_cores = 5)

# Compute permutation test p-value
CLQ_perm_pval <- CLQ_permutation_test(CLQ_obs = CLQ_obs, CLQ_perm = CLQ_perm)
```

Note that the CLQ permutation test can be time consuming. Here we load
the precomputed results and use the Benjamini-Hochberg (BH) adjusted
p-values to prioritize cell type triads where all constituent pairs
exhibit significant CLQs (FDR \< 0.05) in at least 40 tissue sections.

``` r
load(paste0(data_path, dataset_name, "_CLQ_perm_pval_500iter.RData"))
head(CLQ_perm_pval, 5)
#> [38;5;246m# A tibble: 5 × 5[39m
#> [38;5;246m# Groups:   Celltype, Celltype_neighbor [1][39m
#>   Celltype Celltype_neighbor ImageID  p_perm p_perm.adj
#>   [3m[38;5;246m<chr>[39m[23m    [3m[38;5;246m<chr>[39m[23m             [3m[38;5;246m<chr>[39m[23m     [3m[38;5;246m<dbl>[39m[23m      [3m[38;5;246m<dbl>[39m[23m
#> [38;5;250m1[39m _cell    _lood             175A_1    0         0     
#> [38;5;250m2[39m _cell    _lood             175A_10   1         1     
#> [38;5;250m3[39m _cell    _lood             175A_100  0.002     0.028[4m1[24m
#> [38;5;250m4[39m _cell    _lood             175A_102  0         0     
#> [38;5;250m5[39m _cell    _lood             175A_103  0.422     1
```

``` r
# Exclude scarce (< 500 cells in total) & uninformative cell types ("Other") 
ncell_by_celltype <- table(NSCLC_IMC_dat_in$Celltype)
celltype_to_exclude <- c(names(which(ncell_by_celltype < 500)), "Other")
print(paste("Cell type(s) to exclude:", celltype_to_exclude))
#> [1] "Cell type(s) to exclude: Other"

# All combinations of the rest of the cell types
celltype_combs <- combn(setdiff(NSCLC_IMC_celltype_levs,
                                celltype_to_exclude), 3) %>%
  t() %>% as.data.frame()

# Match CLQ permutation test p-values
celltype_combs_w_pval <- celltype_combs %>%
  left_join(CLQ_perm_pval, by = c("V1" = "Celltype", "V2" = "Celltype_neighbor"), 
            relationship = "many-to-many") %>%
  left_join(CLQ_perm_pval, by = c("V1" = "Celltype_neighbor", "V3" = "Celltype", "ImageID"),
            suffix = c(".11", "")) %>%
  left_join(CLQ_perm_pval, by = c("V1" = "Celltype", "V3" = "Celltype_neighbor", "ImageID"),
            suffix = c(".12", "")) %>%
  left_join(CLQ_perm_pval, by = c("V1" = "Celltype_neighbor", "V3" = "Celltype", "ImageID"),
            suffix = c(".21", "")) %>%
  left_join(CLQ_perm_pval, by = c("V2" = "Celltype", "V3" = "Celltype_neighbor", "ImageID"),
            suffix = c(".22", ""))  %>%
  left_join(CLQ_perm_pval, by = c("V2" = "Celltype_neighbor", "V3" = "Celltype", "ImageID"),
            suffix = c(".31", ".32")) 

# Filter CTs where cell type triads where all constituent pairs exhibit significant CLQs (FDR < 0.05)
combination_table_filtered = celltype_combs_w_pval %>% ungroup %>%
  filter((p_perm.adj.11 < 0.05 | p_perm.adj.12 < 0.05) & 
         (p_perm.adj.21 < 0.05 | p_perm.adj.22 < 0.05) & 
         (p_perm.adj.31 < 0.05 | p_perm.adj.32 < 0.05)) %>%
  group_by(V1, V2, V3) %>%
  summarise(n = n()) %>%
  filter(n > 40) %>%
  select(-n)

print(paste("No of remaining CTs:", nrow(combination_table_filtered)))
#> [1] "No of remaining CTs: 0"

# save(combination_table_filtered, 
#     file = paste0(data_path, dataset_name, "_combination_table_filtered.RData"))
```

``` r
celltype_to_exclude <- cell_neighbor_table_maxdist[, unique(unlist(combination_table_filtered))] %>%
  apply(2, function(x) {sum(x > 0.3, na.rm = TRUE)}) %>% sort() %>%
  as.data.frame %>%
  `colnames<-`("n") %>%
  filter(n <= 50) %>%
  rownames()

if(length(celltype_to_exclude) > 1) {
  comb_to_run <- combination_table_filtered %>%
    mutate(combs = paste(V1, V2, V3, sep = "_")) %>%
    filter(!str_detect(combs, paste(celltype_to_exclude, collapse = "|"))) %>%
    #  filter(!combs %in% str_replace_all(str_remove(list.files(paste0(path_ecotype, "_mgcv_K_2_min.prop.3")), ".RData$"), "\\.", "/")) %>%
    select(-combs)
} else {
  comb_to_run <- combination_table_filtered
}
```

An alternative way to prioritize cell type triads is to use the index.

## 4. Cell type triad niche (CTN) detection using SCOPE

Next, we apply SCOPE to the
$`\mbox{B cell/CD4}^+\mbox{ T cell/CD8}^+\mbox{ T cell}`$ triad using
default configurations for the following parameters:

- `min.prop`: The minimum neighborhood cell-type proportion threshold.
  Proportions below this value are set to 0 (default: 0.3).
- `mgcv_df`: The dimension of the basis used for thin-plate regression
  spline smoothing of the driver cell-type proportions (default: 15; see
  `mgcv::s()` for details).

Here we set the number of clusters `clusternum` at 2. When setting
`save_results = TRUE`, the clustering results of each CT will be saved
individual into `output_dir`.

``` r
# Load prefiltered CTs
load(paste0(data_path, dataset_name, "_combination_table_filtered.RData"))

celltype_to_exclude <- cell_neighbor_table_maxdist[, unique(unlist(combination_table_filtered))] %>%
  apply(2, function(x) {sum(x > 0.3, na.rm = TRUE)}) %>% sort() %>%
  as.data.frame %>%
  `colnames<-`("n") %>%
  filter(n <= 50) %>%
  rownames()


run_SCOPE(
  combination_table = combination_table_filtered,
  cell_neighbor_table = cell_neighbor_table_maxdist,
  dat_in = NSCLC_IMC_dat_in,
  save_results = TRUE,
  output_dir = paste0(data_path, "CTN_K_2"),
  clusternum = 2, min.prop = 0.3,
  linear = FALSE, mgcv_df = 15, iter.max = 100)
```

``` r
sessionInfo()
#> R version 4.4.1 (2024-06-14)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS 26.5.1
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
#>  [1] cowplot_1.2.0       ggtext_0.1.2        lubridate_1.9.5    
#>  [4] forcats_1.0.1       stringr_1.6.0       dplyr_1.2.0        
#>  [7] purrr_1.2.1         readr_2.2.0         tidyr_1.3.2        
#> [10] tibble_3.3.1        ggplot2_4.0.2.9000  tidyverse_2.0.0    
#> [13] data.table_1.18.2.1 SCOPE.CTN_0.0.1    
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6       xfun_0.57          bslib_0.10.0       htmlwidgets_1.6.4 
#>  [5] lattice_0.22-6     tzdb_0.5.0         vctrs_0.7.2        tools_4.4.1       
#>  [9] generics_0.1.4     parallel_4.4.1     pbmcapply_1.5.1    pkgconfig_2.0.3   
#> [13] Matrix_1.7-0       RColorBrewer_1.1-3 S7_0.2.1           desc_1.4.3        
#> [17] assertthat_0.2.1   lifecycle_1.0.5    compiler_4.4.1     farver_2.1.2      
#> [21] FNN_1.1.4.1        textshaping_1.0.5  htmltools_0.5.9    sass_0.4.10       
#> [25] yaml_2.3.12        gmp_0.7-5.1        pillar_1.11.1      pkgdown_2.2.0     
#> [29] jquerylib_0.1.4    cachem_1.1.0       nlme_3.1-164       tidyselect_1.2.1  
#> [33] digest_0.6.39      stringi_1.8.7      labeling_0.4.3     splines_4.4.1     
#> [37] fastmap_1.2.0      grid_4.4.1         cli_3.6.5          magrittr_2.0.4    
#> [41] utf8_1.2.6         withr_3.0.2        scales_1.4.0       timechange_0.4.0  
#> [45] rmarkdown_2.30     otel_0.2.0         ragg_1.5.1         hms_1.1.4         
#> [49] evaluate_1.0.5     knitr_1.51         irlba_2.3.7        mgcv_1.9-3        
#> [53] rlang_1.1.7        gridtext_0.1.5     Rcpp_1.1.1         glue_1.8.0        
#> [57] xml2_1.5.2         rstudioapi_0.18.0  jsonlite_2.0.0     R6_2.6.1          
#> [61] ClusterR_1.3.3     systemfonts_1.3.1  fs_1.6.7
```
