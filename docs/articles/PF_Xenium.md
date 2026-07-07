# Detecting TLS-like CTNs in PF Xenium data

This tutorial shows an example of detecting the tertiary lymphoid
structure (TLS)-like niches characterized by the B cell/CD4\\^+\\ T
cell/CD8\\^+\\ T cell triad in the pulmonary fibrosis (PF) Xenium data
using SCOPE. The dataset comes from [Vannan et al. *Nat Genet*
2025](https://www.nature.com/articles/s41588-025-02080-x) and is
composed of 35 PF and 10 unaffected samples.

First, we load the packages that will be used for the analysis and
visualization.

``` r

library(SCOPE.CTN)
library(data.table)
library(tidyverse)
library(ggtext)
library(ggrepel)
library(cowplot)
library(RColorBrewer)

dataset_name <- "PF_Xenium"
data_path <- paste0("vignette_data/", dataset_name, "/")   # !!! Change the directory
source(paste0(data_path, dataset_name, "_celltype_palette.R")) # cell type color palette
```

## 1. Preparing input data

The `rds` file containing a Seurat object of this dataset can be
downloaded from Gene Expression Omnibus (GEO) under accesssion number
[GSE250346](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250346).
We downloaded the
`GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds` file
and saved it under `<data_path>/`.

``` r

library(Seurat)

PF_Xenium_obj <- readRDS(paste0(data_path, "GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds"))
clinical_data <- distinct(PF_Xenium_obj@meta.data, sample, sample_affect)

# save(clinical_data, file = paste0(data_path, dataset_name, "_clinical_data.RData")))
```

We will also use the `cells_partitioned_by_annotation.csv` file in
`GSE250346_HE_annotations.tar.gz` which contains annotations of PF
histopathological features (including TLSs). The samples in this file
need to renamed following [the code from the original
paper](https://github.com/Banovich-Lab/Spatial_PF/blob/828db997a122531e6b907b792f58ae1bbd997662/figures.R#L260):

``` r

load(paste0(data_path, dataset_name, "_clinical_data.RData"))
PF_Xenium_annot <- read.csv(paste0(data_path, "cells_partitioned_by_annotation.csv")) %>%
  ## keeping only the relevant columns to reduce data size
  select(CellID = cell_id, ImageID = sample, Annotation_Type, annotation_type_instance) %>%
  mutate(patient = str_extract(ImageID, "(THD|TILD|VUHD|VUILD)[0-9]+")) %>%
  left_join(clinical_data[, c("patient", "sample_affect")], by = "patient") %>%
  ## Change sample name based on the code from the original paper (
  ## Initial rename
  mutate(ImageID_new = case_when(
    sample_affect == "Less Affected" ~ paste0(patient, "LA"),
    sample_affect == "More Affected"  ~ paste0(patient, "MA"),
    TRUE ~ ImageID)) %>%
  ## Fix additional names
  mutate(ImageID_new = case_when(
    ImageID == "TILD117MF" ~ "TILD117MA1",
    ImageID == "TILD117MFB" ~ "TILD117MA2",
    ImageID == "VUILD104LF" ~ "VUILD104MA1",
    ImageID == "VUILD104MF" ~ "VUILD104MA2",
    ImageID == "VUILD105LF" ~ "VUILD105MA1",
    ImageID == "VUILD105MF" ~ "VUILD105MA2",
    ImageID == "VUILD48LF" ~ "VUILD48LA1",
    ImageID == "VUILD48MF" ~ "VUILD48LA2",
    str_detect(ImageID_new, 'LF$') ~ str_replace(ImageID_new, "LF$", "LA"),
    str_detect(ImageID_new, 'MF$') ~ str_replace(ImageID_new, "MF$", "MA"),
    TRUE ~ ImageID_new)) %>%
  mutate(CellID = paste(ImageID_new, CellID, sep = "_")) %>%
  select(-ImageID) %>%
  rename(ImageID = ImageID_new)

# save(PF_Xenium_annot, file = paste0(data_path, dataset_name, "_annotations.RData"))
```

The input data frame for SCOPE should contain the following columns:

- `CellID`: cell id
- `ImageID`: sample id
- `X`: X coordinate
- `Y`: Y coordinate
- `Celltype`: cell type labels

``` r

PF_Xenium_dat_in <- PF_Xenium_obj@meta.data %>%
  select(CellID = full_cell_id, ImageID = sample, 
         X = x_centroid, Y = y_centroid, Celltype = final_CT)
## save(PF_Xenium_dat_in, file = paste0(data_path, dataset_name, "_dat_in.RData"))
```

Below are the four samples containing TLSs:

``` r

TLS_img_list <- unique(subset(PF_Xenium_annot, Annotation_Type == "TLS")$ImageID)

## Images containing TLSs, colored by cell type
PF_Xenium_dat_in %>%
  filter(ImageID %in% c(TLS_img_list)) %>%
  ggplot(aes(x = X, y = Y, color = Celltype)) +
  geom_point(shape = 16, size= 0.1) +
  scale_color_manual(values = PF_Xenium_celltype_palette, 
                     breaks = names(PF_Xenium_celltype_palette)) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  facet_wrap(~ImageID, scales = "free", nrow = 1) +
  theme_void(base_size = 14) +
  coord_cartesian(expand = FALSE) +
  theme(legend.key.size = unit(2.5, "mm"),
        panel.background = element_blank(),
        aspect.ratio = 1,
        legend.key.spacing.y = unit(0.5, "mm"),
        legend.position = "bottom")

## Images containing TLSs, colored by annotations
PF_Xenium_dat_in %>%
  left_join(PF_Xenium_annot[, c("CellID", "Annotation_Type")], by = "CellID") %>%
  filter(ImageID %in% c(TLS_img_list)) %>%
  ggplot(aes(x = X, y = Y, color = Annotation_Type)) +
  geom_point(shape = 16, size= 0.1) +
  scale_color_manual(breaks = names(annotation_palette), na.value = "grey90",
                     values = annotation_palette, name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  facet_wrap(~ImageID, scales = "free", nrow = 1) +
  theme_void(base_size = 14) +
  coord_cartesian(expand = FALSE) +
  theme(legend.key.size = unit(2.5, "mm"),
        panel.background = element_blank(),
        aspect.ratio = 1,
        legend.key.spacing.y = unit(0.5, "mm"),
        legend.position = "bottom")
```

![](PF_Xenium_files/figure-html/Visualize%20celltype-1.png)![](PF_Xenium_files/figure-html/Visualize%20celltype-2.png)

## 2. Create neighborhood cell proportion matrix

Next, we use the k-nearest neighbors algorithm to define the spatial
neighborhood and compute the neighborhood cell type composition matrix.
The default number of nearest neighbors \\k\\ is 20, however we can
increase \\k\\ for larger tissues. Since this dataset contains 3-5mm
lung tissue cores, here we set \\k=50\\.

In addition, we enforce a distance boundary for the nearest neighbors
and exclude the cells falling outside this radius. The boundary is
initially defined as the 99th percentile of all distances from focal
cells to their nearest neighbors. If it falls outside a predefined
range, we reset it to either the minimum (`min.radius`) or maximum
threshold (`max.radius`).

``` r

slide_list <- unique(PF_Xenium_dat_in$ImageID)

cell_neighbor_results_maxdist <- pbmcapply::pbmclapply(
  slide_list, function(s) {
    Build_cell_neighbor_maxdist(
      dat_slide=subset(PF_Xenium_dat_in, ImageID == s),
      knn_num = 50, cell_types_available = unique(PF_Xenium_dat_in$Celltype),
      min.radius = 100,       # lower limit for the boundary
      max.radius = 500,       # upper limit for the boundary
      min.n.neighbors = 10,   # minimum no. of neighbors,
      return.nn.ids = TRUE,   # Whether to return nearest neighbor Cell IDs
      prob.list = seq(0, 1, 0.1))},
  mc.cores = 4) # set mc.cores=1 on Windows machine

# Extract neighborhood cell count and ids and save it into a separate file
cell_neighbor_counts_nnids <-  lapply(
  cell_neighbor_results_maxdist, function(x) {x[c("count.df", "neighbor.cell.ids")]})
# save(cell_neighbor_counts_nnids,
#      file = paste0(data_path, dataset_name, "_cell_neighbor_counts_nnids.RData"))
# Extract neighborhood cell type proportion and 
# quantiles of distance to nearest neighbors
cell_neighbor_props_dist.quantiles <-  lapply(
  cell_neighbor_results_maxdist, function(x) {x[c("prop.df", "dist.quantiles")]})
# save(cell_neighbor_props_dist.quantiles, 
#      file = paste0(data_path, dataset_name, "_cell_neighbor_props_dist.quantiles.RData"))

rm(cell_neighbor_results_maxdist)
```

We can examine the distribution of the 99th percentile of all distance
to the nearest neighbors to select the threshold parameters min.radius
and max.radius accordingly. Here we set `max.radius` and `min.radius` to
their default values (500 and 100 respectively).

``` r

## 99th percentile of the distance to nearest neighbors
dist_99th <- sapply(cell_neighbor_props_dist.quantiles, "[[", "dist.quantiles") %>%
  `colnames<-`(slide_list) %>%
  as.data.frame() %>% rownames_to_column("perc") %>%
  mutate(perc = as.numeric(str_remove(perc, "%"))) %>%
  pivot_longer(2:ncol(.), names_to = "ImageID", values_to = "value") %>%
  filter(perc == 99)

## Distribution of 99th percentile
dist_99th %>%
  ggplot(aes(x = value)) +
  geom_histogram(color = "white", fill = "skyblue", bins = 20) +
  geom_boxplot() +
  geom_text(data = subset(dist_99th, value > 150), 
            aes(label = ImageID, y = 2), 
            hjust = 0, size = 3.5, angle = 90) +
  labs(y = "99th percentile of Dist.\nto Nearest Neighbors", x = NULL) +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45))
```

![](PF_Xenium_files/figure-html/Examine%20distance%20distribution-1.png)

We can see that the two outlier samples represent sparsely populated
regions of the tissue.

``` r

img_outliers <- arrange(subset(dist_99th, value > 150), value)$ImageID
PF_Xenium_dat_in %>%
  filter(ImageID %in% img_outliers) %>%
  ggplot(aes(x = X, y = Y, color = Celltype)) +
  geom_point(shape = 16, size= 0.1) +
  scale_color_manual(values = PF_Xenium_celltype_palette, 
                     breaks = names(PF_Xenium_celltype_palette),
                     labels = rename_celltype, name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 5)) +
  facet_wrap(~ factor(ImageID, img_outliers), scales = "free") +
  theme_void(base_size = 14) +
  coord_cartesian(expand = FALSE) +
  theme(legend.key.size = unit(2.5, "mm"),
        panel.background = element_blank(),
        aspect.ratio = 1,
        legend.key.spacing.y = unit(0.2, "mm"),
        legend.position = "bottom")
```

![](PF_Xenium_files/figure-html/Sparse%20samples-1.png)

``` r

# Combine neighborhood proportion table
cell_neighbor_table_maxdist <- lapply(
  cell_neighbor_props_dist.quantiles, "[[", "prop.df") %>%
  do.call("rbind", .)
```

For cells retaining fewer than `min.n.neighbors` neighbors within the
distance boundary, the neighborhood cell-type proportions will all be
set to `NaN` and they will be excluded from subsequent niche clustering.
Here we set `min.n.neighbors` to 10.

``` r

cells_to_exclude = apply(cell_neighbor_table_maxdist, 1, function(x) {any(is.na(x))}) %>%
  as.data.frame() %>%
  `colnames<-`("is_na") %>%
  rownames_to_column("CellID") %>%
  left_join(PF_Xenium_dat_in[, c("CellID", "X", "Y", "ImageID")], by = "CellID")

# Number of cells to exclude
table(cells_to_exclude$is_na)
#> 
#>   FALSE    TRUE 
#> 1628906    1413

# Cells to be removed from the two sparse samples
cells_to_exclude %>%
  filter(ImageID %in% img_outliers) %>%
  ggplot(aes(x = X, y = Y, color = is_na, size = is_na)) +
  geom_point(shape = 16) +
  scale_color_manual(values = c("brown2", "grey90"), 
                     breaks = c(TRUE, FALSE), name = "Cells To\nExclude") +
  scale_size_manual(values = c(0.1, 1)) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1), size = "none") +
  facet_wrap(~factor(ImageID, img_outliers), scales = "free") +
  theme_void() +
  coord_cartesian(expand = FALSE) +
  theme(legend.key.size = unit(2.5, "mm"),
        aspect.ratio = 1,
        panel.background = element_blank())
```

![](PF_Xenium_files/figure-html/Examine%20loose%20cells-1.png)

``` r


# Cells to be removed from the samples containing TLSs
cells_to_exclude %>%
  filter(ImageID %in% TLS_img_list) %>%
  ggplot(aes(x = X, y = Y, color = is_na, size = is_na)) +
  geom_point(shape = 16) +
  scale_color_manual(values = c("brown2", "grey90"), 
                     breaks = c(TRUE, FALSE), name = "Cells To\nExclude") +
  scale_size_manual(values = c(0.1, 0.5)) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1), size = "none") +
  facet_wrap(~ImageID, scales = "free", nrow = 1) +
  theme_void() +
  coord_cartesian(expand = FALSE) +
  theme(legend.key.size = unit(2.5, "mm"),
        aspect.ratio = 1,
        panel.background = element_blank())
```

![](PF_Xenium_files/figure-html/Examine%20loose%20cells-2.png)

## 3. Cell type triad niche (CTN) detection using SCOPE

Next, we apply SCOPE to the B cell/CD4\\^+\\ T cell/CD8\\^+\\ T cell
triad using default configurations for the following parameters:

- `min.prop`: The minimum neighborhood cell-type proportion threshold.
  Proportions below this value are set to 0 (default: 0.3).
- `mgcv_df`: The dimension of the basis used for thin-plate regression
  spline smoothing of the core cell-type proportions (default: 15; see
  [`mgcv::s()`](https://rdrr.io/pkg/mgcv/man/s.html) for details).

We set the number of clusters `clusternum` at 3 to separate the B cell
and T cell zones.

``` r

core_celltypes <- c("B cells", "CD4+ T-cells", "CD8+ T-cells")

CTN_df <- run_SCOPE(
  combination_table = as.data.frame(matrix(core_celltypes, nrow = 1)),
  cell_neighbor_table = cell_neighbor_table_maxdist,
  dat_in = PF_Xenium_dat_in,
  save_results = FALSE, 
  clusternum = 3, min.prop = 0.3,
  linear = FALSE, mgcv_df = 15, iter.max = 100)
# save(CTN_df, file = paste0(data_path, dataset_name, "_B_CD4T_CD8T.RData"))
```

Here we can see that the B cell/CD4\\^+\\ T cell/CD8\\^+\\ T cell CTNs
closely aligns with the lymphocyte aggregates in the CRC samples, where
CTN1 resembles the B cell zone and CTN2 resembles the surrounding T cell
zone.

``` r

# Number of cells per CTN
table(CTN_df$label)
#> 
#>       CTN1       CTN2 Unassigned 
#>      42777      44617    1541512
```

``` r

# CTN label with colored dots separating the cell types
CTN_label <- color_CTN_names(
  paste(core_celltypes, collapse = "_"), PF_Xenium_celltype_palette) 

# Visualizing B cell/CD4+ T cell/CD8+ T cell CTN
CTN_df %>%
  left_join(PF_Xenium_dat_in, by = "CellID") %>%
  filter(ImageID %in% TLS_img_list) %>%
  ggplot(aes(x = X, y = Y, color = label)) +
  geom_point(shape = 16, size= 0.1) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  facet_wrap(~ImageID, scales = "free", nrow = 1) +
  theme_void(base_size = 14) +
  labs(title = CTN_label) +
  scale_color_manual(values = c("CTN1" = "#E69F00", "CTN2" = "#56B4E9", 
                                "Unassigned" = "grey90"), name = NULL) +
  coord_cartesian(expand = FALSE) +
  theme(legend.key.size = unit(2.5, "mm"),
        plot.title = element_markdown(hjust = 0.5),
        panel.background = element_blank(),
        aspect.ratio = 1)

# Show core cell types only
PF_Xenium_dat_in %>%
  filter(ImageID %in% TLS_img_list) %>%
  filter(CellID %in% CTN_df$CellID) %>% # dropping cells with few neighbors
  mutate(Celltype = ifelse(Celltype %in% core_celltypes, as.character(Celltype), "Other")) %>%
  ggplot(aes(x = X, y = Y, color = Celltype)) +
  geom_point(shape = 16, size= 0.1) +
  scale_color_manual(values = c(PF_Xenium_celltype_palette, "Other" = "grey90")) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  facet_wrap(~ImageID, scales = "free", nrow = 1) +
  theme_void(base_size = 14) +
  labs(title = "core cell type triad") +
  coord_cartesian(expand = FALSE) +
  theme(legend.key.size = unit(2.5, "mm"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_blank())
```

![](PF_Xenium_files/figure-html/Visualize%20results-1.png)![](PF_Xenium_files/figure-html/Visualize%20results-2.png)

``` r

CTN_df %>%
  left_join(PF_Xenium_dat_in, by = "CellID") %>%
  count(Celltype, label) %>%
  ggplot(aes(x = n, y = label, fill = factor(Celltype, names(PF_Xenium_celltype_palette)))) +
  geom_col(position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = PF_Xenium_celltype_palette, name = NULL,
                    labels = rename_celltype) +
  labs(title = CTN_label) +
  scale_y_discrete(limits = rev) +
  coord_cartesian(expand = FALSE, clip = "off") +
  guides(fill = guide_legend(ncol = 4)) +
  theme_cowplot(font_size = 12) +
  labs(x = "Proportion of Cells", y = NULL) +
  theme(plot.title = element_markdown(hjust = 0.5),
        legend.key.size = unit(3, "mm"),
        legend.key.spacing.y = unit(0.2, "mm"),
        legend.position = "bottom",
        legend.location = "plot",
        legend.justification = "right")
```

![](PF_Xenium_files/figure-html/proportion%20of%20cells-1.png)

We can then compute the evaluation metrics comparing CTN1 with ground
truth TLS annotations

``` r

library(aricode); library(caret) # for clustering evaluation

cluster_eval_metrics <- CTN_df %>% 
  mutate(label_binary = ifelse(label == "CTN1", "TLS", "Other")) %>%
  left_join(PF_Xenium_annot, by = "CellID") %>%
  filter(ImageID %in% TLS_img_list) %>%
  mutate(annot_binary = ifelse(Annotation_Type == "TLS", "TLS", "Other")) %>%
  split(f = .$ImageID) %>%
  lapply(function(df_sub) {
    # ARI and NMI
    ARI <- ARI(df_sub$label_binary, df_sub$annot_binary)
    NMI <- NMI(df_sub$label_binary, df_sub$annot_binary)
    
    # F-1 score
    conf_mat <- confusionMatrix(
      data = factor(df_sub$label_binary, c("TLS", "Other")),
      reference = factor(df_sub$annot_binary, c("TLS", "Other")))
    F1 <- conf_mat$byClass["F1"] #for multiclass classification problems
    
    data.frame(ARI = ARI, NMI = NMI, F1 = F1)
  }) %>%
  rbindlist(idcol = "ImageID")

cluster_eval_metrics
#>        ImageID       ARI       NMI        F1
#>         <char>     <num>     <num>     <num>
#> 1:   TILD028LA 0.7036145 0.5196893 0.7826087
#> 2:   TILD315MA 0.5556865 0.3321899 0.6670288
#> 3: VUILD104MA1 0.6636726 0.4296721 0.7218789
#> 4:  VUILD110LA 0.3804649 0.1776062 0.4194452
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
#>  [1] caret_7.0-1         lattice_0.22-6      aricode_1.0.3      
#>  [4] RColorBrewer_1.1-3  cowplot_1.2.0       ggrepel_0.9.8      
#>  [7] ggtext_0.1.2        lubridate_1.9.5     forcats_1.0.1      
#> [10] stringr_1.6.0       dplyr_1.2.0         purrr_1.2.1        
#> [13] readr_2.2.0         tidyr_1.3.2         tibble_3.3.1       
#> [16] ggplot2_4.0.2.9000  tidyverse_2.0.0     data.table_1.18.2.1
#> [19] SCOPE.CTN_0.0.1    
#> 
#> loaded via a namespace (and not attached):
#>  [1] pROC_1.19.0.1        rlang_1.1.7          magrittr_2.0.4      
#>  [4] otel_0.2.0           e1071_1.7-17         compiler_4.4.1      
#>  [7] mgcv_1.9-3           systemfonts_1.3.1    pbmcapply_1.5.1     
#> [10] vctrs_0.7.2          reshape2_1.4.5       pkgconfig_2.0.3     
#> [13] fastmap_1.2.0        labeling_0.4.3       rmarkdown_2.30      
#> [16] prodlim_2026.03.11   markdown_2.0         tzdb_0.5.0          
#> [19] ragg_1.5.1           xfun_0.57            cachem_1.1.0        
#> [22] litedown_0.9         jsonlite_2.0.0       recipes_1.3.1       
#> [25] gmp_0.7-5.1          irlba_2.3.7          parallel_4.4.1      
#> [28] R6_2.6.1             bslib_0.10.0         stringi_1.8.7       
#> [31] parallelly_1.46.1    ClusterR_1.3.3       rpart_4.1.23        
#> [34] jquerylib_0.1.4      Rcpp_1.1.1           assertthat_0.2.1    
#> [37] iterators_1.0.14     knitr_1.51           future.apply_1.20.2 
#> [40] FNN_1.1.4.1          Matrix_1.7-0         splines_4.4.1       
#> [43] nnet_7.3-19          timechange_0.4.0     tidyselect_1.2.1    
#> [46] rstudioapi_0.18.0    yaml_2.3.12          timeDate_4052.112   
#> [49] codetools_0.2-20     listenv_0.10.1       plyr_1.8.9          
#> [52] withr_3.0.2          S7_0.2.1             evaluate_1.0.5      
#> [55] future_1.70.0        desc_1.4.3           survival_3.6-4      
#> [58] proxy_0.4-29         xml2_1.5.2           pillar_1.11.1       
#> [61] stats4_4.4.1         foreach_1.5.2        generics_0.1.4      
#> [64] hms_1.1.4            commonmark_2.0.0     scales_1.4.0        
#> [67] globals_0.19.1       class_7.3-22         glue_1.8.0          
#> [70] tools_4.4.1          ModelMetrics_1.2.2.2 gower_1.0.2         
#> [73] fs_1.6.7             grid_4.4.1           ipred_0.9-15        
#> [76] nlme_3.1-164         cli_3.6.5            textshaping_1.0.5   
#> [79] lava_1.8.2           gtable_0.3.6         sass_0.4.10         
#> [82] digest_0.6.39        htmlwidgets_1.6.4    farver_2.1.2        
#> [85] htmltools_0.5.9      pkgdown_2.2.0        lifecycle_1.0.5     
#> [88] hardhat_1.4.2        gridtext_0.1.5       MASS_7.3-60.2
```
