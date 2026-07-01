# Detecting lymphocyte aggregate-like niches in CRC VisiumHD data

This tutorial shows an example of detecting the **lymphocyte aggregate
(LA)-like niches** characterized by the
$`\mbox{B cell/CD4}^+\mbox{ T cell/CD8}^+\mbox{ T cell}`$ triad in the
colorectal cancer (CRC) VisiumHD data using SCOPE. The dataset comes
from [Oliveira et al. *Nat Genet*
2025](https://www.nature.com/articles/s41588-025-02193-3) and is
composed of 3 CRC and 2 normal adjacent tissue samples. For this
tutorial, we will focus on the 3 CRC samples collected from different
patients (P1, P2, P5).

First, we load the packages that will be used for the analysis and
visualization.

``` r

library(SCOPE.CTN)
library(data.table)
library(tidyverse)
library(ggtext)
library(cowplot)
dataset_name <- "CRC_VisiumHD"
data_path <- paste0("vignette_data/", dataset_name, "/")   # !!! Change the directory
source(paste0(data_path, dataset_name, "_celltype_palette.R")) # cell type color palette
```

## 1. Preparing input data

The metadata parquet files for this dataset can be downloaded from
[Zenodo](https://zenodo.org/records/15042463) or
[GitHub](https://github.com/10XGenomics/HumanColonCancer_VisiumHD). We
will simplify the `DeconvolutionLabel1` column and use them as the input
cell type labels. It contains the cell types of the 8-$`\mu`$m binned
data predicted by deconvolution using the single-cell RNA-seq reference.

``` r

library(arrow) # for loading parquet files

CRC_VisiumHD_metadata <- sapply(c("P1", "P2", "P5"), function(x) {
  read_parquet(paste0("https://github.com/10XGenomics/HumanColonCancer_VisiumHD/",
                      "raw/refs/heads/main/MetaData/", x, "CRC_Metadata.parquet"))
  }, simplify = FALSE, USE.NAMES = TRUE) %>%
  data.table::rbindlist(idcol = "ImageID") %>%
  # Make cell ids unique
  mutate(barcode = paste(barcode, ImageID, sep = "_")) %>%
  # Simplifying cell type annotations
  mutate(Celltype = case_when(
    DeconvolutionLabel1 %in% c("Mature B", "Memory B") ~ "B cell",  
    DeconvolutionLabel1 %in% c("Proliferating Macrophages", 
                               "Proliferating Immune II") ~ "Prolif_Immune",
    DeconvolutionLabel1 %in% c("Fibroblast", "Vascular Fibroblast", 
                               "Proliferating Fibroblast") ~ "Fibroblast", 
    DeconvolutionLabel1 %in% c("Enteric Glial", "Enterocyte", "Neuroendocrine") ~ "Neuronal", 
    str_detect(DeconvolutionLabel1, "Tumor") ~ "Tumor",
    DeconvolutionLabel1 %in% c("Epithelial", "Tuft", "Goblet") ~ "Epithelial", 
    DeconvolutionLabel1 %in% c("Smooth Muscle", "SM Stress Response", "vSM") ~ "Smooth Muscle",  
    DeconvolutionLabel1 %in% c("Unknown III (SM)", "Adipocyte") ~ "Other",
    TRUE ~ DeconvolutionLabel1
  )) 
```

The input data frame for SCOPE should contain the following columns:

- `CellID`: cell id
- `ImageID`: sample id
- `X`: X coordinate
- `Y`: Y coordinate
- `Celltype`: cell type labels

We will remove the bins that were assigned as “doublets” or “rejected”
by the deconvolution algorithm and only focus on the singlet bins.

``` r

CRC_VisiumHD_dat_in <- CRC_VisiumHD_metadata %>%
  # Keeping only the singlet bins
  filter(DeconvolutionClass == "singlet") %>%
  select(CellID = barcode, ImageID, X = X, Y = Y, Celltype = Celltype) 
# save(CRC_VisiumHD_dat_in, file = paste0(data_path, dataset_name, "_dat_in.RData"))
```

``` r

head(CRC_VisiumHD_dat_in) # 810,330 bins
#>                      CellID ImageID        X         Y Celltype
#>                      <char>  <char>    <num>     <num>   <char>
#> 1: s_008um_00000_00289-1_P1      P1 69.20113 -309.4240    Tumor
#> 2: s_008um_00000_00290-1_P1      P1 69.44765 -309.4263    Tumor
#> 3: s_008um_00000_00291-1_P1      P1 69.69418 -309.4286    Tumor
#> 4: s_008um_00000_00292-1_P1      P1 69.94070 -309.4309    Tumor
#> 5: s_008um_00000_00293-1_P1      P1 70.18723 -309.4332    Tumor
#> 6: s_008um_00000_00294-1_P1      P1 70.43375 -309.4354    Tumor
```

``` r

CRC_VisiumHD_dat_in %>%
  ggplot(aes(x = X, y = Y, color = Celltype)) +
  geom_point(shape = 16, size= 0.1) +
  scale_color_manual(values = celltype_colors, 
                     breaks = names(celltype_colors)) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  facet_wrap(~ImageID, scales = "free") +
  theme_void(base_size = 14) +
  coord_cartesian(expand = FALSE) +
  theme(legend.key.size = unit(2.5, "mm"),
        panel.background = element_blank())
```

![](CRC_VisiumHD_files/figure-html/Visualize%20celltype-1.png)

## 2. Create neighborhood cell proportion matrix

Next, we use the k-nearest neighbors algorithm to define the spatial
neighborhood and compute the neighborhood cell type composition matrix.
The default number of nearest neighbors $`k`$ is 20, however we can
increase $`k`$ for larger tissue slices. Here we set $`k`$=50.

In addition, we enforce a distance boundary for the nearest neighbors
and exclude the cells falling outside this radius. The boundary is
initially defined as the 99th percentile of all distances from focal
cells to their nearest neighbors. If it falls outside a predefined
range, we reset it to either the minimum (`min.radius`) or maximum
threshold (`max.radius`).

``` r

slide_list <- unique(CRC_VisiumHD_dat_in$ImageID)

cell_neighbor_results_maxdist <- pbmcapply::pbmclapply(
  slide_list, function(s) {
    Build_cell_neighbor_maxdist(
      dat_slide=subset(CRC_VisiumHD_dat_in, ImageID == s),
      knn_num = 50, cell_types_available = unique(CRC_VisiumHD_dat_in$Celltype),
      min.radius = 1,       # lower limit for the boundary
      max.radius = 2,       # upper limit for the boundary
      min.n.neighbors = 15, # minimum no. of neighbors
      prob.list = seq(0, 1, 0.1))},
  mc.cores = 3) # set mc.cores to 1 on Windows
# save(cell_neighbor_results_maxdist, 
#      file = paste0(data_path, dataset_name, "_cell_neighbor_results.RData"))
```

We can examine the distribution of 99th percentile of all distance to
the nearest neighbors to select the threshold parameters min.radius and
max.radius accordingly. The 99th percentiles are 2.466, 1.679, and 1.939
for P1, P2, P5 respectively, while the other percentiles remain nearly
identical across the three samples. Therefore we set `max.radius` to 2
and `min.radius` to 1.

``` r

length(cell_neighbor_results_maxdist)
#> [1] 3
sapply(cell_neighbor_results_maxdist, "[[", "dist.quantiles") %>%
  as.data.frame() %>% 
  rownames_to_column("perc") %>%
  mutate(perc = as.numeric(str_remove(perc, "%"))) %>%
  pivot_longer(2:ncol(.), names_to = "ImageID", values_to = "value") %>%
  ggplot(aes(x = perc, y = value)) +
  geom_point(aes(color = ImageID)) +
  scale_color_discrete(labels = slide_list) +
  scale_x_continuous(labels = function(x) {paste0(x, "%")}) +
  labs(y = "Distance to\nNearest Neighbors", x = "Percentile") +
  theme_cowplot(font_size = 12)
```

![](CRC_VisiumHD_files/figure-html/Examine%20distance%20distribution-1.png)

``` r

# Combine neighborhood proportion table from the three samples
cell_neighbor_table_maxdist <- lapply(
  cell_neighbor_results_maxdist, "[[", "prop.df") %>%
  do.call("rbind", .)
```

Cells that retain fewer than `min.n.neighbors` neighbors within the
distance boundary will have their neighborhood cell-type proportions set
to `NaN` and will be excluded from subsequent niche clustering. Here we
set `min.n.neighbors` to 10.

``` r

cells_to_exclude = apply(cell_neighbor_table_maxdist, 1, function(x) {any(is.na(x))}) %>%
  as.data.frame() %>%
  `colnames<-`("is_na") %>%
  rownames_to_column("CellID") %>%
  left_join(CRC_VisiumHD_dat_in[, c("CellID", "X", "Y", "ImageID")], by = "CellID")

# Number of cells to exclude
table(cells_to_exclude$is_na)
#> 
#>  FALSE   TRUE 
#> 809156   1174

cells_to_exclude %>%
  ggplot(aes(x = X, y = Y, color = is_na)) +
  geom_point(shape = 16, size= 0.1) +
  scale_color_manual(values = c("brown2", "grey90"), 
                     breaks = c(TRUE, FALSE), name = "Cells To\nExclude") +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  facet_wrap(~ImageID, scales = "free") +
  theme_void() +
  coord_cartesian(expand = FALSE) +
  theme(legend.key.size = unit(2.5, "mm"),
        panel.background = element_blank())
```

![](CRC_VisiumHD_files/figure-html/Examine%20loose%20cells-1.png)

## 3. Cell type triad niche (CTN) detection using SCOPE

Next, we apply SCOPE to the
$`\mbox{B cell/CD4}^+\mbox{ T cell/CD8}^+\mbox{ T cell}`$ triad using
default configurations for the following parameters:

- `min.prop`: The minimum neighborhood cell-type proportion threshold.
  Proportions below this value are set to 0 (default: 0.3).
- `mgcv_df`: The dimension of the basis used for thin-plate regression
  spline smoothing of the core cell-type proportions (default: 15; see
  [`mgcv::s()`](https://rdrr.io/pkg/mgcv/man/s.html) for details).

We set the number of clusters `clusternum` at 3 to separate the B cell
and T cell zones.

``` r

core_celltypes <- c("B cell", "CD4 T cell", "CD8 T cell")

CTN_df <- run_SCOPE(
  combination_table = as.data.frame(matrix(core_celltypes, nrow = 1)),
  cell_neighbor_table = cell_neighbor_table_maxdist,
  dat_in = CRC_VisiumHD_dat_in,
  save_results = FALSE, 
  clusternum = 3, min.prop = 0.3,
  linear = FALSE, mgcv_df = 15, iter.max = 100)
# save(CTN_df, file = paste0(data_path, dataset_name, "_B_CD4T_CD8T.RData"))
```

Here we can see that the
$`\mbox{B cell/CD4}^+\mbox{ T cell/CD8}^+\mbox{ T cell}`$ CTNs closely
aligns with the lymphocyte aggregates in the CRC samples, where CTN1
resembles the B cell zone and CTN2 resembles the surrounding T cell
zone.

``` r

# Number of cells per CTN
table(CTN_df$label)
#> 
#>       CTN1       CTN2 Unassigned 
#>       9121       4242     795793
```

``` r

# CTN label with colored dots separating the cell types
CTN_label <- color_CTN_names(
  paste(core_celltypes, collapse = "_"), celltype_colors) 

# Visualizing B cell/CD4+ T cell/CD8+ T cell CTN
CTN_df %>%
  left_join(CRC_VisiumHD_dat_in, by = "CellID") %>%
  ggplot(aes(x = X, y = Y, color = label)) +
  geom_point(shape = 16, size= 0.1) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  facet_wrap(~ImageID, scales = "free") +
  theme_void(base_size = 14) +
  labs(title = CTN_label) +
  scale_color_manual(values = c("CTN1" = "#E69F00", "CTN2" = "#56B4E9", 
                                "Unassigned" = "grey90"), name = NULL) +
  coord_cartesian(expand = FALSE) +
  theme(legend.key.size = unit(2.5, "mm"),
        plot.title = element_markdown(hjust = 0.5),
        panel.background = element_blank())

# Show core cell types only
CRC_VisiumHD_dat_in %>%
  mutate(Celltype = ifelse(Celltype %in% core_celltypes, Celltype, "Other")) %>%
  ggplot(aes(x = X, y = Y, color = Celltype)) +
  geom_point(shape = 16, size= 0.1) +
  scale_color_manual(values = celltype_colors, 
                     breaks = names(celltype_colors)) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  facet_wrap(~ImageID, scales = "free") +
  theme_void(base_size = 14) +
  labs(title = "core cell type triad") +
  coord_cartesian(expand = FALSE) +
  theme(legend.key.size = unit(2.5, "mm"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_blank())
```

![](CRC_VisiumHD_files/figure-html/Visualize%20results-1.png)![](CRC_VisiumHD_files/figure-html/Visualize%20results-2.png)

``` r

CTN_df %>%
  left_join(CRC_VisiumHD_dat_in, by = "CellID") %>%
  count(Celltype, label) %>%
  ggplot(aes(x = n, y = label, fill = factor(Celltype, names(celltype_colors)))) +
  geom_col(position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = celltype_colors, name = NULL) +
  labs(title = CTN_label) +
  scale_y_discrete(limits = rev) +
  coord_cartesian(expand = FALSE, clip = "off") +
  guides(fill = guide_legend(ncol = 3)) +
  theme_cowplot(font_size = 12) +
  labs(x = "Proportion of Cells", y = NULL) +
  theme(plot.title = element_markdown(hjust = 0.5),
        legend.key.size = unit(3, "mm"),
        legend.key.spacing.y = unit(0.5, "mm"),
        legend.position = "bottom",
        legend.location = "plot",
        legend.justification = "right")
```

![](CRC_VisiumHD_files/figure-html/proportion%20of%20cells-1.png)

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
#>  [1] tidyselect_1.2.1      farver_2.1.2          S7_0.2.1             
#>  [4] fastmap_1.2.0         digest_0.6.39         timechange_0.4.0     
#>  [7] lifecycle_1.0.5       cluster_2.1.6         ClusterR_1.3.3       
#> [10] magrittr_2.0.4        compiler_4.4.1        rlang_1.1.7          
#> [13] sass_0.4.10           tools_4.4.1           yaml_2.3.12          
#> [16] knitr_1.51            FNN_1.1.4.1           labeling_0.4.3       
#> [19] htmlwidgets_1.6.4     xml2_1.5.2            RColorBrewer_1.1-3   
#> [22] withr_3.0.2           BiocGenerics_0.52.0   desc_1.4.3           
#> [25] grid_4.4.1            stats4_4.4.1          colorspace_2.1-2     
#> [28] scales_1.4.0          iterators_1.0.14      MASS_7.3-60.2        
#> [31] cli_3.6.5             rmarkdown_2.30        crayon_1.5.3         
#> [34] ragg_1.5.1            generics_0.1.4        otel_0.2.0           
#> [37] rstudioapi_0.18.0     tzdb_0.5.0            rjson_0.2.23         
#> [40] commonmark_2.0.0      cachem_1.1.0          splines_4.4.1        
#> [43] assertthat_0.2.1      parallel_4.4.1        matrixStats_1.5.0    
#> [46] vctrs_0.7.2           Matrix_1.7-0          jsonlite_2.0.0       
#> [49] litedown_0.9          IRanges_2.40.1        GetoptLong_1.1.0     
#> [52] hms_1.1.4             S4Vectors_0.44.0      irlba_2.3.7          
#> [55] clue_0.3-66           pbmcapply_1.5.1       systemfonts_1.3.1    
#> [58] foreach_1.5.2         jquerylib_0.1.4       ggdendro_0.2.0       
#> [61] glue_1.8.0            pkgdown_2.2.0         codetools_0.2-20     
#> [64] stringi_1.8.7         shape_1.4.6.1         gtable_0.3.6         
#> [67] gmp_0.7-5.1           ComplexHeatmap_2.22.0 pillar_1.11.1        
#> [70] htmltools_0.5.9       circlize_0.4.17       R6_2.6.1             
#> [73] textshaping_1.0.5     doParallel_1.0.17     evaluate_1.0.5       
#> [76] lattice_0.22-6        markdown_2.0          png_0.1-9            
#> [79] gridtext_0.1.5        bslib_0.10.0          Rcpp_1.1.1           
#> [82] nlme_3.1-164          mgcv_1.9-3            xfun_0.57            
#> [85] fs_1.6.7              pkgconfig_2.0.3       GlobalOptions_0.1.3
```
