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
library(cowplot)
library(ggdendro)
library(ggtext)
library(patchwork)
library(ComplexHeatmap)
dataset_name <- "NSCLC_IMC"
data_path <- paste0("vignette_data/", dataset_name, "/")   # !!! Change the directory
source(paste0(data_path, dataset_name, "_celltype_palette.R")) # cell type color palette
```

## 1. Prepare input data

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
Here the number of nearest neighbors \\k\\ is set at default value 20.

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
      dat_slide=NSCLC_IMC_dat_in[NSCLC_IMC_dat_in$ImageID==slide_list[i], ], knn_num = 20, 
      min.radius = 100, max.radius = 120, min.n.neighbors = 10,
      cell_types_available = levels(NSCLC_IMC_dat_in$Celltype),
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
increase exponentially. The NSCLC IMC dataset has Collagen_CAF,
hypoxic_CAF, mCAF, SMA_CAF, tpCAF, iCAF, vCAF, dCAF, hypoxic_tpCAF,
IDO_CAF, PDPN_CAF, HEV, Blood, Lymphatic, Myeloid, Neutrophil, Bcell,
CD8, CD4, normal, hypoxic, Other available cell types, which can form
1540 triads. Therefore, we utilize a colocation quotient (CLQ)
permutation test to prioritize cell-type triads (CTs) that exhibit
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
#> # A tibble: 5 × 5
#> # Groups:   Celltype, Celltype_neighbor [1]
#>   Celltype Celltype_neighbor ImageID  p_perm p_perm.adj
#>   <chr>    <chr>             <chr>     <dbl>      <dbl>
#> 1 _cell    _lood             175A_1    0         0     
#> 2 _cell    _lood             175A_10   1         1     
#> 3 _cell    _lood             175A_100  0.002     0.0281
#> 4 _cell    _lood             175A_102  0         0     
#> 5 _cell    _lood             175A_103  0.422     1
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

## 4. Cell type triad niche (CTN) detection using SCOPE

Next, we apply SCOPE for the systematic screening of CTNs. We utilize
the default configurations for most parameters, including:

- `min.prop`: The minimum neighborhood cell-type proportion threshold.
  Proportions below this value are set to 0 (default: 0.3).
- `mgcv_df`: The dimension of the basis used for thin-plate regression
  spline smoothing of the driver cell-type proportions (default: 15; see
  [`mgcv::s()`](https://rdrr.io/pkg/mgcv/man/s.html) for details).

Here we set the number of clusters `clusternum` at 2. When setting
`save_results = TRUE`, the clustering results of each CT will be saved
as individual `.RData` files under `output_dir`.

``` r

# Load prefiltered CTs
load(paste0(data_path, dataset_name, "_combination_table_filtered.RData"))

run_SCOPE(
  combination_table = combination_table_filtered,
  cell_neighbor_table = cell_neighbor_table_maxdist,
  dat_in = NSCLC_IMC_dat_in,
  save_results = TRUE,
  output_dir = paste0(data_path, "CTN_K_2"),
  clusternum = 2, min.prop = 0.3,
  linear = FALSE, mgcv_df = 15, iter.max = 100)
```

We compile the individual outputs and merge them into a single data
frame. The first five columns match the original input data, and the
subsequent columns contain the assigned CTN labels.

``` r

Full_CTN_table_mgcv <- data.frame(dir = list.files(
  paste0(data_path, "_CTN_K_2"), ".RData", full.names = TRUE)) %>%
  mutate(group = row_number() %/% 30) %>%
  split(f = .$group) %>%
  lapply(pull, dir) %>%
  lapply(function(dirs) {
    pbmclapply(dirs, function(dir) {
      load(dir)
      CTN_table %>% select(CellID, CTN, label) %>%
        mutate(CTN = str_replace_all(CTN, "_CAF", ".CAF")) %>%
        mutate(CTN = str_replace_all(CTN, "hypoxic_tpCAF", "hypoxic.tpCAF"))
    }, mc.cores = 10) %>%
      do.call(rbind, .)  %>% `rownames<-`(NULL) %>%
      pivot_wider(names_from = "comb", values_from = "label")
})

Full_CTN_table_mgcv <- lapply(Full_CTN_table_mgcv, function(df) {df %>% select(-CellID)}) %>%
  `names<-`(NULL) %>%
  do.call(cbind, .) %>%
  cbind(data.frame(CellID = Full_CTN_table_mgcv[[1]]$CellID), .)

Full_CTN_table_mgcv <- temp %>%
  left_join(NSCLC_IMC_dat_in, by = "CellID") %>%
  mutate(Celltype = as.character(Celltype)) %>%
  mutate(Celltype = str_replace(as.character(Celltype), "_CAF", ".CAF")) %>%
  mutate(Celltype = str_replace(Celltype, "hypoxic_tpCAF", "hypoxic.tpCAF")) 
  select(all_of(colnames(NSCLC_IMC_dat_in)), everything()) %>%
  
save(Full_CTN_table_mgcv, file = paste0(
  data_path, dataset_name, "_Full_CTN_table_mgcv.RData"))
rm(temp)
```

## 5. Merging highly overlapped CTNs

Sometimes CTNs driven by two different triads can be highly similar. To
merge highly overlapping CTNs, we calculate Jaccard distance between
every pair of CTNs and perform hierarchical clustering to reduce
redundancy.

``` r

load(paste0(data_path, dataset_name, "_Full_CTN_table_mgcv.RData"))

# All cell type triads
CTN_list <- setdiff(colnames(Full_CTN_table_mgcv), colnames(NSCLC_IMC_dat_in))

CTN_jaccard_mat_hclust <- jaccard_dist_hclust(
  Full_CTN_table = Full_CTN_table_mgcv, 
  CTN_ls = CTN_list, num_cores = 5)

save(CTN_jaccard_mat_hclust, file = paste0(
     data_path, dataset_name, "_CTN_jaccard_mat_hclust.RData"))
```

Silhouette scores are utilized to determine the optimal number of
clusters. The final cluster count is set to 77, balancing both the mean
and median Silhouette scores.

``` r

n_clusters = 77 # Final number of hierarchical clusters
draw_silhoutte_score(CTN_jaccard_mat_hclust$hclust_silhouette,
                     num_clusters = n_clusters)
```

![](NSCLC_IMC_files/figure-html/Silouhette%20score-1.png)

Below is a heatmap visualization of the Jaccard similarity indices
(1-Jaccard index) and hierarchical clustering results. Black boxes
outline final CTN groupings. Asterisks (\*) and dots (.) indicate highly
similar (Jaccard index \> 0.95) and moderately similar (0.9 \\\leq\\
Jaccard index \\\leq\\ 0.95) CTN pairs, respectively.

``` r

htmap_jaccard_index <- draw_htmap_jaccard_index(
  jaccard_mat = CTN_jaccard_mat_hclust$jaccard_mat_mgcv,
  cl = CTN_jaccard_mat_hclust$cl,
  num_clusters = n_clusters)  
```

![](NSCLC_IMC_files/figure-html/Jaccard%20similarity%20matrix%20heatmap-1.png)

``` r

# Hierarchical clustering results
htmap_row_order <- unlist(row_order(htmap_jaccard_index))
CTN_merged_label <- as.data.frame(cutree(CTN_jaccard_mat_hclust$cl, k = n_clusters)) %>%
  rownames_to_column("CTN") %>%
  `colnames<-`(c("CTN", "cluster")) %>%
  mutate(CTN = factor(CTN, rownames(
    CTN_jaccard_mat_hclust$jaccard_mat_mgcv)[htmap_row_order])) %>% 
  arrange(CTN) %>%
  mutate(cluster = factor(cluster, unique(.$cluster))) %>%
  mutate(cluster = factor(as.numeric(cluster))) %>%
  mutate(CTN = as.character(CTN))

save(CTN_merged_label, file= paste0(data_path, dataset_name, "_CTN_merged_label.RData"))
```

After merging highly overlapping initial CTNs using hierarchical
clustering, the next step is to assign new names to these consolidated
niches. The package uses the
[`CTN_annotate()`](../reference/CTN_annotate.md) function to rename
consolidated CTNs based on the most frequent core cell types within
their original triads. If a tie occurs between cell types, the function
automatically resolves it by identifying the cell type with the highest
overall proportion within the consolidated CTN.

First, we compute the cell type proportion within the newly consolidated
CTNs, which will be used as one of the inputs to the
[`CTN_annotate()`](../reference/CTN_annotate.md) function.

``` r

CTN_merged_celltype_proportion <- CTN_merged_label %>%
  split(f = .$cluster) %>%
  lapply(function(CTN_merged_label_sub) {
  Full_CTN_table_mgcv[, c("Celltype", CTN_merged_label_sub$CTN)] %>% 
      filter(if_any(everything(), ~ .x == "CTN")) %>%
      count(Celltype) }) %>%
  rbindlist(idcol = "cluster") %>%
  group_by(cluster) %>%
  mutate(prop = n/sum(n)) %>% select(-n) %>%
  ungroup


save(CTN_merged_celltype_proportion, file = paste0(
     data_path, dataset_name, "_CTN_merged_celltype_proportion.RData"))
```

Additional inputs include a data frame containing the original CTN names
and their hierarchical cluster labels, and a character vector defining
the levels or ordering of the cell types. The `CTN_annotate` function
returns a list of two data frames:

- `annotation`: A mapping data frame that matches the original CTN names
  and their hierarchical cluster labels to the final consolidated names
- `core_celltypes_summ`: A summary table focusing on clusters that
  contain at least 2 initial CTNs. It displays how many times each core
  cell type appeared in the original CTN names and tracks their relative
  proportions within the merged CTN.

``` r

CTN_merged_label_annotated <- CTN_annotate(
  CTN_cluster_label = CTN_merged_label, 
  CTN_merged_celltype_prop = CTN_merged_celltype_proportion, 
  celltype_levs = names(NSCLC_IMC_palette))
```

If the top two most frequent cell types only appear twice in the
original CTN names, the function flags the consolidated CTN with an
asterisk (\*). We suggest double-checking the flagged CTNs and manually
renaming them if necessary.

Here we also manually rename a few CTNs with high blood (“Blood”),
lymphatic endothelial (“Lymphatic”), and vCAF proportions with “Vessel”.

``` r

head(CTN_merged_label_annotated$core_celltypes_summary, n = 15)

CTN_annotation_detailed <- CTN_merged_label_annotated$annotation %>%
  mutate(annotation_final = case_when(
    cluster == 12 ~ "Bcell_CD4_Vessel",
    cluster == 14 ~ "Bcell_Vessel_vCAF",
    cluster == 18 ~ "CD8_Vessel_vCAF",
    cluster == 33 ~ "T_Myeloid_vCAF",
    cluster == 36 ~ "CD8_Neutrophil_Vessel",
    cluster == 41 ~ "Neutrophil_Vessel_vCAF",
    cluster == 52 ~ "Neutrophil_Vessel_SMA.CAF",
    cluster == 61 ~ "Vessel_SMA.CAF_vCAF",
    TRUE ~ annotation)) 

# save(CTN_annotation_detailed, file = paste0(
#  data_path, dataset_name, "_CTN_annotation_detailed.RData"))
```

``` r

Full_CTN_table_mgcv_merged <- CTN_annotation_detailed %>%
  split(f = .$cluster) %>%
  sapply(function(annot_df) {
    Full_CTN_table_mgcv[, c("CellID", annot_df$CTN)] %>%
      column_to_rownames("CellID") %>%
      apply(1, function(x) {ifelse(any(x == "CTN"), "CTN", "Unassigned")}) %>%
      as.data.frame() %>%
      `colnames<-`(annot_df$annotation_final[1])
  }, simplify = FALSE, USE.NAMES = FALSE) %>%
  do.call(cbind, .) %>%
  rownames_to_column("CellID")

Full_CTN_table_mgcv_merged <- Full_CTN_table_mgcv_merged %>%
  inner_join(NSCLC_IMC_dat_in, by = "CellID") %>%
  select(all_of(colnames(NSCLC_IMC_dat_in)), everything())

save(Full_CTN_table_mgcv_merged, 
     file = paste0(data_path, dataset_name, "_Full_CTN_table_mgcv_merged.RData"))
```

We can further perform hierarchical clustering on the Jaccard distance
between the consolidated CTNs and grouping them into major categories.

``` r

load(paste0(data_path, dataset_name, "_Full_CTN_table_mgcv_merged.RData"))
CTN_merged_list <- setdiff(colnames(Full_CTN_table_mgcv_merged), colnames(NSCLC_IMC_dat_in))
CTN_merged_jaccard_mat_hclust <- jaccard_dist_hclust(
  Full_CTN_table = Full_CTN_table_mgcv_merged, 
  CTN_ls = CTN_merged_list, num_cores = 5)
save(CTN_merged_jaccard_mat_hclust, file = paste0(
  data_path, dataset_name, "_CTN_merged_jaccard_mat_hclust.RData"))
```

These categories can be annotated by examining the cell type abundance
of the constituent CTNs.

``` r

CTN_annotation_merged <- CTN_annotation_detailed %>%
    # Grouping CTNs based on cell type composition  
    mutate(CTN_group = case_when(cluster %in% 1:3 ~ "iCAF-rich",
    cluster %in% 4:6 ~ "vCAF-rich",
    cluster %in% 7:8 ~ "CD4.T-rich",
    cluster %in% 9:10 ~ "Bcell/Neutrophil",
    cluster %in% 11:14 ~ "Bcell-rich",
    cluster %in% 15:16 ~ "Tcell-rich",
    cluster %in% 17:20 ~ "CD8.T-rich",
    cluster %in% 21:23 ~ "Myeloid/mCAF",
    cluster %in% 24    ~ "Bcell/Myeloid", 
    cluster %in% 25:33 ~ "Myeloid-rich",
    cluster %in% 34:41 ~ "Neutrophil-rich",
    cluster %in% 42:43 ~ "Bcell/SMA.CAF",
    cluster %in% 44:46 ~ "Myeloid/SMA.CAF",
    cluster %in% 47:50 ~ "SMA.CAF/mCAF",
    cluster %in% 51:55 ~ "Neutrophil/SMA.CAF",
    cluster %in% 56:62 ~ "SMA.CAF-rich",
    cluster %in% 63:66 ~ "Bcell/mCAF",
    cluster %in% 67:70 ~ "Neutrophil/mCAF",
    cluster %in% 71:77 ~ "mCAF-rich",
    TRUE ~ "Mixed"))  %>%
  distinct(cluster, annotation_final, CTN_group)%>%
  rename("annotation" = "annotation_final")


CTN_merged_celltype_proportion_w_annot <- CTN_merged_celltype_proportion %>%
  left_join(CTN_annotation_merged[, 1:2], by = "cluster") 

CTN_merged_dend <- draw_CTN_dendro(
  cl = CTN_merged_jaccard_mat_hclust$cl, h = 0.6,
  dot_fontsize = 11, y_lim = c(NA, 1),
  x_expand = c(0.007, 0.007),
  celltype_palette = c(NSCLC_IMC_palette, "Vessel" = "#393B79", 
                       "CAF" = "#780116", "T" = "royalblue"),
  CTN_annotation = CTN_annotation_merged,
  CTN_group_var = "CTN_group",
  CTN_group_palette = NSCLC_IMC_CTN_group_palette, 
  rename_celltype = NSCLC_IMC_rename_celltype, 
  show_text = FALSE)


p_CTN_merged_celltype_prop <- draw_CTN_celltype_prop(
  CTN_merged_celltype_prop = CTN_merged_celltype_proportion_w_annot, 
  CTN_dend = CTN_merged_dend, 
  celltype_palette = NSCLC_IMC_palette,
  rename_celltype = NSCLC_IMC_rename_celltype,
  facet_var = "CTN_group_label", x_expand = c(0.05, 0.05))

p_CTN_merged_celltype_prop + CTN_merged_dend$p + 
  plot_layout(guides = "collect", widths = c(1, 0.4))
```

![](NSCLC_IMC_files/figure-html/CTN%20dendrogram-1.png)

We can also visualize similarities in the cell type compositions of the
CTNs using Uniform Manifold Approximation and Projection (UMAP).

``` r

draw_CTN_umap(Full_CTN_table = Full_CTN_table_mgcv_merged,
              CTN_merged_celltype_prop = CTN_merged_celltype_proportion,
              CTN_ls = CTN_merged_list, CTN_group_var = "CTN_group",
              CTN_annotation = CTN_annotation_merged,
              CTN_group_palette = NSCLC_IMC_CTN_group_palette,
              umap_n_neighbors = 10, umap_seed = 42) +
  guides(fill = guide_legend(ncol = 2, override.aes = list(size = 3)),
         radius = guide_legend(ncol = 1)) +
  theme(legend.position = "bottom",
        legend.key.spacing.y = unit(0.2, "mm"),
        legend.direction = "vertical",
        legend.title.position = "top")
```

![](NSCLC_IMC_files/figure-html/CTN%20cell%20type%20composition%20umap-1.png)

## 6. Identify clinically-relevant CTNs

After obtaining the final CTN labels, we will next examine their
prognostic value using the detailed patient clinical data provided by
the original data source. The analysis is restricted to treatment-naive
samples from patients with lung adenocarcinoma (LUAD) or lung squamous
cell carcinoma (LUSC).

``` r

# Image-level clinical data
clinical_data_LUAD_LUSC <- fread(paste0(data_path, "clinical_data_ROI.csv")) %>%
  mutate(DX_group = case_when(
    DX.name == "Adenocarcinoma" ~ "LUAD",
    DX.name == "Squamous cell carcinoma" ~ "LUSC",
    TRUE ~ "Other")) %>%
  mutate(Ev.DFS = ifelse(OS > DFS, 1, 0)) %>%
  rename("ImageID" = "Tma_ac") %>%
  filter(DX_group %in% c("LUAD", "LUSC")) %>%
  filter(NeoAdj == "NoNeoAdjuvantTherapy") 


table(clinical_data_LUAD_LUSC$DX_group)
#> 
#> LUAD LUSC 
#> 1026  681
```

To test whether the CTNs predict patient survival, we fit a Cox
proportional hazards regression model using patient-level CTN presence
as the primary predictor. For a given patient, a CTN is considered
present if it contains 100 or more cells within a single tissue in at
least one tissue section.

While we conduct a simple univariate Cox regression here for
demonstration purposes, users can easily extend this to a multivariate
analysis by specifying confounding variables via the `confounder`
argument.

``` r

# Overall survival
coxph_res_OS <- CTN_presence_coxph(
    Full_CTN_table = Full_CTN_table_mgcv_merged, 
    CTN_ls = CTN_merged_list, 
    clinical_df = clinical_data_LUAD_LUSC,
    time = "OS", event = "Ev.O",
    patient_id = "Patient_ID", confounder = NULL,
    min.n.patients = 3, min.n.cells = 100)

# Disease-free survival
coxph_res_DFS <- CTN_presence_coxph(
    Full_CTN_table = Full_CTN_table_mgcv_merged, 
    CTN_ls = CTN_merged_list, 
    clinical_df = clinical_data_LUAD_LUSC,
    time = "DFS", event = "Ev.DFS",
    patient_id = "Patient_ID", confounder = NULL,
    min.n.patients = 3, min.n.cells = 100)

save(coxph_res_OS, coxph_res_DFS, file = paste0(
     data_path, dataset_name, "_CTN_merged_coxph_res.RData"))
```

The following plots display the CTNs associated with the overall
survival and disease-free survival of NSCLC patients, respectively.

``` r

load(paste0(data_path, dataset_name, "_CTN_merged_coxph_res.RData"))
# Overall survival
draw_CTN_coxph_logHR(
  coxph_res = coxph_res_OS,
  CTN_merged_celltype_prop = CTN_merged_celltype_proportion_w_annot,
  alpha = 0.1, top = 10, 
  celltype_levs = names(NSCLC_IMC_palette), dot_fontsize = 14,
  celltype_palette = c(NSCLC_IMC_palette, "Vessel" = "#393B79", 
                       "CAF" = "#780116", "T" = "royalblue"),
  rename_celltype = NSCLC_IMC_rename_celltype,
  line_width = 0.4, y_text_size = NULL,
  width_ratio = c(0.45, 0.55)) 
```

![](NSCLC_IMC_files/figure-html/Visualize%20cox%20regression%20results%20on%20OS-1.png)

``` r

# Disease-free survival
draw_CTN_coxph_logHR(
  coxph_res = coxph_res_DFS,
  CTN_merged_celltype_prop = CTN_merged_celltype_proportion_w_annot,
  alpha = 0.1, top = 10, 
  celltype_levs = names(NSCLC_IMC_palette), dot_fontsize = 14,
  celltype_palette = c(NSCLC_IMC_palette, "Vessel" = "#393B79", 
                       "CAF" = "#780116", "T" = "royalblue"),
  rename_celltype = NSCLC_IMC_rename_celltype,
  line_width = 0.6, y_text_size = NULL,
  width_ratio = c(0.45, 0.55))
```

![](NSCLC_IMC_files/figure-html/Visualize%20cox%20regression%20results%20on%20DFS-1.png)

Below are two example CTNs associated with the disease-free survival.
One is \\\mbox{B cell/CD4}^+\mbox{ T cell/CD8}^+\mbox{ T cell}\\, which
resembles the lymphocyte aggregates.

``` r

p_CTN_celltype <- draw_CTN_celltype(
  Full_CTN_table = Full_CTN_table_mgcv_merged,
  CTN = "Bcell_CD4_CD8", 
  img_list = c("175B_117", "88A_119", "176C_1", "88C_17", "87B_135", "176C_8"),
  celltype_palette = NSCLC_IMC_palette, 
  rename_celltype = NSCLC_IMC_rename_celltype,
  scales = "free", aspect.ratio = 1, legend.ncol = 2)
p_CTN_celltype$p_celltype / p_CTN_celltype$p_CTN +
  plot_layout(guides = "collect") 
```

![](NSCLC_IMC_files/figure-html/Visualize%20Bcell_CD4_CD8-1.png) The
other is blood endothelial \\\mbox{(BEC)/CD4}^+\mbox{ T
cell/CD8}^+\mbox{ T cell}\\, which preferentialy co-localized with
normoxic tumor cells compared to hypoxic tumor cells.

``` r

p_CTN_celltype <- draw_CTN_celltype(
  Full_CTN_table = Full_CTN_table_mgcv_merged,
  CTN = "CD4_CD8_Blood", 
  img_list = c("178C_85", "88B_125", "86A_38", "175A_27", "175A_31", "87B_73"),
  celltype_palette = NSCLC_IMC_palette, 
  rename_celltype = NSCLC_IMC_rename_celltype,
  scales = "free", aspect.ratio = 1, legend.ncol = 2)
p_CTN_celltype$p_celltype / p_CTN_celltype$p_CTN +
  plot_layout(guides = "collect") 
```

![](NSCLC_IMC_files/figure-html/Visualize%20CD4_CD8_Blood-1.png)

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
#> [1] grid      stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] ComplexHeatmap_2.22.0 patchwork_1.3.2       ggtext_0.1.2         
#>  [4] ggdendro_0.2.0        cowplot_1.2.0         lubridate_1.9.5      
#>  [7] forcats_1.0.1         stringr_1.6.0         dplyr_1.2.0          
#> [10] purrr_1.2.1           readr_2.2.0           tidyr_1.3.2          
#> [13] tibble_3.3.1          ggplot2_4.0.2.9000    tidyverse_2.0.0      
#> [16] data.table_1.18.2.1   SCOPE.CTN_0.0.1      
#> 
#> loaded via a namespace (and not attached):
#>  [1] tidyselect_1.2.1    farver_2.1.2        S7_0.2.1           
#>  [4] fastmap_1.2.0       digest_0.6.39       timechange_0.4.0   
#>  [7] lifecycle_1.0.5     cluster_2.1.6       survival_3.6-4     
#> [10] ClusterR_1.3.3      magrittr_2.0.4      compiler_4.4.1     
#> [13] rlang_1.1.7         sass_0.4.10         tools_4.4.1        
#> [16] utf8_1.2.6          yaml_2.3.12         knitr_1.51         
#> [19] FNN_1.1.4.1         labeling_0.4.3      htmlwidgets_1.6.4  
#> [22] xml2_1.5.2          RColorBrewer_1.1-3  withr_3.0.2        
#> [25] BiocGenerics_0.52.0 desc_1.4.3          stats4_4.4.1       
#> [28] colorspace_2.1-2    scales_1.4.0        iterators_1.0.14   
#> [31] MASS_7.3-60.2       cli_3.6.5           rmarkdown_2.30     
#> [34] crayon_1.5.3        ragg_1.5.1          generics_0.1.4     
#> [37] otel_0.2.0          RSpectra_0.16-2     rstudioapi_0.18.0  
#> [40] tzdb_0.5.0          rjson_0.2.23        commonmark_2.0.0   
#> [43] cachem_1.1.0        splines_4.4.1       assertthat_0.2.1   
#> [46] parallel_4.4.1      matrixStats_1.5.0   vctrs_0.7.2        
#> [49] Matrix_1.7-0        jsonlite_2.0.0      litedown_0.9       
#> [52] S4Vectors_0.44.0    IRanges_2.40.1      hms_1.1.4          
#> [55] GetoptLong_1.1.0    ggrepel_0.9.8       clue_0.3-66        
#> [58] irlba_2.3.7         pbmcapply_1.5.1     systemfonts_1.3.1  
#> [61] foreach_1.5.2       jquerylib_0.1.4     glue_1.8.0         
#> [64] pkgdown_2.2.0       codetools_0.2-20    uwot_0.2.4         
#> [67] shape_1.4.6.1       stringi_1.8.7       gtable_0.3.6       
#> [70] gmp_0.7-5.1         pillar_1.11.1       htmltools_0.5.9    
#> [73] circlize_0.4.17     R6_2.6.1            textshaping_1.0.5  
#> [76] doParallel_1.0.17   evaluate_1.0.5      lattice_0.22-6     
#> [79] markdown_2.0        png_0.1-9           gridtext_0.1.5     
#> [82] bslib_0.10.0        Rcpp_1.1.1          nlme_3.1-164       
#> [85] mgcv_1.9-3          xfun_0.57           fs_1.6.7           
#> [88] pkgconfig_2.0.3     GlobalOptions_0.1.3
```
