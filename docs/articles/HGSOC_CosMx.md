# Systematic CTN screening in HGSOC CosMx data

This tutorial demonstrates how to conduct a systematically screening of
Celltype Triad Niches (CTNs) in a high-grade serous tubo-ovarian Cancer
(HGSOC) CosMx Spatial Molecular Imaging (SMI) dataset using SCOPE. The
dataset is from the Discovery cohort of [Yeh et al. *Nature Immunology*
2024](https://www.nature.com/articles/s41590-024-01943-5) and consists
of 100 tissue samples from 58 patients.

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
dataset_name <- "HGSOC_CosMx"
data_path <- paste0("vignette_data/", dataset_name, "/")   # !!! Change the directory
source(paste0(data_path, dataset_name, "_celltype_palette.R")) # cell type color palette
```

## 1. Prepare input data

The Seurat object of the Discovery cohort (`ST_Discovery_so.rds`) can be
downloaded from [Zenodo](https://zenodo.org/records/12613839). Metadata
of the tissue sections are included in [Supplementary Table
2b](https://static-content.springer.com/esm/art%3A10.1038%2Fs41590-024-01943-5/MediaObjects/41590_2024_1943_MOESM3_ESM.xlsx)
(`41590_2024_1943_MOESM3_ESM.xlsx`).

``` r

library(Seurat)
ST_Discovery_so <- readRDS(paste0(data_path, "ST_Discovery_so.rds"))
```

The input data frame for SCOPE should contain the following columns:

- `CellID`: cell id
- `ImageID`: sample id
- `X`: X coordinate
- `Y`: Y coordinate
- `Celltype`: cell type labels

Some cell types in this dataset are labeled with a \_LC suffix,
indicating “Low Confidence”. In this analysis, we will merge these with
their high-confidence counterparts. Additionally, we will combine all
the CD8\\^-\\ T cell population excluding the regulatory T cells (i.e.,
the double-negative T cells `CD4.T.cell.DN` and the CD4\\^+\\ T cells
`CD4.T.cell`) into one category.

``` r

ST_Discovery_so$cell.subtypes_combined <- str_remove(ST_Discovery_so$cell.subtypes, "_LC$")
ST_Discovery_so$cell.subtypes_combined <- str_remove(ST_Discovery_so$cell.subtypes_combined, ".DN$")

ST_Discovery_so$CellID <- Cells(ST_Discovery_so)

HGSOC_CosMx_dat_in = ST_Discovery_so@meta.data[
  ,c("samples","CellID","cell.subtypes_combined","x","y")]
colnames(HGSOC_CosMx_dat_in) <- c("ImageID","CellID","Celltype","X","Y")

HGSOC_CosMx_dat_in$Celltype <- factor(HGSOC_CosMx_dat_in$Celltype, names(HGSOC_CosMx_palette))
```

Since all 100 tissue sections exceed 500 cells, no samples are excluded
from subsequent downstream analyses.

``` r

ncell_by_slide <- table(HGSOC_CosMx_dat_in$ImageID)
all(ncell_by_slide > 500)

slide_list <- unique(HGSOC_CosMx_dat_in$ImageID) # 100 images left
save(HGSOC_CosMx_dat_in, file = paste0(data_path, dataset_name, "_dat_in.RData"))
```

## 2. Create neighborhood cell proportion matrix

Next, we define the spatial neighborhood using use the k-nearest
neighbors and compute the neighborhood cell type composition matrix.
Here the number of nearest neighbors \\k\\ is set at the default value
20.

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
      dat_slide=HGSOC_CosMx_dat_in[HGSOC_CosMx_dat_in$ImageID==slide_list[i], ], 
      knn_num = 20, min.radius = 100, max.radius = 500, min.n.neighbors = 10,
      cell_types_available = levels(HGSOC_CosMx_dat_in$Celltype))},
  mc.cores = 5)

## cell type proportion within the neighborhood
cell_neighbor_table_maxdist = lapply(
  cell_neighbor_results_maxdist, "[[", "prop.df") %>% 
  do.call(rbind, .)
dist.quantiles = lapply(
  cell_neighbor_results_maxdist, "[[", "dist.quantiles") %>% 
  do.call(rbind, .)

save(cell_neighbor_table_maxdist, file = paste0(
  data_path, dataset_name, "_cell_neighbor_table_maxdist.RData"))
save(dist.quantiles, file = paste0(
  data_path, dataset_name, "_dist.quantiles.RData"))
```

We can examine the distribution of 99th percentile of all distance to
the nearest neighbors to select the threshold parameters min.radius and
max.radius accordingly.

``` r

dist.quantiles %>%
  as.data.frame() %>% 
  ggplot(aes(x = `99%`)) +
  geom_histogram(color = "white", fill = "skyblue", bins = 20) +
  geom_vline(xintercept = 500, linetype = 2, color = "black") +
  labs(y = "99th percentile of Dist.\nto Nearest Neighbors", x = NULL) +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45))
```

![](HGSOC_CosMx_files/figure-html/Examine%20distance%20distribution-1.png)

Cells that retain fewer than `min.n.neighbors` neighbors within the
distance boundary will have their neighborhood cell-type proportions set
to `NaN` and will be excluded from subsequent niche clustering. Here we
set `min.n.neighbors` to 10.

``` r

cells_to_exclude = apply(cell_neighbor_table_maxdist, 1, function(x) {any(is.na(x))}) %>%
  as.data.frame() %>%
  `colnames<-`("is_na") %>%
  rownames_to_column("CellID") %>%
  left_join(HGSOC_CosMx_dat_in[, c("CellID", "X", "Y", "ImageID")], by = "CellID")

# Number of cells to exclude
table(cells_to_exclude$is_na)
#> 
#>  FALSE   TRUE 
#> 489257   2535

cells_to_exclude %>%
  group_by(ImageID) %>%
  mutate(n_NaN = sum(is_na)) %>%
  filter(n_NaN > 10) %>% ungroup %>%
  ggplot(aes(x = X, y = Y, color = is_na)) +
  geom_point(shape = 16, size= 0.1) +
  scale_color_manual(values = c("brown2", "grey90"), 
                     breaks = c(TRUE, FALSE), name = "Cells To\nExclude") +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  facet_wrap(~str_remove(ImageID, "SMI_"), ncol = 15) +
  theme_void() +
  coord_fixed(expand = FALSE) +
  theme(legend.key.size = unit(2.5, "mm"),
        panel.background = element_blank(),
        strip.clip = "off")
```

![](HGSOC_CosMx_files/figure-html/Examine%20loose%20cells-1.png)

## 3. Cell type triad niche (CTN) detection using SCOPE

The HGSOC CosMx dataset has Malignant, Fibroblast, Monocyte, CD8.T.cell,
B.cell, Endothelial, CD4.T.cell, NK.cell, Treg, T.cell.DP, T.B.cell,
Mast.cell available cell types, which can form only 220 possible
cell-type triads. Given the relatively small number of combinations, we
will skip the colocation quotient (CLQ) permutation test step and
instead focus solely on filtering out cell types that are less frequent
within local neighborhoods.

Specifically, we set the minimum neighborhood cell type proportion
threshold at 0.15. The double-positive T cells (`T.cell.DP`) and mast
cells (`Mast.cell`) have fewer than 100 cells that pass this threshold.
We exclude the 10 triads that contain both of these cell types, leaving
210 cell-type triads to screen.

``` r

celltype_to_exclude <- apply(cell_neighbor_table_maxdist, 2, 
                             function(x) {sum(x > 0.15, na.rm = TRUE)})
celltype_to_exclude <- names(celltype_to_exclude)[celltype_to_exclude < 100]
celltype_to_exclude
#> [1] "T.cell.DP" "Mast.cell"

combination_table_all <- as.data.frame(t(combn(levels(HGSOC_CosMx_dat_in$Celltype), 3)))
has_all_targets <- rowSums(do.call(cbind, lapply(celltype_to_exclude, function(x) {
    rowSums(combination_table_all == x) > 0}))) == length(celltype_to_exclude)

combination_table_filtered <- combination_table_all[!has_all_targets, ]
paste("Number of cell-type triads to screen:", nrow(combination_table_filtered))
#> [1] "Number of cell-type triads to screen: 210"
```

Next, we apply SCOPE for the systematic screening of CTNs. We set the
the minimum neighborhood cell-type proportion threshold `min.prop` at
0.15 and utilize the default configurations for the rest of the
parameters. For this analysis, we set the number of clusters
(`clusternum`) to 2. By setting `save_results = TRUE`, the clustering
results for each CTN are automatically saved as individual `.RData`
files in the specified `output_dir.`

``` r

run_SCOPE(
  combination_table = combination_table_filtered[1:2, ],
  cell_neighbor_table = cell_neighbor_table_maxdist,
  dat_in = HGSOC_CosMx_dat_in,
  save_results = TRUE,
  output_dir = paste0(data_path, "CTN_K_2"),
  clusternum = 2, min.prop = 0.15,
  linear = FALSE, mgcv_df = 15, iter.max = 100)
```

We compile the individual outputs and merge them into a single data
frame. The first five columns match the original input data, and the
subsequent columns contain the assigned CTN labels.

``` r

Full_CTN_table_mgcv <- list.files(
  paste0(data_path, "CTN_K_2"), ".RData", full.names = TRUE) %>%
  pbmclapply(function(dir) {
    load(dir)
    return(CTN_table) 
  }, mc.cores = 5) %>%
  do.call(rbind, .)  %>%
  select(-label_new) %>%
  `rownames<-`(NULL) %>%
  pivot_wider(names_from = "CTN", values_from = "label") %>%
  left_join(HGSOC_CosMx_dat_in, by = "CellID") %>%
  select(all_of(colnames(HGSOC_CosMx_dat_in)), everything())

save(Full_CTN_table_mgcv, file = paste0(data_path, dataset_name, "_Full_CTN_table_mgcv.RData"))
```

## 4. Merging highly overlapped CTNs

Sometimes CTNs driven by two different triads can be highly similar. To
merge highly overlapping CTNs, we calculate Jaccard distance between
every pair of CTNs and perform hierarchical clustering to reduce
redundancy.

``` r

load(paste0(data_path, dataset_name, "_Full_CTN_table_mgcv.RData"))

# All cell type triads
CTN_list <- setdiff(colnames(Full_CTN_table_mgcv), colnames(HGSOC_CosMx_dat_in))

CTN_jaccard_mat_hclust <- jaccard_dist_hclust(
  Full_CTN_table = Full_CTN_table_mgcv, 
  CTN_ls = CTN_list, num_cores = 5)

save(CTN_jaccard_mat_hclust, file = paste0(
     data_path, dataset_name, "_CTN_jaccard_mat_hclust.RData"))
```

Silhouette scores are utilized to determine the optimal number of
clusters. The final cluster count is set to 65, balancing both the mean
and median Silhouette scores. As shown in the figure below, while 65 is
not the optimal choice, it provides a trade-off to keep the total number
of clusters manageable while still preserving as much cluster separation
as possible.

``` r

n_clusters = 65 # Final number of hierarchical clusters
draw_silhoutte_score(CTN_jaccard_mat_hclust$hclust_silhouette,
                     num_clusters = n_clusters)
```

![](HGSOC_CosMx_files/figure-html/Silouhette%20score-1.png)

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

![](HGSOC_CosMx_files/figure-html/Jaccard%20similarity%20matrix%20heatmap-1.png)

``` r

# Hierarchical clustering results
htmap_row_order <- unlist(row_order(htmap_jaccard_index))
CTN_merged_label <- as.data.frame(cutree(CTN_jaccard_mat_hclust$cl, k = n_clusters)) %>%
  rownames_to_column("CTN") %>%
  `colnames<-`(c("CTN", "cluster")) %>%
  # keep the CTNs in the same order as the heatmap
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
[`CTN_annotate()`](https://github.com/YufanWu147/SCOPE.CTN/reference/CTN_annotate.md)
function to rename consolidated CTNs based on the most frequent core
cell types within their original triads. If a tie occurs between cell
types, the function automatically resolves it by identifying the cell
type with the highest overall proportion within the consolidated CTN.

First, we compute the cell type proportion within the newly consolidated
CTNs, which will be used as one of the inputs to the
[`CTN_annotate()`](https://github.com/YufanWu147/SCOPE.CTN/reference/CTN_annotate.md)
function.

``` r

CTN_merged_celltype_proportion <- CTN_merged_label %>%
  split(f = .$cluster) %>%
  lapply(function(CTN_merged_label_sub) {
  Full_CTN_table_mgcv[, c("Celltype", CTN_merged_label_sub$CTN)] %>% 
      filter(if_any(everything(), ~ .x == "CTN")) %>%
      count(Celltype) }) %>%
  rbindlist(idcol = "cluster") %>%
  group_by(cluster) %>%
  mutate(prop = n/sum(n)) %>% select(-n)


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
  celltype_levs = names(HGSOC_CosMx_palette))
```

If the top two most frequent cell types only appear twice in the
original CTN names, the function flags the consolidated CTN with an
asterisk (\*). Duplicated new CTN names are marked with a plus sign (+).
We suggest double-checking the flagged CTNs and manually renaming them
if necessary.

Here we also manually adjust a few merged CTNs where the constituent
original CTNs are formed by a combination of four cell types.

``` r

head(CTN_merged_label_annotated$core_celltypes_summary, n = 15)

CTN_annotation_detailed <- CTN_merged_label_annotated$annotation %>%
    mutate(annotation_final = case_when(
    cluster == 5 ~ "B.cell_NK.cell_Malignant_Endothelial",
    cluster == 6 ~ "B.cell_T.cell_NK.cell_Malignant",
    cluster == 9 ~ "B.cell_T.cell_Mast.cell_Malignant",
    cluster == 35 ~ "B.cell_T.cell_NK.cell_Monocyte",
    cluster == 46 ~ "T.cell_NK.cell_Mast.cell_Fibroblast",
    TRUE ~ annotation
  ))

save(CTN_annotation_detailed, file = paste0(
     data_path, dataset_name, "_CTN_annotation_detailed.RData"))
```

Given a newly consolidated CTN, any cell belonging to at least one of
its constituent CTNs will be labeled as “CTN”, while the remaining cells
will be labeled as “Unassigned”.

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
  inner_join(HGSOC_CosMx_dat_in, by = "CellID") %>%
  select(all_of(colnames(HGSOC_CosMx_dat_in)), everything())

save(Full_CTN_table_mgcv_merged, 
     file = paste0(data_path, dataset_name, "_Full_CTN_table_mgcv_merged.RData"))
```

We can further perform hierarchical clustering on the Jaccard distance
between the consolidated CTNs and grouping them into major categories.

``` r

load(paste0(data_path, dataset_name, "_Full_CTN_table_mgcv_merged.RData"))
CTN_merged_list <- setdiff(colnames(Full_CTN_table_mgcv_merged),
                           colnames(HGSOC_CosMx_dat_in))
CTN_merged_jaccard_mat_hclust <- jaccard_dist_hclust(
  Full_CTN_table = Full_CTN_table_mgcv_merged, 
  CTN_ls = CTN_merged_list, num_cores = 5)
save(CTN_merged_jaccard_mat_hclust, file = paste0(
  data_path, dataset_name, "_CTN_merged_jaccard_mat_hclust.RData"))
```

These categories can be annotated by examining the cell type composition
of the constituent CTNs.

``` r

CTN_annotation_merged <- CTN_annotation_detailed %>%
  # Grouping CTNs based on cell type composition  
  mutate(CTN_group = case_when(
    cluster %in% 1:11 ~ "Malignant",
    cluster %in% 12:18 ~ "Immune/Endothelial",
    cluster == 19 ~ "Bcell/Endothelial",
    cluster %in% 20:25 ~ "Bcell-rich",
    cluster %in% 26:28 ~ "Monocyte/Fibroblast",
    cluster %in% 29:37 ~ "Immune/Malignant",
    cluster %in% 38:39 ~ "Mixed",
    cluster %in% 40:52 ~ "Immune/Fibroblast",
    cluster %in% 53:65 ~ "Immune")) %>%
  distinct(cluster, annotation_final, CTN_group)%>%
  rename("annotation" = "annotation_final") 
  

CTN_merged_celltype_proportion_w_annot <- CTN_merged_celltype_proportion %>%
  left_join(CTN_annotation_merged[, 1:2], by = "cluster") 

CTN_merged_dend <- draw_CTN_dendro(
  cl = CTN_merged_jaccard_mat_hclust$cl, h = 0.7,
  annot_var = "annotation",
  rename_celltype = HGSOC_CosMx_rename_celltype, 
  CTN_annotation = CTN_annotation_merged,
  CTN_group_var = "CTN_group",
  CTN_group_palette = HGSOC_CosMx_CTN_group_palette,
  celltype_palette = c(HGSOC_CosMx_palette, "T.cell" = "dodgerblue3"),
  dot_fontsize = 11, y_lim = c(NA, 1),
  x_expand = c(0.003, 0.003),
  show_text = FALSE)


p_CTN_merged_celltype_prop <- draw_CTN_celltype_prop(
  CTN_merged_celltype_prop = CTN_merged_celltype_proportion_w_annot, 
  CTN_dend = CTN_merged_dend, 
  celltype_palette = HGSOC_CosMx_palette,
  rename_celltype = HGSOC_CosMx_rename_celltype,
  facet_var = "CTN_group_label", 
  x_expand = c(0.05, 0.05),
  radius_range = c(1, 6))

p_CTN_merged_celltype_prop + CTN_merged_dend$p + 
  plot_layout(guides = "collect", widths = c(1, 0.4))
```

![](HGSOC_CosMx_files/figure-html/CTN%20dendrogram-1.png)

We can also visualize similarities in the cell type compositions of the
CTNs using Uniform Manifold Approximation and Projection (UMAP).

``` r

draw_CTN_umap(Full_CTN_table = Full_CTN_table_mgcv_merged,
              CTN_merged_celltype_prop = CTN_merged_celltype_proportion,
              CTN_ls = CTN_merged_list, CTN_group_var = "CTN_group",
              CTN_annotation = CTN_annotation_merged,
              CTN_group_palette = HGSOC_CosMx_CTN_group_palette,
              umap_n_neighbors = 10, umap_seed = 42) +
  guides(fill = guide_legend(ncol = 2, override.aes = list(size = 3)),
         radius = guide_legend(ncol = 1)) +
  theme(legend.position = "bottom",
        legend.key.spacing.y = unit(0.2, "mm"),
        legend.direction = "vertical",
        legend.title.position = "top")
```

![](HGSOC_CosMx_files/figure-html/CTN%20cell%20type%20composition%20umap-1.png)

## 5. Identify clinically-relevant CTNs

After obtaining the final CTN labels, we will next examine their
prognostic value using the detailed patient clinical data provided by
the original data source. One patient with unknown cause of death will
be excluded from survival analysis.

``` r

# Image-level clinical data
clinical_data <- readxl::read_xlsx(
  paste0(data_path, "41590_2024_1943_MOESM3_ESM.xlsx"),
  sheet = "Table 2b", skip = 1) %>%
  filter(dataset == "Discovery") %>%
  filter(outcome != "Dead (other)") %>%
  mutate(OS = as.numeric(fu_time1),               # Follow-up time post diagnosis
         Ev.O = case_when(outcome == "Alive" ~ 0, # Outcome (0 = alive, 1 = dead)
                          str_detect(outcome, "[Dd]ead") ~ 1,
                          is.na(outcome) ~ NA)) %>%
  mutate(age = as.numeric(age)) %>%
  mutate(ICB = ifelse(str_to_lower(immunotherapy) %in% c(
    "not reported (lost to follow-up)", 
    "on atezolizumab  vs placebo, stable",
    "no", "na"), "No", "Yes")) %>%
  select(ImageID = profile, Patient_ID = patients, sites_binary,
         OS, Ev.O, treatment, age, stage, ICB)

head(clinical_data)
#> # A tibble: 6 × 9
#>   ImageID      Patient_ID sites_binary    OS  Ev.O treatment   age stage ICB  
#>   <chr>        <chr>      <chr>        <dbl> <dbl> <chr>     <dbl> <chr> <chr>
#> 1 SMI_T10_F001 HGSC1      Adnexa         910     0 Untreated    58 III   No   
#> 2 SMI_T10_F002 HGSC7      Adnexa        1153     0 Untreated    58 III   No   
#> 3 SMI_T10_F003 HGSC8      Adnexa        1056     0 Treated      65 IV    Yes  
#> 4 SMI_T10_F004 HGSC13     Omentum       1005     1 Treated      67 III   No   
#> 5 SMI_T10_F005 HGSC15     Adnexa        1043     0 Treated      33 IV    No   
#> 6 SMI_T10_F006 HGSC2      Adnexa        1492     0 Untreated    36 IV    No
```

To test whether the CTNs predict patient survival, we fit a Cox
proportional hazards regression model using patient-level CTN presence
as the primary predictor. For a given patient, a CTN is considered
present if it contains 100 or more cells within a single tissue in at
least one tissue section.

We conduct both a univariate and a multivariate Cox regression analysis
adjusting for age at diagnosis, tumor stage, and treatment history.

``` r

# Overall survival, univariate Cox regression
coxph_res_OS <- CTN_presence_coxph(
    Full_CTN_table = Full_CTN_table_mgcv_merged, 
    CTN_ls = CTN_merged_list, 
    clinical_df = clinical_data,
    time = "OS", event = "Ev.O",
    patient_id = "Patient_ID", confounder = NULL,
    min.n.patients = 3, min.n.cells = 100)

# Overall survival, multivariate Cox regression
coxph_res_OS_adjusted <- CTN_presence_coxph(
    Full_CTN_table = Full_CTN_table_mgcv_merged, 
    CTN_ls = CTN_merged_list, 
    clinical_df = clinical_data,
    time = "OS", event = "Ev.O",
    patient_id = "Patient_ID", 
    confounder = c("treatment", "age", "stage", "ICB"),
    min.n.patients = 3, min.n.cells = 100)
```

The following plots show the CTNs associated with the overall survival
of HGSOC patients.

``` r

# Overall survival, univariate
draw_CTN_coxph_logHR(
  coxph_res = coxph_res_OS,
  CTN_merged_celltype_prop = CTN_merged_celltype_proportion_w_annot,
  alpha = 1, top = 10, 
  celltype_levs = names(HGSOC_CosMx_palette), dot_fontsize = 14,
  celltype_palette = c(HGSOC_CosMx_palette, "T" = "dodgerblue3"),
  rename_celltype = HGSOC_CosMx_rename_celltype,
  line_width = 0.4, y_text_size = NULL,
  width_ratio = c(0.45, 0.55)) 
```

![](HGSOC_CosMx_files/figure-html/Univariate%20cox%20regression%20on%20OS-1.png)

``` r

# Overall survival, multivariate
draw_CTN_coxph_logHR(
  coxph_res = coxph_res_OS_adjusted,
  CTN_merged_celltype_prop = CTN_merged_celltype_proportion_w_annot,
  alpha = 1, top = 10, 
  celltype_levs = names(HGSOC_CosMx_palette), dot_fontsize = 14,
  celltype_palette = c(HGSOC_CosMx_palette, "T" = "dodgerblue3"),
  rename_celltype = HGSOC_CosMx_rename_celltype,
  line_width = 0.4, y_text_size = NULL,
  width_ratio = c(0.45, 0.55)) 
```

![](HGSOC_CosMx_files/figure-html/Multivariate%20cox%20regression%20on%20OS-1.png)

Below is the lymphocyte aggregate-like B cell/CD4\\^+\\ T cell/CD8\\^+\\
T cell CTN that is associated with the overall survival.

``` r

p_CTN_celltype <- draw_CTN_celltype(
  Full_CTN_table = Full_CTN_table_mgcv_merged,
  CTN = "B.cell_CD4.T.cell_CD8.T.cell", 
  img_list = c("SMI_T14_F011", "SMI_T12_F001", "SMI_T11_F005", "SMI_T14_F006"),
  celltype_palette = HGSOC_CosMx_palette, 
  rename_celltype = HGSOC_CosMx_rename_celltype,
  scales = "fixed", legend.ncol = 1)
p_CTN_celltype$p_celltype / p_CTN_celltype$p_CTN +
  plot_layout(guides = "collect")
```

![](HGSOC_CosMx_files/figure-html/Visualize%20Bcell_CD4_CD8-1.png)

Below shows the CD8\\^+\\T cell/Fibroblast/Malignant CTN which separate
the stromal and malignant compartments of the tissue.

``` r

p_CTN_celltype <- draw_CTN_celltype(
  Full_CTN_table = Full_CTN_table_mgcv_merged,
  CTN = "CD8.T.cell_Fibroblast_Malignant", 
  img_list = c('SMI_T10_F001','SMI_T12_F001',
               'SMI_T12_F009','SMI_T12_F016'),
  celltype_palette = HGSOC_CosMx_palette, 
  rename_celltype = HGSOC_CosMx_rename_celltype,
  scales = "fixed", legend.ncol = 1)
p_CTN_celltype$p_celltype / p_CTN_celltype$p_CTN +
  plot_layout(guides = "collect")
```

![](HGSOC_CosMx_files/figure-html/Visualize%20CD8T_Fibroblast_Malignant-1.png)

## 6. Subclustering CTNs

Once the clincially-associated CTNs are identified, we can refit SCOPE
to generate subclusters and obtain finer-grained CTN substructures.

**Note**: To use the
[`run_SCOPE_subcluster()`](https://github.com/YufanWu147/SCOPE.CTN/reference/run_SCOPE_subcluster.md)
function, the `save_results` parameter must be set to `TRUE` when
running
[`run_SCOPE()`](https://github.com/YufanWu147/SCOPE.CTN/reference/run_SCOPE.md).
Additionally, please ensure that both functions share the same
`output_dir` path.

``` r

Full_CTN_table_mgcv_merged_subclustered <- run_SCOPE_subcluster(
    Full_CTN_table = Full_CTN_table_mgcv,  # Full CTN table before merging
    cell_neighbor_table = cell_neighbor_table_maxdist,
    CTN_annotation = CTN_annotation_detailed,
    annot_var = "annotation_final",
    subclusternum = 2, 
    verbose = FALSE,
    output_dir = paste0(data_path, "CTN_K_2"))

save(Full_CTN_table_mgcv_merged_subclustered, file = paste0(
       data_path, dataset_name, "_Full_CTN_table_mgcv_merged_subclustered.RData"))
```

Subclusters CTN1 and CTN2 of the he B cell/CD4\\^+\\ T cell/CD8\\^+\\ T
cell CTN resemble B cell and T cell aggregates, respectively.

``` r

load(paste0(data_path, dataset_name, "_Full_CTN_table_mgcv_merged_subclustered.RData"))
p_CTN_celltype <- draw_CTN_celltype(
  Full_CTN_table = Full_CTN_table_mgcv_merged_subclustered$Full_CTN_table_reclustered,
  CTN = "B.cell_CD4.T.cell_CD8.T.cell", 
  img_list = c("SMI_T14_F011", "SMI_T12_F001", "SMI_T11_F005", "SMI_T14_F006"),
  celltype_palette = HGSOC_CosMx_palette, 
  CTN_palette = c("CTN1" = "#E69F00", "CTN2" = "#56B4E9", "Unassigned" = "grey90"),
  rename_celltype = HGSOC_CosMx_rename_celltype,
  scales = "fixed", legend.ncol = 1)
p_CTN_celltype$p_celltype / p_CTN_celltype$p_CTN +
  plot_layout(guides = "collect")
```

![](HGSOC_CosMx_files/figure-html/Visualize%20Bcell_CD4_CD8%20subclusters-1.png)

The CD8\\^+\\ T cell/Fibroblast/Malignant CTN subclusters correspond to
the tumor core (“CTN1”), the tumor-stromal boundary (“CTN2”), and the
stromal compartment (“Unassigned”).

``` r

p_CTN_celltype <- draw_CTN_celltype(
  Full_CTN_table = Full_CTN_table_mgcv_merged_subclustered$Full_CTN_table_reclustered,
  CTN = "CD8.T.cell_Fibroblast_Malignant", 
  img_list = c('SMI_T10_F001','SMI_T12_F001',
               'SMI_T12_F009','SMI_T12_F016'),
  celltype_palette = HGSOC_CosMx_palette, 
  CTN_palette = c("CTN1" = "#E69F00", "CTN2" = "#56B4E9", "Unassigned" = "grey90"),
  rename_celltype = HGSOC_CosMx_rename_celltype,
  scales = "fixed", legend.ncol = 1)
p_CTN_celltype$p_celltype / p_CTN_celltype$p_CTN +
  plot_layout(guides = "collect")
```

![](HGSOC_CosMx_files/figure-html/Visualize%20CD8T_Fibroblast_Malignant%20subclusters-1.png)

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
#>  [1] readxl_1.4.5        rlang_1.1.7         magrittr_2.0.4     
#>  [4] clue_0.3-66         GetoptLong_1.1.0    otel_0.2.0         
#>  [7] matrixStats_1.5.0   compiler_4.4.1      mgcv_1.9-3         
#> [10] png_0.1-9           systemfonts_1.3.1   pbmcapply_1.5.1    
#> [13] vctrs_0.7.2         pkgconfig_2.0.3     shape_1.4.6.1      
#> [16] crayon_1.5.3        fastmap_1.2.0       labeling_0.4.3     
#> [19] utf8_1.2.6          rmarkdown_2.30      markdown_2.0       
#> [22] tzdb_0.5.0          ragg_1.5.1          xfun_0.57          
#> [25] cachem_1.1.0        litedown_0.9        jsonlite_2.0.0     
#> [28] gmp_0.7-5.1         irlba_2.3.7         parallel_4.4.1     
#> [31] cluster_2.1.6       R6_2.6.1            bslib_0.10.0       
#> [34] stringi_1.8.7       RColorBrewer_1.1-3  ClusterR_1.3.3     
#> [37] jquerylib_0.1.4     cellranger_1.1.0    Rcpp_1.1.1         
#> [40] assertthat_0.2.1    iterators_1.0.14    knitr_1.51         
#> [43] IRanges_2.40.1      FNN_1.1.4.1         Matrix_1.7-0       
#> [46] splines_4.4.1       timechange_0.4.0    tidyselect_1.2.1   
#> [49] rstudioapi_0.18.0   yaml_2.3.12         doParallel_1.0.17  
#> [52] codetools_0.2-20    lattice_0.22-6      withr_3.0.2        
#> [55] S7_0.2.1            evaluate_1.0.5      desc_1.4.3         
#> [58] survival_3.6-4      xml2_1.5.2          circlize_0.4.17    
#> [61] pillar_1.11.1       foreach_1.5.2       stats4_4.4.1       
#> [64] generics_0.1.4      hms_1.1.4           S4Vectors_0.44.0   
#> [67] commonmark_2.0.0    scales_1.4.0        glue_1.8.0         
#> [70] tools_4.4.1         fs_1.6.7            colorspace_2.1-2   
#> [73] nlme_3.1-164        cli_3.6.5           textshaping_1.0.5  
#> [76] uwot_0.2.4          gtable_0.3.6        sass_0.4.10        
#> [79] digest_0.6.39       BiocGenerics_0.52.0 ggrepel_0.9.8      
#> [82] rjson_0.2.23        htmlwidgets_1.6.4   farver_2.1.2       
#> [85] htmltools_0.5.9     pkgdown_2.2.0       lifecycle_1.0.5    
#> [88] GlobalOptions_0.1.3 gridtext_0.1.5      MASS_7.3-60.2
```
