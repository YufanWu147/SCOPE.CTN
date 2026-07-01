#' Separating cell types in CTN labels with colored dots
#'
#' Adding colored dots before each cell type in cell type triad (CTN) labels. Outputs can be used in titles and labels \code{ggplot} with markdown text enabled using [ggtext::element_markdown].
#'
#' @param CTN CTN label; e.g. B_CD4T_CD8T.
#' @param celltype_palette A named vector with colors as values and cell type labels as names. Must contain all cell types.
#' @param dot_fontsize Size of the dots. Default is 16pt.
#' @param sep A character string separating the cell types in CTN. Note that it can't appear in cell type labels.
#' @import stringr ggtext
#' @export
#' @examples
#' library(ggplot2)
#' library(ggtext)
#'
#' CTN_colored <- color_CTN_names(
#'   CTN = "B_CD4T_CD8T",
#'   celltype_palette = c("B" = "#fee440", "CD4T" = "#3a86ff", "CD8T" = "#8338ec"))
#'
#' # ggplot() +
#' # labs(title = CTN_colored) +
#' #  theme(plot.title = element_markdown(hjust = 0.5))
color_CTN_names <- function(CTN, celltype_palette = c("B" =  "#fee440", "CD4T" = "#3a86ff", "CD8T" = "#8338ec"),
                            dot_fontsize = 16, sep = "_") {
  core_celltypes <- str_split(CTN, pattern = sep)
  CTN_colored <- lapply(core_celltypes, function(celltype) {
    colors = celltype_palette[celltype]
    CTN_new <- paste0("<span style='color:", colors, ";font-size:",
                      dot_fontsize, "pt'>\U25CF</span>", celltype)
    CTN_new <- paste(CTN_new, collapse = " ")
    return(CTN_new)
  }) %>% unlist()

  return(CTN_colored)
}

#' Visualize silhouette score
#'
#' Evaluating hierarchical clustering performance for selecting optimal cluster number
#' @param hclust_silhouette A data frame containing the average (\code{avg_si}) and median Silhouette score (\code{median_si})
#' under cluster number \code{K}. Output of the [jaccard_dist_hclust()] function.
#' @param num_clusters Optimal number of clusters.
#' @import ggplot2
#' @importFrom cowplot theme_cowplot
#' @export
draw_silhoutte_score <- function(hclust_silhouette, num_clusters) {
  hclust_silhouette %>%
    pivot_longer(1:2, names_to = "metric", values_to = "si") %>%
    ggplot(aes(x = K, y = si, color = metric)) +
    geom_point(size = 0.5) +
    geom_line() +
    geom_vline(xintercept = num_clusters,
               linetype = 2, color = "black") +
    scale_color_brewer(palette = "Set1") +
    scale_x_continuous(breaks = c(seq(0, 200, 50), num_clusters)) +
    scale_color_discrete(labels = c("Average Silhouette Score",
                                    "Median Silhouette Score"),
                         name = NULL) + labs(y = NULL) +
    theme_cowplot() +
    theme(legend.position = "bottom",
          legend.direction = "vertical")
}


#' Jaccard similarity index heatmap
#'
#' Heatmap visualization of the Jaccard similarity indices between all pairs of cell-type triad
#' Niches (CTNs) using an ordered heatmap based on hierarchical clustering.
#' @param jaccard_mat A numeric matrix containing the Jaccard similarity indices between all pairs of CTNs. Output of the [jaccard_dist_hclust()] function.
#' @param cl An object of class "hclust" which describes the tree produced by the clustering process. Output of the [jaccard_dist_hclust()] function.
#' @param num_clusters An integer specifying the target number of clusters to slice or highlight on the heatmap.
#' @return A [ComplexHeatmap::Heatmap-class] object.
#' @import ComplexHeatmap grid
#' @export
draw_htmap_jaccard_index <- function(jaccard_mat, cl, num_clusters) {
  p <- Heatmap(
    jaccard_mat,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 7),
    name = "Jaccard\nSimilarity",
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    cluster_rows = cl,
    cluster_columns = cl,
    row_split = num_clusters,
    column_split = num_clusters,
    cluster_row_slices = FALSE,
    cluster_column_slices = FALSE,
    row_gap = unit(0, "mm"),
    column_gap = unit(0, "mm"),
    row_title = NULL, column_title = NULL,
    cell_fun = function(j, i, x, y, width, height, fill) {
      if(jaccard_mat[i, j] > 0.95) {
        grid.text("*", x, y, gp = gpar(fontsize = 8), vjust = 0.85)
      }
      if(between(jaccard_mat[i, j], 0.9, 0.95)) {
        grid.text(".", x, y, gp = gpar(fontsize = 8), vjust = 0)
      }
    }
  )

  p <- draw(p)

  for(i in 1:num_clusters) {
    decorate_heatmap_body("Jaccard\nSimilarity", row_slice = i, column_slice = i, {
      grid.rect(unit(0, "npc"), unit(1, "npc"), # top left
                gp = gpar(fill = "transparent", lwd = 1.5),
                just = c("left", "top")
      )}) }

  return(p)
}

#' Prepare input data CTN dendrogram
#'
#' Internal function for preparing input data to [draw_CTN_dendro()] function
#' @param cl An object of class "hclust" which describes the tree produced by the clustering process. Output of the [jaccard_dist_hclust()] function.
#' @param k An integer scalar or vector with the desired number of groups. See [stats::cutree()] for details.
#' @param h Numeric scalar or vector with heights where the tree should be cut. [stats::cutree()] for details.
#' At least one of \code{k} or \code{h} must be specified, \code{k} overrides \code{h} if both are given.
#' @return A list with components (See [ggdendro::dendro_data()] for details):
#' \describe{
#'   \item{segments}{Line segment data}
#'   \item{labels}{Label data}
#'   }
#' @import ggdendro
#' @importFrom stats cutree
#' @keywords internal
#' @references
#' Atrebas. Dendrograms in R, a lightweight approach. \url{https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/}
prepare_dendro_data <- function(cl, h, k = NULL) {

  dend_data    <-  ggdendro::dendro_data(cl, type = "rectangle")
  seg       <-  dend_data$segments
  if(is.null(k)) {
    labclust  <-  cutree(cl, h = h)[cl$order]
    k         <-  length(unique(unname(labclust)))
  }
  else {
    labclust  <-  cutree(cl, k = k)[cl$order]
  }

  heights   <-  sort(cl$height, decreasing = TRUE)
  height    <-  mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)
  segclust  <-  rep(0L, nrow(seg))

  for (i in 1:k) {
    xi      <-  dend_data$labels$x[labclust == i]
    idx1    <-  seg$x    >= min(xi) & seg$x    <= max(xi)
    idx2    <-  seg$xend >= min(xi) & seg$xend <= max(xi)
    idx3    <-  seg$yend < height
    idx     <-  idx1 & idx2 & idx3
    segclust[idx] <- i
  }

  idx                    <-  which(segclust == 0L)
  segclust[idx]          <-  segclust[idx + 1L]
  dend_data$segments$clust  <-  segclust
  dend_data$segments$line   <-  as.integer(segclust < 1L)
  dend_data$labels$clust    <-  labclust

  return(dend_data)
}

#' Draw CTN dendrogram
#' Dendrogram visualization for the hierarchical clustering results on the pairwise Jaccard distance between CTNs.
#' @param cl An object of class "hclust" which describes the tree produced by the clustering process. Output of the [jaccard_dist_hclust()] function.
#' @param k An integer scalar or vector with the desired number of groups. See [stats::cutree()] for details.
#' @param h Numeric scalar or vector with heights where the tree should be cut. [stats::cutree()] for details. At least one of \code{k} or \code{h} must be specified, \code{k} overrides \code{h} if both are given.
#' @param CTN_annotation (Optional) A data frame with the original CTN names (\code{CTN}), the corresponding cluster label (\code{cluster}), and the new name for the consolidated CTN (\code{annotation_final}) after merging.
#' @param celltype_palette Input of the [color_CTN_names()] function. A named vector with colors as values and cell type labels as names. Must contain all cell types.
#' @param dot_fontsize Input of the [color_CTN_names()] function. Size of the dots. Default is 16pt.
#' @param CTN_group_palette (Optional) A named vector with colors as values and CTN categories as names. Only used when \code{CTN_annotation} is provided.
#' @param show_text Whether to display CTN names alongside the dendrogram (default: \code{TRUE}).
#' @param y_lim A numeric vector of length two providing limits of the scale. Use NA to refer to the existing minimum or maximum.
#' @param x_expand A numeric vector of length two of the multiplicative range expansion factors for the x axis. See [ggplot2::expansion()] for details.
#' @param rename_celltype A function for renaming cell type labels.
#' @return A list with components:
#' \describe{
#'   \item{p}{Dendrogram}
#'   \item{label_df}{Label data. See [ggdendro::dendro_data()] for details.}
#'   }
#' @import stringr dplyr ggplot2 ggdendro ggtext
#' @export
draw_CTN_dendro <- function(cl, h, k = NULL, CTN_annotation = NULL,
                            celltype_palette = NULL, dot_fontsize = 12,
                            CTN_group_palette = NULL, show_text = TRUE,
                            x_expand = c(0.001, 0.001), y_lim = c(-2, 1),
                            rename_celltype = NULL) {
  dend_data <- prepare_dendro_data(cl = cl, h = h, k = k)
  label_df <- dend_data$labels %>%
    mutate(order = row_number()) %>%
    mutate(CTN_colored = color_CTN_names(
      label, celltype_palette = celltype_palette,
      sep = "_", dot_fontsize = dot_fontsize)) %>%
    mutate(clust = factor(clust, unique(.$clust)))
  if(!is.null(rename_celltype)) {
    label_df$CTN_colored <- rename_celltype(label_df$CTN_colored)
  }
  label_df$CTN_colored <- factor(label_df$CTN_colored, unique(label_df$CTN_colored))
  if(is.null(CTN_annotation)) {
    segment_color_var <- 'clust'
    dend_data$segments$clust <- as.factor(dend_data$segments$clust)
  }
  else {
    segment_color_var <- "CTN_group"
    label_df_CTN_group <- label_df %>%
      left_join(distinct(CTN_annotation, cluster, annotation_final, CTN_group),
                by = c("label" = "annotation_final")) %>%
      group_by(CTN_group) %>%
      summarise(x = mean(x)) %>%
      mutate(CTN_group_label = str_replace(CTN_group, "\\/", "/<br/>"))
    dend_data$segments <- dend_data$segments %>%
      left_join(distinct(label_df_CTN_group, clust, CTN_group),
                by = "clust", relationship = "many-to-many")
  }

  p <- dend_data$segments %>%
    ggplot() +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend,
                     color = .data[[segment_color_var]]))+
    scale_x_reverse(expand = expansion(mult = x_expand)) +
    # scale_y_discrete(expand = expansion(mult = c(0.005, 0.005))) +
    geom_hline(yintercept = 0.6, linetype = 2) +
    coord_flip(ylim = y_lim, clip = FALSE) +
    theme_void() +
    theme(legend.position = "none")
  if(show_text) {
    p <- p + geom_richtext(
      data = label_df, aes(x, y, label = CTN_colored),
      fill = NA, label.color = NA, hjust = 1, angle = 0, size = 3)}
  if(!is.null(CTN_annotation)) {
    p <- p + geom_richtext(
      data = label_df_CTN_group,
      aes(x, y = 0.98, label = CTN_group_label),
      fill = NA, label.color = NA, hjust = 1,
      angle = 0, size = 3, lineheight = 0.5) +
      scale_color_manual(values = CTN_group_palette)
  }
  return(list(p = p, label_df = label_df))
}


#' CTN cell type proportion dotplot
#' Dotplot visualization of cell type abundance within each cell-type triad niche (CTN).
#' @param CTN_merged_celltype_prop A data frame representing the
#'   cell type proportion for each CTN with the following columns: \code{CTN} (the CTN names),
#'   \code{Celltype} (cell type labels), and \code{prop} (cell type proportion of the CTN).
#' @param CTN_dend Output of the [draw_CTN_dendro()] function
#' @param CTN_annotation (Optional) A data frame with the original CTN names (\code{CTN}), the corresponding cluster label (\code{cluster}), and the new name for the
#' consolidated CTN (specified using \code{annot_col}) after merging.
#' @param celltype_palette A named vector with colors as values and cell type labels as names. Must contain all cell types.
#' @param annot_col Column name of the newly consolidated CTN in \code{CTN_annotation}.
#' @param facet_var Faceting variables. Default to \code{"clust"} from \code{CTN_dend$label_df}
#' @param x_expand A numeric vector of length two of the multiplicative range expansion factors for the x axis. See [ggplot2::expansion()] for details.
#' @param rename_celltype A function for renaming cell type labels.
#' @details
#' Dotplot visualization of the cell type abundance of cell-type triad niches (CTN). Bubbles are colored by cell type labels, and their sizes represent relative cell type proportions with each CTN.
#' @import ggplot2 ggtext cowplot
#' @export
draw_CTN_celltype_prop <- function(CTN_merged_celltype_prop, CTN_dend, CTN_annotation,
                                   annot_col = "annotation_final",
                                   celltype_palette,
                                   rename_celltype = waiver(),
                                   facet_var = "clust",
                                   x_expand = c(0.05, 0.05)) {
  CTN_merged_celltype_prop %>%
    left_join(distinct(CTN_annotation, cluster, .data[[annot_col]]), by = "cluster") %>%
    left_join(CTN_dend$label_df, by = c("annotation_final" = "label")) %>%
    ggplot(aes(x = Celltype, y = CTN_colored, fill = Celltype)) +
    geom_point(aes(size = prop), shape = 21) +
    scale_fill_manual(values = celltype_palette, name = "Cell Type",
                      breaks = names(celltype_palette),
                      labels = rename_celltype) +
    scale_x_discrete(expand = expansion(mult = x_expand)) +
    scale_y_discrete(limits = rev, position = "right") +
    scale_radius(range = c(0.5, 8), name = "Proportion\nof Cells") +
    labs(x = NULL, y = NULL) +
    theme_cowplot(font_size = 10) +
    coord_cartesian(clip = "off") +
    guides(fill = guide_legend(ncol = 1, override.aes = list(size = 3))) +
    facet_grid(.data[[facet_var]] ~., scales = "free_y", space = "free_y", switch = "y") +
    theme(panel.spacing.y = unit(0, "mm"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_markdown(angle = 45, hjust = 1),
          axis.text.y.right = element_markdown(),
          strip.text.y.right = element_markdown(angle = 0, face = 2),
          strip.background.y = element_rect(fill = "white", color = "grey30"),
          legend.position = "right",
          legend.key.size = unit(3, "mm"))
}





