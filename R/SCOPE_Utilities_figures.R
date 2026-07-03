#' Separating cell types in CTN labels with colored dots
#'
#' Adding colored dots before each cell type in cell type triad (CTN) labels. Outputs can be used in titles and labels \code{ggplot} with markdown text enabled using [ggtext::element_markdown].
#'
#' @param CTN CTN label; e.g. B_CD4T_CD8T.
#' @param celltype_palette A named vector with colors as values and cell type labels as names. Must contain all cell types.
#' @param dot_fontsize Size of the dots. Default is 16pt.
#' @param sep A character string separating the cell types in CTN. Note that it can't appear in cell type labels.
#' @import stringr
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
#' #  labs(title = CTN_colored) +
#' #  theme(plot.title = ggtext::element_markdown(hjust = 0.5))
color_CTN_names <- function(CTN, celltype_palette = c("B" =  "#fee440", "CD4T" = "#3a86ff", "CD8T" = "#8338ec"),
                            dot_fontsize = 16, sep = "_") {
  if (!requireNamespace("ggtext", quietly = TRUE)) {
    stop(paste("Package ggtext is required but not installed."))
  }
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
#' @importFrom rlang .data
#' @export
#' @description
#' \code{cowplot} package is required for the function.
draw_silhoutte_score <- function(hclust_silhouette, num_clusters) {
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop(paste("Package cowplot is required but not installed."))
  }
  hclust_silhouette %>%
    pivot_longer(1:2, names_to = "metric", values_to = "si") %>%
    ggplot(aes(x = .data$K, y = .data$si, color = .data$metric)) +
    geom_point(size = 0.5) +
    geom_line() +
    geom_vline(xintercept = num_clusters,
               linetype = 2, color = "black") +
    scale_color_brewer(palette = "Set1") +
    scale_x_continuous(breaks = c(seq(0, 200, 50), num_clusters)) +
    scale_color_discrete(labels = c("Average Silhouette Score",
                                    "Median Silhouette Score"),
                         name = NULL) + labs(y = NULL) +
    cowplot::theme_cowplot() +
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
#' @description
#' \code{ComplexHeatmap} and \code{grid} packages are required for the function.
#' @export
draw_htmap_jaccard_index <- function(jaccard_mat, cl, num_clusters) {
  packages <- c("ComplexHeatmap", "grid")
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  }
  p <- ComplexHeatmap::Heatmap(
    jaccard_mat,
    row_names_gp = grid::gpar(fontsize = 7),
    column_names_gp = grid::gpar(fontsize = 7),
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
        grid::grid.text("*", x, y, gp = grid::gpar(fontsize = 8), vjust = 0.85)
      }
      if(between(jaccard_mat[i, j], 0.9, 0.95)) {
        grid::grid.text(".", x, y, gp = grid::gpar(fontsize = 8), vjust = 0)
      }
    }
  )

  p <- ComplexHeatmap::draw(p)

  for(i in 1:num_clusters) {
    ComplexHeatmap::decorate_heatmap_body(
      "Jaccard\nSimilarity", row_slice = i, column_slice = i, {
      grid::grid.rect(
        grid::unit(0, "npc"),
        grid::unit(1, "npc"), # top left
        gp = grid::gpar(fill = "transparent", lwd = 1.5),
        just = c("left", "top"))
    }) }

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
#' @importFrom stats cutree
#' @keywords internal
#' @description
#' \code{ggdendro} package is required for the function.
#' @references
#' Atrebas. Dendrograms in R, a lightweight approach. \url{https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/}
prepare_dendro_data <- function(cl, h, k = NULL) {
  if (!requireNamespace("ggdendro", quietly = TRUE)) {
    stop(paste("Package ggdendro is required but not installed."))
  }

  dend_data <-  ggdendro::dendro_data(cl, type = "rectangle")
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
  dend_data$segments$clust    <- factor(segclust, unique(segclust))
  dend_data$segments$line   <- as.integer(segclust < 1L)
  dend_data$labels$clust    <- factor(labclust, unique(labclust))
  dend_data$labels$order    <- 1:nrow(dend_data$labels)

  return(dend_data)
}

#' Draw CTN dendrogram
#' Dendrogram visualization for the hierarchical clustering results on the pairwise Jaccard distance between CTNs.
#' @param cl An object of class "hclust" which describes the tree produced by the clustering process. Output of the [jaccard_dist_hclust()] function.
#' @param k An integer scalar or vector with the desired number of groups. See [stats::cutree()] for details.
#' @param h Numeric scalar or vector with heights where the tree should be cut. [stats::cutree()] for details. At least one of \code{k} or \code{h} must be specified, \code{k} overrides \code{h} if both are given.
#' @param CTN_annotation (Optional) A data frame with the original CTN names (\code{CTN}), the corresponding cluster label (\code{cluster}), and the new name for the consolidated CTN (specified in \code{annot_var}) after merging.
#' @param CTN_group_var Faceting variables (CTN category). Default to \code{"clust"}. Only used when \code{CTN_annotation} is provided.
#' @param CTN_group_palette (Optional) A named vector with colors as values and CTN categories as names. Only used when \code{CTN_annotation} is provided.
#' @param annot_var Column name of the newly consolidated CTN in \code{CTN_annotation}.
#' @param celltype_palette Input of the [color_CTN_names()] function. A named vector with colors as values and cell type labels as names. Must contain all cell types.
#' @param dot_fontsize Input of the [color_CTN_names()] function. Size of the dots. Default is 16pt.
#' @param show_text Whether to display CTN names alongside the dendrogram (default: \code{TRUE}).
#' @param y_lim A numeric vector of length two providing limits of the scale. Use NA to refer to the existing minimum or maximum.
#' @param x_expand A numeric vector of length two of the multiplicative range expansion factors for the x axis. See [ggplot2::expansion()] for details.
#' @param rename_celltype A function for renaming cell type labels.
#' @return A list with components:
#' \describe{
#'   \item{p}{Dendrogram}
#'   \item{label_df}{Label data. See [ggdendro::dendro_data()] for details.}
#'   }
#' @import stringr dplyr ggplot2
#' @importFrom rlang .data
#' @description
#' \code{forcats}, \code{ggdendro}, \code{ggtext} packages are required for the function.
#' @export
draw_CTN_dendro <- function(cl, h, k = NULL, CTN_annotation = NULL,
                            CTN_group_var = "clust", annot_var = "annotation",
                            CTN_group_palette = NULL,
                            celltype_palette = NULL, dot_fontsize = 12, show_text = TRUE,
                            x_expand = c(0.001, 0.001), y_lim = c(-2, 1),
                            rename_celltype = NULL) {
  packages <- c("forcats", "ggdendro", "ggtext")
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  }

  dend_data <- prepare_dendro_data(cl = cl, h = h, k = k)
  label_df <- dend_data$labels
  segment_df <- dend_data$segments
  if(!is.null(CTN_annotation)) {
    label_df <- dend_data$labels %>%
      left_join(distinct(CTN_annotation, .data$cluster, .data[[annot_var]], .data[[CTN_group_var]]),
                by = c("label" = annot_var)) %>%
      mutate(CTN_group = forcats::fct_reorder(.data[[CTN_group_var]], .data$order)) %>%
      mutate(CTN_group_label = str_replace(.data[[CTN_group_var]], "\\/", "/<br/>")) %>%
      mutate(CTN_group_label = forcats::fct_reorder(.data$CTN_group_label, .data$order))
    label_df_CTN_group <- label_df %>%
      group_by(.data$CTN_group, .data$CTN_group_label) %>%
      summarise(x = mean(.data$x)) %>%
      drop_na()
    segment_df <- segment_df %>%
      left_join(distinct(label_df, .data$clust, .data$CTN_group),
                by = "clust", relationship = "many-to-many")
  }

  label_df <- label_df %>%
    mutate(CTN_colored = color_CTN_names(
      .data$label, celltype_palette = celltype_palette,
      sep = "_", dot_fontsize = dot_fontsize))

  if(!is.null(rename_celltype)) {
    label_df$CTN_colored <- rename_celltype(label_df$CTN_colored)
  }
  label_df$CTN_colored <- factor(label_df$CTN_colored, unique(label_df$CTN_colored))


  p <- segment_df %>%
    ggplot() +
    geom_segment(aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
                     color = .data[[CTN_group_var]]))+
    scale_x_reverse(expand = expansion(mult = x_expand)) +
    geom_hline(yintercept = 0.6, linetype = 2) +
    coord_flip(ylim = y_lim, clip = FALSE) +
    theme_void() +
    theme(legend.position = "none")
  if(show_text) {
    p <- p + ggtext::geom_richtext(
      data = label_df, aes(.data$x, .data$y, label = .data$CTN_colored),
      fill = NA, label.color = NA, hjust = 1, angle = 0, size = 3)
    }
  if(!is.null(CTN_annotation)) {
    p <- p + ggtext::geom_richtext(
      data = label_df_CTN_group,
      aes(x = .data$x, y = 0.98, label = .data$CTN_group_label),
      fill = NA, label.color = NA, hjust = 1,
      angle = 0, size = 3, lineheight = 0.5) +
      scale_color_manual(values = CTN_group_palette)
  }
  return(list(p = p, label_df = label_df))
}


#' CTN cell type proportion dotplot
#'
#' Dotplot visualization of cell type abundance within each cell-type triad niche (CTN).
#' @param CTN_merged_celltype_prop A data frame representing the
#'   cell type proportion for each CTN with the following columns: \code{CTN} (the CTN names),
#'   \code{Celltype} (cell type labels), and \code{prop} (cell type proportion of the CTN).
#' @param CTN_dend Output of the [draw_CTN_dendro()] function
#' @param facet_var Faceting variables, usually CTN category. Default to \code{"clust"} from \code{CTN_dend$label_df}
#' @param celltype_palette A named vector with colors as values and cell type labels as names. Must contain all cell types.
#' @param rename_celltype A function for renaming cell type labels.
#' @param x_expand A numeric vector of length two of the multiplicative range expansion factors for the x axis. See [ggplot2::expansion()] for details.
#' @param radius_range A numeric vector of length two that specifies the minimum and maximum size of the bubbles. See [ggplot2::scale_radius()] for details.
#' @details
#' Dotplot visualization of the cell type abundance of cell-type triad niches (CTN). Bubbles are colored by cell type labels, and their sizes represent relative cell type proportions with each CTN.
#' @description
#' \code{cowplot} and \code{ggtext} packages are required for the function.
#' @import ggplot2
#' @importFrom rlang .data
#' @export
draw_CTN_celltype_prop <- function(CTN_merged_celltype_prop, CTN_dend,
                                   facet_var = "clust",
                                   celltype_palette, radius_range = c(0.5, 8),
                                   rename_celltype = waiver(),
                                   x_expand = c(0.05, 0.05)) {
  packages <- c("ggtext", "cowplot")
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  }

  CTN_merged_celltype_prop %>%
    left_join(CTN_dend$label_df, by = c("annotation" = "label")) %>%
    ggplot(aes(x = .data$Celltype, y = .data$CTN_colored, fill = .data$Celltype)) +
    geom_point(aes(size = .data$prop), shape = 21) +
    scale_fill_manual(values = celltype_palette, name = "Cell Type",
                      breaks = names(celltype_palette),
                      labels = rename_celltype) +
    scale_x_discrete(expand = expansion(mult = x_expand)) +
    scale_y_discrete(limits = rev) +
    scale_radius(range = radius_range, name = "Proportion\nof Cells") +
    labs(x = NULL, y = NULL) +
    cowplot::theme_cowplot(font_size = 10) +
    coord_cartesian(clip = "off") +
    guides(fill = guide_legend(ncol = 1, override.aes = list(size = 3))) +
    facet_grid(.data[[facet_var]] ~., scales = "free_y", space = "free_y", switch = "y") +
    theme(panel.spacing.y = unit(0, "mm"),
          panel.background = element_rect(fill = NA, color = "grey90"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1),
          axis.text.y.left = ggtext::element_markdown(),
          strip.text.y.left = ggtext::element_markdown(angle = 0, face = 2),
          strip.background.y = element_rect(fill = NA, color = "grey"),
          strip.clip = "off",
          legend.position = "right",
          legend.key.size = unit(3, "mm"))
}

#' Plot Cox regression results on patient-level CTN presence
#'
#' Visualizes the prognostic value and cellular composition of Cell-type Triad Niches (CTNs).
#' The left panel features a lollipop plot where the x-axis represents the log hazard ratio (HR) from Cox regression models,
#' with lines colored by P-values adjusted via the Benjamini-Hochberg method (red: FDR < 0.05; black: 0.05 \eqn{\le} FDR < 0.1).
#' The right panel features an aligned bar plot illustrating the cell type composition for each CTN,
#' with core cell types are highlighted with black boxes.
#' @param coxph_res A data frame representing the Cox proportional hazards regression results
#' on patient-level CTN presence. Output of the [CTN_presence_coxph()] function
#' @param CTN_merged_celltype_prop A data frame representing the
#'   cell type proportion for each CTN with the following columns: \code{CTN} (the CTN names),
#'   \code{Celltype} (cell type labels), and \code{prop} (cell type proportion of the CTN).
#' @param top An integer specifying the maximum number of top-ranked CTNs ordered by \eqn{|logHR|} to display on the plot.
#' @param alpha A numeric value specifying the False Discovery Rate (FDR) significance threshold for filtering the CTNs.
#' @param celltype_levs A character vector defining the levels or ordering of the cell types.
#' @param celltype_palette Input of the [color_CTN_names()] function. A named vector with colors as values and cell type labels as names. Must contain all cell types.
#' @param dot_fontsize Input of the [color_CTN_names()] function. Size of the dots. Default is 16pt.
#' @param rename_celltype A function for renaming cell type labels.
#' @param line_width A numeric value specifying the line width for the boxes highlighting the core cell types.
#' @param y_text_size A numeric value specifying the font size of the y-axis labels.
#' @param width_ratio A numeric vector of length 2 specifying the width ratios between the log harzard ratio lollipop plot and cell type proportion bar plot.
#' @return A combined [ggplot2::ggplot] object containing both the survival lollipop plot and the cell type composition bar plot.
#' @import dplyr ggplot2
#' @importFrom rlang .data
#' @description
#' \code{ggtext}, \code{forcats}, \code{cowplot}, and \code{patchwork} packages are required for the function.
#' @export
draw_CTN_coxph_logHR <- function(
    coxph_res, CTN_merged_celltype_prop,
    top = NULL, alpha = 0.1, celltype_levs,
    celltype_palette, dot_fontsize = 14,
    rename_celltype = NULL,
    line_width = 0.4, y_text_size = NULL,
    width_ratio = c(0.5, 0.5)) {
  packages <- c("ggtext", "forcats", "cowplot", "patchwork")
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  }
  df_tmp <- as.data.frame(coxph_res) %>%
    group_by(sign(.data$logHR)) %>%
    arrange(desc(abs(.data$logHR))) %>%
    mutate(rank = row_number()) %>% ungroup %>%
    filter(.data$FDR < alpha) %>%
    mutate(CTN_colored = color_CTN_names(
      .data$CTN, celltype_palette = celltype_palette,
      sep = "_", dot_fontsize = dot_fontsize)) %>%
    mutate(signif = case_when(
      .data$FDR < 0.05 ~ "FDR < 0.05",
      .data$FDR < 0.1 ~ "FDR < 0.1",
      TRUE ~ "n.s."))

  if(!is.null(rename_celltype)) {
    df_tmp$CTN_colored <- rename_celltype(df_tmp$CTN_colored)
  }

  df_tmp <- df_tmp %>%
    mutate(CTN_colored = forcats::fct_reorder(
      .data$CTN_colored, .data$logHR, .desc = TRUE))

  if(!is.null(top)) {
    df_tmp <- df_tmp %>%
      filter(.data$rank <= top)
  }

  p_logHR <- df_tmp %>%
    ggplot(aes(x = .data$logHR, y = .data$CTN_colored, color = .data$signif)) +
    geom_segment(aes(xend = 0)) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_color_manual(values = c("FDR < 0.05"="red", "FDR < 0.1"= "black", "n.s." = "grey80"), name = NULL) +
    cowplot::theme_cowplot() +
    labs(y = NULL, title = NULL) +
    theme(axis.text.y = ggtext::element_markdown(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.location = "plot",
          legend.justification = "right",
          legend.margin = margin(t = -5, l = -10))


  # Cell type proportion barplot of top CTNs
    p_celltype_prop <- df_tmp[, c("CTN", "CTN_colored", "logHR")] %>%
      inner_join(CTN_merged_celltype_prop, by = c("CTN" = "annotation")) %>%
      mutate(color = ifelse(str_detect(.data$CTN, as.character(.data$Celltype)),
                            "black", "#FFFFFF00")) %>%
      arrange(.data$color) %>%
      ggplot(aes(x = .data$prop, y = .data$CTN_colored,
                 fill = forcats::fct_rev(.data$Celltype))) +
      geom_col(aes(color = .data$color), linewidth = line_width) +
      labs(x = "Proportion of Cells", y = NULL) +
      cowplot::theme_cowplot() +
      guides(fill = guide_legend(nrow = 4)) +
      scale_radius() +
      coord_cartesian(expand = FALSE) +
      scale_color_identity(breaks = c("black", "#FFFFFF00"))+
      scale_fill_manual(values = celltype_palette, labels = rename_celltype,
                        breaks = celltype_levs) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(hjust = 1, angle = 45),
            axis.text.y = ggtext::element_markdown(hjust = 1, size = y_text_size),
            legend.text = ggtext::element_markdown(),
            legend.position = "bottom",
            legend.key.width = unit(4, "mm"),
            legend.key.height = unit(3, "mm"),
            legend.key.spacing.y = unit(0, "mm"),
            legend.location = "plot",
            legend.justification = "right",
            strip.clip = "off",
            strip.text.y.left = element_text(angle = 0))

    p_celltype_prop <- p_celltype_prop +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_text(angle = 0, hjust = 0.5),
            axis.line.y = element_blank(),
            legend.position = "none")
    p_logHR + p_celltype_prop + patchwork::plot_layout(widths = width_ratio)
}


#' Spatial Visualization of CTNs and Cell Types
#'
#' Generates plots to visualize the spatial distribution of Celltype Triad Niches (CTNs) and individual cell types across tissue sections.
#' @param Full_CTN_table A data frame containing spatial coordinates (\code{X}, \code{Y}), image IDs (\code{ImageID}), cell type
#'   labels (\code{Celltype}), and niche labels for all cell-type triads
#' @param CTN CTN label; e.g. B_CD4T_CD8T.
#' @param img_list A list of images to be visualized.
#' @param celltype_palette A named vector with colors as values and cell type labels as names. Must contain all cell types.
#' @param CTN_palette A named vector with colors as values and CTN labels (e.g. CTN vs. Unassigned) as names.
#' @param scales \code{scales} parameter [ggplot2::facet_wrap()].
#' @param ncol Number of columns passed to [ggplot2::facet_wrap()].
#' @param legend.ncol Number of columns of legend keys, passed to [ggplot2::guide_legend()].
#' @param rename_celltype A function for renaming cell type labels.
#' @param aspect.ratio Aspect ratio of the panel, passed directly to [ggplot2::theme()].
#' @returns A list of [ggplot2::ggplot] objects:
#' \describe{
#'   \item{p_celltype}{A spatial plot showing individual cell type distributions.}
#'   \item{p_CTN}{A spatial plot showing the distribution of the specified CTN.}
#'   }
#' @import ggplot2
#' @importFrom rlang .data
#' @description
#' \code{ggtext} package are required for the function.
#' @export
draw_CTN_celltype <- function(Full_CTN_table, CTN = "Bcell_CD4T_CD8T", img_list,
                              celltype_palette, CTN_palette = c("CTN" = "#E69F00", "Unassigned" = "grey90"),
                              rename_celltype = NULL,
                              scales = "fixed", ncol = length(img_list),
                              legend.ncol = NULL, aspect.ratio = NULL) {
  if (!requireNamespace("ggtext", quietly = TRUE)) {
    stop(paste("Package ggtext is required but not installed."))
  }
  p_celltype <- Full_CTN_table %>%
    filter(.data$ImageID %in% img_list) %>%
    mutate(ImageID = factor(.data$ImageID, img_list)) %>%
    ggplot(aes(x = .data$X, y = .data$Y)) +
    geom_point(aes(color = .data$Celltype), size = 0.1, shape = 16) +
    scale_color_manual(values = celltype_palette, name = "Cell Type",
                       labels = rename_celltype) +
    facet_wrap(~ .data$ImageID, ncol = ncol, scales = scales) +
    labs(title = "Cell Type") +
    guides(color = guide_legend(ncol = legend.ncol, override.aes = list(size = 3))) +
    theme_void() +
    theme(legend.key.size = unit(3, "mm"),
          legend.text = ggtext::element_markdown(),
          plot.title = element_text(size = 12, face = 2))

  p_CTN <- Full_CTN_table %>%
    filter(.data$ImageID %in% img_list) %>%
    mutate(ImageID = factor(.data$ImageID, img_list)) %>%
    ggplot(aes(x = .data$X, y = .data$Y)) +
    geom_point(aes(color = .data[[CTN]]), size = 0.1, shape = 16) +
    scale_color_manual(values = CTN_palette, name = "CTN label") +
    facet_wrap( ~ .data$ImageID, ncol = ncol, scales = scales) +
    guides(color = guide_legend(ncol = legend.ncol, override.aes = list(size = 3))) +
    theme_void() +
    labs(title = rename_celltype(color_CTN_names(
      CTN, celltype_palette = celltype_palette))) +
    theme(legend.key.size = unit(3, "mm"),
          legend.text = ggtext::element_markdown(),
          plot.title = ggtext::element_markdown(size = 12, face = 2))

  if(is.null(aspect.ratio)) {
    p_celltype <- p_celltype + coord_fixed(expand = FALSE)
    p_CTN <- p_CTN + coord_fixed(expand = FALSE)
  }
  else {
    p_celltype <- p_celltype + coord_cartesian(expand = FALSE) +
      theme(aspect.ratio = aspect.ratio)
    p_CTN <- p_CTN + coord_cartesian(expand = FALSE) +
      theme(aspect.ratio = aspect.ratio)
  }

  return(list(p_celltype = p_celltype, p_CTN = p_CTN))
}

