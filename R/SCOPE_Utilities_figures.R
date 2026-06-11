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
#' ggplot() +
#'   labs(title = CTN_colored) +
#'   theme(plot.title = element_markdown(hjust = 0.5))

color_CTN_names <- function(CTN, celltype_palette = c("B" =  "#fee440", "CD4T" = "#3a86ff", "CD8T" = "#8338ec"),
                            dot_fontsize = 16, sep = "_") {
  driver_celltypes <- str_split_1(CTN, sep)
  colors <- celltype_palette[driver_celltypes]
  CTN_colored <- paste0("<span style='color:", colors, ";font-size:",
                        dot_fontsize, "pt'>\U25CF</span>", driver_celltypes)
  CTN_colored <- paste(CTN_colored, collapse = " ")
  return(CTN_colored)
}
