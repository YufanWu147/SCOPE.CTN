#' Computing colocation quotient (CLQ)
#'
#' Computes colocation quotient (CLQ) between available cell type pairs.
#' @param interact_matrix Input
#' @param dat_in Input data frame with \eqn{N} rows. Columns must contain \code{CellID} and \code{Celltype}.
#' @details The colocalization of a pair of cell types \eqn{a} and \eqn{b} in tissue section \eqn{l}
#' can be evaluated with the CLQ metric
#'
#' \deqn{CLQ_{b\rightarrow a} =
#' \left\{\begin{array}{ll}
#' \frac{C_{b\rightarrow a} / N_a}{N_b / (N - 1)} & N_a > 5 \mbox{and} \ N_b > 5 \\
#' 0, & \text{otherwise}
#' \end{array}\right.}
#'
#' where \eqn{C_{b\rightarrow a}} is the number of cells of cell type \eqn{b} within the neighborhood of cell type \eqn{a},
#' and \eqn{N_a}, \eqn{N_b}, \eqn{N} are the number of cells for cell type \eqn{a},
#' \eqn{b} and the total number of cells in tissue section \eqn{l}, respectively.
#'
#' @return A \eqn{T\times T} matrix containing CLQs between every cell type pair, where \eqn{T} is the number of available cell types.
#' The rows of the matrix are the focal cell types and the columns are the neighboring cell types.
#' @keywords internal
#' @references
#' Bouchard, G. et al. A quantitative spatial cell-cell colocalizations framework enabling comparisons between in vitro assembloids and pathological specimens.
#' Nat Commun 16, 1392 (2025). \url{https://doi.org/10.1038/s41467-024-55129-6}
compute_CLQ <- function(interact_matrix, dat_in) {
  Ncell_by_celltype <- table(dat_in$Celltype)[rownames(interact_matrix)] # No. of cells by cell type
  Ncell_total <- nrow(dat_in)                 # No. of cells in total
  n_celltype <- length(Ncell_by_celltype)     # No. of cell types
  Ncell_a_mat <- matrix(rep(Ncell_by_celltype, n_celltype), n_celltype, byrow = FALSE)
  Ncell_b_mat <- matrix(rep(Ncell_by_celltype, n_celltype), n_celltype, byrow = TRUE)
  diag(Ncell_b_mat) <- Ncell_by_celltype - 1
  CLQ = interact_matrix / Ncell_a_mat / Ncell_b_mat * (Ncell_total - 1)
  # set CLQ as 0 if either the center or the neighboring cell type has <= 5 cells
  CLQ[Ncell_a_mat <= 5 | Ncell_b_mat <= 5] <- 0
  return(CLQ)
}


#' Computing observed colocation quotient (CLQ)
#'
#' Computes the observed colocation quotient (CLQ) between available cell type pairs in each tissue section in a spatial omics dataset.
#' @param cell_neighbor_counts \eqn{N \times T} neighborhood cell type count matrix obtained using [Build_cell_neighbor_maxdist()].
#' @param dat_in Input data frame with \eqn{N} rows. Columns must contain \code{CellID}, \code{ImageID} and \code{Celltype}.
#' @param min.ncell Minimum cell number threshold. Tissue sections with fewer than \code{min.ncell} cells will be excluded from the analysis.
#' @param num_cores The number of cores to use in [pbmcapply::pbmclapply()], i.e. at most how many child processes will be run simultaneously. Can only be set at 1 on Windows.
#' @details
#' The colocalization of a pair of cell types \eqn{a} and \eqn{b} in tissue section \eqn{l}
#' can be evaluated with the observed CLQ
#'
#' \deqn{CLQ_{b\rightarrow a}^{(l), obs} =
#' \left\{\begin{array}{ll}
#' 0, & N_a^{(l)} \leq 5 \mbox{ or} \ N_b^{(l)} \leq 5 \\
#' \frac{C_{b\rightarrow a}^{(l), obs} / N_a^{(l)}}{N_b^{(l)} / (N^{(l)} - 1)}, & \text{otherwise}
#' \end{array}\right.}
#'
#' where \eqn{C_{b\rightarrow a}^{(l), obs}} is the number of cells of cell type \eqn{b} within the neighborhood of cell type \eqn{a},
#' and \eqn{N_a^{(l)}}, \eqn{N_b^{(l)}}, \eqn{N^{(l)}} are the number of cells for cell type \eqn{a},
#' \eqn{b} and the total number of cells in tissue section \eqn{l}, respectively.
#'
#' @return A data frame of \eqn{T\times T\times L} rows containing observed CLQs between focal (\code{Celltype})
#' and neighboring cell types (\code{Celltype_neighbor}) in each tissue section (\code{ImageID}),
#' where \eqn{T} is the number of available cell types and \eqn{L} is the number of tissue sections.
#' @import dplyr tibble tidyr data.table pbmcapply
#' @importFrom rlang .data
#' @export
#' @references
#' Bouchard, G. et al. A quantitative spatial cell-cell colocalizations framework enabling comparisons between in vitro assembloids and pathological specimens.
#' Nat Commun 16, 1392 (2025). \url{https://doi.org/10.1038/s41467-024-55129-6}
compute_CLQ_observed <- function(cell_neighbor_counts, dat_in, min.ncell = 100, num_cores = 1) {
  ImageID <- NULL
  rownames(dat_in) <- dat_in$CellID
  cell_neighbor_counts <- drop_na(as.data.frame(cell_neighbor_counts))
  slide_list = unique(names(table(dat_in$ImageID))[which(table(dat_in$ImageID)>min.ncell)])
  celltypes_available <- colnames(cell_neighbor_counts)

  CLQ_obs <- pbmcapply::pbmclapply(slide_list, function(img) {
    dat_in_sub <- subset(dat_in, ImageID == img)[, "Celltype", drop = FALSE]
    cells_sub <- intersect(rownames(dat_in_sub), rownames(cell_neighbor_counts))
    if(length(cells_sub) < 5) {return(NULL)}

    interaction_mat_sub <- bind_cols(dat_in_sub[cells_sub, , drop = FALSE],
                                     cell_neighbor_counts[cells_sub, ]) %>%
      group_by(.data$Celltype) %>%
      summarise_at(celltypes_available, sum, na.rm = TRUE) %>%
      mutate(Celltype = factor(.data$Celltype, celltypes_available)) %>%
      complete(.data$Celltype) %>%
      mutate_at(celltypes_available, replace_na, 0) %>%
      arrange(.data$Celltype) %>%
      column_to_rownames("Celltype")

    if(any(rownames(interaction_mat_sub) != colnames(interaction_mat_sub))) {
      warning(paste0("Image ", img, ": Non-symmetrical interaction matrix"))
      return(NULL)
    }

    CLQ_obs_sub <- compute_CLQ(interaction_mat_sub, dat_in_sub) %>%
      rownames_to_column("Celltype") %>%
      pivot_longer(2:ncol(.data), values_to = "CLQ", names_to = "Celltype_neighbor") %>%
      cbind(ImageID = img, .data)

    return(CLQ_obs_sub)
  }, mc.cores = num_cores)
  CLQ_obs <- do.call(rbind, CLQ_obs)
  return(CLQ_obs)
}

#' Computing colocation quotient (CLQ) on permutated data
#'
#' Randomly shuffles the cell type labels within each tissue section and computes permutated CLQs between every cell type pair,
#' creating a null distribution that is compared with the observed CLQ to compute a p-value.
#' @param cell_neighbor_ids Cell IDs of the nearest neighbors obtained using [Build_cell_neighbor_maxdist()].
#' @param dat_in Input data frame with \eqn{N} rows. Columns must contain \code{CellID}, \code{ImageID} and \code{Celltype}.
#' @param min.ncell Minimum cell number threshold. Tissue sections with fewer than \code{min.ncell} cells will be excluded from the analysis.
#' @param num_cores The number of cores to use in [pbmcapply::pbmclapply()], i.e. at most how many child processes will be run simultaneously. Can only be set at 1 on Windows.
#' @param iter.num Number of permutations (default: 500).
#' @param seed Random seed passed to [base::set.seed()].
#' @details
#' In each tissue section, the cell type labels are randomly shuffled for \code{iter.num} (default: 500) times to generate permutated CLQ values
#'
#' \deqn{CLQ_{b\rightarrow a}^{(l), k} =
#' \left\{\begin{array}{ll}
#' 0, & N_a^{(l)} \leq 5 \mbox{ or} \ N_b^{(l)} \leq 5 \\
#' \frac{C_{b\rightarrow a}^{(l), k} / N_a^{(l)}}{N_b^{(l)} / (N^{(l)} - 1)}, & \text{otherwise}
#' \end{array}\right.}
#'
#' where \eqn{C_{b\rightarrow a}^{(l), k}} is the number of cells of cell type \eqn{b} within the neighborhood of cell type \eqn{a} in iteration \eqn{k, k=}1,2,...,\code{iter.num},
#' and \eqn{N_a^{(l)}}, \eqn{N_b^{(l)}}, \eqn{N^{(l)}} are the number of cells for cell type \eqn{a},
#' \eqn{b} and the total number of cells in tissue section \eqn{l}, respectively.
#'
#'
#' @return A data frame of \eqn{T\times T\times L\times \mbox{iter.num}} rows containing the permutated CLQs between focal (\code{Celltype})
#' and neighboring cell types (\code{Celltype_neighbor}) in each tissue section (\code{ImageID}),
#' where \eqn{T} is the number of available cell types and \eqn{L} is the number of tissue sections.
#' @import dplyr tibble tidyr data.table
#' @importFrom rlang .data
#' @export
#' @references
#' Bouchard, G. et al. A quantitative spatial cell-cell colocalizations framework enabling comparisons between in vitro assembloids and pathological specimens.
#' Nat Commun 16, 1392 (2025). \url{https://doi.org/10.1038/s41467-024-55129-6}
compute_CLQ_permutated <- function(cell_neighbor_ids, dat_in, min.ncell = 100, num_cores = 5,
                                   iter.num = 500, seed = 42) {
  set.seed(seed)
  ImageID <- CellID <- Celltype <- N <- Celltype_neighbor <- NULL
  dat_in_temp <- data.table(dat_in[, c("CellID", "ImageID", "Celltype")], key = "CellID")

  slide_list = unique(names(table(dat_in$ImageID))[which(table(dat_in$ImageID)>min.ncell)])

  # Convert cell neighbor id dataframe to data.table to speed up joining
  cell_neighbor_ids_dt <- data.table(cell_neighbor_ids, key = "CellID_neighbor")

  # Shuffling cell type labels to generate a null distribution of CLQ
  CLQ_perm <- pbmcapply::pbmclapply(1:iter.num, function(i) {

    sapply(slide_list, function(img) {
      dat_in_sub <- subset(dat_in_temp, ImageID == img)
      cell_neighbor_ids_dt_sub <- subset(cell_neighbor_ids_dt, CellID %in% dat_in_sub$CellID)
      dat_in_sub$Celltype <- sample(dat_in_sub$Celltype, size = nrow(dat_in_sub), replace = FALSE)
      interact_matrix_temp = merge(cell_neighbor_ids_dt_sub, dat_in_temp, by.x = "CellID_neighbor", by.y = "CellID", sort = FALSE)
      interact_matrix_temp = interact_matrix_temp[, .N, by = list(CellID, Celltype)]
      setkey(interact_matrix_temp, "CellID")
      interact_matrix_temp = merge(interact_matrix_temp, dat_in_sub, by = "CellID", sort = FALSE, suffixes = c("_neighbor", ""))
      interact_matrix_temp = interact_matrix_temp[, sum(N, na.rm = TRUE),
                                                  by = list(Celltype, Celltype_neighbor)]
      interact_matrix_temp = interact_matrix_temp %>%
        full_join(expand.grid(Celltype = levels(dat_in$Celltype),
                              Celltype_neighbor = levels(dat_in$Celltype)),
                  by = c("Celltype", "Celltype_neighbor")) %>%
        mutate(V1 = ifelse(is.na(.data$V1), 0, .data$V1)) %>%
        arrange(.data$Celltype_neighbor) %>%
        pivot_wider(names_from = "Celltype_neighbor", values_from = "V1", values_fill = 0) %>%
        arrange(.data$Celltype) %>%
        column_to_rownames("Celltype")

      CLQ_df <- compute_CLQ(interact_matrix_temp, dat_in_sub) %>%
        rownames_to_column("Celltype") %>%
        pivot_longer(2:ncol(.data), values_to = "CLQ",
                     names_to = "Celltype_neighbor") %>%
        cbind(iter = i)
      return(CLQ_df) }, simplify = FALSE, USE.NAMES = TRUE) %>%
      data.table::rbindlist(idcol = "ImageID")
    }, mc.cores = num_cores
  )
  CLQ_perm <- do.call(rbind, CLQ_perm)

  return(CLQ_perm)
}

#' Computing colocation quotient (CLQ) on permutated data
#'
#' Randomly shuffles the cell type labels within each tissue section and computes permutated CLQs between every cell type pair,
#' creating a null distribution that is compared with the observed CLQ to compute a permutation test p-value.
#' @param CLQ_obs Observed CLQs computed using [compute_CLQ_observed()].
#' @param CLQ_perm Permutated CLQs computed using [compute_CLQ_permutated()].
#' @details
#' The permutation test p-value of tissue section \eqn{l} is defined as
#'
#' \deqn{p_{b\rightarrow a}^{(l)} = \sum_{k=1}^{\mbox{iter.num}}I(CLQ_{b\rightarrow a}^{(l), k} > CLQ_{b\rightarrow a}^{(l), obs}),}
#'
#' where \eqn{CLQ_{b\rightarrow a}^{(l), k}} is the permutated CLQ computed using [compute_CLQ_permutated()] on the randomly shuffled cell type labels,
#' and \eqn{CLQ_{b\rightarrow a}^{(l), obs}} is the observed CLQ computed using [compute_CLQ_observed()] on the true cell type labels.
#'
#'
#' @return A data frame containing the nominal and Benjamini-Hochberg-adjusted CLQ permutation test p-values.
#' @import dplyr tibble tidyr data.table
#' @importFrom stats p.adjust
#' @importFrom rlang .data
#' @export
#' @references
#' Bouchard, G. et al. A quantitative spatial cell-cell colocalizations framework enabling comparisons between in vitro assembloids and pathological specimens.
#' Nat Commun 16, 1392 (2025). \url{https://doi.org/10.1038/s41467-024-55129-6}
CLQ_permutation_test <- function(CLQ_obs, CLQ_perm) {
  CLQ_pval <- CLQ_perm %>%
    left_join(CLQ_obs, by = c("Celltype", "Celltype_neighbor", "ImageID"),
              suffix = c("", ".obs"))  %>%
    filter(.data$Celltype != .data$Celltype_neighbor) %>%
    group_by(.data$Celltype, .data$Celltype_neighbor, .data$ImageID) %>%
    # right-tailed test p-value
    summarise(p_perm = sum(.data$CLQ >= .data$CLQ.obs)/n())

  CLQ_pval$p_perm.adj <- p.adjust(CLQ_pval$p_perm, method = "BH")
  return(CLQ_pval)
}

