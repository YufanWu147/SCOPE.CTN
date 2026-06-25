#' Constructing neighborhood cell type count and proportion matrices
#'
#' Define spatial neighborhood using [FNN::get.knn()] and compute neighborhood cell counts and proportions by all available cell types.
#' @param dat_slide Input data frame for a single tissue section. Columns must contain cell ids (\code{CellID}), spatial coordinates (\code{X}, \code{Y}), and cell type labels (\code{Celltype}).
#' @param knn_num The number of nearest neighbors. The default value is 20.
#' @param cell_types_available A vector of length \eqn{T} containing all available cell type labels.
#' @param min.radius Lower threshold for neighborhood distance boundary. The default value is 100.
#' @param max.radius Upper threshold for neighborhood distance boundary. The default value is 500.
#' @param min.n.neighbors Minimum number of nearest neighbors within the neighborhood distance boundary. Cells with fewer than \code{min.n.neighbors} (default: 10) neighbors will be excluded from downstream analysis.
#' @param return.nn.ids Whether to return cell ids of the nearest neighbors of each cell.
#' @param prop.list Probabilities between \eqn{[0, 1]} to produce quantiles for neighborhood distance distribution
#' @details
#' The neighborhood distance boundary \code{max.dist} is initially set as the 99th percentile of all distances from focal cells to their the nearest neighbors. If it falls outside the predefined range (default: \eqn{[100, 500]}), it will be reset to the minimum (\code{min.radius}) or maximum threshold (\code{max.radius}).
#'
#' For cells with fewer than \code{min.n.neighbors} neighbors within the boundary, their corresponding neighborhood cell type counts and proportions will all be set to \code{NaN} and they will be subsequently excluded from niche clustering.
#'
#' @return A list contains:
#' \describe{
#'   \item{count.df}{An \eqn{N \times T} neighborhood cell type count matrix. Rows with all \code{NaN} values are cells with fewer than fewer than \code{min.n.neighbors} neighbors within the distance boundary.}
#'   \item{prop.df}{An \eqn{N \times T} neighborhood cell type proportion matrix. Rows with all \code{NaN} values are cells with fewer than fewer than \code{min.n.neighbors} neighbors within the distance boundary.}
#'   \item{dist.quantiles}{Quantiles (including the 99th percentile) of the of the distances to all nearest neighbors corresponding to the probabilities in \code{prop.list}. Can be used to examine the distribution of the distances.}
#'   \item{neighbor.cell.ids}{A data frame containing cell ids of the nearest neighbors of each focal cell.}
#' }
#' @import dplyr
#' @importFrom data.table rbindlist
#' @importFrom FNN get.knn
#' @importFrom stats quantile
#' @export
Build_cell_neighbor_maxdist <- function(
    dat_slide, knn_num, cell_types_available,
    min.radius = 100, max.radius = 500, min.n.neighbors = 10,
    return.nn.ids = FALSE, prob_list = seq(0, 1, 0.25)){
  nn_result <- FNN::get.knn(data = dat_slide[, c("X", "Y")], k = knn_num + 1)
  neighbor_indices <- nn_result$nn.index[, -1]
  neighbor_dists <- nn_result$nn.dist[, -1]


  # Quantiles of the distances to all nearest neighbors
  dist.quantiles <- quantile(neighbor_dists, probs = prob_list)

  # 99th percentile of the distances to all nearest neighbors
  dist.99th.perc <- quantile(neighbor_dists, 0.99)

  # Set neighborhood boundary as dist.perc.99 or min.radius, whichever is greater
  max.dist = max(dist.99th.perc, min.radius)

  # Set an upper limit for neighborhood boundary
  max.dist = min(max.dist, max.radius)

  count_df_neighbor_ids <- sapply(1:nrow(dat_slide), function(i) {
    # Only calculate cell type proportion when there are more than
    # min.n.neighbors cells within the radius
    neighbors_within_radius <- which(neighbor_dists[i, ] <= max.dist)

    if(length(neighbors_within_radius) <= min.n.neighbors) {
      return(list(count_table = rep(NaN, length(cell_types_available)),
                  neighbor.cell.ids = NULL))
    }
    else {
      neighbor_cell_types <- dat_slide$Celltype[neighbor_indices[i, neighbors_within_radius]]
      neighbor.cell.ids <- data.frame(CellID_neighbor = dat_slide$CellID[neighbor_indices[i, neighbors_within_radius]])
      count_table <- table(factor(neighbor_cell_types, levels = cell_types_available))
      return(list(count_table = as.numeric(count_table),
                  neighbor.cell.ids = neighbor.cell.ids))
    }
  }, simplify = FALSE)
  count_df <- lapply(count_df_neighbor_ids, "[[", "count_table")
  count_df <- do.call(rbind, count_df)
  rownames(count_df) <- dat_slide$CellID
  colnames(count_df) <- cell_types_available
  prop_df <- prop.table(count_df, margin = 1)

  out <- list(count.df = count_df,
              prop.df = prop_df,
              dist.quantiles = sort(c(dist.quantiles, dist.99th.perc)))

  if(return.nn.ids) {
    # save neighbor cell ids to speed up permutation test
    neighbor.cell.ids <- lapply(count_df_neighbor_ids, "[[", "neighbor.cell.ids") %>%
      `names<-`(dat_slide$CellID) %>%
      data.table::rbindlist(idcol = "CellID")
    out <- c(out, list(neighbor.cell.ids = neighbor.cell.ids))
    }

  return(out)
}

