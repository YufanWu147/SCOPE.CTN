#' Constructing neighborhood cell type count and proportion matrices
#'
#' Define spatial neighborhood using [FNN::get.knn()] and compute neighborhood cell counts and proportions by all available cell types.
#' @param dat_slide Input data frame for a single tissue section. Must contain columns including cell ids (\code{CellID}), spatial coordinates (\code{X}, \code{Y}), and cell type labels (\code{Celltype}).
#' @param knn_num The number of nearest neighbors. The default value is 20.
#' @param cell_types_available A vector of length \eqn{T} containing all available cell type labels.
#' @param min.radius Lower threshold for neighborhood distance boundary. The default value is 100.
#' @param max.radius Upper threshold for neighborhood distance boundary. The default value is 500.
#' @param min.n.neighbors Minimum number of nearest neighbors within the neighborhood distance boundary. Cells with fewer than \code{min.n.neighbors} (default: 10) neighbors will be excluded from downstream analysis.
#' @param return.nn.ids Whether to return cell ids of the nearest neighbors of each cell.
#' @param prob.list Probabilities between \eqn{[0, 1]} to produce quantiles for neighborhood distance distribution
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
    return.nn.ids = FALSE, prob.list = seq(0, 1, 0.25)){
  nn_result <- FNN::get.knn(data = dat_slide[, c("X", "Y")], k = knn_num + 1)
  neighbor_indices <- nn_result$nn.index[, -1]
  neighbor_dists <- nn_result$nn.dist[, -1]


  # Quantiles of the distances to all nearest neighbors
  dist.quantiles <- quantile(neighbor_dists, probs = prob.list)

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

#' Compute Jaccard similarity index
#'
#' Compute Jaccard similarity index between two sets of niche labels.
#' @param a,b Cell IDs of cells assigned to CTNs
#' @return Jaccard similarity index
#' @export
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

#' Perform hierarchical clustering on Jaccard distance
#'
#' Computes the Jaccard distance matrix between all pairs of cell-type triad niches (CTNs) and performs hierarchical clustering analysis on the result.
#' @param Full_CTN_table A data frame containing clustering labels for all cell-type triads
#' @param CTN_ls A character vector containing the names of the CTNs.
#' @param num_cores The number of cores to use in [pbmcapply::pbmclapply()], i.e. at most how many child processes will be run simultaneously. Can only be set at 1 on Windows.
#' @return A list contains:
#' \describe{
#'   \item{jaccard_mat_mgcv}{A matrix containing the Jaccard similarity indices (1-Jaccard distance) between every pair of CTNs}.
#'   \item{cl}{An object of class "hclust" which describes the tree produced by the clustering process. See [stats::hclust()] for details.}
#'   \item{hclust_silhouette}{A dataframe containing the average and median silhouette scores for evaluating cluster quality. See [cluster::silhouette()] for details.}
#'  }
#' @import dplyr
#' @importFrom pbmcapply pbmclapply
#' @importFrom stats hclust
#' @importFrom stats as.dist
#' @importFrom utils combn
#' @importFrom rlang .data
#' @export
jaccard_dist_hclust <- function(Full_CTN_table, CTN_ls, num_cores = 1) {
  cell_id_ls <- sapply(CTN_ls, function(ctn) {
    Full_CTN_table$CellID[Full_CTN_table[[ctn]] == "CTN"]
  }, simplify = FALSE, USE.NAMES = TRUE)

  if (!requireNamespace("survival", quietly = TRUE)) {
    stop(paste("Package survival is required but not installed."))
  }
  # All pairwise combinations of CTNs
  CTN_combs <- as.data.frame(t(combn(CTN_ls, 2)))

  # Compute jaccard similarity index between every pair of CTNs
  jaccard_vec_mgcv <- pbmcapply::pbmclapply(1:nrow(CTN_combs), function(i) {
    jaccard(a = cell_id_ls[[CTN_combs$V1[i]]],
            b = cell_id_ls[[CTN_combs$V2[i]]])
  }, mc.cores = num_cores) %>% unlist()

  jaccard_vec_mgcv <- cbind(CTN_combs, jaccard = jaccard_vec_mgcv)

  # Convert data frame with Jaccard indices from long to wide format
  jaccard_mat_mgcv <- jaccard_vec_mgcv %>%
    rbind(data.frame(V1 = CTN_ls, V2 = CTN_ls, jaccard = 1), # set diagonal as 1
          data.frame(V1 = jaccard_vec_mgcv$V2, V2 = jaccard_vec_mgcv$V1,
                     jaccard = jaccard_vec_mgcv$jaccard)) %>%
    arrange(.data$V2) %>%
    pivot_wider(names_from = "V2", values_from = "jaccard") %>%
    arrange(.data$V1) %>%
    column_to_rownames("V1") %>%
    as.matrix

  if(!isSymmetric(jaccard_mat_mgcv)) {warning("Non-symmetric jaccard distance matrix")}

  # Hierarchical clustering on Jaccard distance (1 - Jaccard index)
  cl = hclust(as.dist(1-jaccard_mat_mgcv), method = "complete", members = NULL)

  # Compute silhouette score
  hclust_silhouette <- sapply(2:(nrow(jaccard_mat_mgcv) - 1), function(k) {
    si = cluster::silhouette(cutree(cl, k = k), as.dist(1 - jaccard_mat_mgcv))
    data.frame(avg_si = summary(si)$avg.width,
               median_si = summary(si)$si.summary[["Median"]])}, simplify = FALSE)
  hclust_silhouette <- do.call(rbind, hclust_silhouette) %>%
    cbind(K = 2:(nrow(jaccard_mat_mgcv) - 1))

  return(list(jaccard_mat_mgcv = jaccard_mat_mgcv, cl = cl,
              hclust_silhouette = hclust_silhouette))
}

#' Rename Consolidated CTNs Based on Dominant Cell Types
#'
#' Renaming the consolidated Cell-type Triads Niches (CTNs) using the most frequent core
#' cell types within the original triads. If a tie occurs, the function resolves
#' it by identifying the cell type with the highest overall proportion.
#' @param CTN_cluster_label A data frame containing hierarchical clustering results for merging the highly overlapped CTNs.
#'   It must include the following two columns: \code{CTN} (the CTN names) and \code{cluster} (the corresponding cluster label).
#' @param CTN_merged_celltype_prop A data frame representing the
#'   cell type proportion for each CTN with the following columns: \code{CTN} (the CTN names),
#'   \code{Celltype} (cell type labels), and \code{prop} (cell type proportion of the CTN).
#' @param celltype_levs A character vector defining the levels or ordering of the cell types.
#' @return A list contains:
#' \describe{
#'   \item{annotation}{A data frame with the original CTN names (\code{CTN}),
#'   the corresponding cluster label (\code{cluster}), and the new name for the consolidated CTN (\code{annotation}) after merging.}
#'   \item{core_celltypes_summary}{A data frame summarizing clusters with at least 2 CTNs.
#'   It includes the number of appearances of each core cell type in the original CTN names (\code{n_appearance}) and their relative proportion within the
#'   consolidated CTNs (\code{prop}).}
#' }
#' @details
#' If the top two most frequent cell types only appear twice in the original CTN names, the function flags the consolidated CTN with an asterisk (*).
#' Duplicated new CTN names are flagged with an plus sign ("+").
#' We suggest double-checking the flagged CTNs and manually renaming them if necessary.
#' @import dplyr
#' @importFrom data.table rbindlist
#' @importFrom rlang .data
#' @export
CTN_annotate <- function(CTN_cluster_label, CTN_merged_celltype_prop, celltype_levs) {
  res <- CTN_cluster_label %>%
    split(f = .$cluster) %>%
    lapply(function(df_sub) {
      if(nrow(df_sub) == 1) {
        annotation = str_split(df_sub$CTN[1], "_")[[1]] %>%
          factor(levels = celltype_levs) %>% sort() %>%
          paste(collapse = "_")
        core_celltypes_count = NULL}
      else {
        core_celltypes <- str_split(df_sub$CTN, "_") %>%
          lapply(as.data.frame) %>%
          `names<-`(df_sub$CTN) %>%
          data.table::rbindlist(idcol = "CTN") %>%
          `colnames<-`(c("CTN", "CoreCelltypes"))

        core_celltypes_count <- core_celltypes %>%
          dplyr::count(.data$CoreCelltypes, name = "n_appearance") %>%
          inner_join(subset(CTN_merged_celltype_prop, cluster == df_sub$cluster[1]),
                     by = c("CoreCelltypes" = "Celltype")) %>%
          arrange(desc(.data$n_appearance), desc(.data$prop)) %>%
          select(all_of(c("cluster", "CoreCelltypes", "n_appearance", "prop"))) %>%
          distinct()

        # Flag new CTN name with * if the most second most frequent cell types only appeared twice
        core_celltypes_count$Note <- ifelse(
           any(core_celltypes_count$n_appearance[1:2] == 2), "*", "")

        annotation <- core_celltypes_count %>%
          dplyr::slice(1:3) %>%
          mutate(CoreCelltypes = factor(.data$CoreCelltypes, celltype_levs)) %>%
          arrange(.data$CoreCelltypes) %>%
          pull(.data$CoreCelltypes) %>% sort() %>%
          paste(collapse = "_")
      }

      return(list(annotation = cbind(df_sub, annotation),
                  core_celltypes_summary = core_celltypes_count))
      })

  core_celltypes_summary <- do.call(rbind, lapply(res, "[[", "core_celltypes_summary")) %>%
    `rownames<-`(NULL) %>% distinct()
  annotation <- do.call(rbind, lapply(res, "[[", "annotation")) %>%
    `rownames<-`(NULL) %>%
    left_join(distinct(core_celltypes_summary, .data$cluster, .data$Note),
              by = "cluster") %>%
    mutate(Note = ifelse(is.na(.data$Note), "", .data$Note)) %>%
    # Flag new CTN name with '+' if duplicated
    group_by(annotation) %>%
    mutate(Note = paste(c(.data$Note,
      ifelse(length(unique(.data$cluster)) > 1, "+", NULL)),
      sep = ";"))

  return(list(annotation = annotation,
              core_celltypes_summary = core_celltypes_summary))
}


#' Cox regression on patient survival based on CTN presence
#'
#' Evaluates the prognostic ability of Celltype Triad Niches (CTNs) by fitting a
#' Cox proportional hazards regression model using patient-level CTN presence as the main predictor.
#' @param Full_CTN_table A data frame containing clustering labels for all cell-type triads
#' @param CTN_ls A character vector containing the names of the CTNs to evaluate.
#' @param clinical_df A data frame containing image-level clinical data. It must
#'   include columns for patient ID (\code{patient_id}), survival time (\code{time}),
#'   status event (\code{event}), and any specified confounders (\code{confounder}).
#' @param time A character string specifying the survival time column name. See [survival::Surv()] for details.
#' @param event A character string specifying the survival status indicator column name, normally 0=alive, 1=dead. See [survival::Surv()] for details.
#' @param patient_id A character string specifying the patient ID column name.
#' @param confounder A character string or vector specifying column names of
#'   confounders to be adjusted for in the model. Defaults to \code{NULL}.
#' @param min.n.cells An integer specifying the minimum cell count threshold
#'   required to define a CTN as present on a tissue section. Defaults to 100.
#' @param min.n.patients An integer specifying the minimum patient threshold.
#'   Only CTNs that are both present present and absent in \eqn{\ge} \code{min.n.patients} will be evaluated. Defaults to 3.
#' @returns A data frame containing the survival analysis results for each CTN:
#' \describe{
#'  \item{CTN}{The name of the CTN.}
#'  \item{logHR}{Log harzard ratio.}
#'  \item{HR}{Hazard ratio.}
#'  \item{pval}{P-value}
#'  \item{c_index}{Concordance-index (c-index)}
#'  \item{FDR}{The Benjamini-Hochberg adjusted p-value}
#' }
#' @import survival dplyr
#' @importFrom stats as.formula
#' @importFrom rlang .data
#' @export
#' @details
#' Patient-level presence of a CTN is defined as having \eqn{\ge} \code{min.n.cells}
#' cells within that niche in at least one tissue section.
CTN_presence_coxph <- function(
    Full_CTN_table, CTN_ls, clinical_df,
    time = "OS", event = "Ev.O",
    patient_id = "Patient_ID", confounder = NULL,
    min.n.patients = 3, min.n.cells = 100) {

  clinical_df_by_patient <- clinical_df %>%
    select(all_of(c(patient_id, time, event, confounder))) %>%
    distinct()

  CTN_presence_by_patient <- Full_CTN_table %>%
    group_by(.data$ImageID) %>%
    summarise(across(all_of(CTN_ls), ~ ifelse(
      sum(. == "CTN") >= min.n.cells, "CTN", "Unassigned"))) %>%
    inner_join(clinical_df, by = "ImageID") %>%
    group_by(.data[[patient_id]]) %>%
    summarise(across(all_of(CTN_ls), ~ ifelse(any(. == "CTN"), "CTN", "Unassigned"))) %>%
    ungroup %>%
    mutate(across(all_of(CTN_ls), ~ factor(., levels = c("Unassigned", "CTN")))) %>%
    left_join(clinical_df_by_patient, by = patient_id)

  surv_object = Surv(time  = CTN_presence_by_patient[[time]],
                     event = CTN_presence_by_patient[[event]])

  CTN_keep <- CTN_presence_by_patient %>%
    select(all_of(CTN_ls)) %>%
    apply(2, table)
  CTN_keep <- as.data.frame(do.call(rbind, CTN_keep)) %>%
    filter(.data$CTN >= min.n.patients,
           .data$Unassigned >= min.n.patients) %>%
    rownames()

  if(length(CTN_keep) == 0) {return(NULL)}

  df_cox <- sapply(CTN_keep, function(col) {
    form_confounder <- ifelse(is.null(confounder), "",
                              paste(c("", confounder), collapse = "+"))
    form = as.formula(paste("surv_object ~", col, form_confounder))
    fit = coxph(form, data = CTN_presence_by_patient)
    sfit = summary(fit)
    c_index = sfit$concordance[1]
    coef_mat = sfit$coefficients

    res = data.frame(
      CTN = col,
      logHR   = coef_mat[, "coef"],
      HR      = coef_mat[, "exp(coef)"],
      pval    = coef_mat[, "Pr(>|z|)"],
      c_index = c_index,
      stringsAsFactors = FALSE
    )[1, ]

    return(res)
  }, simplify = FALSE, USE.NAMES = TRUE)

  df_cox <- do.call(rbind, df_cox)
  df_cox$FDR <- p.adjust(df_cox$pval, method = "BH")
  return(df_cox)
}


