#' Multi-view cell type triad niche screening
#'
#' Cell type triad niche (CTN) detection using SCOPE.
#'
#' @param combination_table Data frame containing all cell type combinations for CTN screening. Each row is a combination.
#' @param cell_neighbor_table \eqn{N\times T} neighborhood cell type proportion table. Column names must contain all cell type labels in \code{combination_table}.
#' @param dat_in Input data frame with \eqn{N} rows. Columns must contain \code{CellID} and \code{Celltype}.
#' @param clusternum The number of clusters \eqn{K}. The default value is 2.
#' @param min.prop Minimum neighborhood cell type proportion threshold. Proportions below this value will all be set to 0. The default value is 0.3.
#' @param linear Whether to run linear version of SCOPE (i.e. not applying thin plate regression spline smoothing). Default is \code{FALSE}.
#' @param mgcv_df The dimension of the basis used to represent the thin plate regression spline smooth term. The default is 15. See [mgcv::s] for details. Disabled when \code{linear=TRUE}.
#' @param num_cores The number of cores to use in [pbmcapply::pbmclapply], i.e. at most how many child processes will be run simultaneously. Can only be set at 1 on Windows.
#' @param iter.max The maximum number of iterations allowed. The default is 100.
#' @param epsilon Convergence threshold. The default is \eqn{10^{-4}}.
#' @param alpha Decay rate \eqn{\alpha} used for computing the exponential moving average of weight vector \eqn{\theta}. In interation \eqn{l}, \eqn{{\theta^{(l)}}^{EMA} = \alpha \theta^{(l)} + (1-\alpha) {\theta^{(l-1)}}^{EMA}, \alpha\in[0, 1]}. The default value is 0.5.
#' @param rename_CTN Whether to rename niches based on the proportion of core cell types. The cluster with the highest proportion will be named as \code{CTN} when \code{clusternum = 2} and \code{CTN1} otherwise. The cluster with the lowest core cell type proportion will be renamed as \code{Unassigned}. Other clusters will be named as \code{CTN2, ...} accordingly when \code{clusternum > 2}.
#' @param save_results Whether to save the niche screening results to the output directory. If \code{TRUE} (default) then the results for each CTN will be saved into a separate file. If \code{FALSE} then will all CTN screening results will be combined into a single data frame and returned as function output.
#' @param output_dir Output directory. Will automatically create a new directory when the target directory doesn't exist.
#' @param suffix File name suffix.
#' @return
#' \describe{
#'   \item{res}{CTN screening results. Only returns when \code{save_results} is \code{FALSE}.}
#' }
#' @import mgcv dplyr stringr tidyr pbmcapply
#' @export

run_SCOPE <- function(
    combination_table, cell_neighbor_table, dat_in,
    clusternum = 2, min.prop = 0.3,
    linear = FALSE, mgcv_df = 15,
    iter.max = 100, epsilon = 1e-04, alpha = 0.5, num_cores = 1,
    rename_CTN = TRUE, save_results = TRUE, output_dir = getwd(), suffix = "") {

  if(!dir.exists(output_dir) & save_results) { dir.create(output_dir) }
  if(!dir.exists(file.path(output_dir, "mkkcEst")) & save_results) {
    dir.create(file.path(output_dir, "mkkcEst")) }

  dat_in <- as.data.frame(dat_in)
  rownames(dat_in) <- dat_in$CellID

  CTN_res <- pbmclapply(1:nrow(combination_table), function(i) {
    CTN_name <- paste0(combination_table[i,], collapse="_")
    core_celltypes <- unname(unlist(combination_table[i, ]))
    print(sprintf("Processing %d: %s", i, CTN_name))

    # Prepare matrices for clustering
    test_mat <- drop_na(as.data.frame(cell_neighbor_table))
    test_mat[test_mat < min.prop] <- 0
    CTN_index <- match(combination_table[i,], colnames(test_mat))

    dat1 <- test_mat[,CTN_index]
    dat2 <- test_mat[,CTN_index[c(1,2)]]
    dat3 <- test_mat[,CTN_index[c(2,3)]]
    dat4 <- test_mat[,CTN_index[c(1,3)]]

    if(any(colSums(dat1) == 0)) {return(NULL)}

    X_m_list <- list(dat1, dat2, dat3, dat4)

    if(!linear) {
      V1 <- V2 <- V3 <- NULL
      X_m_list <- lapply(X_m_list, function(dat) {
        dat <- as.data.frame(dat)
        colnames(dat) <- paste0("V", seq_len(ncol(dat)))
        if(ncol(dat) == 3) {
          sm <- smoothCon(s(V1, V2, V3, k=mgcv_df, fx=TRUE, bs='tp'), data=dat)[[1]]}
        else if(ncol(dat) == 2) {
          sm <- smoothCon(s(V1, V2, k=mgcv_df, fx=TRUE, bs='tp'), data=dat)[[1]]
        }
        else {message(sprintf("Incorrect dimensions"))}
        basis_matrix <- as.matrix(sm$X)
        rownames(basis_matrix) <- rownames(dat)
        return(basis_matrix)
      })
    }

    # Clustering
    res <- mkkcEst_linear(X_m_list, centers=clusternum, iter.max=iter.max,
                          epsilon = epsilon, alpha = alpha)
    cluster_label <- res$cluster
    names(cluster_label) <- rownames(test_mat)
    CTN_table <- cbind(dat_in[rownames(test_mat), c("CellID", "Celltype"), drop = FALSE],
                       CTN = CTN_name, label_new = cluster_label)
    if(rename_CTN) {
      CTN_table <- label_CTN(CTN_table, core_celltypes = core_celltypes,
                             clusternum = clusternum) }

    if(save_results) {
      save(res, file = file.path(output_dir, paste0("mkkcEst/mkkcEst_", str_replace_all(CTN_name, "\\/", "."), suffix, ".RData")))
      save(CTN_table, file = file.path(output_dir, paste0(str_replace_all(CTN_name, "\\/", "."), suffix, ".RData")))
    }
    else {
      return(CTN_table)
    }

   }, mc.cores = num_cores)

  if(!save_results) {
    CTN_res <- do.call(rbind, CTN_res)
    return(CTN_res)
  }

}


#' Renaming CTN cluster labels
#'
#' Internal function for renaming CTN cluster labels based on the proportion of core cell type triads
#' @keywords internal
#' @import dplyr
#' @importFrom rlang .data
label_CTN <- function(CTN_table, core_celltypes = NULL, clusternum = 2) {
  if(clusternum == 2) {cluster_levs <- c("CTN", "Unassigned")}
  else {cluster_levs <- c(paste0("CTN", 1:(clusternum-1)), "Unassigned")}

  CTN_label_df <- CTN_table %>%
    group_by(.data$label_new) %>%
    summarise(prop_corecelltype = sum(.data$Celltype %in% core_celltypes)/n()) %>%
    ungroup %>%
    arrange(desc(.data$prop_corecelltype)) %>%
    mutate(label = case_when(
      row_number() == clusternum ~ "Unassigned",
      clusternum == 2 ~ "CTN",
      TRUE ~ paste0("CTN", row_number()))) %>%
    select(-all_of("prop_corecelltype")) %>%
    mutate(label = factor(.data$label, cluster_levs))

  CTN_table_relabeled <- CTN_table %>%
    left_join(CTN_label_df, by = c("label_new")) %>%
    select(-all_of("Celltype"))

  return(CTN_table_relabeled)
}




