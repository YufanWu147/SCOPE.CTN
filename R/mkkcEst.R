#' Interal function for standardizing each data matrix X_m
#' @keywords internal
StandardizeData <- function(X, center = TRUE, scale = TRUE) {
  assert_that(all(is.finite(X)), msg="X contains Inf/-Inf")
  assert_that(all(is.numeric(X)), msg="X must be numeric")
  assert_that(noNA(X), msg="X contains NA/NaN")
  assert_that(ncol(X) > 0, msg = "X must have at least one column")

  # Center
  if (center) {
    X <- scale(X, center = TRUE, scale = FALSE)
  }

  # Scale so that ||X||_F^2 = 1
  if (scale) {
    normX <- sqrt(sum(X^2))
    if (normX > 0) {
      X <- X / normX
    }
  }
  return(X)
}

#' Compute top eigenvectors (H) of the combined kernel using linear kernels
#' @keywords internal
compute_H_linear <- function(X_m_list, theta, centers) {
  X_scaled_list <- lapply(1:length(X_m_list), function(m) sqrt(theta[m]) * X_m_list[[m]])
  X_merged <- do.call(cbind, X_scaled_list)

  svd_res <- irlba(X_merged, nv = centers)
  H <- svd_res$u
}

#' Convert cluster labels to a binary cluster indicator matrix
#' @keywords internal
LabelToBinaryMat <- function(labels) {
  n <- length(labels)
  k <- length(unique(labels))
  Z <- matrix(0, n, k)
  for (i in 1:n) {
    Z[i, labels[i]] <- 1
  }
  return(Z)
}

#' Computing the within-cluster sum of squares for the linear kernels
#'
#' Computes the within-cluster sum of squares for the linear kernel associated with X
#' @keywords internal
WithinClusterSS_linear <- function(X, H) {
  # X: data matrix (n x d)
  # H: cluster indicator or embedding matrix (n x k), often the eigenvectors or cluster assignments

  # Compute ||X||_F^2
  normX2 <- sum(X^2)

  # Compute X^T H (d x k)
  XTH <- t(X) %*% H

  # Compute ||X^T H||_F^2
  normXTH2 <- sum(XTH^2)

  # Within cluster sum of squares
  wcss <- normX2 - normXTH2
  return(wcss)
}


#' Computing the between-cluster sum of squares for the linear kernels
#'
#' Computes the between-cluster sum of squares for the linear kernel associated with X
#' @keywords internal
BtwClusterSS_linear <- function(X, H) {
  # X: data matrix (n x d)
  # H: cluster indicator or embedding matrix (n x k)

  # Compute X^T H (d x k)
  XTH <- t(X) %*% H

  # Between cluster sum of squares is ||X^T H||_F^2
  bcss <- sum(XTH^2)

  return(bcss)
}

#' @keywords internal
#' Make the linear kernel matrix H orthonormal
make_orthonormal_basis <- function(H) {
  cluster_sizes <- colSums(H)
  H_orthonormal <- sweep(H, 2, sqrt(cluster_sizes), FUN = "/")
  return(H_orthonormal)
}


#' Robust multiple kernel K-means clustering on a multiview data
#'
#' An internal function adapted from MKKC::mkkcEst for performing robust multiple kernel K-means clustering on a multiview data.
#'
#' @param X_m_list A list of length \eqn{P} (default: 4) containing \eqn{P} kernel matrices.
#' @param centers The number of clusters \eqn{K}.
#' @param iter.max The maximum number of iterations allowed. The default is 100.
#' @param epsilon Convergence threshold. The default is \eqn{10^{-4}}.
#' @param theta Initial values for kernel coefficients. The default is 1/P for all views.
#' @param alpha Decay rate used for computing the exponential moving average of \eqn{\theta}. In interation \eqn{\scriptl}, \eqn{{\theta^{(\scriptl)}}^{EMA} = \alpha \theta^{(\scriptl)} + (1-\alpha) {\theta^{(\scriptl-1)}}^{EMA}, \alpha\in[0, 1]}. The default value is 0.5.
#' @return A list contains:
#' \describe{
#'   \item{cluster}{A vector of integers (from \code{1:K}) indicating the cluster to which each point is allocated.}
#'   \item{totss}{The total sum of squares.}
#'   \item{withinss}{Matrix of within-cluster sum of squares by cluster, one row per view.}
#'   \item{withinsscluster}{Vector of within-cluster sum of squares, one component per cluster.}
#'   \item{withinssview}{Vector of within-cluster sum of squares, one component per view.}
#'   \item{tot.withinss}{Total within-cluster sum of squares, i.e. \code{sum(withinsscluster)}.}
#'   \item{betweenssview}{Vector of  between-cluster sum of squares, one component per view.}
#'   \item{tot.betweenss}{The between-cluster sum of squares, i.e. \code{totss-tot.withinss}.}
#'   \item{clustercount}{The number of clusters \code{K}.}
#'   \item{coefficients}{The kernel coefficients}
#'   \item{H}{The continuous clustering assignment}
#'   \item{size}{The number of points, one component per cluster.}
#'   \item{iter}{The number of iterations.}
#'   \item{MbatchKm_WCSS_per_cluster}{Within-cluster sum of squares of each cluster when recovering cluster labels from the final H using minibatch K-means.}
#'   \item{MbatchKm_centroids}{Final cluster centroids.}
#' }
#'
#' @keywords internal
#' @import assertthat ClusterR irlba

mkkcEst_linear <- function(X_m_list, centers, iter.max = 100, epsilon = 1e-04, theta = rep(1/length(X_m_list), length(X_m_list)), alpha = 0.5) {
  # alpha: parameter for smoothing
  assert_that(centers > 0, round(centers)==centers, msg="centers should be a positive integer.")
  P <- length(X_m_list)
  assert_that(P > 0, msg="X_m_list should have at least one view.")

  # Standardize each X_m
  X_m_list <- lapply(X_m_list, StandardizeData, center = TRUE, scale = TRUE)

  # Initial H
  H <- compute_H_linear(X_m_list, theta, centers)

  iter = 0

  for (iter in 1:iter.max) {

    if (iter == 1) {
      theta0 <- theta
    } else {
      theta0 <- theta

      # Update theta using WithinClusterSS adapted for linear kernels
      wcss <- sapply(1:P, function(m) WithinClusterSS_linear(X = X_m_list[[m]], H = H))
      theta_new <- wcss / sqrt(sum(wcss * wcss))

      # add smoothing on theta:
      theta <- alpha * theta_new + (1-alpha) * theta0  # Smoothing

    }

    cat("iter", iter, "... theta", theta, "\n")

    # Update H based on the new theta
    H <- compute_H_linear(X_m_list, theta, centers)

    # Check convergence
    if (norm(theta0 - theta, "2") < epsilon & iter > 1) {
      break
    }

    if (iter == iter.max) {
      message(paste0("did not converge in ", iter, " iterations"))
    }
  }

  # Recover clusters
  Hnorm <- H / matrix(sqrt(rowSums(H^2)), nrow(H), centers, byrow = FALSE)

  km_model <- MiniBatchKmeans(
    Hnorm, clusters = centers, batch_size = 200, num_init = 5, max_iters = 500,
    init_fraction = 0.2, initializer = 'kmeans++', early_stop_iter = 10,
    verbose = TRUE)

  cluster = predict_MBatchKMeans(Hnorm, km_model$centroids, fuzzy = FALSE)

  #cluster <- res$cluster
  Z <- make_orthonormal_basis(LabelToBinaryMat(cluster))

  # Compute statistics
  # Again, WithinClusterSS and BtwClusterSS must be adapted to handle linear kernels from subsets of X
  # For total kernel: K_theta = X_merged X_merged^T (implicitly)
  # Just recompute them using the appropriate functions:

  # Within-cluster SS of combined kernel
  withinSStotal = WithinClusterSS_linear(X = do.call(cbind, lapply(1:P, function(m) sqrt(theta[m]) * X_m_list[[m]])), H = Z)

  # Within-cluster SS per view
  withinSSviews = sapply(1:P, function(m) WithinClusterSS_linear(X = X_m_list[[m]], H = Z))
  names(withinSSviews) = paste0("view", 1:P)

  # Within-cluster SS per cluster for the combined kernel
  withinSScluster = sapply(1:centers, function(cl) {
    idx <- which(Z[,cl] > 0)
    WithinClusterSS_linear(X = do.call(cbind, lapply(1:P, function(m) sqrt(theta[m]) * X_m_list[[m]][idx,,drop=FALSE])),
                    H = Z[idx, cl, drop=FALSE])
  })
  names(withinSScluster) = paste0("cluster", 1:centers)

  # Within-cluster SS per view and per cluster
  withinSS = sapply(1:P, function(m) sapply(1:centers, function(cl) {
    idx <- which(Z[,cl]>0)
    WithinClusterSS_linear(X = X_m_list[[m]][idx,,drop=FALSE], H = Z[idx,cl,drop=FALSE])
  }))
  withinSS = t(withinSS)
  colnames(withinSS) = paste0("cluster", 1:centers)
  rownames(withinSS) = paste0("view", 1:P)

  # Between-cluster SS for combined kernel
  btwSStotal = BtwClusterSS_linear(X = do.call(cbind, lapply(1:P, function(m) sqrt(theta[m]) * X_m_list[[m]])), H = Z)

  # Between-cluster SS per view
  btwSSviews = sapply(1:P, function(m) BtwClusterSS_linear(X = X_m_list[[m]], H = Z))
  names(btwSSviews) = paste0("view", 1:P)

  # Summary
  state <- list()
  state$cluster <- cluster
  state$totss <- withinSStotal + btwSStotal
  state$withinss <- withinSS
  state$withinsscluster <- withinSScluster
  state$withinssview <- withinSSviews
  state$tot.withinss <- withinSStotal
  state$betweenssview <- btwSSviews
  state$tot.betweenss <- btwSStotal
  state$clustercount <- centers
  state$coefficients <- theta
  state$H <- H
  state$iter <- iter
  state$size <- table(state$cluster)
  state$MbatchKm_WCSS_per_cluster <- km_model$WCSS_per_cluster
  state$MbatchKm_centroids <- km_model$centroids
  return(state)
}
