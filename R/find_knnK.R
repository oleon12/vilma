#' Select optimal k for kNN network regionalization
#'
#' Automatically selects the neighborhood size (\code{k}) for building a k-nearest
#' neighbor (kNN) similarity graph from a dissimilarity matrix and running
#' community detection (Louvain or Leiden). The selection favors values of
#' \code{k} that yield a mostly connected graph and stable community assignments
#' across consecutive k values.
#'
#' @param pd.dis A dissimilarity object of class \code{\link[stats]{dist}}.
#'   Must represent valid pairwise distances among cells.
#' @param kernel Character. Default = \code{"exp"}.
#'   Similarity kernel used to transform distances into weights:
#'   \code{"exp"} uses \eqn{w = exp(-d / sigma)} with \code{sigma = median(d[d>0])};
#'   \code{"linear"} uses \eqn{w = 1 - d / max(d)}.
#' @param algorithm Character. Default = \code{"louvain"}.
#'   Community detection algorithm: \code{"louvain"} or \code{"leiden"}.
#' @param connect_target Numeric. Default = \code{0.95}.
#'   Target fraction of nodes contained in the largest connected component (LCC).
#'   Values of \code{k} that achieve \code{lcc_frac >= connect_target} are preferred.
#'
#' @details
#' The function converts \code{pd.dis} to a dense distance matrix and builds a
#' similarity matrix \code{W} using the selected \code{kernel}. For each candidate
#' \code{k} in the default grid \code{5:25} (restricted to \code{2 <= k <= n-1}),
#' it:
#' \enumerate{
#'   \item builds a symmetric kNN adjacency matrix (weighted by \code{W}),
#'   \item runs community detection (\code{algorithm}),
#'   \item computes graph connectivity diagnostics (number of components and LCC fraction),
#'   \item records modularity of the resulting partition.
#' }
#'
#' Stability is estimated via the Adjusted Rand Index (ARI) between community
#' memberships for consecutive k values.
#'
#' Selection rule:
#' \enumerate{
#'   \item Prefer k with \code{lcc_frac >= connect_target}.
#'   \item Among those, choose the highest \code{ARI_prev} (stability).
#'   \item Break ties by higher \code{modularity}, then smaller k.
#' }
#' If no k reaches \code{connect_target}, it falls back to k values with the
#' highest \code{lcc_frac} and then applies the same stability/modularity ranking.
#'
#' Note: This procedure requires at least 10 cells (\code{n >= 10}) to operate,
#' otherwise it stops.
#'
#' @return A list with:
#' \describe{
#'   \item{best_k}{Integer. Selected neighborhood size k.}
#'   \item{diagnostics}{Data frame with one row per evaluated k, containing:
#'     \code{k}, \code{n_components}, \code{lcc_frac}, \code{modularity},
#'     and \code{ARI_prev} (ARI to the previous k).}
#' }
#'
#' @seealso \code{\link[igraph]{graph_from_adjacency_matrix}},
#'   \code{\link[igraph]{cluster_louvain}},
#'   \code{\link[igraph]{cluster_leiden}},
#'   \code{\link[mclust]{adjustedRandIndex}}
#'
#' @examples
#' \dontrun{
#' # pd.dis must be a 'dist' object
#' out <- find_knnK(pd.dis)
#' out$best_k
#' head(out$diagnostics)
#' }
#'
#' @export

find_knnK <- function(pd.dis, kernel = c("exp","linear"), algorithm = c("louvain","leiden"), connect_target = 0.95){
  
  D <- as.matrix(pd.dis)
  n <- nrow(D)
  
  if (n < 10){
    stop("Need at least 10 cells for automated k selection.")
  }
  
  diag(D) <- 0
  
  if (kernel == "linear") {
    W <- 1 - (D / max(D))
  } else {
    sigma <- median(D[D > 0])
    W <- exp(-D / sigma)
  }
  diag(W) <- 0
  
  #########
  
  # helper: kNN adjacency
  make_knn <- function(W, k) {
    n <- nrow(W)
    A <- matrix(0, n, n) 
    for (i in 1:n) {
      ord <- order(W[i, ], decreasing = TRUE)
      nbr <- ord[ord != i][1:min(k, n - 1)]
      A[i, nbr] <- W[i, nbr]
    }
    # symmetrize to undirected
    A <- pmax(A, t(A))
    A
  }
  
  # helper: run communities
  run_comm <- function(A, algorithm) {
    g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
    comps <- components(g)
    lcc_frac <- max(comps$csize) / sum(comps$csize)
    
    if (algorithm == "louvain") {
      comm <- cluster_louvain(g, weights = E(g)$weight)
      modularity_comm <- modularity(g, membership(comm), weights = E(g)$weight)
    }
    
    if (algorithm == "leiden") {
      comm <- cluster_leiden(g, objective_function = "modularity", weights = E(g)$weight)
      modularity_comm <- modularity(g, membership(comm), weights = E(g)$weight)
    }
    
    list(
      membership = membership(comm),
      modularity = modularity_comm,
      n_components = comps$no,
      lcc_frac = lcc_frac
    )
  }
  
  # evaluate each k
  k_grid <- 5:25
  k_grid <- k_grid[k_grid >= 2 & k_grid <= (n - 1)]
  res <- lapply(k_grid, function(k) {
    A <- make_knn(W, k)
    out <- run_comm(A, algorithm = algorithm)
    data.frame(
      k = k,
      n_components = out$n_components,
      lcc_frac = out$lcc_frac,
      modularity = out$modularity,
      stringsAsFactors = FALSE
    )
  })
  
  score <- do.call(rbind, res)
  
  # compute stability via ARI between consecutive k's
  memb_list <- lapply(k_grid, function(k) run_comm(make_knn(W, k), algorithm = algorithm)$membership)
  ari <- rep(NA_real_, length(k_grid))
  if (length(k_grid) >= 2) {
    for (i in 2:length(k_grid)) {
      ari[i] <- mclust::adjustedRandIndex(memb_list[[i - 1]], memb_list[[i]])
    }
  }
  score$ARI_prev <- ari
  
  # selection rule:
  # 1) prefer k where lcc_frac >= connectivity_target (mostly connected)
  # 2) among those, choose highest ARI_prev (stability), break ties by higher modularity
  ok <- which(score$lcc_frac >= connect_target)
  
  if (length(ok) == 0) {
    # fallback: choose k with best lcc_frac, then stability/modularity
    best_lcc <- max(score$lcc_frac, na.rm = TRUE)
    ok <- which(score$lcc_frac == best_lcc)
  }
  
  # rank within ok
  sub <- score[ok, , drop = FALSE]
  # Replace NA ARI with -Inf so it doesn't win
  sub$ARI_prev2 <- ifelse(is.na(sub$ARI_prev), -Inf, sub$ARI_prev)
  sub <- sub[order(-sub$ARI_prev2, -sub$modularity, sub$k), , drop = FALSE]
  
  best_k <- sub$k[1]
  
  out <- list(
    best_k = best_k,
    diagnostics = score
  )
  return(out)
}
