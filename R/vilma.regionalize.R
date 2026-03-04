#' Regionalize cells into phylogenetic regions
#'
#' Creates spatial/phylogenetic regions by clustering the cell-by-cell
#' beta-diversity (dissimilarity) matrix stored in a \code{vilma.beta} object.
#' Supports hierarchical clustering, PAM (k-medoids), and a kNN graph + community
#' detection approach. Additionally, it computes the mean within-group
#' dissimilarity (mean pairwise PD distance) for each region and returns both
#' a region raster and a mean-within-region raster.
#'
#' @param beta A \code{vilma.beta} object produced by \code{beta.mntd},
#'   \code{beta.mpd}, \code{phylo.beta}, \code{phylosor.calc},
#'   \code{rao.beta}, or \code{unifrac.calc}.
#'
#' @param phyloBeta Character. Default = \code{"total"}.
#'   Only used when \code{beta$algorithm == "PhyloBeta"}.
#'   Which component to regionalize: \code{"total"}, \code{"turnover"},
#'   or \code{"nestedness"}.
#'
#' @param k Integer. Default = \code{3}. Number of regions (clusters).
#'   If \code{optimize = TRUE}, this value is used only as an initial/default
#'   and will be replaced by the optimized \code{k}.
#'
#' @param method Character. Default = \code{"hclust"}.
#'   Regionalization method: \code{"hclust"}, \code{"pam"}, or \code{"network"}.
#'
#' @param hmethod Character. Default = \code{"ward.D2"}.
#'   Linkage method for hierarchical clustering when
#'   \code{method = "hclust"}. One of \code{"ward.D2"}, \code{"average"},
#'   or \code{"complete"}.
#'
#' @param optimize Logical. Default = \code{TRUE}.
#'   If \code{TRUE}, attempt to select an optimal \code{k}
#'   using the selected optimization routine (method-dependent).
#'
#' @param opt.method Character. Default = \code{"sil"}.
#'   Optimization criterion used for \code{method="hclust"}.
#'   One of \code{"sil"} (silhouette) or \code{"diff"} (dissimilarity-based).
#'
#' @param k_thr Character. Default = \code{"2diff"}.
#'   Threshold rule for \code{opt.method="diff"} (passed to
#'   \code{find_kdiss}). One of \code{"2diff"} or \code{"thr"}.
#'
#' @param thr Numeric. Default = \code{0.1}.
#'   Threshold value used when \code{k_thr = "thr"} for \code{find_kdiss}.
#'
#' @param net_kernel Character. Default = \code{"exp"}.
#'   Kernel to transform distances into similarities
#'   for \code{method="network"}. One of \code{"exp"} or \code{"linear"}.
#'
#' @param net_algorithm Character. Default = \code{"louvain"}.
#'   Community detection algorithm for \code{method="network"}.
#'   One of \code{"louvain"} or \code{"leiden"}.
#'
#' @param connect_target Numeric. Default = \code{0.95}.
#'   Target connectivity used by the kNN optimizer \code{find_knnK}
#'   when \code{method="network"} and \code{optimize=TRUE}.
#'
#' @details
#' The function extracts a dissimilarity object from \code{beta} depending on
#' \code{beta$algorithm}. For \code{"PhyloBeta"}, the user must specify which
#' component to use via \code{phyloBeta}. The extracted dissimilarity matrix is
#' coerced to a \code{\link[stats]{dist}} object and used for clustering.
#'
#' The network-based regionalization approach follows the conceptual framework
#' of phylogenetic regionalization described by Daru et al. (2017, 2020), in which
#' pairwise phylogenetic dissimilarities are transformed into similarity graphs
#' and spatial regions are identified as graph communities. However, the present
#' implementation in \code{vilma.regionalize} is an independent implementation
#' and extends that framework by:
#' \itemize{
#'   \item explicitly selecting the k-nearest-neighbor (kNN) parameter using
#'   connectivity and stability diagnostics,
#'   \item incorporating modularity-based community detection (Louvain or Leiden),
#'   \item providing optional silhouette- and dissimilarity-based optimization
#'   for hierarchical clustering,
#'   \item allowing alternative clustering paradigms (hierarchical, PAM, and network).
#' }
#'
#' Thus, while the method is conceptually aligned with the phylogenetic
#' regionalization framework implemented in the \code{phyloregion} package,
#' it constitutes a modular and independently developed implementation within
#' the Vilma framework.
#'
#' After clustering, the mean within-region dissimilarity is computed as the mean
#' of the upper triangle of the within-region submatrix. Regions containing a
#' single cell have undefined within-region mean and are stored as \code{NA} in
#' the tabular outputs. For visualization in the mean raster, \code{NA} values are
#' replaced by \code{-1} (with a message) so they can be distinguished on maps.
#'
#' The output rasters use \code{beta$rasters[[1]]} as the template and are filled
#' in-place using \code{\link[terra]{set.values}} with cell indices.
#'
#' @return An object of class \code{vilma.region}, a list with:
#' \describe{
#'   \item{cell.info}{Data frame with \code{cells} (cell index), \code{group}
#'     (region assignment), and \code{mean.group.pd} (mean within-region dissimilarity).}
#'   \item{group.info}{Data frame with \code{groups}, \code{size} (number of cells),
#'     and \code{mean.pd} (mean within-region dissimilarity for the group).}
#'   \item{pd.dis}{The \code{dist} object used for clustering.}
#'   \item{raster}{List with \code{group.raster} (integer region IDs) and
#'     \code{mean.group.raster} (mean within-region dissimilarity; singleton
#'     groups shown as \code{-1} for visualization).}
#'   \item{pd.algorithm}{The beta-diversity algorithm used in the input
#'     \code{vilma.beta} object.}
#'   \item{cluster.obj}{The clustering/community object returned by the chosen method.}
#'   \item{cluster.param}{Character string describing the method/parameters used.}
#' }
#'
#' @references
#' Daru, B.H., Karunarathne, P., & Schliep, K. (2020).
#' \emph{phyloregion: R package for biogeographic regionalization and macroecology}.
#' Methods in Ecology and Evolution, 11, 1483–1491.
#' \doi{10.1111/2041-210X.13478}
#'
#' Daru, B.H., Farooq, H., Antonelli, A., & Faurby, S. (2020).
#' \emph{Endemism patterns are scale dependent}.
#' Nature Communications, 11, 2115.
#' \doi{10.1038/s41467-020-15921-6}
#'
#' Daru, B.H., Elliott, T.L., Park, D.S., & Davies, T.J. (2017).
#' \emph{Understanding the processes underpinning patterns of phylogenetic regionalization}.
#' Trends in Ecology & Evolution, 32, 845–860.
#' \doi{10.1016/j.tree.2017.08.013}
#'
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @seealso \code{\link[stats]{hclust}}, \code{\link[cluster]{pam}},
#'   \code{\link[igraph]{cluster_louvain}}, \code{\link[igraph]{cluster_leiden}}
#'
#' @examples
#' \dontrun{
#' reg1 <- vilma.regionalize(beta)
#' plot(reg1$raster$group.raster)
#'
#' reg2 <- vilma.regionalize(beta, method = "network")
#' plot(reg2$raster$group.raster)
#'
#' reg_pb <- vilma.regionalize(phylobeta, phyloBeta = "turnover")
#' }
#'
#' @export

vilma.regionalize <- function(beta, phyloBeta = c("total","turnover","nestedness"), k = 3, method = c("hclust", "pam", "network"), 
                              hmethod = c("ward.D2","average","complete"), optimize = TRUE, 
                              opt.method = c("sil","diff"), k_thr = c("2diff","thr"), thr = 0.1,
                              net_kernel = c("exp","linear"), net_algorithm = c("louvain","leiden"), connect_target = 0.95){
  
  ###################
  #  Safety checks  #
  ###################
  
  if(!inherits(beta, "vilma.beta")){
    stop("beta is not a 'vima.beta' object. Please, see functions: 'beta.mntd', 'beta.mpd', 'phylo.beta', 'phylosor.calc', 'rao.beta', 'unifrac.calc'")
  }
  
  phyloBeta <- match.arg(phyloBeta) # total defaul
  method <- match.arg(method) #hclust default
  hmethod <- match.arg(hmethod) #ward.D2 default
  opt.method <- match.arg(opt.method) #sil default
  k_thr <- match.arg(k_thr) #2diff default
  net_kernel <- match.arg(net_kernel) # exp default
  net_algorithm  <- match.arg(net_algorithm) # louvain defaul
  
  #########################
  #  Algorithm selection  #
  #########################
  
  if(beta$algorithm == "beta.MNTD"){
    pd.dis <- beta$bMNTD
  }
  
  if(beta$algorithm == "beta.MPD"){
    pd.dis <- beta$bMPD
  }
  
  if(beta$algorithm == "PhyloSor"){
    pd.dis <- beta$dissimilarity
  }
  
  if(beta$algorithm == "RaoBeta"){
    pd.dis <- beta$Rao.Beta
  }
  
  if(beta$algorithm == "UniFrac"){
    pd.dis <- beta$UniFrac
  }
  
  if(beta$algorithm == "PhyloBeta"){
    if(phyloBeta == "total"){
      pd.dis <- beta$total.dissimilarity
    }
    if(phyloBeta == "turnover"){
      pd.dis <- beta$turnover
    }
    if(phyloBeta == "nestedness"){
      pd.dis <- beta$nestedness
    }
    
  }
  
  pd.dis <- as.dist(pd.dis)
  
  ## Check ##
  
  N <- attr(pd.dis, "Size")
  if (optimize == FALSE) {
    if (!is.finite(k) || k < 2) stop("k must be >= 2.")
    if (k > N) stop("k cannot be larger than the number of cells.")
  }
  
  #####################
  #  Regionalization  #
  #####################
  
  if(method == "hclust"){

    # K-optimization #
    if(optimize == TRUE){
      if(opt.method == "sil"){
        K_eval <- find_ksil(pd.dis = pd.dis, method = hmethod)
      }
      if(opt.method == "diff"){
        K_eval <- find_kdiss(pd.dis = pd.dis, method = hmethod, k_thr = k_thr, thr = thr)
      }
      k <- K_eval$best_k
    }
    # Clustering #
    clust <- hclust(pd.dis, method = hmethod)
    grp <- cutree(clust, k = k)
    
    param.clust <- paste0("Method: ", method, ", with: ", hmethod, 
                          ", optimize = ", optimize, ", opt.method = ", opt.method)
    
  }
  
  if(method == "pam"){
    
    # optimization #
    
    if(optimize == TRUE){
      K_eval <- find_kpam(pd.dis = pd.dis)
      k <- K_eval$best_k
    }
    
    clust <- pam(x = pd.dis, k = k, diss = TRUE)
    grp <- clust$clustering
    
    param.clust <- paste0("Method: ", method, ", with: ", "optimize = ", optimize)
  }
  
  if(method == "network"){
    
    D <- as.matrix(pd.dis)
    n <- nrow(D)
    
    if(n < 2){
      stop("network regionalization needs at least 2 cells.")
    }
    
    
    # Knn optimization #
    if(optimize == TRUE){
      K_eval <- find_knnK(pd.dis = pd.dis, kernel = net_kernel, algorithm = net_algorithm, connect_target = connect_target)
      k <- K_eval$best_k
    }
    if(net_kernel == "linear"){
      
      mx <- max(D)
      
      if(!is.finite(mx) || mx <= 0){
        stop("All distances are 0 or non-finite; network kernel cannot be computed.")
      }
      
      W <- 1 - (D / max(D))
      diag(W) <- 0
    }
    if(net_kernel == "exp"){
      
      sigma <- median(D[D > 0])
      if(!is.finite(sigma) || sigma <= 0){
        stop("Cannot compute exp kernel: no positive distances.")
      }
      
      W <- exp(-D / sigma)
      diag(W) <- 0
    }
    
    k <- k
    
    k_eff <- min(k, n - 1)
    
    W_knn <- matrix(0, nrow(W), ncol(W))
    for (i in seq_len(nrow(W))) {
      ord <- order(W[i, ], decreasing = TRUE)
      nbr <- ord[ord != i][1:k_eff]
      W_knn[i, nbr] <- W[i, nbr]
    }

    # make it undirected (mutual or symmetric)
    W_knn <- pmax(W_knn, t(W_knn))
    
    ####
    g <- graph_from_adjacency_matrix(W_knn, mode = "undirected", weighted = TRUE, diag = FALSE)
    
    if(net_algorithm == "louvain"){
      clust <- cluster_louvain(g, weights = E(g)$weight)
    }
    if(net_algorithm == "leiden"){
      clust <- cluster_leiden(g, weights = E(g)$weight)
    }
    grp <- membership(clust)
    names(grp) <- rownames(D)
    
    param.clust <- paste0("Method: ", method, ", with: kernel = ", net_kernel,
                          ", algorithm = ", net_algorithm, ", optimize = ", optimize)
  }
  
  ## Optimization check ##
  
  if(optimize == TRUE){
    if (!is.finite(k) || k < 2) stop("Optimized k is < 2 (this should not happen).")
    if (k > N) stop("Optimized k is larger than the number of cells (this should not happen).")
  }

  ## Sanity Check ##
  lab <- labels(pd.dis)
  if (is.null(names(grp)) || !all(names(grp) %in% lab)) {
    stop("Cluster labels do not match distance labels; cannot compute group mean distances.")
  }
  
  uniq.grp <- unique(grp)
  
  mean.grp <- numeric(length(uniq.grp))
  names(mean.grp) <- uniq.grp
  
  mean.grpC <- rep(NA_real_, length(grp))
  names(mean.grpC) <- names(grp)
  
  sd.grp <- numeric(length(uniq.grp))
  names(sd.grp) <- uniq.grp
  
  sd.grpC <- rep(NA_real_, length(grp))
  names(sd.grpC) <- names(grp)
  
  Dmat <- as.matrix(pd.dis)
  
  for (i in seq_along(uniq.grp)) {
    
    cells_i <- names(grp)[grp == uniq.grp[i]]
    n_i <- length(cells_i)
    
    if (n_i < 2) {
      # singleton group: no pairwise distances exist
      mean_i <- NA_real_        
    } else {
      pd_grp <- Dmat[cells_i, cells_i, drop = FALSE]
      vals <- pd_grp[upper.tri(pd_grp, diag = FALSE)]
      mean_i <- mean(vals, na.rm = TRUE)
      sd_i <- sd(vals, na.rm = TRUE)
    }
    
    mean.grp[i] <- mean_i
    mean.grpC[names(grp) %in% cells_i] <- mean_i
    
    sd.grp[i] <- sd_i
    sd.grpC[names(grp) %in% cells_i] <- mean_i
    
  }
  
  out_dt <- data.frame(
    cells = as.numeric(names(grp)),
    group = as.numeric(grp),
    mean.group.pd = as.numeric(mean.grpC),
    sd.group.pd = as.numeric(sd.grpC)
  )
  
  ##############################################################
  #                       RASTER CREATION                      #
  ##############################################################
  
  grp.raster <- beta$rasters[[1]]
  values(grp.raster) <- NA
  
  set.values(
    grp.raster,
    as.numeric(out_dt$cells),
    as.numeric(out_dt$group)
  )
  
  grp.mean.raster <- beta$rasters[[1]]
  values(grp.mean.raster) <- NA
  
  mean_v <- out_dt$mean.group.pd
  na_val <- which(is.na(mean_v))
  
  if(length(na_val)>0){
    mean_v[na_val] <- c(-1)
    message(paste0("Group(s): ",out_dt$group[na_val]," with NA values. ","This means that these are groups containing only one cell; therefore, these values were changed to '-1' for visualization purposes. Please note that these are not real PD values."))
  }
  
  set.values(
    grp.mean.raster,
    as.numeric(out_dt$cells),
    as.numeric(mean_v)
  )
  
  grp_levels <- sort(unique(out_dt$group))
  group.info <- data.frame(
    groups = grp_levels,
    size = as.integer(table(out_dt$group)[as.character(grp_levels)]),
    mean.pd = as.numeric(mean.grp[as.character(grp_levels)]),
    sd.pd = as.numeric(sd.grp[as.character(grp_levels)])
  )
  ##############################################################
  #                        VILMA OBJECT                        #
  ##############################################################
  
  structure(
    list(
      cell.info = out_dt,
      group.info = group.info,
      pd.dis = pd.dis,
      raster = list(group.raster = grp.raster,
                    mean.group.raster = grp.mean.raster),
      pd.algorithm = beta$algorithm,
      cluster.obj = clust,
      cluster.param = param.clust
      ),
    class = "vilma.region"
  )

}  
