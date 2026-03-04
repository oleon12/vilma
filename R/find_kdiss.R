#' Select optimal K using explained dissimilarity
#'
#' Computes hierarchical clustering on a dissimilarity matrix and evaluates
#' explained dissimilarity (R2) across a range of K values. Explained dissimilarity
#' is computed as \code{R2 = 1 - (W/T)}, where \code{W} is the sum of within-cluster
#' pairwise dissimilarities and \code{T} is the sum of all pairwise dissimilarities.
#' The function selects an optimal K using either a second-difference ("2diff")
#' knee rule or a relative threshold ("thr") based on the incremental gain in R2.
#'
#' @param pd.dis A \code{dist} object representing pairwise dissimilarities
#'   among cells.
#' @param method Character string indicating the agglomeration method to be
#'   used in \code{hclust}. Options are \code{"ward.D2"}, \code{"average"},
#'   or \code{"complete"}. Default is \code{"ward.D2"}.
#' @param k_thr Character string indicating the K-selection rule. Options are
#'   \code{"2diff"} (knee via second difference) or \code{"thr"} (relative threshold).
#'   Default is \code{"2diff"}.
#' @param thr Numeric value used only when \code{k_thr = "thr"}. Represents the
#'   proportion of the maximum incremental gain in R2 used as a cutoff.
#'   The function selects the smallest K such that
#'   \code{dR2 < thr * max(dR2)}. For example, \code{thr = 0.1} selects the
#'   first K where the gain in explained dissimilarity is less than 10% of
#'   the maximum observed gain.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{best_k} The selected K value.
#'     \item \code{R_vals} A data.frame with tested K values, \code{R2.vals}, and \code{dR2}.
#'   }
#'
#' @examples
#' \dontrun{
#' d <- dist(matrix(rnorm(100), 20, 5))
#' res <- find_kdiss(d, method = "ward.D2", k_thr = "2diff")
#' res$best_k
#' head(res$R_vals)
#' }
#'
#' @export

find_kdiss <- function(pd.dis, method = c("ward.D2","average","complete"), 
                       k_thr = c("2diff","thr"), thr = 0.1){

  if(!inherits(pd.dis,"dist")){
    stop("pd.dist must be a 'dist' object.")
  }
  
  method <- match.arg(method)
  k_thr <- match.arg(k_thr)
  
  hc <- hclust(pd.dis, method = method)
  k_grid <- 2:15
  
  d <- as.matrix(pd.dis)
  
  W_out <- c()
  R2_out <- c()
  
  for(i in seq_along(k_grid)){
  
    grp <- cutree(hc, k_grid[i])
    W <- 0
  
    for(j in unique(grp)){
  
      idx <- which(grp==j)
      
      if(length(idx) > 1){
        subd <- d[idx,idx]
        W <- W + sum(subd[upper.tri(subd)])
      }
  
    }
    
   T <- sum(d[upper.tri(d)])
   R2 <- 1 - (W/T)  
   
   W_out <-  c(W_out, W)
   R2_out <- c(R2_out, R2) 
  
  }
  
  out_m <- data.frame(K.vals = k_grid,
                      R2.vals = R2_out,
                      dR2 =c(NA, diff(R2_out)))
                      
  #return(out_m)
  
  # K_val selection
  
  if(k_thr == "thr"){
  
    cut <- thr * max(out_m$dR2, na.rm = TRUE)
    idx_thr <- which(out_m$dR2 < cut)
    best_k <- if(length(idx_thr)) out_m$K.vals[idx_thr[1]] else NA_integer_
  
  }
  
  if(k_thr == "2diff"){
  
    ddR2 <- diff(out_m$dR2)
    best_k <- out_m$K.vals[which.min(ddR2)+1]
  
  }
  
  out <- list(best_k = best_k, R_vals = out_m)
  return(out)
  
  
}
