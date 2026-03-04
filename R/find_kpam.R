#' Select optimal k for PAM clustering using silhouette width
#'
#' Evaluates a range of cluster numbers (k) for Partitioning Around Medoids
#' (PAM) clustering based on average silhouette width and returns the
#' optimal number of clusters.
#'
#' @param pd.dis A dissimilarity object of class \code{\link[stats]{dist}}.
#'   This must represent a valid pairwise distance matrix.
#'
#' @details
#' The function evaluates cluster numbers in the range \code{2:15}
#' (default fixed range) using \code{\link[cluster]{pam}} with
#' \code{diss = TRUE}. For each k, the mean silhouette width across clusters
#' is computed.
#'
#' The optimal k is selected as the value that maximizes the average
#' silhouette width.
#'
#' @return A list with:
#' \describe{
#'   \item{best_k}{Integer. The number of clusters that maximizes
#'   the average silhouette width.}
#'   \item{sil_val}{Data frame with columns:
#'     \code{K} (evaluated cluster numbers) and
#'     \code{avg_sil} (mean silhouette width for each k).}
#' }
#'
#' @seealso \code{\link[cluster]{pam}}, \code{\link[cluster]{silhouette}}
#'
#' @examples
#' \dontrun{
#' # pd.dis must be a 'dist' object
#' res <- find_kpam(pd.dis)
#' res$best_k
#' plot(res$sil_val$K, res$sil_val$avg_sil, type = "b")
#' }
#'
#' @export

find_kpam <- function(pd.dis){

  if(!inherits(pd.dis,"dist")){
    stop("pd.dist must be a 'dist' object.")
  }
  
  k_range <- 2:15
  
  sil <- sapply(k_range, function(k){
                fit <- pam(pd.dis, k = k, diss = TRUE)
                mean(fit$silinfo$clus.avg.widths)
                })
  
  out <- list(
           best_k = k_range[which.max(sil)],
           sil_val = data.frame(K=k_range, avg_sil = sil)
           )
        
  return(out)

}
