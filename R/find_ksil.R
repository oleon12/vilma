#' Select optimal K using silhouette width
#'
#' Computes hierarchical clustering on a phylogenetic dissimilarity matrix
#' and evaluates the mean silhouette width for a range of K values.
#' The function returns the K value that maximizes the average silhouette width.
#'
#' @param pd.dis A \code{dist} object representing pairwise dissimilarities
#'   among cells.
#' @param method Character string indicating the agglomeration method to be
#'   used in \code{hclust}. Options are \code{"ward.D2"}, \code{"average"},
#'   or \code{"complete"}. Default is \code{"ward.D2"}.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{best_k} The K value with the highest mean silhouette width.
#'     \item \code{sil_values} A data.frame with tested K values and their
#'     corresponding average silhouette widths.
#'   }
#'
#' @examples
#' \dontrun{
#' d <- dist(matrix(rnorm(100), 20, 5))
#' res <- find_ksil(d)
#' res$best_k
#' }
#'
#' @export

find_ksil <- function(pd.dis, method = c("ward.D2", "average", "complete")) {

  if (!inherits(pd.dis, "dist")) {
    stop("pd.dis must be a 'dist' object.")
  }

  method <- match.arg(method)

  hc <- hclust(pd.dis, method = method)

  k_grid <- 2:15

  sil_mean <- sapply(k_grid, function(k) {

    grp <- cutree(hc, k = k)
    mean(silhouette(grp, pd.dis)[, "sil_width"])

  })

  best_k <- k_grid[which.max(sil_mean)]

  sil_values <- data.frame(
    K.values = k_grid,
    avg.Sil = sil_mean
  )

  outSil <- list(
    best_k = best_k,
    sil_values = sil_values
  )

  return(outSil)
}

