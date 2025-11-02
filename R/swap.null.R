#' Randomly Swap Cells in a Presence-Absence Matrix
#' @description
#' Performs null-model swaps on a presence-absence matrix by randomly exchanging 
#' species occurrences between cells. Only swaps that preserve row and column totals 
#' (checkerboard swaps) are allowed.
#'
#' @param mat A numeric matrix of 0s and 1s representing species presence-absence.
#' @param n_swaps Number of attempted swaps to perform. 
#'   Defaults to \code{5 * nrow(mat) * ncol(mat)}.
#'
#' @details
#' This function is commonly used in null model analyses for community ecology.
#' It randomly selects two rows and two columns, checks for a checkerboard pattern,
#' and swaps the occurrences to preserve species richness and cell totals.  
#' Only 2x2 submatrices matching the checkerboard patterns
#' \code{[1, 0; 0, 1]} or \code{[0, 1; 1, 0]} are swapped.
#'
#' @return A matrix of the same dimensions as \code{mat}, with swapped entries.
#'
#' @author 
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @export
swap.null <- function(mat, n_swaps = NULL) {
  mat <- as.matrix(mat)
  nrow <- nrow(mat)
  ncol <- ncol(mat)
  
  if (nrow < 2 || ncol < 2) {
    stop("Matrix must have at least 2 rows and 2 columns for swaps.")
  }
  
  if (is.null(n_swaps)) {
    n_swaps <- 5 * nrow * ncol
  }
  
  for (i in seq_len(n_swaps)) {
    rows <- sample(seq_len(nrow), 2, replace = FALSE)
    cols <- sample(seq_len(ncol), 2, replace = FALSE)
    
    sub <- mat[rows, cols]
    
    if (all(sub == c(1,0,0,1)) || all(sub == c(0,1,1,0))) {
      mat[rows, cols] <- 1 - sub
    }
  }
  return(mat)
}

