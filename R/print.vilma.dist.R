#' Print Summary of a VILMA Distribution Object
#' @description
#' Provides a concise textual summary of a \code{vilma.dist} object,
#' including the number of species, number of cells, total records, and basic species abundance statistics.
#'
#' @param x An object of class \code{vilma.dist}.
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' This function is designed to give a quick overview of the distribution data stored in a VILMA object.
#' It prints the number of unique taxa, spatial cells, total records, and summary statistics for species abundance.
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @author Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @export
#' @method print vilma.dist

print.vilma.dist <- function(x, ...) {
  
  # 1. Input validation
  if (!inherits(x, "vilma.dist")) {
    stop("Object must be of class 'vilma.dist'")
  }
  
  if (is.null(x$distribution)) {
    stop("'distribution' component is missing from this vilma.dist object")
  }
  
  # 2. Extract information for cleaner code
  dist_matrix <- x$distribution
  n_species <- length(unique(dist_matrix[, 1]))
  n_cells <- length(unique(dist_matrix[, 4]))
  abundance_table <- table(dist_matrix[, 1])
  
  # 3. Formatted output
  cat("\n")
  cat("vilma.dist Object Summary")
  cat("\n")
  cat("-------------------------")
  cat("\n\n")
  
  cat("Distribution Matrix:\n")
  cat(sprintf("  - Species: %d unique taxa\n", n_species))
  cat(sprintf("  - Cells:   %d spatial units\n", n_cells))
  cat(sprintf("  - Records: %d total observations\n", nrow(dist_matrix)))
  cat("\n")
  
  cat("Species Abundance Summary:\n")
  if (length(abundance_table) > 0) {
    cat(sprintf("  - Range:    %d - %d observations per species\n",
                min(abundance_table, na.rm = TRUE), max(abundance_table, na.rm = TRUE)))
    cat(sprintf("  - Total:    %d observations\n", sum(abundance_table)))
    cat(sprintf("  - Mean:     %.1f +/- %.1f (SD) observations per species\n",
                mean(abundance_table, na.rm = TRUE), sd(abundance_table, na.rm = TRUE)))
  } else {
    cat("  No abundance data available\n")
  }
  cat("\n")
  
  # 5. Return object invisibly
  invisible(x)
}

