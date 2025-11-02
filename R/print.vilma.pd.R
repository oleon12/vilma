#' Print Summary of a VILMA Phylogenetic Diversity Object
#' @description
#' Provides a concise textual summary of a \code{vilma.pd} object, including:
#' the number of species, number of spatial cells, total records, species abundance statistics,
#' and summary statistics of phylogenetic diversity (PD) metrics.
#'
#' @param x An object of class \code{vilma.pd}.
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' This function prints a quick overview of the phylogenetic diversity stored in a VILMA object.
#' It summarizes the distribution data and PD values per cell for quick inspection.
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @author Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @export
#' @method print vilma.pd

print.vilma.pd <- function(x, ...) {
  
  # 1. Input validation
  if (!inherits(x, "vilma.pd")) {
    stop("Object must be of class 'vilma.pd'")
  }
  
  if (is.null(x$distribution)) {
    stop("'distribution' component is missing from this vilma.pd object")
  }
  
  if (is.null(x$pd.table)) {
    stop("'pd table' component is missing from this vilma.pd object")
  }
  
  # 2. Extract information
  dist_matrix <- x$distribution
  n_species <- length(unique(dist_matrix[, 1]))
  n_cells <- length(unique(dist_matrix[, 4]))
  abundance_table <- table(dist_matrix[, 1])
  
  # 3. Formatted output
  cat("\n")
  cat("vilma.pd Object Summary")
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
                min(abundance_table), max(abundance_table)))
    cat(sprintf("  - Total:    %d observations\n", sum(abundance_table)))
    cat(sprintf("  - Mean:     %.1f +/- %.1f (SD) observations per species\n",
                mean(abundance_table), sd(abundance_table)))
  } else {
    cat("  No abundance data available\n")
  }
  cat("\n")
  
  for (i in 2:length(colnames(x$pd.table))) {
    cat(paste(colnames(x$pd.table)[i], "Summary:", collapse = ","), "\n")
    cat(sprintf("  - Range:    %.2f - %.2f \n",
                min(x$pd.table[, i], na.rm = TRUE),
                max(x$pd.table[, i], na.rm = TRUE)))
    cat(sprintf("  - Mean:     %.2f +/- %.2f (SD) \n",
                mean(x$pd.table[, i], na.rm = TRUE),
                sd(x$pd.table[, i], na.rm = TRUE)))
    cat("\n")
  }
  
  cat(paste0("  - Index: ", x$index, "\n"))
  cat("\n")
}

