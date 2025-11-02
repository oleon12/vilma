#' Example points for demonstration
#' @description
#' Loads an example points (CSV format) included in the package.
#'
#' @return An object of class \code{data.frame} or \code{matrix}.
#' @examples
#' dist <- example_dist()
#' head(dist)
#' @export
   example_dist <- function() {
     read.csv(system.file("extdata", "Sp_localities.csv", package = "vilma"))
}
