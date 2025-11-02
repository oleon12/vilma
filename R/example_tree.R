#' Example phylogenetic tree for demonstration
#' @description
#' Loads an example tree (Newick format) included in the package.
#'
#' @return An object of class \code{phylo}.
#' @examples
#' tree <- example_tree()
#' plot(tree)
#' @export
example_tree <- function() {
  read.tree(system.file("extdata", "eg_tree3.tre", package = "vilma"))
}
