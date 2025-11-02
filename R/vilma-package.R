#' vilma: Spatial Phylogenetic Diversity
#' @description
#' Functions for PD (Faith), MPD, MNTD, PE, Rao, beta (UniFrac, PhyloSor, phylo.beta),
#' and null models. Includes a Shiny app in inst/app.
#'
#' @keywords internal
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  op <- options()
  op.vilma <- list(
    vilma.n.cores = 1L,
    vilma.seed    = 123L
  )
  toset <- !(names(op.vilma) %in% names(op))
  if (any(toset)) options(op.vilma[toset])
  invisible()
}

