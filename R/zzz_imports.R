#' Internal: central imports for vilma
#'
#' @name zzz_imports
#' @keywords internal
#'
#' @importFrom utils head read.csv write.csv capture.output setTxtProgressBar txtProgressBar
#' @importFrom stats sd setNames as.dist na.omit
#' @importFrom graphics par text legend hist abline
#' @importFrom grDevices png dev.off
#'
#' @importFrom leaflet leaflet leafletOptions addProviderTiles setView addRasterImage addLegend addLayersControl colorNumeric hideGroup
#' @importFrom leafem addImageQuery
#' @importFrom viridisLite viridis
#' @importFrom shiny runApp
#'
#' @importFrom ape read.tree node.depth.edgelength pcoa getMRCA keep.tip drop.tip cophenetic.phylo is.rooted
#' @importFrom phytools getDescendants
#' @importFrom vegan monoMDS
#' @importFrom terra values 'values<-' ncell adjacent xyFromCell writeRaster set.values
#' @import methods
#' @importFrom magrittr %>%
#' @importFrom htmltools htmlEscape
NULL

utils::globalVariables(c("method"))
utils::globalVariables(c("vilma"))
