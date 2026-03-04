#' Internal: central imports for vilma
#'
#' @name zzz_imports
#' @keywords internal
#'
#' @import methods
#'
#' @importFrom utils head read.csv write.csv capture.output setTxtProgressBar txtProgressBar installed.packages
#' @importFrom stats AIC Gamma as.dist as.formula binomial complete.cases cutree gaussian glm hclust inverse.gaussian lm median na.omit sd setNames predict
#' @importFrom graphics par text legend hist abline plot
#' @importFrom grDevices png dev.off
#'
#' @importFrom magrittr %>%
#'
#' @importFrom ape read.tree node.depth.edgelength pcoa getMRCA keep.tip drop.tip cophenetic.phylo is.rooted
#' @importFrom phytools getDescendants
#' @importFrom vegan monoMDS
#'
#' @importFrom terra values `values<-` ncell adjacent xyFromCell project writeRaster set.values
#' @importFrom raster extract getValues raster stack crop mask extend extent res crs xmin xmax ymin ymax rasterToPoints
#'
#' @importFrom sf st_buffer st_make_valid
#' @importFrom sp SpatialPoints CRS proj4string
#'
#' @importFrom leaflet leaflet leafletOptions addProviderTiles setView addRasterImage addLegend addLayersControl colorNumeric hideGroup addControl
#' @importFrom leafem addImageQuery
#' @importFrom viridisLite viridis
#' @importFrom htmltools htmlEscape
#'
#' @importFrom maxnet maxnet maxnet.formula
#' @importFrom dismo randomPoints evaluate threshold
#' @importFrom randomForest randomForest
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom mgcv gam s
#' @importFrom usdm vifstep exclude
#' @importFrom adehabitatHR mcp
#'
#' @importFrom runApp 
#'
#' @importFrom cluster pam silhouette
#' @importFrom igraph graph_from_adjacency_matrix components cluster_louvain cluster_leiden modularity membership E
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ
NULL

utils::globalVariables(c(".data", "sr.vals", "pd.vals"))

utils::globalVariables(c("method"))
utils::globalVariables(c("vilma"))
