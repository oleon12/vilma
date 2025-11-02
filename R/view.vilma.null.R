#' Interactive viewer for vilma null-model results (SES raster)
#' @description
#' Renders a \strong{leaflet} map to explore the SpatialRaster of
#' standardized effect sizes (SES) produced by a \code{vilma.null} object
#' (e.g., from \code{\link{rao.calc.null}}). The map centers on the raster
#' extent, uses a \code{viridis} palette, shows a live value readout on mouse
#' move, and provides base-layer and legend controls.
#'
#' @param null A \code{vilma.null} object containing a SES raster in
#'   \code{null$Raster}. The function will error if \code{null$Raster} is
#'   missing or \code{NULL}.
#'
#' @details
#' \itemize{
#'   \item Centers the map at the midpoint of \code{terra::ext(null$Raster)}.
#'   \item Uses \code{viridisLite::viridis(256, option = "D")} via \code{leaflet::colorNumeric}
#'         for continuous SES coloring (\code{na.color = NA}).
#'   \item Adds \code{"Esri"} and \code{"CartoDB"} provider tiles as base layers.
#'   \item Displays a legend titled \emph{"SES values"} (bottom-left).
#'   \item Displays an on-hover pixel value panel using
#'         \code{leafem::addImageQuery} (top-right).
#'   \item Map starts at \code{minZoom = 3}; initial zoom is set to 3.
#' }
#'
#' @return A \code{leaflet} \code{htmlwidget} map for interactive viewing.
#'
#' @seealso \code{\link{rao.calc.null}}, \code{\link{rao.calc}},
#'   \code{\link{points_to_raster}}
#'
#' @examples
#' \dontrun{
#' tree <- example_tree()
#' dist <- example_dist()
#' dist <- points_to_raster(dist)
#' out.c <- pe.calc.null(tree, dist)
#' view.vilma.null(out.c)
#' }
#'
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @export
view.vilma.null <- function(null){
  
  ################################################################
  #                         Verification                         #
  ################################################################
  if (!inherits(null, "vilma.null")) {
    stop("null is not a vilma.null object")
  }
  if (is.null(null$Raster)) {
    stop("No raster found ")
  }
  
  ################################################################
  #                     Preparing Raster                         #
  ################################################################
  
  ############## Helpers ###############################
  .finvals <- function(r) {
    v <- try(terra::values(r), silent = TRUE)
    if (inherits(v, "try-error") || is.null(v)) return(numeric(0))
    v[is.finite(v)]
  }
  
  ex_try <- try(terra::ext(null$Raster), silent = TRUE)
  if (!inherits(ex_try, "SpatExtent")) stop("Invalid extent in null$Raster")
  
  ############## Find raster center ######################
  ex <- terra::ext(null$Raster)
  center.x <- as.vector((ex$xmin + ex$xmax) / 2)
  center.y <- as.vector((ex$ymin + ex$ymax) / 2)
  
  ############### Set palette ###########################
  vals <- .finvals(null$Raster)
  pal  <- if (length(vals)) leaflet::colorNumeric(
    palette = viridisLite::viridis(256, option = "D"),
    domain  = vals, na.color = NA, reverse = FALSE
  ) else NULL
  
  ################# Base Map ############################
  map <- leaflet::leaflet(options = leaflet::leafletOptions(minZoom = 3)) %>%
    leaflet::addProviderTiles("Esri", group = "Esri") %>%
    leaflet::addProviderTiles("CartoDB", group = "Carto") %>%
    leaflet::setView(lng = center.x ,lat = center.y ,zoom = 3)
  
  ################## Add Raster #########################
  if (!is.null(pal)) {
    map <- map %>%
      leaflet::addRasterImage(x = null$Raster, group = "Null model", colors = pal) %>%
      leaflet::addLegend(pal = pal, values = vals,
                         position = "bottomleft", group = "Null model",
                         opacity = 1, title = "SES values") %>%
      leafem::addImageQuery(x = null$Raster, layerId = "Null model",
                            type = "mousemove", project = TRUE, digits = 2,
                            position = "topright", prefix = "Value: ") %>%
      leaflet::addLayersControl(baseGroups = c("Esri","Carto"),
                                overlayGroups = "Null model")
  } else {
    map <- map %>%
      leaflet::addControl(
        html = "<div style='background:#fff;padding:8px;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.2)'>No finite SES values to display (all NA)</div>",
        position = "topright"
      ) %>%
      leaflet::addLayersControl(baseGroups = c("Esri","Carto"))
  }
  
  map
}

