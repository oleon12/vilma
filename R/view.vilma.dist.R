#' Interactive viewer for vilma distributions (richness & abundance)
#' @description
#' Renders a \strong{leaflet} map to explore the richness and abundance rasters
#' stored in a \code{vilma.dist} object (e.g., from \code{\link{points_to_raster}}).
#' The map centers on the raster extent, uses \code{viridis} palettes, shows
#' on-hover cell values, and provides base-layer and legend controls. By default,
#' the \emph{Richness} layer is shown and \emph{Abundance} is hidden.
#'
#' @param dist A \code{vilma.dist} object containing \code{dist$r.raster} (richness)
#'   and \code{dist$ab.raster} (abundance). The function will error if either raster
#'   is missing or \code{NULL}.
#'
#' @details
#' \itemize{
#'   \item Centers the map at the midpoint of \code{terra::ext(dist$r.raster)}.
#'   \item Builds continuous color scales with \code{viridisLite::viridis(256, option = "D")}
#'         via \code{leaflet::colorNumeric} (\code{na.color = NA}).
#'   \item Adds \code{"Esri"} and \code{"CartoDB"} provider tiles as base layers.
#'   \item Adds legends titled \emph{"richness"} and \emph{"abundance"} (bottom-left).
#'   \item Shows an on-hover value panel for each raster using
#'         \code{leafem::addImageQuery} (top-right, \code{digits = 2}).
#'   \item Uses \code{minZoom = 3}; initial zoom is 3. \emph{Abundance} layer
#'         is hidden initially.
#' }
#'
#' @return A \code{leaflet} \code{htmlwidget} map with richness and abundance overlays.
#'
#' @seealso \code{\link{points_to_raster}} to create \code{vilma.dist} objects,
#'   \code{\link{view.vilma.beta}}, \code{\link{view.vilma.null}}
#'
#' @examples
#' \dontrun{
#' # Assuming 'd' is a vilma.dist produced by points_to_raster(...)
#' view.vilma.dist(d)
#' }
#'
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @export
view.vilma.dist <- function(dist){
  
  ################################################################
  #                         Verification                         #
  ################################################################
  if (!inherits(dist, "vilma.dist")) {
    stop("dist is not a vilma.dist object")
  }
  if (is.null(dist$r.raster)) {
    stop("No richness raster found")
  }
  if (is.null(dist$ab.raster)) {
    stop("No abundance raster found")
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
  
  # Choose a valid extent for center (prefer richness, fallback to abundance)
  r_center <- NULL
  for (cand in list(dist$r.raster, dist$ab.raster)) {
    ex_try <- try(terra::ext(cand), silent = TRUE)
    if (inherits(ex_try, "SpatExtent")) { r_center <- cand; break }
  }
  if (is.null(r_center)) stop("No valid raster extent found in dist")
  
  ############## Find raster center ######################
  ex <- terra::ext(r_center)
  center.x <- as.vector((ex$xmin + ex$xmax) / 2)
  center.y <- as.vector((ex$ymin + ex$ymax) / 2)
  
  ############### Set palettes by Raster #################
  vals_r <- .finvals(dist$r.raster)
  vals_a <- .finvals(dist$ab.raster)
  
  pal_r <- if (length(vals_r)) leaflet::colorNumeric(
    palette = viridisLite::viridis(256, option = "D"),
    domain  = vals_r, na.color = NA, reverse = FALSE
  ) else NULL
  
  pal_a <- if (length(vals_a)) leaflet::colorNumeric(
    palette = viridisLite::viridis(256, option = "D"),
    domain  = vals_a, na.color = NA, reverse = FALSE
  ) else NULL
  
  ################# Base Map ############################
  map <- leaflet::leaflet(options = leaflet::leafletOptions(minZoom = 3)) %>%
    leaflet::addProviderTiles("Esri", group = "Esri") %>%
    leaflet::addProviderTiles("CartoDB", group = "Carto") %>%
    leaflet::setView(lng = center.x ,lat = center.y ,zoom = 3)
  
  ################## Add Rasters ########################
  shown_groups <- character(0)
  
  if (!is.null(pal_r)) {
    map <- map %>%
      leaflet::addRasterImage(x = dist$r.raster, group = "Richness", colors = pal_r) %>%
      leaflet::addLegend(pal = pal_r, values = vals_r,
                         position = "bottomleft", group = "Richness",
                         opacity = 1, title = "richness") %>%
      leafem::addImageQuery(x = dist$r.raster, layerId = "Richness",
                            type = "mousemove", project = TRUE, digits = 2,
                            position = "topright", prefix = "Value: ")
    shown_groups <- c(shown_groups, "Richness")
  } else {
    map <- map %>%
      leaflet::addControl(
        html = "<div style='background:#fff;padding:6px;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.2)'><b>Richness</b>: no finite values (all NA)</div>",
        position = "topright"
      )
  }
  
  if (!is.null(pal_a)) {
    map <- map %>%
      leaflet::addRasterImage(x = dist$ab.raster, group = "Abundance", colors = pal_a) %>%
      leaflet::addLegend(pal = pal_a, values = vals_a,
                         position = "bottomleft", group = "Abundance",
                         opacity = 1, title = "abundance") %>%
      leafem::addImageQuery(x = dist$ab.raster, layerId = "Abundance",
                            type = "mousemove", project = TRUE, digits = 2,
                            position = "topright", prefix = "Value: ")
    shown_groups <- c(shown_groups, "Abundance")
  } else {
    map <- map %>%
      leaflet::addControl(
        html = "<div style='background:#fff;padding:6px;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.2)'><b>Abundance</b>: no finite values (all NA)</div>",
        position = "topright"
      )
  }
  
  ################# Control Layers ######################
  if (length(shown_groups)) {
    map <- map %>%
      leaflet::addLayersControl(
        baseGroups = c("Esri","Carto"),
        overlayGroups = shown_groups
      )
    if (length(shown_groups) > 1) {
      map <- map %>% leaflet::hideGroup(setdiff(shown_groups, shown_groups[1]))
    }
  } else {
    map <- map %>%
      leaflet::addControl(
        html = "<div style='background:#fff;padding:8px;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.2)'>No layers with finite values to display.</div>",
        position = "topright"
      )
  }
  
  map
}


