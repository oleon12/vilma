#' View \code{vilma.pd} Rasters in an Interactive Leaflet Map
#' @description
#' Renders all rasters stored in a \code{vilma.pd} object as toggleable
#' overlays in a Leaflet map, using \pkg{viridisLite} color ramps and
#' on-map legends. Pixel values are shown on hover via \pkg{leafem}.
#'
#' @param pd A \code{vilma.pd} object containing a named list \code{pd$rasters}
#'   of \pkg{terra} \code{SpatRaster} layers. The element
#'   \code{pd$rasters$ab.raster} is used to compute the initial map center.
#'
#' @details
#' \strong{Verification:} The function checks that \code{pd} is of class
#' \code{"vilma.pd"} and that \code{pd$rasters} exists and is non-empty.
#'
#' \strong{Centering:} The initial view is centered on the midpoint of the
#' extent of \code{pd$rasters$ab.raster}.
#'
#' \strong{Color mapping:} For each raster, a numeric color function is created
#' with \code{viridis(256, option = "D")}, mapped over the raster's value range
#' (NAs ignored). Legends are added with the raster name as the title.
#'
#' \strong{Interaction:} Raster values are queried on mouse move using
#' \code{leafem::addImageQuery(project = TRUE)}, displaying a live readout
#' (rounded to 2 decimals) in the top-right corner.
#'
#' \strong{Layers:} Two base maps labeled "Esri" and "Carto" are added and a
#' layers control allows toggling of raster overlays. By default, all overlays
#' except the third (if present) are hidden, so the third raster is visible at
#' load.
#'
#' @return A \pkg{leaflet} \code{htmlwidget} map.
#'
#' @section Performance:
#' Very large rasters can be slow in browsers. Consider aggregating
#' (\code{terra::aggregate()}) or tiling before viewing.
#'
#' @examples
#' \dontrun{
#' tree <- example_tree()
#' dist <- example_dist()
#' dist <- points_to_raster(dist)
#' pd <- faith.pd(tree, dist)
#' view.vilma.pd(pd)
#' }
#'
#' @seealso
#' \code{\link[leaflet]{leaflet}},
#' \code{\link[leaflet]{addRasterImage}},
#' \code{\link[leafem]{addImageQuery}},
#' \code{\link[viridisLite]{viridis}},
#' \code{\link[terra]{ext}}
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @export

view.vilma.pd <- function(pd){
  
  ################################################################
  #                         Verification                         #
  ################################################################
  if (!inherits(pd, "vilma.pd")) {
    stop("pd is not a vilma.pd object")
  }
  if (is.null(pd$rasters) || length(pd$rasters) == 0) {
    stop("No raster found ")
  }
  
  ################################################################
  #                     Preparing Rasters                         #
  ################################################################
  
  ############## Helpers ###############################
  .finvals <- function(r) {
    v <- try(terra::values(r), silent = TRUE)
    if (inherits(v, "try-error") || is.null(v)) return(numeric(0))
    v[is.finite(v)]
  }
  
  # Choose center: prefer ab.raster if available and valid, else first valid layer
  r_center <- NULL
  if (!is.null(pd$rasters$ab.raster)) {
    ex_try <- try(terra::ext(pd$rasters$ab.raster), silent = TRUE)
    if (inherits(ex_try, "SpatExtent")) r_center <- pd$rasters$ab.raster
  }
  if (is.null(r_center)) {
    for (r in pd$rasters) {
      ex_try <- try(terra::ext(r), silent = TRUE)
      if (inherits(ex_try, "SpatExtent")) { r_center <- r; break }
    }
  }
  if (is.null(r_center)) stop("No valid raster extent found in pd$rasters")
  
  ############## Find raster center ######################
  ex <- terra::ext(r_center)
  center.x <- as.vector((ex$xmin + ex$xmax) / 2)
  center.y <- as.vector((ex$ymin + ex$ymax) / 2)
  
  ############### Set palettes by Raster #################
  cols <- vector("list", length(pd$rasters))
  vals_list <- vector("list", length(pd$rasters))
  for (i in seq_along(pd$rasters)) {
    vals <- .finvals(pd$rasters[[i]])
    vals_list[[i]] <- vals
    cols[[i]] <- if (length(vals)) leaflet::colorNumeric(
      palette = viridisLite::viridis(256, option = "D"),
      domain  = vals, na.color = NA, reverse = FALSE
    ) else NULL
  }
  
  # Ensure names
  nm <- names(pd$rasters)
  if (is.null(nm) || any(!nzchar(nm))) {
    nm <- paste0("layer_", seq_along(pd$rasters))
    names(pd$rasters) <- nm
  }
  names(cols) <- names(pd$rasters)
  
  ################# Base Map ############################
  map <- leaflet::leaflet(options = leaflet::leafletOptions(minZoom = 3)) %>%
    leaflet::addProviderTiles("Esri",    group = "Esri")  %>%
    leaflet::addProviderTiles("CartoDB", group = "Carto") %>%
    leaflet::setView(lng = center.x, lat = center.y, zoom = 3)
  
  ################## Add Rasters ########################
  shown_groups <- character(0)
  for (i in seq_along(cols)) {
    rname <- names(cols)[i]
    pal   <- cols[[i]]
    if (is.null(pal) || length(vals_list[[i]]) == 0) {
      map <- map %>%
        leaflet::addControl(
          html = sprintf(
            "<div style='background:#fff;padding:6px;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.2)'><b>%s</b>: no finite values (all NA)</div>",
            htmltools::htmlEscape(rname)
          ),
          position = "topright"
        )
      next
    }
    
    map <- map %>%
      leaflet::addRasterImage(x = pd$rasters[[i]], group = rname, colors = pal) %>%
      leaflet::addLegend(pal = pal, values = vals_list[[i]],
                         position = "bottomleft", group = rname,
                         opacity = 1, title = rname) %>%
      leafem::addImageQuery(x = pd$rasters[[i]], layerId = rname,
                            type = "mousemove", project = TRUE, digits = 2,
                            position = "topright", prefix = "Value: ")
    shown_groups <- c(shown_groups, rname)
  }
  
  ################# Control Layers ######################
  if (length(shown_groups)) {
    map <- map %>%
      leaflet::addLayersControl(
        baseGroups    = c("Esri","Carto"),
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

