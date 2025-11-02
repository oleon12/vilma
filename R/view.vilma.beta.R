#' Interactive viewer for vilma β-diversity rasters
#' @description
#' Renders a \strong{leaflet} map to explore raster layers contained in a
#' \code{vilma.beta} object (e.g., mean Rao β dissimilarity, PCoA axes,
#' nMDS axes). The map centers on the extent of the first raster, uses a
#' \code{viridis} palette per layer, shows on-hover pixel values, and
#' provides base-layer and legend controls. All overlay rasters are added
#' as separate toggleable groups; by default only the first group is shown.
#'
#' @param beta A \code{vilma.beta} object with a named list of raster layers in
#'   \code{beta$rasters} (e.g., \code{mean.dissimilarity}, \code{pcoa.1},
#'   \code{pcoa.2}, \code{ndms.1}, \code{ndms.2}). The function will error if
#'   \code{beta$rasters} is \code{NULL}.
#'
#' @details
#' \itemize{
#'   \item Centers the map at the midpoint of \code{terra::ext(beta$rasters[[1]])}.
#'   \item For each raster, builds a continuous color function with
#'         \code{viridis(256, option = "D")} via \code{leaflet::colorNumeric}
#'         (\code{na.color = NA}).
#'   \item Adds \code{"Esri"} and \code{"CartoDB"} provider tiles as base layers.
#'   \item Adds a legend titled with the raster name (bottom-left).
#'   \item Shows an on-hover value panel using \code{leafem::addImageQuery}
#'         (top-right, \code{digits = 2}).
#'   \item Uses \code{minZoom = 3}; initial zoom is 3. All rasters are added as
#'         overlay groups; groups other than the first are hidden initially.
#' }
#'
#' @return A \code{leaflet} \code{htmlwidget} with the β-diversity rasters and controls.
#'
#' @seealso \code{\link{rao.beta}} for computing Rao β dissimilarity,
#'   \code{\link{points_to_raster}} for building \code{vilma.dist}
#'
#' @examples
#' \dontrun{
#' tree <- example_tree()
#' dist <- example_dist()
#' dist <- points_to_raster(dist)
#' b <- bmpd(tree, dist)
#' view.vilma.beta(b)
#' }
#'
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @export

view.vilma.beta <- function(beta){

  ################################################################
  #                         Verification                         #
  ################################################################

  if (!inherits(beta, "vilma.beta")) {
    stop("beta is not a vilma.beta object")
  }

  if (is.null(beta$rasters) || length(beta$rasters) == 0) {
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

  # Find first raster with a valid extent for centering
  idx_ok <- which(vapply(beta$rasters, function(r) {
    inherits(try(terra::ext(r), TRUE), "SpatExtent")
  }, logical(1)))[1]
  if (is.na(idx_ok)) stop("No valid raster extent found in beta$rasters")

  ############## Find raster center ######################
  ex <- terra::ext(beta$rasters[[idx_ok]])
  center.x <- as.vector((ex$xmin + ex$xmax) / 2)
  center.y <- as.vector((ex$ymin + ex$ymax) / 2)

  ############### Set palettes by Raster #################
  cols <- vector("list", length(beta$rasters))
  vals_list <- vector("list", length(beta$rasters))
  for (i in seq_along(beta$rasters)) {
    vals <- .finvals(beta$rasters[[i]])
    vals_list[[i]] <- vals
    if (length(vals)) {
      cols[[i]] <- leaflet::colorNumeric(
        palette = viridisLite::viridis(256, option = "D"),
        domain  = vals,
        na.color = NA,
        reverse  = FALSE
      )
    } else {
      cols[[i]] <- NULL  # mark as empty; we'll skip it later
    }
  }

  # Ensure layer names exist
  nm <- names(beta$rasters)
  if (is.null(nm) || any(!nzchar(nm))) {
    nm <- paste0("layer_", seq_along(beta$rasters))
    names(beta$rasters) <- nm
  }
  names(cols) <- names(beta$rasters)

  ################# Base Map ############################
  map <- leaflet::leaflet(options = leaflet::leafletOptions(minZoom = 3)) |>
    leaflet::addProviderTiles("Esri", group = "Esri") |>
    leaflet::addProviderTiles("CartoDB", group = "Carto") |>
    leaflet::setView(lng = center.x, lat = center.y, zoom = 3)

  ################## Add Rasters ########################
  shown_groups <- character(0)
  for (i in seq_along(cols)) {
    rname <- names(cols)[i]
    pal   <- cols[[i]]
    if (is.null(pal) || length(vals_list[[i]]) == 0) {
      # Informative note for empty layers (keeps UI feedback without breaking)
      map <- map |>
        leaflet::addControl(
          html = sprintf(
            "<div style='background:#fff;padding:6px;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.2)'>
             <b>%s</b>: no finite values (all NA)</div>", htmltools::htmlEscape(rname)
          ),
          position = "topright"
        )
      next
    }

    map <- map |>
      leaflet::addRasterImage(
        x = beta$rasters[[i]],
        group = rname,
        colors = pal
      ) |>
      leaflet::addLegend(
        pal      = pal,
        values   = vals_list[[i]],
        position = "bottomleft",
        group    = rname,
        opacity  = 1,
        title    = rname
      ) |>
      leafem::addImageQuery(
        x       = beta$rasters[[i]],
        layerId = rname,
        type    = "mousemove",
        project = TRUE,
        digits  = 2,
        position= "topright",
        prefix  = "Value: "
      )

    shown_groups <- c(shown_groups, rname)
  }

  ################# Control Layers ######################
  if (length(shown_groups)) {
    map <- map |>
      leaflet::addLayersControl(
        baseGroups    = c("Esri","Carto"),
        overlayGroups = shown_groups
      )
    if (length(shown_groups) > 1) {
      map <- map |> leaflet::hideGroup(shown_groups[-1])
    }
  } else {
    map <- map |>
      leaflet::addControl(
        html =
          "<div style='background:#fff;padding:8px;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.2)'>
           No layers with finite values to display.</div>",
        position = "topright"
      )
  }

  map
}

