# View `vilma.pd` Rasters in an Interactive Leaflet Map

Renders all rasters stored in a `vilma.pd` object as toggleable overlays
in a Leaflet map, using viridisLite color ramps and on-map legends.
Pixel values are shown on hover via leafem.

## Usage

``` r
view.vilma.pd(pd)
```

## Arguments

- pd:

  A `vilma.pd` object containing a named list `pd$rasters` of terra
  `SpatRaster` layers. The element `pd$rasters$ab.raster` is used to
  compute the initial map center.

## Value

A leaflet `htmlwidget` map.

## Details

**Verification:** The function checks that `pd` is of class `"vilma.pd"`
and that `pd$rasters` exists and is non-empty.

**Centering:** The initial view is centered on the midpoint of the
extent of `pd$rasters$ab.raster`.

**Color mapping:** For each raster, a numeric color function is created
with `viridis(256, option = "D")`, mapped over the raster's value range
(NAs ignored). Legends are added with the raster name as the title.

**Interaction:** Raster values are queried on mouse move using
`leafem::addImageQuery(project = TRUE)`, displaying a live readout
(rounded to 2 decimals) in the top-right corner.

**Layers:** Two base maps labeled "Esri" and "Carto" are added and a
layers control allows toggling of raster overlays. By default, all
overlays except the third (if present) are hidden, so the third raster
is visible at load.

## Performance

Very large rasters can be slow in browsers. Consider aggregating
([`terra::aggregate()`](https://rspatial.github.io/terra/reference/aggregate.html))
or tiling before viewing.

## See also

[`leaflet`](https://rstudio.github.io/leaflet/reference/leaflet.html),
[`addRasterImage`](https://rstudio.github.io/leaflet/reference/addRasterImage.html),
[`addImageQuery`](https://r-spatial.github.io/leafem/reference/addImageQuery.html),
[`viridis`](https://sjmgarnier.github.io/viridisLite/reference/viridis.html),
[`ext`](https://rspatial.github.io/terra/reference/ext.html)

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com> J. Angel
Soto-Centeno <https://www.mormoops.com>

## Examples

``` r
if (FALSE) { # \dontrun{
tree <- example_tree()
dist <- example_dist()
dist <- points_to_raster(dist)
pd <- faith.pd(tree, dist)
view.vilma.pd(pd)
} # }
```
