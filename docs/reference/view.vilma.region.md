# View `vilma.region` Rasters in an Interactive Leaflet Map

Renders all rasters stored in a `vilma.region` object as toggleable
overlays in an interactive Leaflet map. Color scales are generated with
viridisLite and per-layer legends are added automatically. Raster values
can be queried on hover via leafem.

## Usage

``` r
view.vilma.region(region)
```

## Arguments

- region:

  A `vilma.region` object containing a named list `region$raster` of
  terra `SpatRaster` layers (e.g., `group.raster` and
  `mean.group.raster`). The function uses `region$raster$group.raster`
  to compute the initial map center when available; otherwise, it falls
  back to the first valid raster layer.

## Value

A leaflet `htmlwidget` map.

## Details

**Verification:** The function checks that `region` is of class
`"vilma.region"` and that `region$raster` exists and is non-empty.

**Centering:** The initial map view is centered on the midpoint of the
extent of `region$raster$group.raster` when available. If that layer is
missing or invalid, the first raster layer with a valid extent is used.

**Color mapping:** For each raster, a numeric color function is created
with `viridisLite::viridis(256, option = "D")` mapped over the finite
value range (NAs ignored). Legends are added with the raster name as the
title.

**Interaction:** Raster values are queried on mouse move using
`leafem::addImageQuery(project = TRUE)`, displaying a live readout
(rounded to 2 decimals) in the top-right corner.

**Layers:** Two base maps labeled "Esri" and "Carto" are added, and a
layers control allows toggling raster overlays. If multiple overlays are
available, the function hides all but the first overlay at load.

## Performance

Very large rasters can be slow in browsers. Consider aggregating
([`terra::aggregate()`](https://rspatial.github.io/terra/reference/aggregate.html))
or reducing resolution prior to viewing.

## See also

[`leaflet`](https://rstudio.github.io/leaflet/reference/leaflet.html),
[`addRasterImage`](https://rstudio.github.io/leaflet/reference/addRasterImage.html),
[`addImageQuery`](https://r-spatial.github.io/leafem/reference/addImageQuery.html),
[`viridis`](https://sjmgarnier.github.io/viridisLite/reference/viridis.html),
[`ext`](https://rspatial.github.io/terra/reference/ext.html),
[`vilma.regionalize`](https://oleon12.github.io/vilma/reference/vilma.regionalize.md)

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com> J. Angel
Soto-Centeno <https://www.mormoops.com>

## Examples

``` r
if (FALSE) { # \dontrun{
reg <- vilma.regionalize(beta)
view.vilma.region(reg)
} # }
```
