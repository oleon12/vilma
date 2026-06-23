# Convert polygon ranges to occurrence records and optional raster layers

Converts polygon range maps (e.g., IUCN-style distributions) into a
standard occurrence table (`Sp`, `Lon`, `Lat`) for use in the *vilma*
workflow. Polygons are rasterized on a global template in geographic
coordinates (EPSG:4326) at a user-selected resolution. Optionally, the
function also returns the per-species rasterized layers used to generate
occurrences.

## Usage

``` r
pol2occ(pols, sp_id, res = c("30s", "1m", "2.5m", "5m"), return_raster = TRUE)
```

## Arguments

- pols:

  An `sf` object containing polygon geometries and an attribute column
  with species identifiers.

- sp_id:

  Character. Name of the column in `pols` that contains species
  names/IDs. This is used to label occurrences and to name any returned
  raster layers.

- res:

  Character. Default = `"1m"`. Output raster resolution. One of: `"30s"`
  (30 arc-seconds), `"1m"` (1 arc-minute), `"2.5m"` (2.5 arc-minutes),
  or `"5m"` (5 arc-minutes). If multiple values are provided, the
  function defaults to `"1m"`.

- return_raster:

  Logical. Default = `TRUE`. If `TRUE`, return both the occurrence table
  and the list of rasterized polygon layers. If `FALSE`, return only the
  occurrence table.

## Value

If `return_raster = TRUE`, a list with:

- occ:

  A data frame with columns `Sp`, `Lon`, `Lat` containing one occurrence
  per occupied raster cell.

- rasters:

  A named list of
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  objects (one per polygon), where list names correspond to species
  identifiers in `pols[[sp_id]]`.

If `return_raster = FALSE`, only the occurrence data frame (`Sp`, `Lon`,
`Lat`) is returned.

## Details

The function builds a global template raster spanning `-180:180`
longitude and `-90:90` latitude with CRS `EPSG:4326`. Each polygon
geometry is cropped to the template extent and rasterized to a binary
layer.

If rasterization yields an all-`NA` raster (which can occur when
polygons are very small relative to the chosen resolution), the function
iteratively increases raster resolution (halving cell size) within the
polygon extent until at least one raster cell is assigned a value, then
resamples back to the original cropped template. A warning is issued
when the resolution must be refined.

Occurrence records are extracted by converting each rasterized polygon
layer to points (one point per occupied raster cell) and combining them
into a single three-column data frame: `Sp` (species identifier), `Lon`,
and `Lat`.

## See also

[`st_geometry`](https://r-spatial.github.io/sf/reference/st_geometry.html),
[`rast`](https://rspatial.github.io/terra/reference/rast.html),
[`crop`](https://rspatial.github.io/terra/reference/crop.html),
[`rasterize`](https://rspatial.github.io/terra/reference/rasterize.html),
[`resample`](https://rspatial.github.io/terra/reference/resample.html),
[`rasterToPoints`](https://rdrr.io/pkg/raster/man/rasterToPoints.html),
[`rast2occ`](https://oleon12.github.io/vilma/reference/rast2occ.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# pols is an sf object with polygons and a species column, e.g. "species"
x <- pol2occ(pols, sp_id = "species", res = "1m", return_raster = TRUE)
head(x$occ)
names(x$rasters)

# Return only occurrences
occ <- pol2occ(pols, sp_id = "species", res = "2.5m", return_raster = FALSE)
head(occ)
} # }
```
