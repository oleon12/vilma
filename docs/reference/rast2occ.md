# Convert binary ENM rasters to occurrence points

Converts binary ensemble niche model (ENM) rasters (values 1 and 0) into
occurrence points by extracting one point for every cell with a value
of 1. The function handles various input types and behaves differently
depending on whether `enmRasters` is a list or a single raster object.

## Usage

``` r
rast2occ(enmRasters, occ, includeSkipped = TRUE)
```

## Arguments

- enmRasters:

  Binary ENM rasters to convert. Can be one of the following:

  - A `vilma.ENM` object (output from
    [`vilma.ENM`](https://oleon12.github.io/vilma/reference/vilma.ENM.md))

  - A `RasterStack` (or `RasterBrick`) containing binary layers

  - A `SpatRaster` (from the terra package)

  - A `list` of raster objects (each element must be a `Raster*` or
    `SpatRaster`)

  **All layers or list elements must have names** – these are used as
  species identifiers. For non‑list inputs, names must match (at least
  partially) those in the `occ` data frame when `includeSkipped = TRUE`.

- occ:

  A data frame of original occurrence records with exactly three
  columns. The function expects either `Sp`, `Lon`, `Lat` or `Sp`,
  `Longitude`, `Latitude`. Coordinates must be in decimal degrees.
  **Note:** This argument is *required* as a formal parameter, but is
  **only used** when `enmRasters` is a single raster object (not a
  list). For list inputs, `occ` is ignored.

- includeSkipped:

  Logical. Default = `TRUE`. If `TRUE`, the function appends original
  occurrence records for species that are present in `occ` but absent
  from `enmRasters` (e.g., species skipped during ENM fitting). This
  option only applies when `enmRasters` is **not** a list; for list
  inputs, this parameter has no effect.

## Value

A data frame containing occurrence points. The column names depend on
the input:

- If `enmRasters` is a **list**: always columns `Sp`, `Lon`, `Lat`.

- If `enmRasters` is a **Other object**: columns are renamed to match
  the column names of `occ` (either `Sp`, `Lon`, `Lat` or `Sp`,
  `Longitude`, `Latitude`).

In both cases, the data frame contains one row per cell with value 1,
plus optionally original records for skipped species (only in the
single‑raster mode).

## Details

The function operates in two distinct modes:

1.  **List input:** When `enmRasters` is a list, the function:

    - Validates that every element is a raster object.

    - Iterates over each element, extracts points from cells with value
      1 using
      [`rasterToPoints`](https://rdrr.io/pkg/raster/man/rasterToPoints.html).

    - Returns a data frame with fixed column names `Sp`, `Lon`, `Lat`
      (these are *not* standardised to the column names of `occ`).

    - Does *not* use `occ` or `includeSkipped` – skipped species are not
      appended.

2.  **Other input:** When `enmRasters` is a `RasterStack`,
    `RasterBrick`, `SpatRaster`, or `vilma.ENM` object (the latter is
    treated as a `RasterStack`), the function:

    - Validates the column names of `occ`.

    - Identifies common species between raster layers and `occ`, and
      species present only in `occ` (if `includeSkipped = TRUE`).

    - Extracts points from each layer (each layer is treated as a
      separate species).

    - Returns a data frame with columns renamed to match the column
      names of `occ` (i.e., `Sp`, `Lon`, `Lat` or `Sp`, `Longitude`,
      `Latitude`).

    - If `includeSkipped = TRUE`, appends the original occurrence
      records for species not found in the rasters.

## Note

- The function assumes that binary rasters contain only values 0 and 1,
  where 1 indicates predicted presence.

- For large rasters with many presence cells, the output can be very
  large.

- Names are required for `enmRasters` – if missing, the function stops.

- When `enmRasters` is a list, elements may have different extents,
  projections, or origins; they are not stacked, so no alignment is
  forced.

- A progress bar is shown during point extraction for both modes.

## See also

[`vilma.ENM`](https://oleon12.github.io/vilma/reference/vilma.ENM.md)
for generating binary ENM rasters,
[`rasterToPoints`](https://rdrr.io/pkg/raster/man/rasterToPoints.html)
for the underlying extraction,
[`xyFromCell`](https://rspatial.github.io/terra/reference/xyCellFrom.html)
for alternative coordinate extraction.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: Convert a RasterStack (non‑list) with skipped species
out <- vilma.ENM(occ = occ_data, envs = envs, ...)
new_occ <- rast2occ(enmRasters = out, 
                    occ = occ_data, 
                    includeSkipped = TRUE)
head(new_occ)  # Column names match occ_data

# Example 2: Convert a list of rasters (occ is ignored)
binary_list <- list(Species_A = raster_a, Species_B = raster_b)
# occ must still be supplied (even though it is not used)
dummy_occ <- data.frame(Sp = "dummy", Lon = 0, Lat = 0)
points_list <- rast2occ(enmRasters = binary_list, 
                        occ = dummy_occ, 
                        includeSkipped = FALSE)
# Output columns are always Sp, Lon, Lat
colnames(points_list)  # "Sp" "Lon" "Lat"

# Example 3: Visualise extracted points
library(sp)
new_points <- rast2occ(out, occ_data)
coordinates(new_points) <- ~ Lon + Lat
proj4string(new_points) <- CRS("+proj=longlat +datum=WGS84")
plot(envs[[1]], main = "Original (red) vs Extracted (blue)")
points(occ_data[, c("Lon", "Lat")], col = "red", cex = 0.7)
points(new_points, col = "blue", cex = 0.3, pch = 3)
} # }
```
