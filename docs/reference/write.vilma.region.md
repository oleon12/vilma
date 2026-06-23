# Write a `vilma.region` Object to Disk

Exports a `vilma.region` object produced by
[`vilma.regionalize`](https://oleon12.github.io/vilma/reference/vilma.regionalize.md)
by writing its tabular outputs (`cell.info` and `group.info`) and its
raster outputs (region IDs and mean within-region dissimilarity) to
disk. File names are generated from the user-provided `file` prefix.

## Usage

``` r
write.vilma.region(
  vilma.region,
  file,
  raster.format = c("tif", "grd", "img"),
  overwrite = TRUE
)
```

## Arguments

- vilma.region:

  A `vilma.region` object (see
  [`vilma.regionalize`](https://oleon12.github.io/vilma/reference/vilma.regionalize.md)).

- file:

  Character string giving the prefix for output file names. The function
  appends appropriate suffixes and extensions automatically.

- raster.format:

  Character; output format for raster layers. One of `"tif"`, `"grd"`,
  or `"img"`.

- overwrite:

  Logical; whether to overwrite existing files. Default = `TRUE`.

## Value

Invisibly returns a named list with normalized file paths:

- cell.info.csv:

  Path to the per-cell CSV file.

- group.info.csv:

  Path to the per-group CSV file.

- group.raster:

  Path to the region-ID raster.

- group.mean:

  Path to the mean within-region dissimilarity raster.

## Details

The function writes four files:

- `<file>cell_info.csv` — per-cell assignments, including the cell
  index, region (group) ID, and mean within-region dissimilarity for
  each cell.

- `<file>group_info.csv` — per-group summary, including group size and
  mean within-region dissimilarity.

- `<file>_groups.<format>` — raster of integer region (group) IDs.

- `<file>_mean_groups.<format>` — raster of mean within-region
  dissimilarity per cell (singleton groups may appear as `-1` if encoded
  that way in `vilma.region` for visualization).

Output paths are printed to the console, and the function invisibly
returns normalized (absolute) paths.

## See also

[`vilma.regionalize`](https://oleon12.github.io/vilma/reference/vilma.regionalize.md)

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com/>

## Examples

``` r
if (FALSE) { # \dontrun{
reg <- vilma.regionalize(beta)
write.vilma.region(reg, file = "results/region_export", raster.format = "tif")
} # }
```
