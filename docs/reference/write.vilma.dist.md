# Write a `vilma.dist` Object to Disk

Exports a `vilma.dist` object to disk by writing its core outputs — the
species distribution table and the corresponding richness and abundance
raster layers — along with a text summary generated from
[`print.vilma.dist()`](https://oleon12.github.io/vilma/reference/print.vilma.dist.md).
The function automatically names and saves these files based on the
user-provided `file` prefix.

## Usage

``` r
write.vilma.dist(
  vilma.dist,
  file,
  raster.format = c("tif", "grd", "img"),
  overwrite = TRUE
)
```

## Arguments

- vilma.dist:

  A `vilma.dist` object (see
  [`points_to_raster`](https://oleon12.github.io/vilma/reference/points_to_raster.md)).

- file:

  Character string giving the prefix for output file names. The function
  appends appropriate extensions automatically.

- raster.format:

  Character; output format for raster layers. One of `"tif"`, `"grd"`,
  or `"img"`. Defaults to `"tif"` if multiple or invalid options are
  provided.

- overwrite:

  Logical; whether to overwrite existing files. Defaults to `TRUE`.

## Value

Invisibly returns a named list with absolute file paths:

- `$csv` — path to the distribution CSV table.

- `$richness` — path to the richness raster file.

- `$abundance` — path to the abundance raster file.

- `$log` — path to the summary text file.

## Details

The function generates four files in the current working directory:

- `<file>.csv` — species-by-cell distribution table.

- `<file>_richness.<format>` — raster map of species richness.

- `<file>_abundance.<format>` — raster map of abundance.

- `<file>_log.txt` — textual summary produced by
  [`print.vilma.dist()`](https://oleon12.github.io/vilma/reference/print.vilma.dist.md).

All files are written to the current working directory; absolute paths
are printed to the console upon completion.

## See also

[`points_to_raster`](https://oleon12.github.io/vilma/reference/points_to_raster.md),
[`print.vilma.dist`](https://oleon12.github.io/vilma/reference/print.vilma.dist.md),
[`faith.pd`](https://oleon12.github.io/vilma/reference/faith.pd.md),
[`write.vilma.pd`](https://oleon12.github.io/vilma/reference/write.vilma.pd.md)

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com>

## Examples

``` r
if (FALSE) { # \dontrun{
data("example_vilma_dist")
write.vilma.dist(example_vilma_dist, file = "results/dist_export")
} # }
```
