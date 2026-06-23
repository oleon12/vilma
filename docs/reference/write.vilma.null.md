# Write a `vilma.null` Object to Disk

Exports a `vilma.null` object to disk. For cell-based null models
(`Method == "cell"`), the function writes a CSV table of cell values, a
SES raster, and a text summary produced by
[`print.vilma.null()`](https://oleon12.github.io/vilma/reference/print.vilma.null.md).
For global null models (`Method != "cell"`), the function saves a PNG
histogram via `plot(vilma.null)` and a text summary.

## Usage

``` r
write.vilma.null(
  vilma.null,
  file,
  raster.format = c("tif", "grd", "img"),
  overwrite = TRUE
)
```

## Arguments

- vilma.null:

  A `vilma.null` object returned by one of the null model functions
  (e.g.,
  [`faith.pd.null`](https://oleon12.github.io/vilma/reference/faith.pd.null.md),
  [`mpd.calc.null`](https://oleon12.github.io/vilma/reference/mpd.calc.null.md),
  [`mntd.calc.null`](https://oleon12.github.io/vilma/reference/mntd.calc.null.md),
  [`pe.calc.null`](https://oleon12.github.io/vilma/reference/pe.calc.null.md),
  [`rao.calc.null`](https://oleon12.github.io/vilma/reference/rao.calc.null.md)).

- file:

  Character string giving the prefix for output file names. The function
  appends appropriate extensions automatically.

- raster.format:

  Character; output format for raster layers (cell mode). One of
  `"tif"`, `"grd"`, or `"img"`. Defaults to `"tif"` if multiple or
  invalid options are provided.

- overwrite:

  Logical; whether to overwrite existing files (cell mode rasters).
  Defaults to `TRUE`.

## Value

Invisibly returns a named list of absolute file paths:

- *Cell method*: `$pd.table`, `$pd.raster`, `$log`.

- *Global method*: `$png.img`, `$log`.

## Details

Output depends on `vilma.null$Method`:

**Cell method** (`Method == "cell"`):

- `<file>.csv` — cell-wise null results (e.g., SES, p-values), from
  `vilma.null$CellValues`.

- `<file>_ses_raster.<format>` — SES raster from `vilma.null$Raster`.

- `<file>_log.txt` — textual summary via
  [`print.vilma.null()`](https://oleon12.github.io/vilma/reference/print.vilma.null.md).

**Global method** (`Method != "cell"`):

- `<file>null_hist.png` — histogram saved from `plot(vilma.null)`.

- `<file>_log.txt` — textual summary via
  [`print.vilma.null()`](https://oleon12.github.io/vilma/reference/print.vilma.null.md).

All files are written to the current working directory, and absolute
paths are displayed in the console upon completion.

## See also

[`faith.pd.null`](https://oleon12.github.io/vilma/reference/faith.pd.null.md),
[`mpd.calc.null`](https://oleon12.github.io/vilma/reference/mpd.calc.null.md),
[`mntd.calc.null`](https://oleon12.github.io/vilma/reference/mntd.calc.null.md),
[`pe.calc.null`](https://oleon12.github.io/vilma/reference/pe.calc.null.md),
[`rao.calc.null`](https://oleon12.github.io/vilma/reference/rao.calc.null.md),
[`print.vilma.null`](https://oleon12.github.io/vilma/reference/print.vilma.null.md),
[`plot.vilma.null`](https://oleon12.github.io/vilma/reference/plot.vilma.null.md)

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# For a cell-based null result:
write.vilma.null(vilma_null_cell, file = "results/mpd_null_cell")

# For a global null result:
write.vilma.null(vilma_null_global, file = "results/mpd_null_global")
} # }
```
