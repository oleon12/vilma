# Write a `vilma.pd` Object to Disk

Exports a `vilma.pd` object to disk by writing its core outputs — the
phylogenetic diversity table and corresponding raster layer — along with
a text summary file generated from
[`print.vilma.pd()`](https://oleon12.github.io/vilma/reference/print.vilma.pd.md).
The function automatically names and saves these files based on the
user-provided `file` prefix.

## Usage

``` r
write.vilma.pd(
  vilma.pd,
  file,
  raster.format = c("tif", "grd", "img"),
  overwrite = TRUE
)
```

## Arguments

- vilma.pd:

  A `vilma.pd` object, typically created using one of the
  alpha-diversity functions such as
  [`faith.pd`](https://oleon12.github.io/vilma/reference/faith.pd.md),
  [`mpd.calc`](https://oleon12.github.io/vilma/reference/mpd.calc.md),
  or
  [`mntd.calc`](https://oleon12.github.io/vilma/reference/mntd.calc.md).

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

Invisibly returns a named list containing the absolute file paths of all
written outputs:

- `$pd.table` — path to the CSV table of PD values.

- `$pd.raster` — path to the raster file of PD.

- `$log` — path to the summary text file.

## Details

The function generates three files in the current working directory:

- `<file>.csv` — cell-wise phylogenetic diversity values.

- `<file>_pd_raster.<format>` — raster map of PD values.

- `<file>_log.txt` — textual summary produced by
  [`print.vilma.pd()`](https://oleon12.github.io/vilma/reference/print.vilma.pd.md).

All files are written to the current working directory, and absolute
paths are displayed in the console upon completion.

## See also

[`faith.pd`](https://oleon12.github.io/vilma/reference/faith.pd.md),
[`mpd.calc`](https://oleon12.github.io/vilma/reference/mpd.calc.md),
[`mntd.calc`](https://oleon12.github.io/vilma/reference/mntd.calc.md),
[`print.vilma.pd`](https://oleon12.github.io/vilma/reference/print.vilma.pd.md)

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com>

## Examples

``` r
if (FALSE) { # \dontrun{
data("example_vilma_pd")
write.vilma.pd(example_vilma_pd, file = "results/pd_export")
} # }

```
