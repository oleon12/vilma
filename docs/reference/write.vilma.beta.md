# Write a `vilma.beta` Object to Disk

Exports a `vilma.beta` object to disk. The set of files written depends
on the beta-diversity algorithm used and stored in
`vilma.beta$algorithm`:

- **PhyloSor**: writes similarity and dissimilarity CSV matrices, a set
  of rasters (one per entry in `vilma.beta$rasters`), and a text summary
  via
  [`print.vilma.beta()`](https://oleon12.github.io/vilma/reference/print.vilma.beta.md).

- **PhyloBeta**: writes total dissimilarity, turnover, and nestedness
  CSV tables, a set of rasters (one per entry in `vilma.beta$rasters`),
  and a text summary.

- **Other algorithms** (e.g., UniFrac, Rao β): writes a single CSV named
  after the algorithm, a set of rasters (one per entry in
  `vilma.beta$rasters`), and a text summary.

## Usage

``` r
write.vilma.beta(
  vilma.beta,
  file,
  raster.format = c("tif", "grd", "img"),
  overwrite = TRUE
)
```

## Arguments

- vilma.beta:

  A `vilma.beta` object returned by one of the beta-diversity functions
  (e.g.,
  [`phylo.beta`](https://oleon12.github.io/vilma/reference/phylo.beta.md),
  [`phylosor.calc`](https://oleon12.github.io/vilma/reference/phylosor.calc.md),
  [`rao.beta`](https://oleon12.github.io/vilma/reference/rao.beta.md),
  [`unifrac.calc`](https://oleon12.github.io/vilma/reference/unifrac.calc.md)).

- file:

  Character string giving the prefix for output file names. The function
  appends appropriate extensions automatically.

- raster.format:

  Character; output format for raster layers. One of `"tif"`, `"grd"`,
  or `"img"`. Defaults to `"tif"` if multiple or invalid options are
  provided.

- overwrite:

  Logical; whether to overwrite existing raster files. Defaults to
  `TRUE`.

## Value

Invisibly returns `NULL`. Files are written to disk.

## Details

The following files are created in the current working directory:

**If `vilma.beta$algorithm == "PhyloSor"`:**

- `<file>_dissimilarity.csv` — dissimilarity matrix.

- `<file>_similarity.csv` — similarity matrix.

- `<file>_<raster_name>.<format>` — one raster per entry in
  `vilma.beta$rasters`.

- `<file>_log.txt` — textual summary from
  [`print.vilma.beta()`](https://oleon12.github.io/vilma/reference/print.vilma.beta.md).

**If `vilma.beta$algorithm == "PhyloBeta"`:**

- `<file>_dissimilarity.csv` — total dissimilarity matrix.

- `<file>_turnover.csv` — turnover component.

- `<file>_nestedness.csv` — nestedness component.

- `<file>_<raster_name>.<format>` — one raster per entry in
  `vilma.beta$rasters`.

- `<file>_log.txt` — textual summary from
  [`print.vilma.beta()`](https://oleon12.github.io/vilma/reference/print.vilma.beta.md).

**Otherwise** (e.g., UniFrac, Rao β):

- `<file>_<algorithm>.csv` — dissimilarity matrix for the algorithm.

- `<file>_<raster_name>.<format>` — one raster per entry in
  `vilma.beta$rasters`.

- `<file>_log.txt` — textual summary from
  [`print.vilma.beta()`](https://oleon12.github.io/vilma/reference/print.vilma.beta.md).

Absolute paths to the written files are printed to the console on
completion.

## See also

[`phylo.beta`](https://oleon12.github.io/vilma/reference/phylo.beta.md),
[`phylosor.calc`](https://oleon12.github.io/vilma/reference/phylosor.calc.md),
[`rao.beta`](https://oleon12.github.io/vilma/reference/rao.beta.md),
[`unifrac.calc`](https://oleon12.github.io/vilma/reference/unifrac.calc.md),
[`print.vilma.beta`](https://oleon12.github.io/vilma/reference/print.vilma.beta.md),
[`view.vilma.beta`](https://oleon12.github.io/vilma/reference/view.vilma.beta.md)

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# PhyloSor output
write.vilma.beta(vilma_beta_phylosor, file = "results/phylosor")

# PhyloBeta output
write.vilma.beta(vilma_beta_phylobeta, file = "results/phylobeta")

# Other algorithm (e.g., UniFrac / Rao beta)
write.vilma.beta(vilma_beta_unifrac, file = "results/unifrac")
} # }

```
