# Mean Pairwise Distance (MPD) Phylogenetic Diversity Index

Computes the Mean Pairwise Distance (MPD) among species within spatial
cells based on a phylogenetic tree and species distribution data. The
function allows different treatments of single-species cells and can
optionally incorporate species abundances.

## Usage

``` r
mpd.calc(tree, dist, method = c("root", "node", "exclude"), abundance = FALSE)
```

## Arguments

- tree:

  A rooted phylogenetic tree of class `phylo` with branch lengths. For
  methods `"root"` and `"node"`, the tree must be rooted.

- dist:

  An object of class `vilma.dist`, typically produced by
  `points.to.raster()`, containing the species distribution data.

- method:

  Character string specifying how to handle single-species cells.
  Options are:

  - `"root"` – For cells with single species, uses distance from root to
    taxon for single-species cells.

  - `"node"` – For cells with single species, uses branch length to the
    nearest ancestral node.

  - `"exclude"` – Excludes single-species cells from calculations
    (default).

- abundance:

  Logical. If `TRUE`, weights pairwise distances by relative species
  abundances within each cell. Default is `FALSE`.

## Value

An object of class `vilma.pd`, which is a list containing:

- `distribution` – Original species distribution.

- `grid` – The reference raster grid.

- `pd.table` – Data frame with cells, species richness, and MPD values.

- `rasters` – List with abundance, richness, and MPD rasters.

- `calculation.method` – Method used for single-species cells.

- `abundance` – Logical, whether abundance-weighting was applied.

- `index` – The string `"mpd.calc"`.

## Details

MPD is calculated as the average phylogenetic distance between all
possible pairs of species in a community (cell). When
`abundance = TRUE`, MPD is weighted by the relative abundance of each
species. If a cell contains only one species, treatment depends on the
`method` argument.

## References

Webb, C.O., Ackerly, D.D., McPeek, M.A., & Donoghue, M.J. (2002).
Phylogenies and community ecology. *Annual Review of Ecology and
Systematics*, 33, 475–505.
[doi:10.1146/annurev.ecolsys.33.010802.150448](https://doi.org/10.1146/annurev.ecolsys.33.010802.150448)

Hardy, O.J. (2008). Testing the spatial phylogenetic structure of local
communities: statistical performances of different null models and test
statistics on a locally neutral community. *Journal of Ecology*, 96(5),
914–926.
[doi:10.1111/j.1365-2745.2008.01421.x](https://doi.org/10.1111/j.1365-2745.2008.01421.x)

## See also

[`faith.pd`](https://oleon12.github.io/vilma/reference/faith.pd.md),
[`pd.taxon`](https://oleon12.github.io/vilma/reference/pd.taxon.md),
[`points_to_raster`](https://oleon12.github.io/vilma/reference/points_to_raster.md)

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com/>

## Examples

``` r
if (FALSE) { # \dontrun{
tree <- examplae_tree()
dist <- example_dist()
dist <- points_to_raster(dist)

# Calculate abundance-unweighted MPD
result <- mpd.calc(tree = tree, dist = dist, method = "exclude")

print(result)
plot(result)
view.vilma(result)

# Calculate abundance-weighted MPD
result_ab <- mpd.calc(tree = tree, dist = dist, 
                      method = "exclude", abundance = TRUE)
} # }
```
