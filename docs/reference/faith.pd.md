# Faith's Phylogenetic Diversity (PD)

Calculates Faith’s Phylogenetic Diversity (PD) for species assemblages
across spatial cells using a rooted phylogenetic tree and species
distributions. The function supports different treatments of
single-species cells through the `method` argument.

## Usage

``` r
faith.pd(tree, dist, method = c("root", "node", "exclude"))
```

## Arguments

- tree:

  A rooted phylogenetic tree of class `phylo` with branch lengths.

- dist:

  An object of class `vilma.dist`, typically produced by
  `points.to.raster()`, containing the species distribution data.

- method:

  Character string specifying how to handle single-species cells.
  Options are:

  - `"root"` – For cells with single species, PD is calculated
    considering the root path to the species.

  - `"node"` – For cells with single species, PD considers the closest
    ancestral node length.

  - `"exclude"` – Single-species cells are excluded from the calculation
    (default).

## Value

An object of class `vilma.pd`, which is a list containing:

- `distribution` – Original species distribution.

- `grid` – The reference raster grid.

- `pd.table` – Data frame with cells, species richness, and PD values.

- `rasters` – List with abundance, richness, and PD rasters.

- `calculation.method` – Method used for single-species cells.

- `index` – The string `"faith.pd"`.

## Details

Faith’s PD is calculated as the sum of the branch lengths of the minimum
spanning path connecting all species in a cell. Species not shared
between `tree` and `dist` are ignored, and informative messages are
returned about missing species.

## References

Faith, D.P. (1992). Conservation evaluation and phylogenetic diversity.
*Biological Conservation*, 61(1), 1–10.
[doi:10.1016/0006-3207(92)91201-3](https://doi.org/10.1016/0006-3207%2892%2991201-3)

Webb, C.O., Ackerly, D.D., McPeek, M.A., & Donoghue, M.J. (2002).
Phylogenies and community ecology. *Annual Review of Ecology and
Systematics*, 33, 475–505.
[doi:10.1146/annurev.ecolsys.33.010802.150448](https://doi.org/10.1146/annurev.ecolsys.33.010802.150448)

## See also

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

# Calculate PD excluding single-species cells
fpd <- faith.pd(tree = tree, dist = dist, method = "exclude")

print(fpd)
plot(fpd)
view.vilma(fpd)
} # }
```
