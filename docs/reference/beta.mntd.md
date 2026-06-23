# Between-community Mean Nearest Taxon Distance (betaMNTD) between communities (cells)

Computes pairwise **betaMNTD** among grid cells and returns: (i) a
betaMNTD distance matrix (larger = more phylogenetically distant at the
tip scale), and (ii) rasters summarizing mean betaMNTD per cell as well
as ordination (PCoA and NMDS) axes mapped to space.

Two variants are provided:

- **Unweighted (presence/absence)** - betaMNTD is the mean of
  nearest-taxon distances in both directions (A to B and B to A).

- **Weighted (abundances)** - betaMNTD is the expected nearest-taxon
  distance using within-cell relative abundances as sampling
  probabilities.

## Usage

``` r
beta.mntd(
  tree,
  dist,
  mntd.method = c("exclude", "root", "node"),
  abundance = FALSE,
  exclude.conspecific = FALSE,
  normalize = FALSE,
  scale01 = FALSE
)
```

## Arguments

- tree:

  A rooted phylogenetic tree of class `phylo` with branch lengths.

- dist:

  An object of class `vilma.dist` (e.g., from `points.to.raster()`),
  containing species occurrences per grid cell and a raster grid.

- mntd.method:

  Character; how within-cell MNTD is handled when `normalize = TRUE`.
  Options:

  - `"exclude"` - single-species cells excluded (Default).

  - `"root"` - singletons use root-to-tip distance.

  - `"node"` - singletons use tip-to-nearest-node branch.

- abundance:

  Logical; if TRUE, compute abundance-weighted betaMNTD.

- exclude.conspecific:

  Logical; if TRUE, conspecific matches across cells are excluded.
  Default FALSE.

- normalize:

  Logical; if TRUE, return \\betaMNTD_rel = betaMNTD / ((MNTD_A +
  MNTD_B)/2)\\.

- scale01:

  Logical; if TRUE, divide betaMNTD matrix by its maximum to rescale
  values to \\\[0,1\]\\. Default FALSE.

## Value

A `vilma.beta` list containing:

- `distribution` - species distribution

- `bMNTD` - distance object (pairwise betaMNTD)

- `rasters` - SpatRaster layers:

  - `mean.bMNTD`

  - `pcoa.1`, `pcoa.2`

  - `nmds.1`, `nmds.2`

- `pcoa.eig` - eigenvalues

- `nmds.stress` - NMDS stress

- `calculation.method`

- `algorithm` - "beta.mntd"

## Details

Let \\D_sp\\ be the patristic distance matrix among species. For two
cells:

**Unweighted betaMNTD** (symmetrized): \$\$ betaMNTD(A,B) = 0.5 \* \[
mean\_{i in S_A} min\_{j in S_B} d(i,j) + mean\_{j in S_B} min\_{i in
S_A} d(j,i) \] \$\$

**Weighted betaMNTD**: \$\$ betaMNTD_w(A,B) = 0.5 \* \[ sum p_A(i) \*
min\_{j in S_B} d(i,j) + sum p_B(j) \* min\_{i in S_A} d(j,i) \] \$\$

Ordinations (PCoA, NMDS) are computed from the resulting distance
matrix.

## Notes

- Tree must be rooted and have branch lengths.

- Only shared species are used.

- `normalize = TRUE` uses `mntd.calc`.

- `scale01 = TRUE` rescales values for visualization.

## References

Webb CO, Ackerly DD, McPeek MA, Donoghue MJ (2002) Phylogenies and
community ecology. Annual Review of Ecology and Systematics 33:475-505.

Miller ET, Farine DR, Trisos CH (2017) Phylogenetic community structure
metrics and null models. Ecology Letters 20:807-818.

## See also

[`mntd.calc`](https://oleon12.github.io/vilma/reference/mntd.calc.md),
[`beta.mpd`](https://oleon12.github.io/vilma/reference/beta.mpd.md),
[`phylosor.calc`](https://oleon12.github.io/vilma/reference/phylosor.calc.md),
[`unifrac.calc`](https://oleon12.github.io/vilma/reference/unifrac.calc.md)

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com>

## Examples

``` r
if (FALSE) { # \dontrun{
tree <- examplae_tree()
dist <- example_dist()
dist <- points_to_raster(dist)
# Unweighted betaMNTD
bn <- beta.mntd(tree, dist, abundance = FALSE)
print(bn)
plot(bn)

# Weighted betaMNTD
bn_w <- beta.mntd(tree, dist, abundance = TRUE)

# Relative betaMNTD
bn_rel <- beta.mntd(tree, dist, abundance = FALSE, normalize = TRUE)

# Scaled to 0-1
bn_scaled <- beta.mntd(tree, dist, abundance = TRUE, scale01 = TRUE)
} # }
```
