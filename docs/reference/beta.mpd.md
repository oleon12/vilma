# Between-community Mean Pairwise Distance (betaMPD) between communities (cells)

Computes pairwise between-community MPD among grid cells and returns:
(i) a betaMPD matrix (larger = more phylogenetically distant), and (ii)
rasters summarizing mean betaMPD per cell as well as ordination (PCoA
and NMDS) axes mapped to space.

Two variants are provided:

- Unweighted (presence/absence) - betaMPD is the mean patristic distance
  over all cross-cell species pairs.

- Weighted (abundances) - betaMPD is the expectation of the patristic
  distance under within-cell relative abundance weights.

Optionally, the result can be normalized by the mean within-cell MPDs
and/or rescaled to the unit interval for visualization.

## Usage

``` r
beta.mpd(
  tree,
  dist,
  mpd.method = c("exclude", "root", "node"),
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
  containing species occurrences per grid cell and a template raster.

- mpd.method:

  Character; how within-cell MPD is handled when `normalize = TRUE`. One
  of:

  - `"exclude"` - single-species cells are excluded (Default).

  - `"root"` - singletons use root-to-tip distance.

  - `"node"` - singletons use tip-to-nearest-node branch.

  Ignored if `normalize = FALSE`.

- abundance:

  Logical; if `TRUE`, compute abundance-weighted betaMPD using relative
  abundances within each cell. Default `FALSE`.

- exclude.conspecific:

  Logical; if `TRUE`, conspecific matches across cells are excluded from
  betaMPD (and abundance weights are renormalized in the weighted case).
  Default `FALSE`.

- normalize:

  Logical; if `TRUE`, return the relative (dimensionless) index
  \\betaMPD_rel = betaMPD / ((MPD_A + MPD_B)/2)\\, where within-cell
  MPDs follow `mpd.method` and `abundance`. Default `FALSE`.

- scale01:

  Logical; if `TRUE`, divide the final betaMPD matrix by its maximum
  (over all finite entries) to rescale values to \\\[0,1\]\\ for
  visualization. This is a display transform; it does not change the
  underlying metric. Default `FALSE`.

## Value

An object of class `vilma.beta`, a list with:

- `distribution` - Original species distribution (from `dist`).

- `bMPD` - Pairwise betaMPD as a `dist` object (labels are cell IDs).

- `rasters` - A named list of `SpatRaster` layers:

  - `mean.bMPD` - Raster of mean betaMPD per cell (diagonal excluded).

  - `pcoa.1`, `pcoa.2` - PCoA axes 1 and 2 mapped to cells (from
    `bMPD`).

  - `nmds.1`, `nmds.2` - NMDS axes 1 and 2 mapped to cells (from
    `bMPD`).

- `pcoa.eig` - PCoA eigen values of the first two axes.

- `nmds.stress` - NMDS stress value.

- `calculation.method` - The resolved method used (see Notes).

- `algorithm` - The string `"beta.mpd"`.

## Details

Let \\D\_{sp}\\ be the patristic (cophenetic) distance matrix among
species. For two cells with species sets \\S_A\\ and \\S_B\\:

Unweighted betaMPD: \$\$betaMPD(A,B) = mean\\ d(i,j) : i \in S_A,\\ j
\in S_B \\.\$\$

Weighted betaMPD (relative abundances \\p_A\\, \\p_B\\):
\$\$betaMPD_w(A,B) = \sum\_{i \in S} \sum\_{j \in S} p_A(i)\\ p_B(j)\\
d(i,j) = \mathbf{p}\_A^\top D\_{sp}\\ \mathbf{p}\_B.\$\$

If `exclude.conspecific = TRUE`, conspecific pairs are removed; in the
weighted case, weights are renormalized within each cell. Comparisons
involving an empty cell return `NA`. The diagonal of the returned
betaMPD matrix is set to 0 (purely between-community distance).

When `normalize = TRUE`, each pairwise value is divided by the mean of
the corresponding within-cell MPDs, computed via `mpd.calc` with the
provided `mpd.method` and `abundance` so that numerator and denominator
are consistent. Note that `normalize` produces a dimensionless ratio
(not bounded to \\\[0,1\]\\). If `scale01 = TRUE`, the final matrix
(after any normalization) is divided by its maximum for plotting.

## Notes

- The tree must be rooted and have branch lengths.

- Only species present in both `tree` and `dist` are used; informative
  messages are printed for mismatches.

- `normalize = TRUE` uses `mpd.calc` internally with `mpd.method` and
  `abundance` to compute within-cell MPDs for the denominator.

- `scale01 = TRUE` rescales the final matrix by its maximum for
  visualization; values then lie in \\\[0,1\]\\ but are not comparable
  across datasets/trees.

## References

Webb CO, Ackerly DD, McPeek MA, Donoghue MJ (2002). Phylogenies and
community ecology. Annual Review of Ecology and Systematics, 33,
475-505.

Miller ET, Farine DR, Trisos CH (2017). Phylogenetic community structure
metrics and null models: a review with new methods. Ecology Letters,
20(7), 807-818.

## See also

[`mpd.calc`](https://oleon12.github.io/vilma/reference/mpd.calc.md),
[`beta.mntd`](https://oleon12.github.io/vilma/reference/beta.mntd.md),
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
# Unweighted betaMPD (raw, in tree units)
bm <- beta.mpd(tree, dist, abundance = FALSE)
print(bm)
plot(bm)
view.vilma(bm)

# Weighted betaMPD (relative abundances)
bm_w <- beta.mpd(tree, dist, abundance = TRUE)

# Relative (normalized) betaMPD
bm_rel <- beta.mpd(tree, dist, abundance = FALSE,
                   normalize = TRUE, mpd.method = "exclude")

# Rescale to 0-1 for visualization
bm_scaled <- beta.mpd(tree, dist, abundance = TRUE,
                      normalize = FALSE, scale01 = TRUE)
} # }
```
