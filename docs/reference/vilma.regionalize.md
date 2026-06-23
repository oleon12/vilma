# Regionalize cells into phylogenetic regions

Creates spatial/phylogenetic regions by clustering the cell-by-cell
beta-diversity (dissimilarity) matrix stored in a `vilma.beta` object.
Supports hierarchical clustering, PAM (k-medoids), and a kNN graph +
community detection approach. Additionally, it computes the mean
within-group dissimilarity (mean pairwise PD distance) for each region
and returns both a region raster and a mean-within-region raster.

## Usage

``` r
vilma.regionalize(
  beta,
  phyloBeta = c("total", "turnover", "nestedness"),
  k = 3,
  method = c("hclust", "pam", "network"),
  hmethod = c("ward.D2", "average", "complete"),
  optimize = TRUE,
  opt.method = c("sil", "diff"),
  k_thr = c("2diff", "thr"),
  thr = 0.1,
  net_kernel = c("exp", "linear"),
  net_algorithm = c("louvain", "leiden"),
  connect_target = 0.95
)
```

## Arguments

- beta:

  A `vilma.beta` object produced by `beta.mntd`, `beta.mpd`,
  `phylo.beta`, `phylosor.calc`, `rao.beta`, or `unifrac.calc`.

- phyloBeta:

  Character. Default = `"total"`. Only used when
  `beta$algorithm == "PhyloBeta"`. Which component to regionalize:
  `"total"`, `"turnover"`, or `"nestedness"`.

- k:

  Integer. Default = `3`. Number of regions (clusters). If
  `optimize = TRUE`, this value is used only as an initial/default and
  will be replaced by the optimized `k`.

- method:

  Character. Default = `"hclust"`. Regionalization method: `"hclust"`,
  `"pam"`, or `"network"`.

- hmethod:

  Character. Default = `"ward.D2"`. Linkage method for hierarchical
  clustering when `method = "hclust"`. One of `"ward.D2"`, `"average"`,
  or `"complete"`.

- optimize:

  Logical. Default = `TRUE`. If `TRUE`, attempt to select an optimal `k`
  using the selected optimization routine (method-dependent).

- opt.method:

  Character. Default = `"sil"`. Optimization criterion used for
  `method="hclust"`. One of `"sil"` (silhouette) or `"diff"`
  (dissimilarity-based).

- k_thr:

  Character. Default = `"2diff"`. Threshold rule for `opt.method="diff"`
  (passed to `find_kdiss`). One of `"2diff"` or `"thr"`.

- thr:

  Numeric. Default = `0.1`. Threshold value used when `k_thr = "thr"`
  for `find_kdiss`.

- net_kernel:

  Character. Default = `"exp"`. Kernel to transform distances into
  similarities for `method="network"`. One of `"exp"` or `"linear"`.

- net_algorithm:

  Character. Default = `"louvain"`. Community detection algorithm for
  `method="network"`. One of `"louvain"` or `"leiden"`.

- connect_target:

  Numeric. Default = `0.95`. Target connectivity used by the kNN
  optimizer `find_knnK` when `method="network"` and `optimize=TRUE`.

## Value

An object of class `vilma.region`, a list with:

- cell.info:

  Data frame with `cells` (cell index), `group` (region assignment), and
  `mean.group.pd` (mean within-region dissimilarity).

- group.info:

  Data frame with `groups`, `size` (number of cells), and `mean.pd`
  (mean within-region dissimilarity for the group).

- pd.dis:

  The `dist` object used for clustering.

- raster:

  List with `group.raster` (integer region IDs) and `mean.group.raster`
  (mean within-region dissimilarity; singleton groups shown as `-1` for
  visualization).

- pd.algorithm:

  The beta-diversity algorithm used in the input `vilma.beta` object.

- cluster.obj:

  The clustering/community object returned by the chosen method.

- cluster.param:

  Character string describing the method/parameters used.

## Details

The function extracts a dissimilarity object from `beta` depending on
`beta$algorithm`. For `"PhyloBeta"`, the user must specify which
component to use via `phyloBeta`. The extracted dissimilarity matrix is
coerced to a [`dist`](https://rdrr.io/r/stats/dist.html) object and used
for clustering.

The network-based regionalization approach follows the conceptual
framework of phylogenetic regionalization described by Daru et al.
(2017, 2020), in which pairwise phylogenetic dissimilarities are
transformed into similarity graphs and spatial regions are identified as
graph communities. However, the present implementation in
`vilma.regionalize` is an independent implementation and extends that
framework by:

- explicitly selecting the k-nearest-neighbor (kNN) parameter using
  connectivity and stability diagnostics,

- incorporating modularity-based community detection (Louvain or
  Leiden),

- providing optional silhouette- and dissimilarity-based optimization
  for hierarchical clustering,

- allowing alternative clustering paradigms (hierarchical, PAM, and
  network).

Thus, while the method is conceptually aligned with the phylogenetic
regionalization framework implemented in the `phyloregion` package, it
constitutes a modular and independently developed implementation within
the Vilma framework.

After clustering, the mean within-region dissimilarity is computed as
the mean of the upper triangle of the within-region submatrix. Regions
containing a single cell have undefined within-region mean and are
stored as `NA` in the tabular outputs. For visualization in the mean
raster, `NA` values are replaced by `-1` (with a message) so they can be
distinguished on maps.

The output rasters use `beta$rasters[[1]]` as the template and are
filled in-place using
[`set.values`](https://rspatial.github.io/terra/reference/inplace.html)
with cell indices.

## References

Daru, B.H., Karunarathne, P., & Schliep, K. (2020). *phyloregion: R
package for biogeographic regionalization and macroecology*. Methods in
Ecology and Evolution, 11, 1483–1491.
[doi:10.1111/2041-210X.13478](https://doi.org/10.1111/2041-210X.13478)

Daru, B.H., Farooq, H., Antonelli, A., & Faurby, S. (2020). *Endemism
patterns are scale dependent*. Nature Communications, 11, 2115.
[doi:10.1038/s41467-020-15921-6](https://doi.org/10.1038/s41467-020-15921-6)

Daru, B.H., Elliott, T.L., Park, D.S., & Davies, T.J. (2017).
*Understanding the processes underpinning patterns of phylogenetic
regionalization*. Trends in Ecology & Evolution, 32, 845–860.
[doi:10.1016/j.tree.2017.08.013](https://doi.org/10.1016/j.tree.2017.08.013)

## See also

[`hclust`](https://rdrr.io/r/stats/hclust.html),
[`pam`](https://rdrr.io/pkg/cluster/man/pam.html),
[`cluster_louvain`](https://r.igraph.org/reference/cluster_louvain.html),
[`cluster_leiden`](https://r.igraph.org/reference/cluster_leiden.html)

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com>

## Examples

``` r
if (FALSE) { # \dontrun{
reg1 <- vilma.regionalize(beta)
plot(reg1$raster$group.raster)

reg2 <- vilma.regionalize(beta, method = "network")
plot(reg2$raster$group.raster)

reg_pb <- vilma.regionalize(phylobeta, phyloBeta = "turnover")
} # }
```
