# Phylogenetic Sørensen Similarity (PhyloSor) between communities (cells)

Computes pairwise PhyloSor similarity among grid cells and returns: (i)
a similarity matrix, (ii) its complementary dissimilarity matrix, and
(iii) rasters summarizing mean similarity/dissimilarity per cell plus
ordination (PCoA and NMDS) axes mapped to space.

Two calculation modes are supported:

- **Unweighted** (`abundance = FALSE`; default): presence–absence
  PhyloSor.

- **Weighted** (`abundance = TRUE`): edge contributions are weighted by
  per-cell species *abundances*. If `normalize = TRUE`, abundances are
  converted to relative abundances within each cell (so each cell sums
  to 1), yielding a normalized, 0–1 scaled similarity.

The PD convention for single-species cells follows `method`: `"root"`
uses root→tip distance; `"node"` uses the terminal branch; `"exclude"`
removes singletons from the analysis. For pairs sharing exactly one
species, behavior is controlled by `singleton_overlap`.

## Usage

``` r
phylosor.calc(
  tree,
  dist,
  method = c("root", "node", "exclude"),
  singleton_overlap = FALSE,
  abundance = FALSE,
  normalize = TRUE
)
```

## Arguments

- tree:

  A rooted phylogenetic tree of class `phylo` with branch lengths.

- dist:

  An object of class `vilma.dist` (e.g., from `points.to.raster()`),
  containing species occurrences per grid cell and a template raster.

- method:

  Character; PD convention for single-species cases:

  - `"root"` – singletons use root→tip; multi-species PD includes
    root→MRCA (default).

  - `"node"` – singletons use terminal branch; multi-species PD is the
    minimal subtree.

  - `"exclude"` – single-species cells are dropped before computing
    pairwise similarity.

- singleton_overlap:

  Logical; when two cells share exactly one species:

  - `FALSE` – count that single shared species as shared PD under
    `method`.

  - `TRUE` – treat that case as no shared PD (similarity = 0).

  Note: `method="exclude"` is incompatible with
  `singleton_overlap = FALSE`.

- abundance:

  Logical; if `TRUE`, use abundance-weighted PhyloSor by partitioning
  shared/unique branch lengths with abundance weights (via
  `branch.partition.weighted()`). When `FALSE` (default),
  presence–absence is used.

- normalize:

  Logical; only relevant when `abundance = TRUE`. If `TRUE` (default),
  per-cell abundances are normalized to relative abundances before
  computing edge fractions. If `FALSE`, raw counts are used.

## Value

An object of class `vilma.beta` with:

- `distribution` – Original species distribution (from `dist`).

- `similarity` – PhyloSor similarity as a `dist` object (labels = cell
  IDs).

- `dissimilarity` – `1 - similarity` as a `dist` object.

- `rasters` – A named list of `SpatRaster` layers: `mean.similarity`,
  `mean.dissimilarity`, `pcoa.1`, `pcoa.2`, `ndms.1`, `ndms.2`.

- `pcoa.eig` – PCoA eigenvalues (first two axes).

- `ndms.stress` – NMDS stress.

- `calculation.method` – The resolved `method`.

- `algorithm` – the string `"PhyloSor"`.

## Details

Let \\PD_A\\, \\PD_B\\ be (convention-consistent) phylogenetic diversity
of communities A and B, and \\PD\_{shared}\\ their shared branch length.
The classic PhyloSor similarity is: \$\$\mathrm{PhyloSor}(A,B) =
\frac{2\\PD\_{shared}}{PD_A + PD_B}.\$\$

In the unweighted mode, presence–absence is used to compute \\PD\\ and
\\PD\_{shared}\\. In the weighted mode, for each tree edge \\e\\ with
length \\L_e\\, the fractions of (relative) abundance descending from
\\e\\ in A and B are computed (\\p_A(e)\\, \\p_B(e)\\). These yield
three edge-weighted branch-length sums: shared \\a=\sum_e L_e
\min\[p_A(e),p_B(e)\]\\, unique-to-A \\b=\sum_e L_e
\max\[p_A(e)-p_B(e),0\]\\, and unique-to-B \\c=\sum_e L_e
\max\[p_B(e)-p_A(e),0\]\\. The weighted PhyloSor then uses
\\PD\_{shared}=a\\ and \\PD_A+PD_B=2a+b+c\\, so
\$\$\mathrm{PhyloSor}\_w(A,B)=\frac{2a}{2a+b+c}.\$\$ With
`normalize = TRUE`, \\p\_\cdot(e)\\ are computed from relative
abundances (summing to 1 within each cell), keeping the index in
\\\[0,1\]\\.

Only species present in both `tree` and `dist` are used; informative
messages are printed for mismatches. The similarity matrix has diagonal
1; the dissimilarity matrix is `1 - similarity` with diagonal 0.

After computing matrices, each cell is summarized by its mean similarity
(and `1 -` that value for mean dissimilarity). Ordinations (PCoA with
Cailliez correction, and NMDS) are done on the dissimilarity matrix; the
first two axes are rasterized to the study grid.

## Notes

- The tree must be rooted and have branch lengths.

- `abundance = TRUE` calls `branch.partition.weighted()` with
  `normalize` passed through; `abundance = FALSE` uses presence–absence
  via `branch.partition()`.

- `method = "exclude"` removes singletons prior to pairwise
  computations; such cells appear as `NA` in rasters.

- Ordinations are applied to dissimilarities; PCoA uses
  `ape::pcoa(..., correction = "cailliez")`.

## References

Bryant JA, Lamanna C, Morlon H, Kerkhoff AJ, Enquist BJ, Green JL (2008)
Microbes on mountainsides: contrasting elevational patterns of bacterial
and plant diversity. *Proceedings of the National Academy of Sciences*
105:11505–11511. (Introduces the PhyloSor concept/formulation.)

Kembel SW, Cowan PD, Helmus MR, Cornwell WK, Morlon H, Ackerly DD,
Blomberg SP, Webb CO (2010) Picante: R tools for integrating phylogenies
and ecology. *Bioinformatics* 26:1463–1464. (Implements `phylosor` and
related functions.)

Chiu C-H, Jost L, Chao A (2014) Phylogenetic beta diversity, similarity,
and differentiation measures based on Hill numbers. *Ecological
Monographs* 84:21–44. (Abundance-based extensions and normalization
considerations.)

Gower JC (1966) Some distance properties of latent root and vector
methods used in multivariate analysis. *Biometrika* 53:325–338.
(Classical basis of principal coordinates analysis.)

Cailliez F (1983) The analytical solution of the additive constant
problem. *Psychometrika* 48:305–312. (PCoA negative-eigenvalue
correction used here.)

Kruskal JB (1964) Nonmetric multidimensional scaling: a numerical
method. *Psychometrika* 29:115–129. (Foundational NMDS procedure.)

## See also

[`faith.pd`](https://oleon12.github.io/vilma/reference/faith.pd.md),
[`unifrac.calc`](https://oleon12.github.io/vilma/reference/unifrac.calc.md),
[`points_to_raster`](https://oleon12.github.io/vilma/reference/points_to_raster.md)

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com>

## Examples

``` r
if (FALSE) { # \dontrun{
tree <- examplae_tree()
dist <- example_dist()
dist <- points_to_raster(dist)

# Presence–absence PhyloSor (default):
ps <- phylosor.calc(tree, dist)
print(ps)
plot(ps)
view.vilma(ps)

# Abundance-weighted PhyloSor with relative abundances:
ps_w <- phylosor.calc(tree, dist, abundance = TRUE, normalize = TRUE)

# Abundance-weighted without normalization (use raw counts):
ps_wc <- phylosor.calc(tree, dist, abundance = TRUE, normalize = FALSE)
} # }

```
