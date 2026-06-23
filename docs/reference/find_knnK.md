# Select optimal k for kNN network regionalization

Automatically selects the neighborhood size (`k`) for building a
k-nearest neighbor (kNN) similarity graph from a dissimilarity matrix
and running community detection (Louvain or Leiden). The selection
favors values of `k` that yield a mostly connected graph and stable
community assignments across consecutive k values.

## Usage

``` r
find_knnK(
  pd.dis,
  kernel = c("exp", "linear"),
  algorithm = c("louvain", "leiden"),
  connect_target = 0.95
)
```

## Arguments

- pd.dis:

  A dissimilarity object of class
  [`dist`](https://rdrr.io/r/stats/dist.html). Must represent valid
  pairwise distances among cells.

- kernel:

  Character. Default = `"exp"`. Similarity kernel used to transform
  distances into weights: `"exp"` uses \\w = exp(-d / sigma)\\ with
  `sigma = median(d[d>0])`; `"linear"` uses \\w = 1 - d / max(d)\\.

- algorithm:

  Character. Default = `"louvain"`. Community detection algorithm:
  `"louvain"` or `"leiden"`.

- connect_target:

  Numeric. Default = `0.95`. Target fraction of nodes contained in the
  largest connected component (LCC). Values of `k` that achieve
  `lcc_frac >= connect_target` are preferred.

## Value

A list with:

- best_k:

  Integer. Selected neighborhood size k.

- diagnostics:

  Data frame with one row per evaluated k, containing: `k`,
  `n_components`, `lcc_frac`, `modularity`, and `ARI_prev` (ARI to the
  previous k).

## Details

The function converts `pd.dis` to a dense distance matrix and builds a
similarity matrix `W` using the selected `kernel`. For each candidate
`k` in the default grid `5:25` (restricted to `2 <= k <= n-1`), it:

1.  builds a symmetric kNN adjacency matrix (weighted by `W`),

2.  runs community detection (`algorithm`),

3.  computes graph connectivity diagnostics (number of components and
    LCC fraction),

4.  records modularity of the resulting partition.

Stability is estimated via the Adjusted Rand Index (ARI) between
community memberships for consecutive k values.

Selection rule:

1.  Prefer k with `lcc_frac >= connect_target`.

2.  Among those, choose the highest `ARI_prev` (stability).

3.  Break ties by higher `modularity`, then smaller k.

If no k reaches `connect_target`, it falls back to k values with the
highest `lcc_frac` and then applies the same stability/modularity
ranking.

Note: This procedure requires at least 10 cells (`n >= 10`) to operate,
otherwise it stops.

## See also

[`graph_from_adjacency_matrix`](https://r.igraph.org/reference/graph_from_adjacency_matrix.html),
[`cluster_louvain`](https://r.igraph.org/reference/cluster_louvain.html),
[`cluster_leiden`](https://r.igraph.org/reference/cluster_leiden.html),
[`adjustedRandIndex`](https://mclust-org.github.io/mclust/reference/adjustedRandIndex.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# pd.dis must be a 'dist' object
out <- find_knnK(pd.dis)
out$best_k
head(out$diagnostics)
} # }
```
