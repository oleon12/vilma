# Select optimal k for PAM clustering using silhouette width

Evaluates a range of cluster numbers (k) for Partitioning Around Medoids
(PAM) clustering based on average silhouette width and returns the
optimal number of clusters.

## Usage

``` r
find_kpam(pd.dis)
```

## Arguments

- pd.dis:

  A dissimilarity object of class
  [`dist`](https://rdrr.io/r/stats/dist.html). This must represent a
  valid pairwise distance matrix.

## Value

A list with:

- best_k:

  Integer. The number of clusters that maximizes the average silhouette
  width.

- sil_val:

  Data frame with columns: `K` (evaluated cluster numbers) and `avg_sil`
  (mean silhouette width for each k).

## Details

The function evaluates cluster numbers in the range `2:15` (default
fixed range) using [`pam`](https://rdrr.io/pkg/cluster/man/pam.html)
with `diss = TRUE`. For each k, the mean silhouette width across
clusters is computed.

The optimal k is selected as the value that maximizes the average
silhouette width.

## See also

[`pam`](https://rdrr.io/pkg/cluster/man/pam.html),
[`silhouette`](https://rdrr.io/pkg/cluster/man/silhouette.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# pd.dis must be a 'dist' object
res <- find_kpam(pd.dis)
res$best_k
plot(res$sil_val$K, res$sil_val$avg_sil, type = "b")
} # }
```
