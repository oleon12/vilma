# Randomly Swap Cells in a Presence-Absence Matrix

Performs null-model swaps on a presence-absence matrix by randomly
exchanging species occurrences between cells. Only swaps that preserve
row and column totals (checkerboard swaps) are allowed.

## Usage

``` r
swap.null(mat, n_swaps = NULL)
```

## Arguments

- mat:

  A numeric matrix of 0s and 1s representing species presence-absence.

- n_swaps:

  Number of attempted swaps to perform. Defaults to
  `5 * nrow(mat) * ncol(mat)`.

## Value

A matrix of the same dimensions as `mat`, with swapped entries.

## Details

This function is commonly used in null model analyses for community
ecology. It randomly selects two rows and two columns, checks for a
checkerboard pattern, and swaps the occurrences to preserve species
richness and cell totals. Only 2x2 submatrices matching the checkerboard
patterns `[1, 0; 0, 1]` or `[0, 1; 1, 0]` are swapped.

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com>
