# Calculate Faith's PD for Individual Taxa

Computes the phylogenetic diversity (PD) associated with one or more
taxa in a phylogenetic tree. Two methods are available:

- **root** – PD from each taxon to the root of the tree.

- **node** – PD as the branch length leading to each taxon.

## Usage

``` r
pd.taxon(tree = NULL, sp = NULL, method = c("root", "node"))
```

## Arguments

- tree:

  A rooted phylogenetic tree of class `phylo`.

- sp:

  A character vector of tip labels to calculate PD for. Defaults to all
  tips in the tree.

- method:

  Character, either `"root"` or `"node"` indicating the PD calculation
  method. Defaults to `"root"`.

## Value

A named numeric vector of PD values for the specified taxa.

## Details

This function calculates Faith's PD for individual taxa. The `"root"`
method sums branch lengths from the taxon to the root of the tree, while
the `"node"` method returns only the branch leading directly to the
taxon.

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com>
