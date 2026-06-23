# Phylogenetic Endemism (PE)

Compute Phylogenetic Endemism (PE) for each grid cell from a rooted
phylogeny and species distributions. Optionally returns Relative PE
(RPE) as PE / Faith's PD per cell.

## Usage

``` r
pe.calc(
  tree,
  dist,
  RPE = c(TRUE, FALSE),
  faith.method = c("node", "root", "exclude")
)
```

## Arguments

- tree:

  A rooted phylogenetic tree of class `phylo` with branch lengths.

- dist:

  An object of class `vilma.dist` (see
  [`points_to_raster`](https://oleon12.github.io/vilma/reference/points_to_raster.md)).

- RPE:

  Logical; if `TRUE` also compute Relative PE (RPE = PE / PD). Default
  `TRUE`.

- faith.method:

  Character string specifying how to handle single-species cells.
  Options are:

  - `"root"` – For cells with single species, PD is calculated
    considering the root path to the species.

  - `"node"` – For cells with single species, PD considers the closest
    ancestral node length.

  - `"exclude"` – Single-species cells are excluded from the calculation
    (default).

## Value

A `vilma.pd` object with:

- `distribution` – original species distribution (input).

- `grid` – reference raster grid (from `dist`).

- `pd.table` – a data frame with:

  - `Cell`:

    Cell identifier.

  - `SR`:

    Species richness per cell.

  - `PE`:

    Phylogenetic Endemism value per cell.

  - `RPE`:

    Relative PE (if `RPE = TRUE`).

- `rasters` – list with abundance, richness, and PE rasters.

- `index` – the string `"pe.calc"`.

## Details

PE measures how much evolutionary history in a cell is geographically
restricted, summing branch lengths weighted by (i) the fraction of
descendants present in the cell and (ii) the inverse of the number of
cells those descendants occupy.

## References

Rosauer D, Laffan SW, Crisp MD, Donnellan SC, Cook LG (2009)
Phylogenetic endemism: a new approach for identifying geographical
concentrations of evolutionary history. *Molecular Ecology*
18(19):4061–4072.
[doi:10.1111/j.1365-294X.2009.04311.x](https://doi.org/10.1111/j.1365-294X.2009.04311.x)

Mishler BD, Knerr N, González-Orozco CE, Thornhill AH, Laffan SW, Miller
JT (2014) Phylogenetic measures of biodiversity and neo- and
paleo-endemism in Australian *Acacia*. *Nature Communications* 5:4473.
[doi:10.1038/ncomms5473](https://doi.org/10.1038/ncomms5473)

Faith DP (1992) Conservation evaluation and phylogenetic diversity.
*Biological Conservation* 61(1):1–10.
[doi:10.1016/0006-3207(92)91201-3](https://doi.org/10.1016/0006-3207%2892%2991201-3)

## See also

[`points_to_raster`](https://oleon12.github.io/vilma/reference/points_to_raster.md)

## Examples

``` r
if (FALSE) { # \dontrun{
tree <- example_tree()
dist <- example_dist()
dist <- points_to_raster(dist)
pe <- pe.calc(tree, dist, RPE = TRUE)
print(pe)
plot(pe)
view.vilma(pe)
} # }
```
