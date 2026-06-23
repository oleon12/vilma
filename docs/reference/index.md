# Package index

## Package overview

General package documentation.

- [`vilma`](https://oleon12.github.io/vilma/reference/vilma-package.md)
  [`vilma-package`](https://oleon12.github.io/vilma/reference/vilma-package.md)
  : vilma: Spatial Phylogenetic Diversity

## Example datasets

Example occurrence, phylogenetic, raster, and polygon datasets included
with vilma.

- [`Artibeus_occ`](https://oleon12.github.io/vilma/reference/Artibeus_occ.md)
  :

  Occurrence records for *Artibeus* (example dataset)

- [`artibeus_tree`](https://oleon12.github.io/vilma/reference/artibeus_tree.md)
  :

  Phylogeny for *Artibeus* species

- [`car_tree`](https://oleon12.github.io/vilma/reference/car_tree.md) :
  Example carnivore tree

- [`data_sa`](https://oleon12.github.io/vilma/reference/data_sa.md) :
  Carnivores occurrences from South Africa

- [`example_dist()`](https://oleon12.github.io/vilma/reference/example_dist.md)
  : Example points for demonstration

- [`example_tree()`](https://oleon12.github.io/vilma/reference/example_tree.md)
  : Example phylogenetic tree for demonstration

- [`sa_pol`](https://oleon12.github.io/vilma/reference/sa_pol.md) :
  South Africa polygon

## Create and convert spatial inputs

Functions for converting occurrence points, rasters, and polygons into
spatial distribution objects.

- [`points_to_raster()`](https://oleon12.github.io/vilma/reference/points_to_raster.md)
  : Convert Species Occurrence Points to Raster Distribution
- [`rast2occ()`](https://oleon12.github.io/vilma/reference/rast2occ.md)
  : Convert binary ENM rasters to occurrence points
- [`pol2occ()`](https://oleon12.github.io/vilma/reference/pol2occ.md) :
  Convert polygon ranges to occurrence records and optional raster
  layers
- [`return.vilma.dist()`](https://oleon12.github.io/vilma/reference/return.vilma.dist.md)
  : Convert a VILMA Distribution Matrix to a Data Frame
- [`return.vilma.dist2()`](https://oleon12.github.io/vilma/reference/return.vilma.dist2.md)
  : Convert a VILMA Distribution Matrix to a Data Frame

## Alpha phylogenetic diversity

Functions for calculating alpha phylogenetic diversity indices.

- [`faith.pd()`](https://oleon12.github.io/vilma/reference/faith.pd.md)
  : Faith's Phylogenetic Diversity (PD)
- [`mpd.calc()`](https://oleon12.github.io/vilma/reference/mpd.calc.md)
  : Mean Pairwise Distance (MPD) Phylogenetic Diversity Index
- [`mntd.calc()`](https://oleon12.github.io/vilma/reference/mntd.calc.md)
  : Mean Nearest Taxon Distance (MNTD) Phylogenetic Diversity Index
- [`pe.calc()`](https://oleon12.github.io/vilma/reference/pe.calc.md) :
  Phylogenetic Endemism (PE)
- [`rao.calc()`](https://oleon12.github.io/vilma/reference/rao.calc.md)
  : Rao's Q (within-community) phylogenetic diversity
- [`nri.calc()`](https://oleon12.github.io/vilma/reference/nri.calc.md)
  : Net Relatedness Index (NRI)
- [`nti.calc()`](https://oleon12.github.io/vilma/reference/nti.calc.md)
  : Nearest Taxon Index (NTI)
- [`pd.taxon()`](https://oleon12.github.io/vilma/reference/pd.taxon.md)
  : Calculate Faith's PD for Individual Taxa

## Spatial null models

Functions for randomizing communities, calculating null expectations,
and evaluating standardized phylogenetic diversity patterns.

- [`faith.pd.null()`](https://oleon12.github.io/vilma/reference/faith.pd.null.md)
  : Null model for Faith's Phylogenetic Diversity (PD)
- [`mpd.calc.null()`](https://oleon12.github.io/vilma/reference/mpd.calc.null.md)
  : Null model for Mean Pairwise Distance (MPD) phylogenetic diversity
- [`mntd.calc.null()`](https://oleon12.github.io/vilma/reference/mntd.calc.null.md)
  : Null model for Mean Nearest Taxon Distance (MNTD) Phylogenetic
  Diversity
- [`pe.calc.null()`](https://oleon12.github.io/vilma/reference/pe.calc.null.md)
  : Null model for Faith's Phylogenetic Diversity (PD)
- [`rao.calc.null()`](https://oleon12.github.io/vilma/reference/rao.calc.null.md)
  : Rao's Q (α) null models for vilma objects
- [`swap.null()`](https://oleon12.github.io/vilma/reference/swap.null.md)
  : Randomly Swap Cells in a Presence-Absence Matrix

## Beta phylogenetic diversity

Functions for calculating phylogenetic dissimilarity, turnover,
nestedness, UniFrac, and beta diversity indices.

- [`beta.mpd()`](https://oleon12.github.io/vilma/reference/beta.mpd.md)
  : Between-community Mean Pairwise Distance (betaMPD) between
  communities (cells)
- [`beta.mntd()`](https://oleon12.github.io/vilma/reference/beta.mntd.md)
  : Between-community Mean Nearest Taxon Distance (betaMNTD) between
  communities (cells)
- [`phylo.beta()`](https://oleon12.github.io/vilma/reference/phylo.beta.md)
  : Phylogenetic Beta Diversity Partitioning (β_sor, β_sim, β_sne)
- [`phylosor.calc()`](https://oleon12.github.io/vilma/reference/phylosor.calc.md)
  : Phylogenetic Sørensen Similarity (PhyloSor) between communities
  (cells)
- [`rao.beta()`](https://oleon12.github.io/vilma/reference/rao.beta.md)
  : Rao Beta (phylogenetic) dissimilarity index
- [`unifrac.calc()`](https://oleon12.github.io/vilma/reference/unifrac.calc.md)
  : UniFrac Phylogenetic Dissimilarity between communities (cells)

## Spatial regionalization

Functions for clustering and regionalizing phylogenetic diversity
patterns.

- [`vilma.regionalize()`](https://oleon12.github.io/vilma/reference/vilma.regionalize.md)
  : Regionalize cells into phylogenetic regions
- [`find_kdiss()`](https://oleon12.github.io/vilma/reference/find_kdiss.md)
  : Select optimal K using explained dissimilarity
- [`find_knnK()`](https://oleon12.github.io/vilma/reference/find_knnK.md)
  : Select optimal k for kNN network regionalization
- [`find_kpam()`](https://oleon12.github.io/vilma/reference/find_kpam.md)
  : Select optimal k for PAM clustering using silhouette width
- [`find_ksil()`](https://oleon12.github.io/vilma/reference/find_ksil.md)
  : Select optimal K using silhouette width

## Environmental niche modeling and spatial correlation

Additional tools for environmental niche modeling and correlation
analyses.

- [`vilma.ENM()`](https://oleon12.github.io/vilma/reference/vilma.ENM.md)
  : Fit Vilma ensemble ENMs for multiple species with optional M-area
  definition

- [`fast.ensemble.enm()`](https://oleon12.github.io/vilma/reference/fast.ensemble.enm.md)
  : Fast ensemble ecological niche model (ENM) with cross-validated
  weighting

- [`vilma.cor()`](https://oleon12.github.io/vilma/reference/vilma.cor.md)
  :

  Correlate phylogenetic diversity with species richness in a `vilma.pd`
  object

## Visualization

Plotting and interactive mapping methods for vilma objects.

- [`plot(`*`<vilma.beta>`*`)`](https://oleon12.github.io/vilma/reference/plot.vilma.beta.md)
  :

  Plot method for `vilma.beta` objects

- [`plot(`*`<vilma.cor>`*`)`](https://oleon12.github.io/vilma/reference/plot.vilma.cor.md)
  : Plot VILMA Correlation Results

- [`plot(`*`<vilma.dist>`*`)`](https://oleon12.github.io/vilma/reference/plot.vilma.dist.md)
  : Plot VILMA Distribution Rasters

- [`plot(`*`<vilma.null>`*`)`](https://oleon12.github.io/vilma/reference/plot.vilma.null.md)
  : Plot Results of a VILMA Null Model

- [`plot(`*`<vilma.pd>`*`)`](https://oleon12.github.io/vilma/reference/plot.vilma.pd.md)
  : Plot Rasters from a VILMA PD Object

- [`plot(`*`<vilma.region>`*`)`](https://oleon12.github.io/vilma/reference/plot.vilma.region.md)
  : Plot VILMA Regionalization Rasters

- [`view.vilma()`](https://oleon12.github.io/vilma/reference/view.vilma.md)
  :

  Unified interactive viewer for *Vilma* spatial phylogenetic outputs

- [`view.vilma.beta()`](https://oleon12.github.io/vilma/reference/view.vilma.beta.md)
  : Interactive viewer for vilma β-diversity rasters

- [`view.vilma.dist()`](https://oleon12.github.io/vilma/reference/view.vilma.dist.md)
  : Interactive viewer for vilma distributions (richness & abundance)

- [`view.vilma.null()`](https://oleon12.github.io/vilma/reference/view.vilma.null.md)
  : Interactive viewer for vilma null-model results (SES raster)

- [`view.vilma.pd()`](https://oleon12.github.io/vilma/reference/view.vilma.pd.md)
  :

  View `vilma.pd` Rasters in an Interactive Leaflet Map

- [`view.vilma.region()`](https://oleon12.github.io/vilma/reference/view.vilma.region.md)
  :

  View `vilma.region` Rasters in an Interactive Leaflet Map

## Printing methods

Console print methods for vilma objects.

- [`print(`*`<vilma.beta>`*`)`](https://oleon12.github.io/vilma/reference/print.vilma.beta.md)
  :

  Print summary for a `vilma.beta` object

- [`print(`*`<vilma.dist>`*`)`](https://oleon12.github.io/vilma/reference/print.vilma.dist.md)
  : Print Summary of a VILMA Distribution Object

- [`print(`*`<vilma.null>`*`)`](https://oleon12.github.io/vilma/reference/print.vilma.null.md)
  : Print Summary of a VILMA Null Model

- [`print(`*`<vilma.pd>`*`)`](https://oleon12.github.io/vilma/reference/print.vilma.pd.md)
  : Print Summary of a VILMA Phylogenetic Diversity Object

## Export results

Functions for writing vilma outputs to disk.

- [`write.vilma.beta()`](https://oleon12.github.io/vilma/reference/write.vilma.beta.md)
  :

  Write a `vilma.beta` Object to Disk

- [`write.vilma.dist()`](https://oleon12.github.io/vilma/reference/write.vilma.dist.md)
  :

  Write a `vilma.dist` Object to Disk

- [`write.vilma.null()`](https://oleon12.github.io/vilma/reference/write.vilma.null.md)
  :

  Write a `vilma.null` Object to Disk

- [`write.vilma.pd()`](https://oleon12.github.io/vilma/reference/write.vilma.pd.md)
  :

  Write a `vilma.pd` Object to Disk

- [`write.vilma.region()`](https://oleon12.github.io/vilma/reference/write.vilma.region.md)
  :

  Write a `vilma.region` Object to Disk

## Shiny application

Launch the bundled graphical user interface.

- [`run.vilma.app()`](https://oleon12.github.io/vilma/reference/run.vilma.app.md)
  : Run the Vilma Shiny App

## Internal or developer utilities

Helper functions mainly used internally by vilma.

- [`descendant.matrix()`](https://oleon12.github.io/vilma/reference/descendant.matrix.md)
  : Build Descendant Matrix (edges × tips) for a Phylogenetic Tree
- [`find.root.node()`](https://oleon12.github.io/vilma/reference/find.root.node.md)
  : Find the Root Node of a Rooted Phylogeny
- [`zzz_imports`](https://oleon12.github.io/vilma/reference/zzz_imports.md)
  : Internal: central imports for vilma
