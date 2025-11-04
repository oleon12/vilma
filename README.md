<!-- badges: start -->
[![R-CMD-check](https://github.com/oleon12/vilma/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/oleon12/vilma/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![GitHub release](https://img.shields.io/github/v/release/oleon12/vilma?display_name=tag)](https://github.com/oleon12/vilma/releases)
[![GitHub stars](https://img.shields.io/github/stars/oleon12/vilma?style=social)](https://github.com/oleon12/vilma)
<!-- badges: end -->



# vilma (v1.0.0) <img src="inst/app/www/Logo.png" align="right" width="120"/>

**Vilma** is an R package to quantify and visualize **spatial phylogenetic diversity** across geographic landscapes.  

It provides a complete workflow to:

- Convert species occurrences into spatial grids
- Calculate **alpha diversity** indices (Faith’s PD, MPD, MNTD, PE, NRI, NTI, RaoQ)
- Compute **beta diversity** metrics (PhyloBeta, PhyloSor, UniFrac, Rao β, βMPD, βMNTD)
- Run **null models** with multiple randomization schemes
- Export results and visualize patterns interactively through a built-in **Shiny app**

Vilma is designed for biogeography, ecology, and conservation applications, with special emphasis on **integrative phylogenomics + spatial analysis**.

---

###### <i>Created in loving memory of <b>Vilma Alvarado</b> whose kindness and love continue to inspire this work. </i>

---
<br>

## Installation 

**Development version (recommended)**
Requires `devtools`:

```r
install.packages("devtools")
devtools::install_github("oleon12/vilma", build_vignettes = FALSE)
```
<br>

## Shiny App

Vilma can be used through a <b>Shiny interface</b>. This interface was built for users with little experience in R or programming in general. However, the user is required to understand every parameter and index in order to perform an appropriate analysis for their research question. Therefore, we always advise users to keep in mind the principle of <b>GIGO (Garbage In, Garbage Out)</b>. Running the app is very easy: simply open an <i>R</i> or <i>RStudio</i> session and run the following lines of code.

```r
library(vilma)

# Launch the analysis app
run.vilma.app()
```
<br>

Vilma is intended to be a pipeline where the user starts by creating a raster file and then uses a phylogeny to calculate &alpha;-diversity, &alpha;-null models, and &beta;-diversity. However, the user can skip some steps and run only the desired analyses. For the Shiny app, please note that the phylogeny (of class <i>phylo</i>) must be uploaded on the &alpha;-diversity page.

<br>

## Interactive view

Although the Shiny app is the most remarkable feature for interactivity, Vilma also provides interactive map visualizations. For every analysis (raster,  &alpha;-diversity, &alpha;-null models, and &beta;), results can be viewed using the <i>view.vilma()</i> function. This function will display a Leaflet map with the raster layers. You can toggle between layers using the Layer Control button <img src="inst/app/www/Layer_Icon.png" width="15"/>. Furthermore, if you move the mouse over any pixel, a pop-up window will appear above the Layer Control button showing the value for that specific pixel.

<br>

## Rasters

The idea of this step is to reduce the time and stress associated with creating and manipulating raster files for spatial analysis. Many other packages aimed at calculating Phylogenetic Diversity indices use community tables. This approach requires several steps to convert rasters to tables and then back to rasters again, which can be an obstacle for many users. In contrast, the <i>points_to_rasters</i> function requires only a single data frame (or table) with three columns: Species, Longitude, and Latitude (decimal degrees). This data is easy to gather, especially from sources like GBIF. Users can also set up the pixel resolution and CRS. The function return a <i>vilma.dist</i> object
<br>

<div align="center">

|Species|Longitude|Latitude|
|:---:|:---:|:---:|
|<i>Artibeus lituratus</i>|-73.45|6.35|
|<i>Artibeus lituratus</i>|-73.55|5.03|
|<i>Artibeus jamaicensis</i>|-70.05|2.05|
|<i>Artibeus obscurus</i>|-71.25|5.93|

<sup>1</sup><i>The column names can be different; the function handles them automatically.</i><br>
<sup>2</sup><i>The columns must follow this exact order: Species, Longitude, and Latitude.</i>

</div>

<br>

```r
points_to_raster(points, crs = 4326, ext = NULL, res = 1, doRast = TRUE, symmetrical = FALSE) 

# Example

dist_ex <- example_dist()
raster_out <- points_to_raster(points = dist_ex, res = 5)
print(raster_out)
plot(raster_out)
view.vilma(raster_out)
```
<br>

## &alpha;-Diversity

<br>

This is the most basic analysis. Vilma provides seven different indices based on different approaches (e.g., minimum spanning tree, distance-based). The goal is to increase the number of available indices over time. Each &alpha;-diversity function requires, at a minimum, a vilma.dist object (created with points_to_raster) and a rooted tree with branch lengths. These functions work with both ultrametric and non-ultrametric trees. The result of each function is a vilma.pd object containing: a table exhibiting the Species Richness, Abundance, and PD values for each cell; a set of rasters (the number of which may vary between indices); and the original distribution matrix.

<br>
```r

# Faith PD
faith.pd(tree, dist, method = c("root", "node", "exclude"))

# MNTD
mntd.calc(tree, dist, method = c("root","node","exclude"), abundance = FALSE)

# MPD
mpd.calc(tree, dist, method = c("root","node","exclude"), abundance = FALSE)

# Phylogenetic Endemisms
pe.calc(tree, dist, RPE = c(TRUE,FALSE), faith.method = c("node","root","exclude"))

# Rao's Q (PD)
rao.calc(tree, dist, abundance = FALSE, scale01 = TRUE)

# NRI
nri.calc <- function(tree, dist, mpd.method = c("root","node","exclude"), abundance = FALSE, iterations = 999, sampling = c("taxa.label","range","neigbor","regional"), n.directions = c("rook","bishop","queen"), regional.weight = c("uniform","frequency","range"))

#NTI
nti.calc <- function(tree, dist, mntd.method = c("root","node","exclude"), abundance = FALSE, iterations = 999, sampling = c("taxa.label","range","neigbor","regional"), n.directions = c("rook","bishop","queen"), regional.weight = c("uniform","frequency","range"))

dist_ex <- example_dist()
raster_out <- points_to_raster(points = dist_ex, res = 5)

tree_ex <- example_tree()

mpd_out <- mpd(tree = tree_ex, dist = raster_out, method = "root", abundance = FALSE)
print(mpd_out)
plot(mpd_out)
view.vilma(mpd_out)
```
<br>

## Null models (&alpha;-Diversity)

<br>

This is one of the main features of Vilma: the implementation of null models for &alpha;-diversity indices. With the exception of NRI and NTI, all indices have their own null models. The null models offer two different approaches: <b>global</b> and <b>cell</b>. When "global" is selected, the function calculates the SES and p-value for the entire area, while "cell" calculates the <i>SES</i> and <i>p-value</i> for each individual cell. Only the "cell" option returns a raster.

The null models feature four different sampling options: taxa.label, range, neighbor, and regional. For neighbor, the algorithm re-samples the points using three methods: "rook", "bishop", and "queen". For regional, the species pool re-sampling is performed with different weights: "uniform", "frequency", and "range".

Null models require a vilma.pd object, a phylogenetic tree, and a vilma.dist object. Each function returns a vilma.null object containing the associated statistics (SES and p-value) in tables, a raster (if the "cell" option was selected), and the original distribution matrix.

<br>
```r
faith.pd.null(pd, tree, dist, iterations = 999, method = c("global","cell"), sampling = c("taxa.label","range","neighbor","regional"), n.directions = c("rook","bishop","queen"),regional.weight = c("uniform","frequency","range"))

dist_ex <- example_dist()
raster_out <- points_to_raster(points = dist_ex, res = 5)

tree_ex <- example_tree()
faith_out <- faith.pd(tree = tree_ex, dist = raster_out, method = "root")

faith_null <- faith.pd.null(pd = faith_out, tree = tree_ex, dist = raster_out, method = "global", sampling = "taxa.label")
print(faith_null)
plot(faith_null)
view.vilma(faith_null)

faith_null <- faith.pd.null(pd = faith_out, tree = tree_ex, dist = raster_out, method = "cell", sampling = "taxa.label")
print(faith_null)
plot(faith_null)
view.vilma(faith_null)

```
<br>
