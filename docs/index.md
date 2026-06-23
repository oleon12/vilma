# vilma: an R package for spatial phylogenetic diversity analyses with a Shiny interface. ![vilma package logo](reference/figures/Logo_0.png)

------------------------------------------------------------------------

## About vilma (1.0.0)

`vilma` is an R-based platform designed to quantify and visualize
**spatial phylogenetic diversity** across geographic landscapes. It
integrates multiple analytical frameworks to examine how evolutionary
history and ecological processes shape biodiversity through space and
time.

------------------------------------------------------------------------

### Vilma provides a complete and reproducible workflow to:

- Compute **alpha diversity** indices such as Faith’s PD, MPD, MNTD, PE,
  Rao’s Q, NRI, and NTI.

- Estimate **beta diversity** metrics, including PhyloSor, UniFrac,
  Phylobeta weighted and unweighted, betaMPD, and betaMNTD.

- Implement and test **null models** at both regional and local
  cell-based levels to assess whether observed spatial patterns deviate
  from random expectations.

- Visualize and export spatially explicit outputs as tables and raster
  maps for further ecological and evolutionary analyses.

------------------------------------------------------------------------

Developed as part of an ongoing research project at the **Department of
Mammalogy, American Museum of Natural History**, and the **Department of
Earth and Environmental Sciences, Rutgers University–Newark**, `vilma`
provides a flexible framework to investigate phylogenetic structure,
endemism, and connectivity in ecological and evolutionary studies.

------------------------------------------------------------------------

  

![vilma](reference/figures/vilma_pixel.png)

*Created in loving memory of **Vilma Alvarado**, whose kindness and love
continue to inspire this work.*

  

------------------------------------------------------------------------

### Installation

**Development version (recommended)** Requires `devtools`:

``` r
install.packages("devtools")
devtools::install_github("oleon12/vilma", build_vignettes = FALSE)
```

  

### Shiny App

Vilma can be used through a **Shiny interface**. This interface was
built for users with little experience in R or programming in general.
However, the user must understand all parameters and indices to perform
an appropriate analysis for their research question. Therefore, we
always advise users to keep in mind the principle of **GIGO (Garbage In,
Garbage Out)**. Running the app is very easy: simply open an *R* or
*RStudio* session and run the following lines of code.

``` r
library(vilma)

# Launch the analysis app
run.vilma.app()
```

  

Vilma is intended to be a pipeline where the user starts by creating a
raster file and then uses a phylogeny to calculate α-diversity, α-null
models, and β-diversity. However, the user can skip some steps and run
only the desired analyses. For the Shiny app, please note that the
phylogeny (of class *phylo*) must be uploaded on the α-diversity page.

  

### Interactive view

Although the Shiny app is the most remarkable feature for interactivity,
Vilma also provides interactive map visualizations. For every analysis
(raster, α-diversity, α-null models, and β), results can be viewed using
the *view.vilma()* function. This function displays a Leaflet map with
raster layers. You can toggle between layers using the Layer Control
button ![Layer control icon](reference/figures/Layer_Icon.png).
Furthermore, if you move the mouse over any pixel, a pop-up window will
appear above the Layer Control button showing the value for that
specific pixel.

  
