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

## Installation </br>

### **Development version (recommended)**
Requires `devtools`:

```r
install.packages("devtools")
devtools::install_github("oleon12/vilma", build_vignettes = FALSE)
```

## Shiny App

Vilma can be used through a <b>Shiny interface</b>. This interface was built for users with little experience in R or programming in general. However, the user is required to understand every parameter and index in order to perform an appropriate analysis for their research question. Therefore, we always advise users to keep in mind the principle of <b>GIGO (Garbage In, Garbage Out)</b>. Running the app is very easy: simply open an <i>R</i> or <i>RStudio</i> session and run the following lines of code.

```r
library(vilma)

# Launch the analysis app
run.vilma.app()
```

