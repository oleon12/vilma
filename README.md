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

## 🚀 Installation

### **Development version (recommended)**
Requires `devtools`:

```r
install.packages("devtools")
devtools::install_github("oleon12/vilma", build_vignettes = FALSE)

library(vilma)

# Launch the analysis app
run.vilma.app()

# Load data
dist <- points_to_raster(points, res = 1, crs = 4326)

# Calculate Faith's PD
pd <- faith.pd(dist, tree)

# Null model example
pd_null <- faith.pd.null(dist, tree, method = "taxa_label", iterations = 100)

help(package = "vilma")

