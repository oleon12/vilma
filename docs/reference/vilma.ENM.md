# Fit Vilma ensemble ENMs for multiple species with optional M-area definition

Fits ensemble ecological niche models (ENMs) for each species in a
multi-species occurrence table using
[`fast.ensemble.enm`](https://oleon12.github.io/vilma/reference/fast.ensemble.enm.md).
For each species, the function can define an accessible area (M; sensu
BAM) from the full study area, from all occurrences pooled, or
separately per taxon using a minimum convex polygon (MCP) or buffer
workflow. Optionally, collinear predictors within M are removed using a
VIF-based stepwise procedure. Species with fewer than 3 unique
occurrences are automatically skipped. Finally, an estimated species
richness raster is computed by summing the binary predictions across all
successfully modeled species.

## Usage

``` r
vilma.ENM(
  occ,
  envs,
  Marea = c("byTaxon", "byArea", "user"),
  background_points = 10000,
  cv_folds = 5,
  n_cores = 1,
  test.VIF = TRUE,
  buffer.dist = 10000,
  maxent_classes = "default",
  maxent_regmult = 1,
  rf_ntree = 500,
  rf_mtry = NULL,
  rf_nodesize = 1,
  glmnet_alpha = 1,
  glmnet_lambda = "lambda.1se",
  weight_method = "auc",
  thresh_stat = c("no_omission", "spec_sens", "equal_sens_spec", "prevalence", "kappa",
    "sensitivity")
)
```

## Arguments

- occ:

  A data frame of occurrences with exactly three columns. The function
  expects either `Sp`, `Lon`, `Lat` or `Sp`, `Longitude`, `Latitude`.
  Coordinates must be in decimal degrees. Columns are internally
  standardized to `Sp`, `Lon`, `Lat`. Duplicate coordinates are handled
  within
  [`fast.ensemble.enm`](https://oleon12.github.io/vilma/reference/fast.ensemble.enm.md).

- envs:

  Environmental predictors. Can be a raster `RasterStack`,
  `RasterBrick`, or a terra `SpatRaster`. If a `SpatRaster` is provided,
  it is coerced to a `RasterStack` via
  [`raster::stack()`](https://rdrr.io/pkg/raster/man/stack.html) for
  compatibility with downstream functions.

- Marea:

  Character. Default = `"byTaxon"`. Defines how the accessible area (M)
  is constructed:

  `"byTaxon"`

  :   Build M separately for each species. For species with `>5`
      occurrences, uses an MCP buffered by `buffer.dist` meters. For
      species with `≤5` occurrences, uses a circular buffer around all
      points (also of radius `buffer.dist`) to avoid unstable MCPs.

  `"byArea"`

  :   Build a single M from all occurrences pooled and apply it to all
      species. Currently implemented as a placeholder; requires
      additional development.

  `"user"`

  :   Use `envs` as provided (no M-area cropping).

- background_points:

  Integer. Default = `10000`. Number of background points passed to
  [`fast.ensemble.enm`](https://oleon12.github.io/vilma/reference/fast.ensemble.enm.md).
  Automatically capped at 30\\ of valid raster cells.

- cv_folds:

  Integer. Default = `5`. Number of cross-validation folds passed to
  [`fast.ensemble.enm`](https://oleon12.github.io/vilma/reference/fast.ensemble.enm.md).
  Automatically reduced for small sample sizes.

- n_cores:

  Integer. Default = `1`. Number of CPU cores for parallel
  cross-validation inside
  [`fast.ensemble.enm`](https://oleon12.github.io/vilma/reference/fast.ensemble.enm.md).

- test.VIF:

  Logical. Default = `TRUE`. If `TRUE`, remove collinear predictors
  inside each M-area. For species with `>5` occurrences, applies
  [`vifstep`](https://rdrr.io/pkg/usdm/man/vif.html) to environmental
  values extracted at occurrence locations, then uses
  [`exclude`](https://rdrr.io/pkg/usdm/man/exclude.html) to remove
  correlated variables. For species with `≤5` occurrences, VIF is
  skipped due to insufficient sample size.

- buffer.dist:

  Numeric. Default = `10000`. Buffer distance in meters for MCP
  expansion or point buffers when constructing M areas. Used in
  `.envs_m_area` (for MCPs) and `.envs_m_buffer` (for point buffers).

- maxent_classes:

  Character. Default = `"default"`. Maxent feature classes passed to
  [`fast.ensemble.enm`](https://oleon12.github.io/vilma/reference/fast.ensemble.enm.md).
  Can be `"default"`, `"l"`, `"lq"`, `"lh"`, `"lqhp"`, etc.

- maxent_regmult:

  Numeric. Default = `1`. Regularization multiplier for Maxent passed to
  [`fast.ensemble.enm`](https://oleon12.github.io/vilma/reference/fast.ensemble.enm.md).

- rf_ntree:

  Integer. Default = `500`. Number of trees in Random Forest passed to
  [`fast.ensemble.enm`](https://oleon12.github.io/vilma/reference/fast.ensemble.enm.md).

- rf_mtry:

  Integer or `NULL`. Default = `NULL`. Number of candidate variables at
  each split. If `NULL`, determined internally by
  [`fast.ensemble.enm`](https://oleon12.github.io/vilma/reference/fast.ensemble.enm.md)
  as `floor(sqrt(p))`.

- rf_nodesize:

  Integer. Default = `1`. Minimum node size in Random Forest passed to
  [`fast.ensemble.enm`](https://oleon12.github.io/vilma/reference/fast.ensemble.enm.md).

- glmnet_alpha:

  Numeric. Default = `1`. Elastic net mixing parameter (`1`=LASSO,
  `0`=ridge) passed to
  [`fast.ensemble.enm`](https://oleon12.github.io/vilma/reference/fast.ensemble.enm.md).

- glmnet_lambda:

  Character. Default = `"lambda.1se"`. Lambda choice for predictions
  from
  [`cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html)
  passed to
  [`fast.ensemble.enm`](https://oleon12.github.io/vilma/reference/fast.ensemble.enm.md).
  For species with `<10` unique presences, a fixed lambda of 0.01 is
  used without cross-validation.

- weight_method:

  Character. Default = `"auc"`. Weighting method used inside
  [`fast.ensemble.enm`](https://oleon12.github.io/vilma/reference/fast.ensemble.enm.md)
  to combine models. Options: `"auc"`, `"tss"`, or any other value
  (falls back to equal weights).

- thresh_stat:

  Character. Threshold statistic passed to
  [`fast.ensemble.enm`](https://oleon12.github.io/vilma/reference/fast.ensemble.enm.md)
  and used to select an optimal threshold for binarizing the ensemble
  prediction. One of `"no_omission"`, `"spec_sens"`,
  `"equal_sens_spec"`, `"prevalence"`, `"kappa"`, or `"sensitivity"`.

## Value

An object of class `"vilma.ENM"` containing a named list with four
elements:

- Prediction:

  A list of continuous suitability rasters (one per successfully modeled
  species) from the ensemble prediction. Names correspond to species
  names.

- Boolean:

  A list of binary rasters (one per successfully modeled species)
  indicating suitability above the optimal threshold (based on
  `thresh_stat`).

- Est_richness:

  A single raster layer of class `RasterLayer` representing estimated
  species richness, calculated as the sum of all binary predictions
  across successfully modeled species, aligned to a common extent. Areas
  with no predicted suitability for any species are set to NA.

- Skipped:

  A character vector with names of species that were skipped due to
  insufficient occurrences (\<3 unique points).

## Details

For each species, `vilma.ENM`:

1.  Subsets `occ` to the focal species.

2.  Checks if the species has at least 3 unique occurrences. If not,
    skips it with a message and continues to the next species.

3.  Constructs the accessible area (M) according to `Marea`:

    - `"byTaxon"`: Uses `.envs_m_area` (MCP + buffer) for species with
      \>5 occurrences, and `.envs_m_buffer` (circular buffers around
      points) for species with ≤5 occurrences to avoid unstable MCPs.

    - `"byArea"`: Placeholder for pooled M area.

    - `"user"`: Uses full `envs` without cropping.

4.  Optionally removes collinear predictors in M using a VIF-based
    filter (if `test.VIF = TRUE`). For species with \>5 occurrences,
    applies `vifstep` to environmental values at occurrence points; for
    ≤5 occurrences, VIF is skipped due to insufficient sample size.

5.  Fits an ensemble ENM using
    [`fast.ensemble.enm`](https://oleon12.github.io/vilma/reference/fast.ensemble.enm.md)
    with all specified parameters.

6.  Produces two outputs: (i) a continuous suitability raster (ensemble
    prediction), and (ii) a binary prediction raster using the model's
    selected optimal threshold.

After processing all species, the function:

- Aligns all binary rasters to a common extent (the maximum extent
  covered by any successfully modeled species).

- Sums the binary rasters to generate an estimated species richness map.

- Sets pixels with zero summed value to NA (indicating areas not
  predicted suitable for any species).

## Note

- The `"byArea"` option is currently a placeholder and requires
  additional implementation to function properly.

- For species with very few occurrences (\<10 unique presences),
  [`fast.ensemble.enm`](https://oleon12.github.io/vilma/reference/fast.ensemble.enm.md)
  issues warnings and adjusts modeling parameters (e.g., uses GLMNET
  without CV). Results for such species should be interpreted with
  caution.

- VIF-based variable selection requires the usdm package. If not
  available, the function will attempt to install it.

- The estimated richness raster is based solely on the binary
  predictions of successfully modeled species; skipped species are not
  included.

- All binary rasters are extended to a common extent before summation,
  ensuring that the richness estimate covers the maximum geographic
  range of any modeled species.

- The function returns an S3 object of class `"vilma.ENM"` to enable
  custom printing and plotting methods in future package versions.

## See also

[`fast.ensemble.enm`](https://oleon12.github.io/vilma/reference/fast.ensemble.enm.md),
[`maxnet`](https://rdrr.io/pkg/maxnet/man/maxnet.html),
[`randomForest`](https://rdrr.io/pkg/randomForest/man/randomForest.html),
[`cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html),
[`vifstep`](https://rdrr.io/pkg/usdm/man/vif.html),
[`exclude`](https://rdrr.io/pkg/usdm/man/exclude.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with per-species M areas
out <- vilma.ENM(
  occ = occ_data,
  envs = envs,
  Marea = "byTaxon",
  background_points = 5000,
  cv_folds = 5,
  n_cores = 2,
  test.VIF = TRUE,
  buffer.dist = 10000,
  maxent_classes = "lq",
  maxent_regmult = 1,
  rf_ntree = 500,
  glmnet_alpha = 1,
  glmnet_lambda = "lambda.1se",
  weight_method = "auc",
  thresh_stat = "spec_sens"
)

# Check object class and structure
class(out)  # "vilma.ENM"
names(out)

# Check results
names(out$Prediction)       # Successfully modeled species
out$Skipped                 # Species skipped due to insufficient data

# Visualize outputs
plot(out$Prediction[[1]])    # Continuous suitability for first species
plot(out$Boolean[[1]])       # Binary prediction for first species
plot(out$Est_richness)       # Estimated species richness map

# Access specific species results
sp_name <- names(out$Prediction)[1]
sp_continuous <- out$Prediction[[sp_name]]
sp_binary <- out$Boolean[[sp_name]]

# Summary of modeling results
cat("Successfully modeled:", length(out$Prediction), "species\n")
cat("Skipped:", length(out$Skipped), "species\n")
} # }
```
