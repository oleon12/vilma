# Fast ensemble ecological niche model (ENM) with cross-validated weighting

Fits an ensemble ecological niche model (ENM) using three algorithms:
Maxent (via maxnet), Random Forest (via randomForest), and regularized
logistic regression (via glmnet). Presence points are calibrated against
background (pseudo-absence) points sampled from valid raster cells (or
provided by the user). Models are evaluated via stratified
cross-validation and combined using performance-based weights (AUC or
TSS). The function returns an `"ensemble_enm"` object containing fitted
models, ensemble weights, an internal prediction function, evaluation
results, and calibration data.

## Usage

``` r
fast.ensemble.enm(
  occurrences,
  predictors,
  background_points = 10000,
  cv_folds = 5,
  n_cores = 1,
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

- occurrences:

  A data frame with occurrence coordinates. Required columns: `Lon` and
  `Lat` (decimal degrees). Rows with missing predictor values after
  extraction are removed.

- predictors:

  A raster predictor object (e.g., `RasterStack` or `RasterBrick`)
  containing environmental predictors used for calibration and
  prediction.

- background_points:

  Background specification. Default = `10000`. If numeric, the number of
  random background points sampled with
  [`randomPoints`](https://rdrr.io/pkg/dismo/man/randomPoints.html)
  (internally capped to 30\\ the first predictor layer). If a data
  frame, it is interpreted as background coordinates and used directly.
  If `NULL`, defaults to sampling 10,000 points (also subject to the
  valid-cell cap).

- cv_folds:

  Integer. Default = `5`. Number of cross-validation folds. Folds are
  assigned separately to presences and background (stratified). If there
  are too few presences or background points, `cv_folds` is reduced
  automatically.

- n_cores:

  Integer. Default = `1`. Number of CPU cores to use for parallel
  cross-validation predictions. If `> 1`, folds are evaluated via
  foreach/doParallel. Model fitting on the full dataset remains
  sequential.

- maxent_classes:

  Character. Default = `"default"`. Feature classes passed to
  [`maxnet.formula`](https://rdrr.io/pkg/maxnet/man/maxnet.html) (e.g.,
  `"l"`, `"lq"`, `"lqh"`, `"default"`).

- maxent_regmult:

  Numeric. Default = `1`. Regularization multiplier for maxnet
  (`regmult`).

- rf_ntree:

  Integer. Default = `500`. Number of trees for Random Forest.

- rf_mtry:

  Integer or `NULL`. Default = `NULL`. Number of variables randomly
  sampled as candidates at each split. If `NULL`, set to
  `floor(sqrt(p))` where `p` is the number of predictors.

- rf_nodesize:

  Integer. Default = `1`. Minimum size of terminal nodes for Random
  Forest.

- glmnet_alpha:

  Numeric. Default = `1`. Elastic net mixing parameter (`1` = LASSO, `0`
  = ridge).

- glmnet_lambda:

  Character. Default = `"lambda.1se"`. Which lambda to use for
  prediction from
  [`cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html)
  (commonly `"lambda.min"` or `"lambda.1se"`).

- weight_method:

  Character. Default = `"auc"`. Method to derive ensemble weights from
  cross-validated predictions. One of `"auc"`, `"tss"`, or any other
  value (falls back to equal weights).

- thresh_stat:

  Character. Default = `"no_omission"`. Threshold statistic used by
  [`threshold`](https://rdrr.io/pkg/dismo/man/threshold.html) for the
  reported optimal threshold. One of `"no_omission"`, `"spec_sens"`,
  `"equal_sens_spec"`, `"prevalence"`, `"kappa"`, or `"sensitivity"`.

## Value

An object of class `"ensemble_enm"` (a list) with:

- models:

  List with fitted models: `$maxent`, `$rf`, `$glmnet`.

- ensemble_weights:

  Named numeric vector of length 3 with weights for `"Maxent"`,
  `"RandomForest"`, and `"GLMNET"`.

- ensemble_predict:

  A function that generates weighted ensemble suitability predictions
  for new data (raster or tabular).

- evaluation:

  List with ensemble evaluation statistics: `auc`, `cor`, `tss`,
  `optimal_threshold`, `threshold_stat`, and `confusion_matrix` (using
  the chosen threshold).

- calibration_data:

  List with calibration inputs: `occ_env`, `bg_env`, `occ_coords`,
  `bg_coords`, and `n_unique_presences`.

- parameters:

  List summarizing key parameter settings for Maxent, RF, GLMNET, and
  cross-validation/parallel options.

- summary_matrix:

  A data frame summarizing the modeling run, including number of unique
  presences, number of background points, AUC, TSS, and optimal
  threshold.

## Details

Workflow:

1.  Extract predictor values for occurrence coordinates and remove rows
    with missing predictor data.

2.  Create background points: if numeric, sample from valid raster cells
    using
    [`randomPoints`](https://rdrr.io/pkg/dismo/man/randomPoints.html);
    if a data frame, use provided background coordinates. The number of
    sampled points is capped at 30\\ in the first predictor layer to
    avoid requesting more points than available.

3.  Fit three models on the full dataset (presence/background labels):

    - Maxent via maxnet with
      `maxnet.formula(..., classes = maxent_classes)` and
      `regmult = maxent_regmult`.

    - Random Forest via randomForest with `ntree`, `mtry`, and
      `nodesize`.

    - Regularized logistic regression via glmnet using `cv.glmnet`
      (binomial family) with `alpha = glmnet_alpha`.

4.  Compute out-of-fold predictions for each model using stratified
    folds. Fold-level model fits are repeated within each fold (and can
    be parallelized).

5.  Derive model weights from cross-validated performance using
    `weight_method` (AUC or TSS). Weights are normalized to sum to 1.

6.  Evaluate the ensemble using cross-validated ensemble predictions and
    compute threshold-based diagnostics using `thresh_stat`.

The returned `ensemble_predict` function supports:

- a raster `Raster*` object as `newdata`, returning a single-layer
  suitability raster (with `NA` where predictors are incomplete), or

- a data frame/matrix of predictor values, returning a numeric vector of
  suitability.

## See also

[`randomPoints`](https://rdrr.io/pkg/dismo/man/randomPoints.html),
[`evaluate`](https://rdrr.io/pkg/dismo/man/evaluate.html),
[`threshold`](https://rdrr.io/pkg/dismo/man/threshold.html),
[`maxnet`](https://rdrr.io/pkg/maxnet/man/maxnet.html),
[`maxnet.formula`](https://rdrr.io/pkg/maxnet/man/maxnet.html),
[`randomForest`](https://rdrr.io/pkg/randomForest/man/randomForest.html),
[`cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fast.ensemble.enm(occurrences, predictors,
                         background_points = 1000,
                         cv_folds = 5,
                         n_cores = 2,
                         maxent_classes = "lq",
                         maxent_regmult = 1,
                         rf_ntree = 500,
                         glmnet_alpha = 1,
                         glmnet_lambda = "lambda.1se",
                         weight_method = "auc",
                         thresh_stat = "spec_sens")

fit$evaluation
print(fit$summary_matrix)

# Predict on a RasterStack/RasterBrick:
suit <- fit$ensemble_predict(predictors)
plot(suit)
} # }
```
