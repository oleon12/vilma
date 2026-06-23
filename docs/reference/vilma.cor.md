# Correlate phylogenetic diversity with species richness in a `vilma.pd` object

Fits simple statistical models relating phylogenetic diversity (PD) to
species richness (SR) using the per-cell summary table stored in a
`vilma.pd` object. The function fits three model families—linear model
(LM), generalized linear model (GLM), and generalized additive model
(GAM)—and selects the best model within the GLM family and within the
GAM family using AIC. It then compares the best LM, best GLM, and best
GAM using AIC and returns a `vilma.cor` object for downstream
printing/plotting.

## Usage

``` r
vilma.cor(vilma.pd)
```

## Arguments

- vilma.pd:

  A `vilma.pd` object containing `vilma.pd$pd.table` with at least two
  numeric columns: `PD` (phylogenetic diversity) and `SR` (species
  richness).

## Value

An object of class `vilma.cor`, a list with:

- Data:

  A data frame with columns `pd.vals` and `sr.vals` used for fitting.

- Models:

  A list with elements `lm`, `glm`, and `gam` containing the fitted
  model objects.

- AIC:

  Named numeric vector with AIC values for `lm`, `glm`, and `gam`.

## Details

The following models are fit:

- **LM:** `lm(PD ~ SR)`.

- **GLM:** Gamma(log), inverse Gaussian(log), and Gaussian(log); the GLM
  with minimum AIC is retained.

- **GAM:** a linear term `gam(PD ~ SR)` and a smooth term
  `gam(PD ~ s(SR))`; the GAM with minimum AIC is retained.

The final AIC comparison among `lm`, selected `glm`, and selected `gam`
is returned in `$AIC`.

## See also

[`lm`](https://rdrr.io/r/stats/lm.html),
[`glm`](https://rdrr.io/r/stats/glm.html),
[`AIC`](https://rdrr.io/r/stats/AIC.html),
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html),
[`plot.vilma.cor`](https://oleon12.github.io/vilma/reference/plot.vilma.cor.md)

## Examples

``` r
if (FALSE) { # \dontrun{
cor_out <- vilma.cor(vilma_pd)
cor_out$AIC
plot(cor_out)
} # }
```
