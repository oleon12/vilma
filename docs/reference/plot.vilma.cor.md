# Plot VILMA Correlation Results

Provides a quick visualization of a `vilma.cor` object by plotting the
relationship between species richness and phylogenetic diversity (PD)
values. The method automatically selects the best-supported model based
on the minimum AIC value stored in `x$AIC`. If ggplot2 is available, a
scatterplot with a fitted trend line is produced; otherwise, the
function falls back to a base R plot using the selected fitted model
object.

## Usage

``` r
# S3 method for class 'vilma.cor'
plot(x, ...)
```

## Arguments

- x:

  An object of class `vilma.cor`.

- ...:

  Additional arguments (currently not used).

## Value

A plot is produced. Invisibly returns `x`.

## Details

Model selection is performed by choosing the model name corresponding to
the smallest value in `x$AIC`. When ggplot2 is installed, the function
plots `x$Data` using `sr.vals` on the x-axis and `pd.vals` on the
y-axis, and overlays a fitted smooth with
`geom_smooth(method = <best model>)`. If ggplot2 is not installed, the
function plots the best model stored in `x$Models[[method]]` using base
R.

## See also

[`vilma.cor`](https://oleon12.github.io/vilma/reference/vilma.cor.md)

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com/>
