# Plot Results of a VILMA Null Model

Provides a quick visualization of null model results from `vilma.null`
objects. Depending on the null model method, it plots either:

- **Global** – A histogram of randomized mean PD values with the
  observed value indicated.

- **Cell** – A raster map of SES values, with numerical SES displayed
  per cell.

## Usage

``` r
# S3 method for class 'vilma.null'
plot(x, ...)
```

## Arguments

- x:

  An object of class `vilma.null`, typically the output of a null model
  function.

- ...:

  Additional arguments passed to plotting functions.

## Value

Produces a plot (histogram or raster map) depending on the method used
to generate the null model. No values are returned.

## Details

This function is designed for quick inspection of null model outputs,
either summarizing global PD deviations or mapping SES values per cell.
For more advanced plotting and customization, users may directly access
the raster or data stored in the `vilma.null` object.

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com>
