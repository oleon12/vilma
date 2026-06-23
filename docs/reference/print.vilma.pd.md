# Print Summary of a VILMA Phylogenetic Diversity Object

Provides a concise textual summary of a `vilma.pd` object, including:
the number of species, number of spatial cells, total records, species
abundance statistics, and summary statistics of phylogenetic diversity
(PD) metrics.

## Usage

``` r
# S3 method for class 'vilma.pd'
print(x, ...)
```

## Arguments

- x:

  An object of class `vilma.pd`.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the input object `x`.

## Details

This function prints a quick overview of the phylogenetic diversity
stored in a VILMA object. It summarizes the distribution data and PD
values per cell for quick inspection.

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com>
