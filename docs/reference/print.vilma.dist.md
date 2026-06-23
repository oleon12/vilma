# Print Summary of a VILMA Distribution Object

Provides a concise textual summary of a `vilma.dist` object, including
the number of species, number of cells, total records, and basic species
abundance statistics.

## Usage

``` r
# S3 method for class 'vilma.dist'
print(x, ...)
```

## Arguments

- x:

  An object of class `vilma.dist`.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the input object `x`.

## Details

This function is designed to give a quick overview of the distribution
data stored in a VILMA object. It prints the number of unique taxa,
spatial cells, total records, and summary statistics for species
abundance.

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com>
