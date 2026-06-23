# Convert a VILMA Distribution Matrix to a Data Frame

Converts a presence-absence matrix (or similar object) into a
long-format data frame compatible with VILMA functions. Each row
represents a species-cell record.

## Usage

``` r
return.vilma.dist(x, ...)
```

## Arguments

- x:

  A matrix or object representing species distributions
  (presence/absence).

- ...:

  Additional arguments (currently unused).

## Value

A data frame with columns:

- Sp:

  Species name or identifier.

- Lon:

  Longitude (default 0).

- Lat:

  Latitude (default 0).

- Cell:

  Cell identifier corresponding to the row in the input matrix.

## Details

This function takes a matrix where rows are cells and columns are
species, and returns a data frame with columns `Sp`, `Lon`, `Lat`, and
`Cell`. The `Lon` and `Lat` columns are filled with zeros by default, as
coordinates are not required for internal VILMA calculations.

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com>
