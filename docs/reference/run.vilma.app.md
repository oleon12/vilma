# Run the Vilma Shiny App

Launches the Shiny application bundled in `inst/app`.

## Usage

``` r
run.vilma.app(
  host = "127.0.0.1",
  port = NULL,
  launch.browser = TRUE,
  display.mode = c("auto", "showcase"),
  ...
)
```

## Arguments

- host:

  Host interface. Defaults to `"127.0.0.1"`.

- port:

  TCP port or `NULL` to auto-pick.

- launch.browser:

  Open a browser automatically? Default `TRUE`.

- display.mode:

  Shiny display mode, `"auto"` or `"showcase"`.

- ...:

  Passed to
  [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html) (e.g.,
  `test.mode = TRUE`).

## Value

Invisibly returns the result of
[`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html).
