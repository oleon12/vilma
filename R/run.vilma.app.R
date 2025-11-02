#' Run the Vilma Shiny App
#' @description
#' Launches the Shiny application bundled in \code{inst/app}.
#' @param host Host interface. Defaults to \code{"127.0.0.1"}.
#' @param port TCP port or \code{NULL} to auto-pick.
#' @param launch.browser Open a browser automatically? Default \code{TRUE}.
#' @param display.mode Shiny display mode, \code{"auto"} or \code{"showcase"}.
#' @param ... Passed to \code{shiny::runApp()} (e.g., \code{test.mode = TRUE}).
#' @return Invisibly returns the result of \code{shiny::runApp()}.
#' @export
run.vilma.app <- function(host = "127.0.0.1",
                          port = NULL,
                          launch.browser = TRUE,
                          display.mode = c("auto", "showcase"),
                          ...) {
  display.mode <- match.arg(display.mode)

  # Ensure shiny is available (keep shiny in Suggests in DESCRIPTION)
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required to run the app. Please install it first.")
  }

  # Locate app directory from the installed package
  app.dir <- system.file("app", package = "vilma")
  if (identical(app.dir, "") || !dir.exists(app.dir)) {
    stop("App directory not found inside the installed package (inst/app). ",
         "Did you install the package correctly?")
  }

  # Optional: quick sanity checks for assets you expect
  # if (!file.exists(file.path(app.dir, "app.R"))) {
  #   stop("app.R not found in inst/app.")
  # }

  shiny::runApp(
    appDir = app.dir,
    host = host,
    port = port,
    launch.browser = launch.browser,
    display.mode = display.mode,
    ...
  )
}

