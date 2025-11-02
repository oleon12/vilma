# R/mod_points_to_raster.R
mod_points_to_raster_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h4(br(), "Build ", em("vilma.dist"), " from occurrence points", br(), br()),
    fluidRow(
      column(
        3,
        div(
          class = "soft-card",
          fileInput(ns("file_points"), "Upload points CSV (Species, Lon, Lat)", accept = ".csv"),
          numericInput(ns("epsg"), "CRS (EPSG code)", value = 4326, min = 2000, step = 1),
          numericInput(ns("res"),  "Resolution (cell size)", value = 1, min = 0.0001, step = 1),
          checkboxInput(ns("sym"), "Force symmetrical pixels (square cells)", value = FALSE),
          actionButton(ns("run"), "Run analysis", class = "btn btn-primary w-100"),
          br(), br(),
          downloadButton(ns("download_rda"), "Download results (.rda)", class = "btn btn-sm btn-primary w-100"),
          br(), br(),
          actionButton(ns("next_btn"), "Next ▶︎", class = "btn btn-secondary w-100")
        )
      ),
      column(
        9,
        div(
          class = "right-pane",
          div(class = "map-wrap", leaflet::leafletOutput(ns("map_preview"), height = "100%")),
          div(
            class = "sum-wrap",
            tags$h4("Summary"),
            verbatimTextOutput(ns("dist_summary"))
          )
        )
      )
    )
  )
}

mod_points_to_raster_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    # tiny local helper: attach code lines to any object
    add_code <- function(obj, ...) {
      prev <- attr(obj, "vilma_code")
      new  <- unique(c(prev, unlist(list(...))))
      attr(obj, "vilma_code") <- new
      obj
    }

    has_fun  <- function(fn) exists(fn, mode = "function", inherits = TRUE)
    dist_obj <- reactiveVal(NULL)
    pts_df   <- reactiveVal(NULL)

    # Upload: read & store only
    observeEvent(input$file_points, {
      req(input$file_points)
      df <- tryCatch(
        readr::read_csv(input$file_points$datapath, show_col_types = FALSE),
        error = function(e) e
      )
      if (inherits(df, "error")) {
        showNotification(paste("Error reading CSV:", df$message), type = "error")
        pts_df(NULL)
        return()
      }
      pts_df(df)
      showNotification("Points loaded. Adjust parameters and click 'Run analysis'.", type = "message")
    })

    # Run analysis (with auto-closing spinner modal)
    observeEvent(input$run, {
      req(pts_df())

      # SHOW BUSY POPUP (auto-closes on exit)
      showModal(modalDialog(
        title = NULL, easyClose = FALSE, footer = NULL, size = "m",
        tagList(
          div(style="display:flex; align-items:center; gap:12px; padding:6px 2px;",
              tags$div(class="spinner-border", role="status", `aria-hidden`="true"),
              tags$strong("Building vilma.dist…")
          ),
          p("This may take a moment depending on the number of points, CRS, and grid resolution.",
            style="margin-top:8px;")
        )
      ))
      on.exit(removeModal(), add = TRUE)

      epsg <- as.integer(input$epsg)
      res  <- as.numeric(input$res)
      sym  <- isTRUE(input$sym)

      built <- tryCatch({
        if (has_fun("points_to_raster")) {
          points_to_raster(
            points      = pts_df(),
            crs         = epsg,
            ext         = NULL,
            res         = res,
            doRast      = TRUE,
            symmetrical = sym
          )
        } else if (has_fun("points.to.raster")) {
          points.to.raster(
            points      = pts_df(),
            crs         = epsg,
            ext         = NULL,
            res         = res,
            doRast      = TRUE,
            symmetrical = sym
          )
        } else if (has_fun("point.to.rast")) {
          point.to.rast(
            points      = pts_df(),
            crs         = epsg,
            ext         = NULL,
            res         = res,
            doRast      = TRUE,
            symmetrical = sym
          )
        } else {
          stop("points_to_raster / points.to.raster / point.to.rast not found in the session.")
        }
      }, error = function(e) e)

      if (inherits(built, "error")) {
        showNotification(paste("Build failed:", built$message), type = "error", duration = 10)
        dist_obj(NULL)
      } else {
        # ---- attach exact code lines for the exporter ----
        # NOTE: POINTS_CSV is a placeholder the exporter writes at the top of the script.
        pts_line <- "points <- readr::read_csv(POINTS_CSV, show_col_types = FALSE)"
        call_line <- if (has_fun("points_to_raster")) {
          sprintf(
            "dist <- points_to_raster(points, crs = %d, ext = NULL, res = %s, doRast = TRUE, symmetrical = %s)",
            epsg, deparse(res), if (sym) "TRUE" else "FALSE"
          )
        } else if (has_fun("points.to.raster")) {
          sprintf(
            "dist <- points.to.raster(points, crs = %d, ext = NULL, res = %s, doRast = TRUE, symmetrical = %s)",
            epsg, deparse(res), if (sym) "TRUE" else "FALSE"
          )
        } else {
          sprintf(
            "dist <- point.to.rast(points, crs = %d, ext = NULL, res = %s, doRast = TRUE, symmetrical = %s)",
            epsg, deparse(res), if (sym) "TRUE" else "FALSE"
          )
        }
        built <- add_code(built, pts_line, call_line)
        # --------------------------------------------------

        dist_obj(built)
        showNotification("vilma.dist built. Map and summary updated.", type = "message")
      }
    })

    # Map (¾)
    output$map_preview <- leaflet::renderLeaflet({
      obj <- dist_obj()
      req(obj)
      if (has_fun("view.vilma.dist")) return(view.vilma.dist(obj))
      leaflet::leaflet() |> leaflet::addTiles()
    })

    # Summary (¼)
    output$dist_summary <- renderPrint({
      obj <- dist_obj()
      req(obj)
      suppressMessages(print(obj))
    })

    # Download .rda
    output$download_rda <- downloadHandler(
      filename = function() "vilma_dist.rda",
      content = function(file) {
        obj <- dist_obj()
        req(obj)
        vilma_dist <- obj
        save(vilma_dist, file = file)
      }
    )

    # Optional guard for Next
    observeEvent(input$next_btn, {
      if (is.null(dist_obj())) {
        showNotification("Please run the analysis first.", type = "warning")
      }
    })

    # Expose downstream
    return(list(
      vilma_dist = reactive(dist_obj()),
      go_next    = reactive(input$next_btn)
    ))
  })
}

