# ======================= Beta Diversity UI =======================
mod_beta_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Beta diversity"),

    fluidRow(
      # Left: inputs
      column(
        3,
        div(class = "soft-card",
          h4("Inputs"),

          # Index selector (independent from Alpha / Nulls)
          selectInput(
            ns("index"), "Index",
            choices = c("Phylo Beta" = "phylo_beta",
                        "Phylosor"   = "phylosor",
                        "Rao β"      = "rao_beta",
                        "UniFrac"    = "unifrac",
                        "bMPD"       = "beta_mpd",
                        "bMNTD"      = "beta_mntd"),
            selected = "rao_beta"
          ),

          # ----- Per-index controls -----

          # phylo.beta
          conditionalPanel(
            condition = sprintf("input['%s'] == 'phylo_beta'", ns("index")),
            selectInput(ns("pb_method"), "Method",
                        choices = c("unweighted","weighted"), selected = "unweighted"),
            checkboxInput(ns("pb_normalize"), "Normalize", value = TRUE)
          ),

          # phylosor
          conditionalPanel(
            condition = sprintf("input['%s'] == 'phylosor'", ns("index")),
            selectInput(ns("ps_method"), "Faith method",
                        choices = c("root","node","exclude"), selected = "root"),
            checkboxInput(ns("ps_singleton_overlap"), "Singleton overlap", value = FALSE),
            checkboxInput(ns("ps_abundance"), "Abundance-weighted", value = FALSE),
            checkboxInput(ns("ps_normalize"), "Normalize", value = TRUE)
          ),

          # rao.beta
          conditionalPanel(
            condition = sprintf("input['%s'] == 'rao_beta'", ns("index")),
            checkboxInput(ns("rb_abundance"), "Abundance-weighted", value = FALSE),
            checkboxInput(ns("rb_scale01"), "Scale to [0,1]", value = TRUE)
          ),

          # unifrac
          conditionalPanel(
            condition = sprintf("input['%s'] == 'unifrac'", ns("index")),
            selectInput(ns("uf_method"), "Method",
                        choices = c("unweighted","weighted"), selected = "unweighted")
          ),

          # bMPD
          conditionalPanel(
            condition = sprintf("input['%s'] == 'beta_mpd'", ns("index")),
            selectInput(ns("bmpd_method"), "MPD method",
                        choices = c("exclude","root","node"), selected = "exclude"),
            checkboxInput(ns("bmpd_abundance"), "Abundance-weighted", value = FALSE),
            checkboxInput(ns("bmpd_exclude_conspecific"), "Exclude conspecific pairs", value = FALSE),
            checkboxInput(ns("bmpd_normalize"), "Normalize", value = FALSE),
            checkboxInput(ns("bmpd_scale01"), "Scale to [0,1]", value = FALSE)
          ),

          # bMNTD
          conditionalPanel(
            condition = sprintf("input['%s'] == 'beta_mntd'", ns("index")),
            selectInput(ns("bmntd_method"), "MNTD method",
                        choices = c("exclude","root","node"), selected = "exclude"),
            checkboxInput(ns("bmntd_abundance"), "Abundance-weighted", value = FALSE),
            checkboxInput(ns("bmntd_exclude_conspecific"), "Exclude conspecific pairs", value = FALSE),
            checkboxInput(ns("bmntd_normalize"), "Normalize", value = FALSE),
            checkboxInput(ns("bmntd_scale01"), "Scale to [0,1]", value = FALSE)
          ),

          br(),
          actionButton(ns("run"), "Run", class = "btn btn-primary w-100"),
          br(), br(),
          downloadButton(ns("download_beta"), "Download results (.rda)",
                         class = "btn btn-sm btn-primary w-100"),
          br(), br(),
          div(class = "d-flex gap-2",
            actionButton(ns("back_btn"), "◀︎ Back", class = "btn btn-secondary flex-fill"),
            actionButton(ns("next_btn"), "Next ▶︎", class = "btn btn-secondary flex-fill")
          )
        )
      ),

      # Right: results (map + summary)
      column(
        9,
        div(class = "right-pane",
          div(class = "map-wrap",
              uiOutput(ns("beta_view"))  # server chooses what to render (map)
          ),
          div(class = "sum-wrap",
            tags$h4("Result summary"),
            verbatimTextOutput(ns("beta_summary"))
          )
        )
      )
    )
  )
}

# ===================== Beta Diversity SERVER =====================
mod_beta_server <- function(id, vilma_dist, beta_tree) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns
    has_fun <- function(fn) exists(fn, mode = "function", inherits = TRUE)

    # ---- helper: attach code lines to an object (read by exporter) ----
    add_code <- function(obj, ...) {
      prev <- attr(obj, "vilma_code")
      new  <- unique(c(prev, unlist(list(...))))
      attr(obj, "vilma_code") <- new
      obj
    }

    # small helpers to format code arguments for readable R
    q  <- function(x) paste0("'", gsub("'", "\\\\'", x), "'")
    vc <- function(v) {
      if (is.null(v)) return("NULL")
      if (is.symbol(v)) return(as.character(v))
      if (is.language(v)) return(paste(deparse(v, width.cutoff = 500), collapse = ""))
      if (is.logical(v)) return(ifelse(v, "TRUE", "FALSE"))
      if (is.numeric(v)) {
        if (length(v) == 1) return(as.character(v))
        return(paste0("c(", paste(v, collapse = ", "), ")"))
      }
      if (is.character(v)) {
        if (length(v) == 1) return(q(v))
        return(paste0("c(", paste(q(v), collapse = ", "), ")"))
      }
      paste(deparse(v, width.cutoff = 500), collapse = "")
    }
    argline <- function(named_list) {
      paste(
        mapply(function(k, v) paste0(k, " = ", vc(v)), names(named_list), named_list, USE.NAMES = FALSE),
        collapse = ", "
      )
    }

    # canonicalize helper names to public API for export
    canon_fn <- function(fn) {
      fn <- sub("^phylosor\\.calc$", "phylosor", fn)
      fn <- sub("^unifrac\\.calc$", "unifrac", fn)
      fn # rao.beta, phylo.beta, beta.mpd, beta.mntd already canonical
    }

    result_obj <- reactiveVal(NULL)

    # Choose a dynamic view container (we'll render a leaflet map if possible)
    output$beta_view <- renderUI({
      leaflet::leafletOutput(ns("beta_map"), height = "100%")
    })

    observeEvent(input$run, {
      req(vilma_dist()); req(beta_tree())

      # === SHOW BUSY POPUP (auto-closes on exit) ===
      showModal(modalDialog(
        title = NULL, easyClose = FALSE, footer = NULL, size = "m",
        tagList(
          div(style="display:flex; align-items:center; gap:12px; padding:6px 2px;",
              tags$div(class="spinner-border", role="status", `aria-hidden`="true"),
              tags$strong("Calculating Beta Diversity…")
          ),
          p("This may take a moment depending on tree size, raster resolution, and index settings.",
            style="margin-top:8px;")
        )
      ))
      on.exit(removeModal(), add = TRUE)

      dist <- vilma_dist()
      tree <- beta_tree()
      idx  <- input$index

      # ---- Build arguments per index ----
      args <- switch(idx,
        phylo_beta = {
          stopifnot(has_fun("phylo.beta"))
          list(tree = tree, dist = dist,
               method = input$pb_method,
               normalize = isTRUE(input$pb_normalize))
        },
        phylosor = {
          stopifnot(has_fun("phylosor.calc"))
          list(tree = tree, dist = dist,
               method = input$ps_method,
               singleton_overlap = isTRUE(input$ps_singleton_overlap),
               abundance = isTRUE(input$ps_abundance),
               normalize = isTRUE(input$ps_normalize))
        },
        rao_beta = {
          stopifnot(has_fun("rao.beta"))
          list(tree = tree, dist = dist,
               abundance = isTRUE(input$rb_abundance),
               scale01   = isTRUE(input$rb_scale01))
        },
        unifrac = {
          stopifnot(has_fun("unifrac.calc"))
          list(tree = tree, dist = dist,
               method = input$uf_method)
        },
        beta_mpd = {
          stopifnot(has_fun("beta.mpd"))
          list(tree = tree, dist = dist,
               mpd.method = input$bmpd_method,
               abundance  = isTRUE(input$bmpd_abundance),
               exclude.conspecific = isTRUE(input$bmpd_exclude_conspecific),
               normalize  = isTRUE(input$bmpd_normalize),
               scale01    = isTRUE(input$bmpd_scale01))
        },
        beta_mntd = {
          stopifnot(has_fun("beta.mntd"))
          list(tree = tree, dist = dist,
               mntd.method = input$bmntd_method,
               abundance   = isTRUE(input$bmntd_abundance),
               exclude.conspecific = isTRUE(input$bmntd_exclude_conspecific),
               normalize   = isTRUE(input$bmntd_normalize),
               scale01     = isTRUE(input$bmntd_scale01))
        },
        stop("Unknown beta index.")
      )

      # ---- Compute ----
      fun <- switch(idx,
        phylo_beta = phylo.beta,
        phylosor   = phylosor.calc,
        rao_beta   = rao.beta,
        unifrac    = unifrac.calc,
        beta_mpd   = beta.mpd,
        beta_mntd  = beta.mntd
      )

      res <- tryCatch(do.call(fun, args), error = function(e) e)

      if (inherits(res, "error")) {
        showNotification(paste("Beta diversity failed:", conditionMessage(res)),
                         type = "error", duration = 10)
        result_obj(NULL)
        return(invisible())
      }

      # -------- attach exact beta code line + structured meta for the exporter --------
      raw_fun_name <- switch(idx,
        phylo_beta = "phylo.beta",
        phylosor   = "phylosor.calc",
        rao_beta   = "rao.beta",
        unifrac    = "unifrac.calc",
        beta_mpd   = "beta.mpd",
        beta_mntd  = "beta.mntd"
      )
      fun_name <- canon_fn(raw_fun_name)

      # readable printed args: show symbols for tree/dist
      print_args <- args
      print_args$tree <- quote(tree)
      print_args$dist <- quote(dist)

      code_lines <- c(
        "# --- Beta Diversity ---",
        paste0("beta <- ", fun_name, "(", argline(print_args), ")")
      )
      res <- add_code(res, code_lines)

      # structured fallback so exporter can reconstruct even if code attr is lost
      meta_args <- args
      meta_args$tree <- "tree"
      meta_args$dist <- "dist"
      attr(res, "vilma_meta") <- list(fn = fun_name, args = meta_args)
      # ----------------------------------------------------------------

      result_obj(res)
      showNotification("Beta diversity computed.", type = "message")
    })

    # ---- MAP: use view.vilma.beta if available; fallback to raster ----
    output$beta_map <- leaflet::renderLeaflet({
      req(result_obj())
      res <- result_obj()

      if (has_fun("view.vilma.beta")) {
        m <- try(view.vilma.beta(res), silent = TRUE)
        if (!inherits(m, "try-error")) return(m)
      }

      # Fallbacks: try common raster fields
      r <- NULL
      if (!is.null(res$r.raster)    && inherits(res$r.raster,    "SpatRaster")) r <- res$r.raster
      if (is.null(r) && !is.null(res$beta.raster) && inherits(res$beta.raster, "SpatRaster")) r <- res$beta.raster

      if (!is.null(r)) {
        ext  <- as.vector(terra::ext(r))
        vals <- try(terra::values(r[, drop = FALSE]), silent = TRUE)
        if (inherits(vals, "try-error")) vals <- try(terra::values(r), silent = TRUE)
        vals <- as.numeric(vals); vals <- vals[is.finite(vals)]
        pal  <- leaflet::colorNumeric("viridis", domain = if (length(vals)) vals else c(0, 1), na.color = NA)

        return(
          leaflet::leaflet() |>
            leaflet::addTiles() |>
            leaflet::fitBounds(ext[1], ext[3], ext[2], ext[4]) |>
            leaflet::addRasterImage(r, project = TRUE, opacity = 0.85, colors = pal) |>
            leaflet::addLegend(pal = pal, values = if (length(vals)) vals else c(0, 1), title = "Beta")
        )
      }

      # If no raster, just show a base map (summary will still print)
      leaflet::leaflet() |> leaflet::addTiles()
    })

    # ---- SUMMARY ----
    output$beta_summary <- renderPrint({
      req(result_obj())
      res <- result_obj()

      # Prefer the print method if available
      ok <- try(suppressMessages(print(res)), silent = TRUE)
      if (!inherits(ok, "try-error")) return(invisible())

      # Fallback: structure dump
      str(res)
    })

    # ---- DOWNLOAD ----
    output$download_beta <- downloadHandler(
      filename = function() {
        idx <- input$index
        paste0("beta_", idx, ".rda")
      },
      content = function(file) {
        res <- result_obj(); req(res)
        vilma_beta <- res
        save(vilma_beta, file = file)
      }
    )

    # Expose nav hooks (optional; mirrors other modules)
    return(list(
      beta_result = reactive(result_obj()),
      go_back     = reactive(input$back_btn),
      go_next     = reactive(input$next_btn)
    ))
  })
}

