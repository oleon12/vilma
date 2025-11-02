# ============================
# Null Models UI
# ============================
mod_nulls_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Null models"),

    fluidRow(
      # Left: inputs
      column(
        3,
        div(class = "soft-card",
          h4("Inputs"),

          # Index (locked by server to Alpha’s choice)
          selectInput(
            ns("index"), "Index (locked to Alpha)",
            choices = c("Faith PD" = "faith_pd",
                        "MPD"      = "mpd",
                        "MNTD"     = "mntd",
                        "PE"       = "pe",
                        "RaoQ"     = "rao",
                        "NRI"      = "nri",
                        "NNI"      = "nti"),
            selected = "mpd"
          ),

          numericInput(ns("iterations"), "Iterations", value = 999, min = 1, step = 1),
          selectInput(ns("method"), "Null method", choices = c("cell","global"), selected = "cell"),
          selectInput(ns("sampling"), "Null sampling",
                      choices = c("taxa.label","range","neighbor","regional"),
                      selected = "taxa.label"),

          # neighbor/regional extras
          conditionalPanel(
            condition = sprintf("input['%s'] == 'neighbor' || input['%s'] == 'regional'",
                                ns("sampling"), ns("sampling")),
            selectInput(ns("n_directions"), "Neighborhood",
                        choices = c("rook","bishop","queen"), selected = "rook")
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'regional'", ns("sampling")),
            selectInput(ns("regional_weight"), "Regional weight",
                        choices = c("uniform","frequency","range"), selected = "uniform")
          ),

          # Rao-only params
          conditionalPanel(
            condition = sprintf("input['%s'] == 'rao'", ns("index")),
            checkboxInput(ns("rao_abundance"), "Abundance-weighted", value = FALSE),
            checkboxInput(ns("rao_scale01"), "Scale to [0,1]", value = TRUE)
          ),

          br(),
          actionButton(ns("run"), "Run", class = "btn btn-primary w-100"),
          br(), br(),
          downloadButton(ns("download_null"), "Download results (.rda)",
                         class = "btn btn-sm btn-primary w-100"),
          br(), br(),
          div(class = "d-flex gap-2",
            actionButton(ns("back_btn"), "◀︎ Back", class = "btn btn-secondary flex-fill"),
            actionButton(ns("next_btn"), "Next ▶︎", class = "btn btn-secondary flex-fill")
          )
        )
      ),

      # Right: results (map OR plot + summary)
      column(
        9,
        div(class = "right-pane",
          conditionalPanel(
            condition = sprintf("input['%s'] == 'cell'", ns("method")),
            leaflet::leafletOutput(ns("null_map"), height = "520px")
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'global'", ns("method")),
            plotOutput(ns("null_plot"), height = "520px")
          ),
          tags$div(style="height: 12px;"),
          div(class = "sum-wrap",
            tags$h4("Result summary"),
            verbatimTextOutput(ns("null_summary"))
          )
        )
      )
    )
  )
}

# ============================
# Null Models SERVER
# ============================
mod_nulls_server <- function(id, vilma_dist, alpha_result, alpha_tree, alpha_index) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns
    has_fun <- function(fn) exists(fn, mode = "function", inherits = TRUE)

    # ---- helper: attach code lines to an object (for exporter) ----
    add_code <- function(obj, ...) {
      prev <- attr(obj, "vilma_code")
      new  <- unique(c(prev, unlist(list(...))))
      attr(obj, "vilma_code") <- new
      obj
    }

    # small helpers to format code arguments
    q  <- function(x) paste0("'", gsub("'", "\\\\'", x), "'")
    vc <- function(v) {
      if (is.null(v)) return("NULL")
      if (is.symbol(v)) return(as.character(v))
      if (is.language(v)) return(paste(deparse(v, width.cutoff = 500), collapse = ""))
      if (is.logical(v)) return(ifelse(v, "TRUE", "FALSE"))
      if (is.numeric(v)) return(paste(v, collapse = ", "))
      if (is.character(v)) return(paste(q(v), collapse = ", "))
      paste(deparse(v, width.cutoff = 500), collapse = "")
    }
    argline <- function(named_list) {
      paste(
        mapply(function(k, v) paste0(k, " = ", vc(v)), names(named_list), named_list, USE.NAMES = FALSE),
        collapse = ", "
      )
    }

    result_obj   <- reactiveVal(NULL)
    go_beta_trig <- reactiveVal(0)

    # Next button behavior for NRI/NTI
    observeEvent(input$next_btn, {
      if (input$index %in% c("nri","nti")) {
        showNotification("NRI/NNI use MPD/MNTD nulls. Moving to Beta.", type = "warning")
        go_beta_trig(isolate(go_beta_trig()) + 1)
      }
    }, ignoreInit = TRUE)

    # Lock index to Alpha’s choice where applicable
    observeEvent(alpha_index(), {
      idx <- alpha_index()
      if (is.null(idx)) return()
      allowed <- c("faith_pd","mpd","mntd","pe","rao")
      if (idx %in% allowed) {
        locked_choices <- setNames(idx, switch(
          idx,
          faith_pd = "Faith PD", mpd = "MPD", mntd = "MNTD", pe = "PE", rao = "RaoQ"
        ))
        updateSelectInput(session, "index", choices = locked_choices, selected = idx)
        showNotification(sprintf("Null model locked to Alpha: %s", names(locked_choices)),
                         type = "message", duration = 3)
      } else {
        showNotification("NRI/NTI have no dedicated null model. Redirecting to Beta.",
                         type = "warning", duration = 5)
        go_beta_trig(isolate(go_beta_trig()) + 1)
      }
    }, ignoreInit = FALSE)

    # =====================
    # Run (with spinner modal like Alpha/Beta)
    # =====================
    observeEvent(input$run, {
      req(vilma_dist()); req(alpha_tree())

      # Spinner popup in the new unified style
      showModal(modalDialog(
        title = NULL, easyClose = FALSE, footer = NULL, size = "m",
        tagList(
          div(style="display:flex; align-items:center; gap:12px; padding:6px 2px;",
              tags$div(class="spinner-border", role="status", `aria-hidden`="true"),
              tags$strong("Calculating Null Models…")
          ),
          p("This may take a moment depending on iterations, method (cell/global), and sampling.",
            style="margin-top:8px;")
        )
      ))
      on.exit(removeModal(), add = TRUE)

      # Do the compute synchronously
      res <- tryCatch({
        dist  <- vilma_dist()
        tree  <- alpha_tree()
        idx   <- input$index
        alpha <- isolate(alpha_result())

        if (idx %in% c("nri","nti")) stop("NRI/NNI have no dedicated nulls. Choose MPD/MNTD here, or go to Beta.")

        ndir <- if (!is.null(input$n_directions) && nzchar(input$n_directions)) input$n_directions else "rook"
        rwt  <- if (!is.null(input$regional_weight) && nzchar(input$regional_weight)) input$regional_weight else "uniform"

        base_args <- list(
          pd                = alpha,
          tree              = tree,
          dist              = dist,
          iterations        = as.integer(input$iterations),
          method            = input$method,
          sampling          = input$sampling,
          `n.directions`    = ndir,
          `regional.weight` = rwt
        )

        fun <- NULL
        args <- base_args
        null_fun_export <- NULL

        if (identical(idx, "faith_pd")) {
          stopifnot(has_fun("faith.pd.null")); fun <- faith.pd.null; null_fun_export <- "faith.pd.null"
        } else if (identical(idx, "mpd")) {
          stopifnot(has_fun("mpd.calc.null")); fun <- mpd.calc.null; null_fun_export <- "mpd.calc.null"
        } else if (identical(idx, "mntd")) {
          stopifnot(has_fun("mntd.calc.null")); fun <- mntd.calc.null; null_fun_export <- "mntd.calc.null"
        } else if (identical(idx, "pe")) {
          stopifnot(has_fun("pe.calc.null"));  fun <- pe.calc.null;  null_fun_export <- "pe.calc.null"
        } else if (identical(idx, "rao")) {
          stopifnot(has_fun("rao.calc.null")); fun <- rao.calc.null; null_fun_export <- "rao.calc.null"
          args$abundance <- isTRUE(input$rao_abundance)
          args$scale01   <- isTRUE(input$rao_scale01)
        } else stop("Unknown index.")

        ans <- do.call(fun, args)

        # Attach reproducible code/meta
        print_args <- args
        print_args$pd   <- quote(alpha)
        print_args$tree <- quote(tree)
        print_args$dist <- quote(dist)

        ans <- add_code(ans, c(
          "# --- Null Models ---",
          paste0("nulls <- ", null_fun_export, "(", argline(print_args), ")")
        ))

        meta_args <- args
        meta_args$pd   <- "alpha"
        meta_args$tree <- "tree"
        meta_args$dist <- "dist"
        attr(ans, "vilma_meta") <- list(fn = null_fun_export, args = meta_args)

        ans
      }, error = function(e) e)

      if (inherits(res, "error")) {
        showNotification(paste("Null model failed:", conditionMessage(res)), type = "error", duration = 10)
        result_obj(NULL)
      } else {
        result_obj(res)
        showNotification("Null model computed.", type = "message")
      }
    })

    # =====================
    # RENDERERS
    # =====================

    # CELL map
    output$null_map <- leaflet::renderLeaflet({
      req(result_obj()); req(identical(input$method, "cell"))
      res <- result_obj()

      if (has_fun("view.vilma.null")) {
        m <- try(view.vilma.null(res), silent = TRUE)
        if (!inherits(m, "try-error")) return(m)
      }

      # Fallback simple raster map
      r <- NULL
      if (!is.null(res$r.raster)  && inherits(res$r.raster,  "SpatRaster")) r <- res$r.raster
      if (is.null(r) && !is.null(res$pd.raster) && inherits(res$pd.raster, "SpatRaster")) r <- res$pd.raster
      validate(need(!is.null(r), "No raster to display."))

      ext  <- as.vector(terra::ext(r))
      vals <- try(terra::values(r[, drop = FALSE]), silent = TRUE)
      if (inherits(vals, "try-error")) vals <- try(terra::values(r), silent = TRUE)
      vals <- as.numeric(vals); vals <- vals[is.finite(vals)]

      pal <- leaflet::colorNumeric("viridis", domain = if (length(vals)) vals else c(0, 1), na.color = NA)

      leaflet::leaflet() |>
        leaflet::addTiles() |>
        leaflet::fitBounds(ext[1], ext[3], ext[2], ext[4]) |>
        leaflet::addRasterImage(r, project = TRUE, opacity = 0.85, colors = pal) |>
        leaflet::addLegend(pal = pal, values = if (length(vals)) vals else c(0, 1), title = "Null (cell)")
    })

    # GLOBAL plot
    output$null_plot <- renderPlot({
      req(result_obj()); req(identical(input$method, "global"))
      res <- result_obj()

      # Prefer your S3 method if present; otherwise safe fallback
      ok <- try(plot(res), silent = TRUE)
      if (inherits(ok, "try-error")) {
        if (!is.null(res$null_values) && is.numeric(res$null_values)) {
          hist(res$null_values, main = "Null distribution", xlab = "", ylab = "Freq")
          if (!is.null(res$observed) && is.numeric(res$observed)) abline(v = res$observed, lwd = 2)
        }
      }
    })

    # SUMMARY
    output$null_summary <- renderPrint({
      req(result_obj())
      x <- result_obj()

      ok <- try(suppressMessages(print(x)), silent = TRUE)
      if (!inherits(ok, "try-error")) return(invisible())

      obs <- tryCatch(x$observed, error = function(e) NULL)
      nv  <- tryCatch(if (!is.null(x$null_values)) x$null_values else x$null, error = function(e) NULL)

      if (is.numeric(nv) && length(nv) > 1) {
        m   <- mean(nv, na.rm = TRUE)
        s   <- stats::sd(nv, na.rm = TRUE)
        ses <- if (is.finite(s) && s > 0) (obs - m) / s else NA_real_
        p_ge <- mean(nv >= obs, na.rm = TRUE)
        p_le <- mean(nv <= obs, na.rm = TRUE)

        cat("Index:", isolate(input$index), "\n",
            "Method:", isolate(input$method), "\n",
            "Sampling:", isolate(input$sampling), "\n",
            "Iterations:", isolate(as.integer(input$iterations)), "\n",
            "Observed:", sprintf("%.6g", obs), "\n",
            "Null mean:", sprintf("%.6g", m), "  SD:", sprintf("%.6g", s), "\n",
            "SES:", sprintf("%.6g", ses), "\n",
            "P(>=obs):", sprintf("%.6g", p_ge), "  P(<=obs):", sprintf("%.6g", p_le), "\n", sep = "")
        return(invisible())
      }

      is_spat <- function(z) inherits(z, "SpatRaster")
      p_r  <- tryCatch(x$p.raster,  error = function(e) NULL)
      sesr <- tryCatch(x$ses.raster, error = function(e) NULL)
      idxr <- tryCatch(if (!is.null(x$pd.raster)) x$pd.raster else x$r.raster, error = function(e) NULL)

      if (is_spat(p_r) || is_spat(sesr) || is_spat(idxr)) {
        g <- function(r, fn) tryCatch(terra::global(r, fn, na.rm = TRUE)[[1]], error = function(e) NA_real_)
        cat("Index:", isolate(input$index), "\n",
            "Method:", isolate(input$method), "\n",
            "Sampling:", isolate(input$sampling), "\n", sep = "")
        if (is_spat(idxr)) {
          cat("--- Index raster ---\n",
              "  mean:", sprintf("%.6g", g(idxr,"mean")),
              " sd:",   sprintf("%.6g", g(idxr,"sd")),
              " min:",  sprintf("%.6g", g(idxr,"min")),
              " max:",  sprintf("%.6g", g(idxr,"max")), "\n")
        }
        if (is_spat(sesr)) {
          cat("--- SES raster ---\n",
              "  mean:", sprintf("%.6g", g(sesr,"mean")),
              " sd:",   sprintf("%.6g", g(sesr,"sd")),
              " min:",  sprintf("%.6g", g(sesr,"min")),
              " max:",  sprintf("%.6g", g(sesr,"max")), "\n")
        }
        if (is_spat(p_r)) {
          pvals <- tryCatch(as.numeric(terra::values(p_r)), error = function(e) NULL)
          if (is.numeric(pvals)) {
            pvals <- pvals[is.finite(pvals)]
            cat("--- P raster ---\n",
                "  mean:", sprintf("%.6g", g(p_r,"mean")),
                "  (<=0.05):", sum(pvals <= 0.05, na.rm = TRUE),
                "  (>=0.95):", sum(pvals >= 0.95, na.rm = TRUE), "\n")
          }
        }
        return(invisible())
      }

      str(x, max.level = 1)
      invisible()
    })

    # Keep outputs active when hidden by conditionalPanel
    outputOptions(output, "null_map",     suspendWhenHidden = FALSE)
    outputOptions(output, "null_plot",    suspendWhenHidden = FALSE)
    outputOptions(output, "null_summary", suspendWhenHidden = FALSE)

    # DOWNLOAD
    output$download_null <- downloadHandler(
      filename = function() paste0("null_", input$index, "_", input$method, ".rda"),
      content = function(file) {
        res <- result_obj(); req(res)
        vilma_null <- res
        save(vilma_null, file = file)
      }
    )

    return(list(
      null_result = reactive(result_obj()),
      go_back     = reactive(input$back_btn),
      go_next     = reactive(input$next_btn),
      go_beta     = eventReactive(go_beta_trig(), TRUE, ignoreInit = TRUE)
    ))
  })
}

