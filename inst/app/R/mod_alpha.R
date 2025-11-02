# ===================== mod_alpha.R =====================
# Safe infix coalesce used by exporter code
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

# ----------------------- UI -----------------------
mod_alpha_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Phylogenetic Diversity"),

    fluidRow(
      # Left: Inputs
      column(
        3,
        div(class = "soft-card",
          h4("Inputs"),

          # Upload tree
          fileInput(ns("tree_file"), "Upload tree (.tre / .tree / Newick)"),

          # Index selector (same keys as server switch)
          selectInput(
            ns("idx"), "Index",
            choices = c("Faith PD" = "faith_pd",
                        "MPD"      = "mpd",
                        "MNTD"     = "mntd",
                        "PE"       = "pe",
                        "RaoQ"     = "rao",
                        "NRI"      = "nri",
                        "NNI"      = "nti"),
            selected = "faith_pd"
          ),

          # ----- Per-index controls (kept minimal & robust) -----

          # Faith PD
          conditionalPanel(
            condition = sprintf("input['%s'] == 'faith_pd'", ns("idx")),
            selectInput(ns("faith_method"), "Faith method",
                        choices = c("root","node","exclude"), selected = "root")
          ),

          # MPD
          conditionalPanel(
            condition = sprintf("input['%s'] == 'mpd'", ns("idx")),
            selectInput(ns("mpd_method"), "MPD method",
                        choices = c("root","node","exclude"), selected = "root"),
            checkboxInput(ns("mpd_abundance"), "Abundance-weighted", value = FALSE)
          ),

          # MNTD
          conditionalPanel(
            condition = sprintf("input['%s'] == 'mntd'", ns("idx")),
            selectInput(ns("mntd_method"), "MNTD method",
                        choices = c("root","node","exclude"), selected = "exclude"),
            checkboxInput(ns("mntd_abundance"), "Abundance-weighted", value = FALSE)
          ),

          # NRI (wrapper for MPD)
          conditionalPanel(
            condition = sprintf("input['%s'] == 'nri'", ns("idx")),
            selectInput(ns("nri_mpd_method"), "MPD method",
                        choices = c("root","node","exclude"), selected = "root"),
            checkboxInput(ns("nri_abundance"), "Abundance-weighted", value = FALSE),
            numericInput(ns("nri_iterations"), "Iterations", value = 999, min = 1, step = 1),
            selectInput(ns("nri_sampling"), "Null sampling",
                        choices = c("taxa.label","range","neighbor","regional"),
                        selected = "taxa.label"),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'neighbor' || input['%s'] == 'regional'",
                                  ns("nri_sampling"), ns("nri_sampling")),
              selectInput(ns("nri_n_directions"), "Neighborhood",
                          choices = c("rook","bishop","queen"), selected = "rook")
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'regional'", ns("nri_sampling")),
              selectInput(ns("nri_regional_weight"), "Regional weight",
                          choices = c("uniform","frequency","range"), selected = "uniform")
            )
          ),

          # NTI (wrapper for MNTD)
          conditionalPanel(
            condition = sprintf("input['%s'] == 'nti'", ns("idx")),
            selectInput(ns("nti_mntd_method"), "MNTD method",
                        choices = c("exclude","root","node"), selected = "exclude"),
            checkboxInput(ns("nti_abundance"), "Abundance-weighted", value = FALSE),
            numericInput(ns("nti_iterations"), "Iterations", value = 999, min = 1, step = 1),
            selectInput(ns("nti_sampling"), "Null sampling",
                        choices = c("taxa.label","range","neighbor","regional"),
                        selected = "taxa.label"),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'neighbor' || input['%s'] == 'regional'",
                                  ns("nti_sampling"), ns("nti_sampling")),
              selectInput(ns("nti_n_directions"), "Neighborhood",
                          choices = c("rook","bishop","queen"), selected = "rook")
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'regional'", ns("nti_sampling")),
              selectInput(ns("nti_regional_weight"), "Regional weight",
                          choices = c("uniform","frequency","range"), selected = "uniform")
            )
          ),

          # PE
          conditionalPanel(
            condition = sprintf("input['%s'] == 'pe'", ns("idx")),
            checkboxInput(ns("pe_RPE"), "Relative PE (RPE)", value = FALSE),
            selectInput(ns("pe_faith_method"), "Faith method",
                        choices = c("root","node","exclude"), selected = "root")
          ),

          # RaoQ
          conditionalPanel(
            condition = sprintf("input['%s'] == 'rao'", ns("idx")),
            checkboxInput(ns("rao_abundance"), "Abundance-weighted", value = FALSE),
            checkboxInput(ns("rao_scale01"), "Scale to [0,1]", value = TRUE)
          ),

          br(),
          actionButton(ns("run"), "Run", class = "btn btn-primary w-100"),
          br(), br(),
          downloadButton(ns("download_alpha"), "Download results (.rda)",
                         class = "btn btn-sm btn-primary w-100"),
          br(), br(),
          div(class = "d-flex gap-2",
            actionButton(ns("back_btn"), "◀︎ Back", class = "btn btn-secondary flex-fill"),
            actionButton(ns("next_btn"), "Next ▶︎", class = "btn btn-secondary flex-fill")
          )
        )
      ),

      # Right: Results (map + summary)
      column(
        9,
        div(class = "right-pane",
          div(class = "map-wrap",
              uiOutput(ns("alpha_view"))  # leafletOutput injected by server
          ),
          div(class = "sum-wrap",
            tags$h4("Result summary"),
            verbatimTextOutput(ns("alpha_summary"))
          )
        )
      )
    )
  )
}

# --------------------- SERVER ---------------------
mod_alpha_server <- function(id, vilma_dist) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns
    has_fun <- function(fn) exists(fn, mode = "function", inherits = TRUE)

    # --- Track state for navigation logic ---
    last_run_idx <- reactiveVal(NULL)  # which index actually ran successfully
    result_obj   <- reactiveVal(NULL)
    tree_obj     <- reactiveVal(NULL)

    # ---- helper: attach code lines to the result (used by exporter) ----
    add_code <- function(obj, ...) {
      prev <- attr(obj, "vilma_code")
      new  <- unique(c(prev, unlist(list(...))))
      attr(obj, "vilma_code") <- new
      obj
    }

    # small helpers to format readable R code
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

    # Canonicalize helper names to public API for export
    canon_fn <- function(fn) {
      switch(fn,
        "mpd.calc"  = "mpd",
        "mntd.calc" = "mntd",
        "nri.calc"  = "nri",
        "nti.calc"  = "nti",
        "pe.calc"   = "pe",
        "rao.calc"  = "rao.alpha",
        fn # faith.pd stays faith.pd
      )
    }

    # --- TREE UPLOAD ---
    observeEvent(input$tree_file, {
      req(input$tree_file)
      tr <- tryCatch(ape::read.tree(input$tree_file$datapath),
                     error = function(e) e)
      if (inherits(tr, "error")) {
        showNotification(paste("Error reading tree:", tr$message), type = "error")
        tree_obj(NULL)
      } else {
        tree_obj(tr)
        showNotification("Tree loaded.", type = "message")
      }
    })

    # --- Dynamic result container (like Nulls) ---
    output$alpha_view <- renderUI({
      leaflet::leafletOutput(ns("alpha_map"), height = "100%")
    })

    # --- RUN ---
    observeEvent(input$run, {
      req(vilma_dist()); req(tree_obj())

      # SHOW BUSY POPUP (auto-closes on exit)
      showModal(modalDialog(
        title = NULL, easyClose = FALSE, footer = NULL, size = "m",
        tagList(
          div(style="display:flex; align-items:center; gap:12px; padding:6px 2px;",
              tags$div(class="spinner-border", role="status", `aria-hidden`="true"),
              tags$strong("Calculating Alpha Index…")
          ),
          p("This may take a moment depending on tree size, raster resolution, and abundance settings.",
            style="margin-top:8px;")
        )
      ))
      on.exit(removeModal(), add = TRUE)

      dist <- vilma_dist()
      tree <- tree_obj()
      idx  <- input$idx

      res <- tryCatch({
        switch(idx,
          faith_pd = {
            stopifnot(has_fun("faith.pd"))
            faith.pd(tree = tree, dist = dist, method = input$faith_method)
          },
          mpd = {
            stopifnot(has_fun("mpd.calc"))
            mpd.calc(tree = tree, dist = dist,
                     method    = input$mpd_method,
                     abundance = isTRUE(input$mpd_abundance))
          },
          mntd = {
            stopifnot(has_fun("mntd.calc"))
            mntd.calc(tree = tree, dist = dist,
                      method    = input$mntd_method,
                      abundance = isTRUE(input$mntd_abundance))
          },
          nri = {
            stopifnot(has_fun("nri.calc"))
            nri.calc(tree = tree, dist = dist,
                     mpd.method        = input$nri_mpd_method,
                     abundance         = isTRUE(input$nri_abundance),
                     iterations        = as.integer(input$nri_iterations),
                     sampling          = input$nri_sampling,
                     `n.directions`    = input$nri_n_directions,
                     `regional.weight` = input$nri_regional_weight)
          },
          nti = {
            stopifnot(has_fun("nti.calc"))
            nti.calc(tree = tree, dist = dist,
                     mntd.method       = input$nti_mntd_method,
                     abundance         = isTRUE(input$nti_abundance),
                     iterations        = as.integer(input$nti_iterations),
                     sampling          = input$nti_sampling,
                     `n.directions`    = input$nti_n_directions,
                     `regional.weight` = input$nti_regional_weight)
          },
          pe = {
            stopifnot(has_fun("pe.calc"))
            pe.calc(tree = tree, dist = dist,
                    RPE          = isTRUE(input$pe_RPE),
                    faith.method = input$pe_faith_method)
          },
          rao = {
            stopifnot(has_fun("rao.calc"))
            rao.calc(tree = tree, dist = dist,
                     abundance = isTRUE(input$rao_abundance),
                     scale01   = isTRUE(input$rao_scale01))
          },
          stop("Unknown index.")
        )
      }, error = function(e) e)

      if (inherits(res, "error")) {
        showNotification(paste("Alpha index failed:", res$message), type = "error")
        result_obj(NULL)
        last_run_idx(NULL)                # clear on failure
      } else {
        # -------- attach exact alpha code line + structured meta for exporter --------
        raw_fun <- switch(idx,
          faith_pd = "faith.pd",
          mpd      = "mpd.calc",
          mntd     = "mntd.calc",
          nri      = "nri.calc",
          nti      = "nti.calc",
          pe       = "pe.calc",
          rao      = "rao.calc"
        )
        fun_name <- canon_fn(raw_fun)

        # build printed args (use symbols for tree/dist)
        args_print <- switch(idx,
          faith_pd = list(tree = quote(tree), dist = quote(dist), method = input$faith_method),
          mpd      = list(tree = quote(tree), dist = quote(dist),
                          method = input$mpd_method, abundance = isTRUE(input$mpd_abundance)),
          mntd     = list(tree = quote(tree), dist = quote(dist),
                          method = input$mntd_method, abundance = isTRUE(input$mntd_abundance)),
          nri      = list(tree = quote(tree), dist = quote(dist),
                          mpd.method = input$nri_mpd_method,
                          abundance = isTRUE(input$nri_abundance),
                          iterations = as.integer(input$nri_iterations),
                          sampling = input$nri_sampling,
                          `n.directions` = input$nri_n_directions,
                          `regional.weight` = input$nri_regional_weight),
          nti      = list(tree = quote(tree), dist = quote(dist),
                          mntd.method = input$nti_mntd_method,
                          abundance = isTRUE(input$nti_abundance),
                          iterations = as.integer(input$nti_iterations),
                          sampling = input$nti_sampling,
                          `n.directions` = input$nti_n_directions,
                          `regional.weight` = input$nti_regional_weight),
          pe       = list(tree = quote(tree), dist = quote(dist),
                          RPE = isTRUE(input$pe_RPE), faith.method = input$pe_faith_method),
          rao      = list(tree = quote(tree), dist = quote(dist),
                          abundance = isTRUE(input$rao_abundance), scale01 = isTRUE(input$rao_scale01))
        )

        # exact code line (canonical function name)
        code_lines <- c(
          "# --- Alpha Diversity ---",
          paste0("alpha <- ", fun_name, "(", argline(args_print), ")")
        )
        res <- add_code(res, code_lines)

        # structured fallback (meta) so exporter can reconstruct even if vilma_code is lost
        meta_args <- args_print
        meta_args$tree <- "tree"
        meta_args$dist <- "dist"
        attr(res, "vilma_meta") <- list(fn = fun_name, args = meta_args)
        # ----------------------------------------------------------------

        result_obj(res)
        last_run_idx(idx)                 # remember which index actually ran
        showNotification("Alpha index computed.", type = "message")
      }
    })

    # --- VIEW: map (leaflet) ---
    output$alpha_map <- leaflet::renderLeaflet({
      req(result_obj())
      res <- result_obj()

      if (has_fun("view.vilma.pd") && inherits(res, "vilma.pd")) {
        m <- try(view.vilma.pd(res), silent = TRUE)
        if (!inherits(m, "try-error")) return(m)
      }

      if (!is.null(res$r.raster) && inherits(res$r.raster, "SpatRaster")) {
        r <- res$r.raster
        ext <- as.vector(terra::ext(r))
        vals <- try(terra::values(r[, drop = FALSE]), silent = TRUE)
        if (inherits(vals, "try-error")) vals <- try(terra::values(r), silent = TRUE)
        vals <- as.numeric(vals); vals <- vals[is.finite(vals)]
        pal  <- leaflet::colorNumeric("viridis", domain = if (length(vals)) vals else c(0,1), na.color = NA)

        leaflet::leaflet() |>
          leaflet::addTiles() |>
          leaflet::fitBounds(ext[1], ext[3], ext[2], ext[4]) |>
          leaflet::addRasterImage(r, project = TRUE, opacity = 0.85, colors = pal) |>
          leaflet::addLegend(pal = pal, values = if (length(vals)) vals else c(0,1),
                    title = input$idx %||% "Alpha")
      } else {
        leaflet::leaflet() |> leaflet::addTiles()
      }
    })

    # --- SUMMARY (print method first; fallback to raster stats) ---
    output$alpha_summary <- renderPrint({
      req(result_obj())
      rs <- result_obj()

      ok <- try(suppressMessages(print(rs)), silent = TRUE)
      if (!inherits(ok, "try-error")) return(invisible())

      if (!is.null(rs$r.raster) && inherits(rs$r.raster, "SpatRaster")) {
        v <- terra::values(rs$r.raster)
        v <- v[is.finite(v)]
        if (length(v)) {
          qs <- stats::quantile(v, probs = c(0, .25, .5, .75, 1), na.rm = TRUE, names = FALSE)
          cat("Cells:", terra::ncell(rs$r.raster), "\n",
              "Mean:", mean(v), " SD:", stats::sd(v), "\n",
              "Min:", qs[1], " Q25:", qs[2], " Median:", qs[3],
              " Q75:", qs[4], " Max:", qs[5], "\n")
        } else {
          cat("[Raster contains no finite values]\n")
        }
      } else {
        str(rs)
      }
    })

    # --- DOWNLOAD ---
    output$download_alpha <- downloadHandler(
      filename = function() paste0("alpha_", input$idx, ".rda"),
      content = function(file) {
        res <- result_obj(); req(res)
        alpha_result <- res
        save(alpha_result, file = file)
      }
    )

    # --- Return for other tabs ---
    list(
      alpha_result        = reactive(result_obj()),
      alpha_tree          = reactive(tree_obj()),
      alpha_index         = reactive(input$idx),            # current selection
      alpha_last_run_idx  = reactive(last_run_idx()),       # index that actually ran
      alpha_has_run       = reactive(!is.null(result_obj())),# TRUE only after successful run
      go_back             = reactive(input$back_btn),
      go_next             = reactive(input$next_btn)
    )
  })
}

