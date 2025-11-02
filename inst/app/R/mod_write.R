# mod_write.R — flat export with simple confirmation + script export
`%||%` <- function(a, b) if (!is.null(a)) a else b

# =====================================================
# USE YOUR WRITER FUNCTIONS (keep them as-is)
# =====================================================
# write.vilma.dist()
# write.vilma.pd()
# write.vilma.beta()
# write.vilma.null()

# =====================================================
# Helpers for formatting & robust detection
# =====================================================

# Symbols: avoid quoting object names in exported calls
.sym <- function(name) structure(name, class = c("vilma_sym", "character"))

.vfmt <- function(x) {
  if (inherits(x, "vilma_sym")) return(as.character(x))      # e.g., tree -> tree
  if (is.null(x)) return("NULL")
  if (is.logical(x)) return(ifelse(x, "TRUE", "FALSE"))
  if (is.numeric(x) || is.integer(x)) return(as.character(x))
  if (is.character(x)) return(paste0("\"", x, "\""))
  if (is.factor(x)) return(paste0("\"", as.character(x), "\""))
  if (is.list(x)) {
    ch <- tryCatch(unlist(x, use.names = FALSE), error = function(e) x)
    if (is.character(ch)) return(paste0("c(", paste(paste0("\"", ch, "\""), collapse = ", "), ")"))
    return(paste0("list(", paste(vapply(x, .vfmt, character(1)), collapse = ", "), ")"))
  }
  deparse(x, width.cutoff = 80)
}

.vcall_line <- function(lhs, fn, args_named_list, comment = NULL) {
  arg_pairs <- paste(sprintf("%s = %s",
                             names(args_named_list),
                             vapply(args_named_list, .vfmt, character(1))),
                     collapse = ", ")
  line <- if (!nzchar(lhs)) sprintf("%s(%s)", fn, arg_pairs)
          else sprintf("%s <- %s(%s)", lhs, fn, arg_pairs)
  if (!is.null(comment) && nzchar(comment)) paste0(line, "  # ", comment) else line
}

# Pull a value from attributes or fields
.g <- function(obj, ...) {
  keys <- c(...)
  for (k in keys) {
    v <- attr(obj, k, exact = TRUE)
    if (!is.null(v)) return(v)
    if (is.list(obj) && !is.null(obj[[k]])) return(obj[[k]])
  }
  NULL
}

# Canonicalize calc helpers to public API names for export
.canon_fn <- function(fn) {
  fn <- gsub("^phylosor\\.calc$", "phylosor", fn)
  fn <- gsub("^unifrac\\.calc$", "unifrac", fn)
  fn <- gsub("^rao\\.beta$", "rao.beta", fn)
  fn <- gsub("^mpd\\.calc$", "mpd", fn)
  fn <- gsub("^mntd\\.calc$", "mntd", fn)
  fn <- gsub("^nri\\.calc$", "nri", fn)
  fn <- gsub("^nti\\.calc$", "nti", fn)
  fn <- gsub("^pe\\.calc$",  "pe",  fn)
  fn <- gsub("^rao\\.calc$", "rao.alpha", fn)
  ## Map old dist builder to new API
  fn <- gsub("^points\\.to\\.raster$", "points_to_raster", fn)
  fn
}

# Extract function token from a "lhs <- fn(...)" line
.fn_from_call_line <- function(call_line) {
  if (is.null(call_line)) return(NULL)
  s <- sub("^[^<]*<-\\s*", "", call_line)
  m <- regexpr("^([A-Za-z0-9_.]+)\\s*\\(", s)
  if (m <= 0) return(NULL)
  sub("\\($", "", sub("\\s*\\($", "", regmatches(s, m)))
}

# Parse a function call string into canon(fn) + args string (no deep parse; just clean)
# Accepts a full assignment like: "alpha <- mpd.calc(tree = tree, ...)"
# or just the RHS "mpd.calc(tree = tree, ...)" — we normalize either case.
.extract_rhs_call <- function(line) {
  rhs <- sub("^[^<]*<-\\s*", "", line)             # remove leading "x <- "
  rhs <- sub("^\\s*", "", rhs)
  m <- regexpr("^([A-Za-z0-9_.]+)\\s*\\(", rhs)
  if (m <= 0) return(rhs)
  fn <- regmatches(rhs, m)
  fn <- sub("\\($", "", sub("\\s*\\($", "", fn))   # strip "("
  canon <- .canon_fn(fn)
  sub(paste0("^", fn), canon, rhs)
}

# ------------- deep recursive search for vilma_code -------------
# returns character vector (the code lines) or NULL
.deep_find_vilma_code <- function(x, depth = 0, max_depth = 4, seen = NULL) {
  if (depth > max_depth) return(NULL)
  # Make a lightweight id to avoid cycles
  id <- paste0(typeof(x), ":", length(attributes(x)), ":", if (is.list(x)) length(x) else 0)
  if (!is.null(seen) && id %in% seen) return(NULL)
  seen <- c(seen, id)

  # 1) attribute on this object
  vc <- attr(x, "vilma_code", exact = TRUE)
  if (is.character(vc) && length(vc)) return(vc)

  # 2) attributes that themselves may carry vilma_code
  at <- attributes(x)
  if (!is.null(at) && length(at)) {
    for (nm in names(at)) {
      v <- at[[nm]]
      res <- .deep_find_vilma_code(v, depth + 1, max_depth, seen)
      if (!is.null(res)) return(res)
    }
  }

  # 3) list elements
  if (is.list(x)) {
    for (i in seq_along(x)) {
      res <- .deep_find_vilma_code(x[[i]], depth + 1, max_depth, seen)
      if (!is.null(res)) return(res)
    }
  }
  NULL
}

# Preferred extractor: read exact call from vilma_code attached by modules
# target_lhs is one of "alpha","beta","nulls","dist"
.extract_from_vilma_code <- function(lhs, obj, target_lhs = c("alpha","beta","nulls","dist")) {
  target_lhs <- match.arg(target_lhs)
  code <- .deep_find_vilma_code(obj)
  if (is.null(code)) return(NULL)

  # Pick a line that assigns the target_lhs if available; otherwise the last call-ish line.
  patt <- paste0("^\\s*", target_lhs, "\\s*<-")
  hits <- grep(patt, code)
  if (!length(hits)) {
    hits <- grep("\\w+\\s*\\(", code)
    if (!length(hits)) return(NULL)
  }
  line <- code[hits[length(hits)]]
  rhs  <- .extract_rhs_call(line)
  if (nzchar(lhs) && !startsWith(rhs, paste0(lhs, " <-"))) paste0(lhs, " <- ", rhs) else rhs
}

# Try to get an explicit call from other object metadata if present
.extract_from_meta <- function(lhs, obj) {
  cs <- attr(obj, "vilma_call_str", exact = TRUE)
  if (!is.null(cs) && is.character(cs) && length(cs) == 1) {
    rhs <- .extract_rhs_call(cs)
    if (nzchar(lhs) && !startsWith(rhs, paste0(lhs, " <-"))) return(paste0(lhs, " <- ", rhs))
    return(rhs)
  }
  cl <- attr(obj, "vilma_call", exact = TRUE) %||% attr(obj, "call", exact = TRUE)
  if (!is.null(cl)) {
    s <- paste(deparse(cl, width.cutoff = 80), collapse = "")
    rhs <- .extract_rhs_call(s)
    if (nzchar(lhs) && !startsWith(rhs, paste0(lhs, " <-"))) return(paste0(lhs, " <- ", rhs))
    return(rhs)
  }
  meta <- attr(obj, "vilma_meta", exact = TRUE)
  if (is.list(meta) && !is.null(meta$fn) && is.list(meta$args)) {
    fn <- .canon_fn(meta$fn)
    args <- lapply(meta$args, function(v) {
      if (is.character(v) && length(v) == 1 && grepl("^[A-Za-z.][A-Za-z0-9_.]*$", v)) .sym(v) else v
    })
    return(.vcall_line(lhs, fn, args))
  }
  NULL
}

# =====================================================
# Builders (Alpha / Beta / Dist / Nulls)
# =====================================================

# ---------- Alpha ----------
.build_alpha_call <- function(alpha_obj, tree_nm = "tree", dist_nm = "dist", lhs = "alpha") {
  from_code <- .extract_from_vilma_code(lhs, alpha_obj, target_lhs = "alpha")
  if (!is.null(from_code)) return(from_code)
  from_meta <- .extract_from_meta(lhs, alpha_obj)
  if (!is.null(from_meta)) return(from_meta)

  .vcall_line(lhs, "faith.pd", list(tree = .sym(tree_nm), dist = .sym(dist_nm)),
              comment = "TODO: supply exact alpha call (no vilma_code/vilma_meta found)")
}

# ---------- Beta ----------
.build_beta_call <- function(beta_obj, tree_nm = "tree", dist_nm = "dist", lhs = "beta") {
  from_code <- .extract_from_vilma_code(lhs, beta_obj, target_lhs = "beta")
  if (!is.null(from_code)) return(from_code)

  from_meta <- .extract_from_meta(lhs, beta_obj)
  if (!is.null(from_meta)) return(from_meta)

  .vcall_line(lhs, "unifrac",
              list(tree = .sym(tree_nm), dist = .sym(dist_nm), weighted = FALSE, normalize = TRUE),
              comment = "TODO: could not infer beta index (no vilma_code/meta)")
}

# ---------- Dist ----------
.build_dist_call <- function(dist_obj, lhs = "dist") {
  from_code <- .extract_from_vilma_code(lhs, dist_obj, target_lhs = "dist")
  if (!is.null(from_code)) return(from_code)
  from_meta <- .extract_from_meta(lhs, dist_obj)
  if (!is.null(from_meta)) return(from_meta)
  .vcall_line(lhs, "points_to_raster", list(points = .sym("points")),
              comment = "TODO: supply the same args used in app")
}

# Map null fn from alpha fn (if we need a fallback)
.map_null_from_alpha <- function(alpha_fn) {
  if (is.null(alpha_fn)) return(NULL)
  alpha_fn <- .canon_fn(alpha_fn)
  switch(alpha_fn,
    "mntd"       = "mntd.calc.null",   # matches your server naming
    "mpd"        = "mpd.calc.null",
    "faith.pd"   = "faith.pd.null",
    "pe"         = "pe.calc.null",
    "rao.alpha"  = "rao.calc.null",
    NULL
  )
}

# ---------- Nulls ----------
.build_nulls_call <- function(null_obj, alpha_nm = "alpha", tree_nm = "tree", dist_nm = "dist", lhs = "nulls", alpha_obj = NULL) {
  from_code <- .extract_from_vilma_code(lhs, null_obj, target_lhs = "nulls")
  if (!is.null(from_code)) return(from_code)
  from_meta <- .extract_from_meta(lhs, null_obj)
  if (!is.null(from_meta)) return(from_meta)

  # If neither code nor meta are present, try to infer from alpha
  alpha_call <- if (!is.null(alpha_obj)) .extract_from_vilma_code("alpha", alpha_obj, target_lhs = "alpha") else NULL
  if (is.null(alpha_call) && !is.null(alpha_obj)) alpha_call <- .extract_from_meta("alpha", alpha_obj)
  alpha_fn <- .fn_from_call_line(alpha_call)
  null_fn  <- .map_null_from_alpha(alpha_fn) %||% (.g(null_obj, "null_fn") %||% "faith.pd.null")

  args <- list(
    pd               = .sym(alpha_nm),
    tree             = .sym(tree_nm),
    dist             = .sym(dist_nm),
    iterations       = .g(null_obj, "iterations") %||% 999,
    method           = .g(null_obj, "method") %||% "cell",
    sampling         = .g(null_obj, "sampling") %||% "taxa.label",
    n.directions     = .g(null_obj, "n.directions") %||% "rook",
    regional.weight  = .g(null_obj, "regional.weight") %||% "uniform"
  )
  .vcall_line(lhs, null_fn, args)
}

# =====================================================
# UI
# =====================================================
mod_write_export_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h3("Write / Export"),
    shiny::fluidRow(
      shiny::column(
        6,
        shiny::div(class = "soft-card",
          shiny::h4("Destination"),
          shiny::textInput(
            ns("out_dir"),
            "Output folder (will be created if missing)",
            value = normalizePath(file.path(getwd(), "vilma_export"), winslash = "/", mustWork = FALSE)
          ),
          shiny::actionButton(ns("make_dir"), "Create / Check folder", class = "btn btn-outline-secondary"),
          shiny::br(), shiny::br(),
          shiny::h4("Select outputs"),
          shiny::uiOutput(ns("which_objs_ui")),
          shiny::helpText("Only analyses already run appear here."),
          shiny::radioButtons(
            ns("raster_fmt"),
            "Raster format",
            choices = c("tif", "grd", "img"),
            selected = "tif", inline = TRUE
          ),
          shiny::checkboxInput(ns("overwrite"), "Overwrite existing files", value = TRUE)
        )
      ),
      shiny::column(
        6,
        shiny::div(class = "soft-card",
          shiny::h4("Run"),
          shiny::actionButton(ns("write_now"), "Write files to folder", class = "btn btn-primary"),
          shiny::br(), shiny::br(),
          shiny::actionButton(ns("export_script"), "Export session code (.R)", class = "btn btn-secondary"),
          shiny::br(), shiny::br(),
          # --- Close app: synchronous JS so the tab actually closes ---
          shiny::actionButton(
            ns("close_app"),
            "Close app",
            class = "btn btn-outline-secondary",
            onclick = sprintf(
              "Shiny.setInputValue('%s', Date.now());
               try { window.close(); } catch(e) {}
               setTimeout(function(){
                 try { window.open('', '_self'); window.close(); } catch(e) {}
               }, 120);",
              ns("close_direct")
            )
          ),
          shiny::tags$hr(),
          shiny::verbatimTextOutput(ns("diag"), placeholder = TRUE),
          shiny::tags$hr(),
          shiny::verbatimTextOutput(ns("status"), placeholder = TRUE)
        )
      )
    )
  )
}

# =====================================================
# SERVER
# =====================================================
mod_write_export_server <- function(id, results) {
  shiny::moduleServer(id, function(input, output, session) {

    # Ensure R process stops when window closes OR button clicked
    session$onSessionEnded(function() {
      try(shiny::stopApp(), silent = TRUE)
      try(invisible(gc()), silent = TRUE)
    })
    # legacy button id still stops app even if tab refuses to close
    shiny::observeEvent(input$close_app, {
      try(shiny::stopApp(), silent = TRUE)
      try(invisible(gc()), silent = TRUE)
    })
    # NEW: user-initiated synchronous click -> close tab + stop app
    shiny::observeEvent(input$close_direct, {
      try(shiny::stopApp(), silent = TRUE)
      try(invisible(gc()), silent = TRUE)
    }, ignoreInit = TRUE)

    # Diagnostics
    output$diag <- shiny::renderText({
      r <- tryCatch(results(), error = function(e) list())
      lines <- c("Detected objects:")
      add <- function(label, obj) {
        if (is.null(obj)) lines <<- c(lines, sprintf("  - %s: <none>", label))
        else lines <<- c(lines, sprintf("  - %s: %s", label, paste(class(obj), collapse = "/")))
      }
      add("dist",  r$dist); add("alpha", r$alpha); add("beta", r$beta); add("null", r$null)
      paste(lines, collapse = "\n")
    })

    # Output selection
    output$which_objs_ui <- shiny::renderUI({
      r <- tryCatch(results(), error = function(e) list())
      opts <- names(Filter(Negate(is.null), r[c("dist","alpha","beta","null")]))
      shiny::checkboxGroupInput(session$ns("which_objs"), NULL,
        choices = opts, selected = opts)
    })

    # Folder check
    shiny::observeEvent(input$make_dir, {
      out_dir <- input$out_dir %||% ""
      if (!nzchar(out_dir)) return(output$status <- shiny::renderText("Provide folder path."))
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      output$status <- shiny::renderText(paste("Folder ready:", normalizePath(out_dir, winslash = "/")))
    })

    # Write everything flat
    shiny::observeEvent(input$write_now, {
      r <- tryCatch(results(), error = function(e) list())
      sel <- input$which_objs %||% character()
      out_dir <- input$out_dir %||% ""
      fmt <- input$raster_fmt %||% "tif"
      ow  <- isTRUE(input$overwrite)
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      if ("dist"  %in% sel && !is.null(r$dist))  write.vilma.dist(r$dist,  file = file.path(out_dir, "dist"),  raster.format = fmt, overwrite = ow)
      if ("alpha" %in% sel && !is.null(r$alpha)) write.vilma.pd(r$alpha,  file = file.path(out_dir, "alpha"), raster.format = fmt, overwrite = ow)
      if ("beta"  %in% sel && !is.null(r$beta))  write.vilma.beta(r$beta, file = file.path(out_dir, "beta"),  raster.format = fmt, overwrite = ow)
      if ("null"  %in% sel && !is.null(r$null))  write.vilma.null(r$null, file = file.path(out_dir, "null"),  raster.format = fmt, overwrite = ow)

      shiny::showNotification("All files exported successfully!", type = "message", duration = 6)
      output$status <- shiny::renderText(
        paste0("✅ Everything exported successfully to:\n", normalizePath(out_dir, winslash = "/"))
      )
    })

    # Export Script — explicit named args, now reading vilma_code deeply
    shiny::observeEvent(input$export_script, {
      r <- tryCatch(results(), error = function(e) list())
      out_dir <- input$out_dir %||% getwd()
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      dist_nm  <- "dist"
      tree_nm  <- "tree"
      alpha_nm <- "alpha"
      beta_nm  <- "beta"
      null_nm  <- "nulls"

      ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      lines <- c(
        "# ==========================================================",
        "#  Vilma Session Export (Reproducible Script)",
        paste0("#  Date: ", ts),
        "# ==========================================================",
        "",
        "library(terra); library(ape); # library(vilma) # if installed as a package",
        "",
        "# ---- User I/O placeholders (adjust to your files) ----",
        "# points <- read.csv('points.csv')",
        "# tree   <- read.tree('tree.tre')",
        ""
      )

      if (!is.null(r$dist)) {
        lines <- c(lines, "# ---- Dist ----",
                   .build_dist_call(r$dist, lhs = dist_nm),
                   "")
      }
      if (!is.null(r$alpha)) {
        lines <- c(lines, "# ---- Alpha (explicit arguments) ----",
                   .build_alpha_call(r$alpha, tree_nm = tree_nm, dist_nm = dist_nm, lhs = alpha_nm),
                   "")
      }
      if (!is.null(r$beta)) {
        lines <- c(lines, "# ---- Beta (explicit arguments) ----",
                   .build_beta_call(r$beta, tree_nm = tree_nm, dist_nm = dist_nm, lhs = beta_nm),
                   "")
      }
      if (!is.null(r$null)) {
        lines <- c(lines, "# ---- Null models ----",
                   .build_nulls_call(r$null, alpha_nm = alpha_nm, tree_nm = tree_nm, dist_nm = dist_nm, lhs = null_nm, alpha_obj = r$alpha),
                   "")
      }

      lines <- c(lines,
        "# ---- Export results ----",
        if (!is.null(r$dist))  sprintf("write.vilma.dist(%s,  file='dist',  raster.format='tif', overwrite=TRUE)", dist_nm) else NULL,
        if (!is.null(r$alpha)) sprintf("write.vilma.pd(%s,   file='alpha', raster.format='tif', overwrite=TRUE)", alpha_nm) else NULL,
        if (!is.null(r$beta))  sprintf("write.vilma.beta(%s, file='beta',  raster.format='tif', overwrite=TRUE)", beta_nm) else NULL,
        if (!is.null(r$null))  sprintf("write.vilma.null(%s, file='null',  raster.format='tif', overwrite=TRUE)", null_nm) else NULL,
        "",
        "sessionInfo()"
      )

      script_path <- file.path(out_dir, "vilma_session_export.R")
      writeLines(lines, script_path)
      shiny::showNotification(paste("Script exported to", normalizePath(script_path, winslash = "/")),
                              type = "message", duration = 6)
      output$status <- shiny::renderText(
        paste0("✅ Script exported successfully to:\n", normalizePath(script_path, winslash = "/"))
      )
    })
  })
}

