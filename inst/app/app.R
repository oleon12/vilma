# ========================= app.R (complete) =========================
library(shiny)
library(bslib)
library(readr)
library(leaflet)
library(leaflet.extras)
library(ape)
library(terra)
library(later)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- Load modules (adapt to your project layout as needed) ----
safely_load_modules <- function() {
  app_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)

  # If running from the *installed* package's inst/app, do not source or hot-reload anything.
  installed_app_dir <- ""
  if (requireNamespace("vilma", quietly = TRUE)) {
    installed_app_dir <- tryCatch(
      normalizePath(system.file("app", package = "vilma"), winslash = "/", mustWork = TRUE),
      error = function(e) ""
    )
  }
  if (nzchar(installed_app_dir) && identical(app_dir, installed_app_dir)) {
    return(invisible(TRUE))
  }

  # Otherwise we're in a dev tree. Do NOT call pkgload::load_all() here.
  guess_root <- normalizePath(file.path(app_dir, "..", ".."), winslash = "/", mustWork = FALSE)

  cand_dirs <- unique(c(
    file.path(app_dir, "R"),
    file.path(app_dir, "inst", "app", "R"),
    if (dir.exists(guess_root)) file.path(guess_root, "R")               else character(),
    if (dir.exists(guess_root)) file.path(guess_root, "inst", "app", "R") else character()
  ))
  cand_dirs <- cand_dirs[dir.exists(cand_dirs)]

  for (d in cand_dirs) {
    rfiles <- list.files(d, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
    rfiles <- rfiles[basename(rfiles) != "app.R"]
    for (f in rfiles) try(source(f, local = FALSE, chdir = TRUE), silent = TRUE)
  }

  required_funs <- c(
    "mod_points_to_raster_ui","mod_points_to_raster_server",
    "mod_alpha_ui","mod_alpha_server",
    "mod_nulls_ui","mod_nulls_server",
    "mod_beta_ui","mod_beta_server",
    "mod_write_export_ui","mod_write_export_server"
  )
  missing <- required_funs[!vapply(required_funs, exists, logical(1), mode = "function", inherits = TRUE)]
  if (length(missing)) stop(paste("Missing module functions:", paste(missing, collapse=", ")), call. = FALSE)

  invisible(TRUE)
}
safely_load_modules()

# =============================== UI ================================
ui <- page_fluid(
  theme = bs_theme(version = 5),
  tags$head(
    tags$style(HTML("
      html, body { height:100%; width:100%; margin:0; padding:0; overflow:hidden; }
      .container-fluid { height:100%; width:100%; }
      .vilma-header {
        background:#FFA4A4; height:85px; display:flex; align-items:center; justify-content:space-between;
        padding:0 20px; box-shadow:0 2px 6px rgba(0,0,0,0.08);
      }
      .vilma-title { font-weight:800; font-size:28px; color:white; letter-spacing:0.5px; font-style:italic; }
      .vilma-logo { height:80px; width:auto; }
      .vilma-page { margin-top:16px; height:calc(100% - 85px); overflow-y:auto; padding:20px; }
      .btn-primary { background-color:#FFA4A4; border-color:#FFA4A4; }
      .btn-primary:hover, .btn-primary:focus { background-color:#ff2c95; border-color:#ff2c95; }
      .btn-secondary, .btn-outline-secondary { background-color:#FFA4A4; border-color:#FFA4A4; color:white; }
      .btn-secondary:hover, .btn-outline-secondary:hover { background-color:#ff2c95; border-color:#ff2c95; color:white; }
      .soft-card { background:#f9f9fb; border-radius:14px; padding:16px; box-shadow:0 1px 5px rgba(0,0,0,0.06); }

      /* How-to styling */
      .howto h2, .howto h3, .howto h4 { color:#444; }
      .howto p, .howto li { line-height:1.6; }
      .howto table { width:100%; border-collapse:collapse; margin: 8px 0 16px; }
      .howto th, .howto td { border: 1px solid #eee; padding: 8px 10px; vertical-align: top; }
      .howto th { background: #ffe6eb; font-weight:700; }
      .divider { border-top: 2px solid #FFD1DC; margin: 28px 0; }

      .leaflet, .leaflet-container { min-height: 520px !important; }
      .tab-pane .leaflet-container { height: 520px !important; }
    ")),
    tags$script(HTML("
      document.addEventListener('shown.bs.tab', function(e) {
        setTimeout(function(){ window.dispatchEvent(new Event('resize')); }, 100);
      }, true);
    ")),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('vilma-close', function() {
        try { window.close(); } catch(e) {}
        setTimeout(function(){
          try { window.open('', '_self'); window.close(); } catch(e) {}
        }, 120);
      });
    "))
  ),

  # ============================ HEADER =============================
  div(class = "vilma-header",
      div(class = "vilma-title", "Vilma: Spatial Phylogenetic Diversity Analyses"),
      tags$img(src = "logo.png", class = "vilma-logo", alt = "Vilma logo")),

  # ============================== PAGES =============================
  div(class = "vilma-page",
      tabsetPanel(
        id = "tabs",

        # ----------------------- Welcome -----------------------
        tabPanel("Welcome",
          fluidRow(
            column(12,
              div(style="text-align:center; margin-top:60px; margin-bottom:30px;",
                h2("Welcome to Vilma", style="font-weight:800; color:#444;"),
                p(br(),
                  "This app was designed to help users without programming experience conduct ",
                  strong("Spatial Phylogenetic Diversity Analyses."), " ",
                  "Although the app makes it easier to perform the analyses, you should be aware of the importance of truly understanding the logic, mathematics, and parameters behind each index.",
                  br(), br(),
                  "Every analysis requires high-quality data and a clearly defined question. ",
                  "Remember: ", strong("GIGO — Garbage In, Garbage Out."), " Use the app carefully.",
                  style="max-width:700px; margin:0 auto; text-align:justify; line-height:1.6;"
                ),
                br(),
                actionButton("continue_button","Continue",
                  class="btn btn-lg",
                  style="background-color:#ff9999; border:none; color:white; font-weight:bold; padding:10px 25px;")
              )
            )
          )
        ),

        # ------------------------ About ------------------------
        tabPanel("About",
          fluidRow(
            column(
              10, offset = 1,
              div(class="soft-card",
                h2(br(), "About Vilma", br(), style="font-weight:800; color:#444;"),
                p(br(),"Vilma is an R-based platform designed to quantify and visualize ",
                  strong("spatial phylogenetic diversity"),
                  " across geographic landscapes. It integrates multiple analytical frameworks ",
                  "to examine how evolutionary history and ecological processes shape biodiversity through space and time.",
                  style="text-align:justify; line-height:1.6;"),
                div(class="divider"),
                h4("Vilma provides a complete and reproducible workflow to:",br()),
                tags$ul(
                  tags$li(HTML("Compute <b>alpha diversity</b> indices such as Faith’s PD, MPD, MNTD, PE, Rao’s Q, NRI, and NTI.")),
                  tags$li(HTML("Estimate <b>beta diversity</b> metrics including PhyloSor, UniFrac, Phylobeta (weighted and unweighted), &beta;MPD, and &beta;MNTD.")),
                  tags$li(HTML("Implement and test <b>null models</b> at both regional and local (cell-based) levels to assess whether observed spatial patterns deviate from random expectations.")),
                  tags$li("Visualize and export spatially explicit outputs as tables and raster maps for further ecological and evolutionary analyses.")
                ),
                div(class="divider"),
                p(HTML(
                  "Developed as part of an ongoing research project at the ",
                  "<b>Department of Mammalogy, American Museum of Natural History</b>, ",
                  "and the <b>Department of Earth and Environmental Sciences, Rutgers University–Newark</b>, ",
                  "Vilma provides a flexible framework to investigate phylogenetic structure, endemism, and connectivity in ecological and evolutionary studies."
                )),
                br(),
                em("Created in loving memory of ", strong("Vilma Alvarado, "),
                   "whose kindness and love continue to inspire this work.", br()),
                br(), br(),
                div(style="text-align:center;",
                  actionButton("continue_about","Continue",
                    class="btn btn-lg",
                    style="background-color:#ff9999; border:none; color:white; font-weight:bold; padding:10px 25px;"))
              )
            )
          )
        ),

        # --------------------- How to Use ----------------------
        tabPanel("How to Use",
          fluidRow(
            column(
              10, offset = 1,
              div(class="soft-card howto",
                h2("How to Use Vilma", style="font-weight:800;"),
                p("New to Phylogenetic Diversity Indices? Please read this easy-to-follow manual. Then, scroll down to the Continue button."),
                div(class="divider"),

                # ---------- General Overview ----------
                h3(strong("General Overview")),
                p(br(), "Vilma is designed as a ", strong("pipeline"), ", where analyses follow a sequential order. ",
                  "You begin by creating a raster object (", code("vilma.dist"), "), and then continue with ",
                  strong("\u03B1-diversity"), ", ", strong("null models"), ", and ", strong("\u03B2-diversity"),
                  ". Finally, you can export all tables and raster layers generated by the different functions."),
                p("You may choose to run only specific analyses (for example, only \u03B1-diversity or \u03B2-diversity), ",
                  "but remember that the ", strong("phylogenetic tree is uploaded only in the \u03B1-diversity page"),
                  ". If you skip this step, later functions will not recognize the tree."),
                div(class="divider"),

                # ---------- Interactive Maps ----------
                h3("Interactive Maps"),
                p(br(), "Each analysis generates an interactive map with multiple layers. ",
                  "Use the ", strong("Layer Control Button"), " ",
                  img(src="Layer_Icon.png", height="22px", style="margin:0 4px; vertical-align:middle;"),
                  "to toggle among the layers."),
                p("Hovering over cells displays their corresponding values. Depending on the analysis, layers may represent:"),
                tags$ul(
                  tags$li("Species richness or abundance"),
                  tags$li("Diversity index values"),
                  tags$li("Standardized effect sizes (SES)"),
                  tags$li("P-values or significance maps")
                ),
                p("These layers help interpret how evolutionary history and ecological processes vary across space."),
                div(class="divider"),

                # ---------- 1) Raster Creation ----------
                h3(br(), strong("1) Raster Creation (Points \u2192 Raster)"), br(), br()),
                p("Input format (required, in this order):", br(), br()),
                tags$table(
                  tags$thead(tags$tr(tags$th("Column"), tags$th("Description"), tags$th("Example"))),
                  tags$tbody(
                    tags$tr(tags$td("Species"),   tags$td("Species name"), tags$td(em("Artibeus_jamaicensis"))),
                    tags$tr(tags$td("Longitude"), tags$td("X (decimal degrees)"), tags$td("\u221274.125")),
                    tags$tr(tags$td("Latitude"),  tags$td("Y (decimal degrees)"), tags$td("4.678"))
                  )
                ),
                p(em("Accepted files: .csv, .txt, .rds (one record per row).")),
                tags$hr(),

                h4("Resolution Units in Vilma", br(),br()),
                p("Vilma uses a resolution defined in ", strong("degrees"),
                  ", with the default projection ", strong("WGS84 (EPSG:4326)"), ". ",
                  "If you want your raster resolution to match other datasets that use ",
                  strong("arc-seconds"), " (e.g., 30 arc-seconds \u2248 1 km), follow the conversion table below. ",
                  "Be mindful that this conversion works only with the projection ", strong("WGS84"),
                  ", and the actual distance per degree varies slightly with latitude.", br(), br()),
                tags$table(
                  tags$thead(
                    tags$tr(
                      tags$th("Resolution (arc-seconds)"),
                      tags$th("Approximate meters"),
                      tags$th("Degrees (use in Vilma)")
                    )
                  ),
                  tags$tbody(
                    tags$tr(tags$td("300"),  tags$td("\u2248 10,000 m"), tags$td("0.083333")),
                    tags$tr(tags$td("150"),  tags$td("\u2248 5,000 m"),  tags$td("0.041667")),
                    tags$tr(tags$td("30"),   tags$td("\u2248 1,000 m"),  tags$td("0.008333")),
                    tags$tr(tags$td("10"),   tags$td("\u2248 300 m"),    tags$td("0.002778")),
                    tags$tr(tags$td("5"),    tags$td("\u2248 150 m"),    tags$td("0.001389")),
                    tags$tr(tags$td("2.5"),  tags$td("\u2248 77 m"),     tags$td("0.000694")),
                    tags$tr(tags$td("1"),    tags$td("\u2248 31 m"),     tags$td("0.000278")),
                    tags$tr(tags$td("0.5"),  tags$td("\u2248 15 m"),     tags$td("0.000139")),
                    tags$tr(tags$td("0.25"), tags$td("\u2248 7.7 m"),    tags$td("0.000069"))
                  )
                ),
                p(strong("Outputs (object: "), code("vilma.dist"), strong(")")),
                tags$ul(
                  tags$li("Table with each species and its respective cell ID"),
                  tags$li("Richness and abundance rasters")
                ),
                div(class="divider"),

                # ---------- 2) Alpha Diversity ----------
                h3(strong("2) \u03B1-Diversity (within-cell)"), br(), br()),
                h4("Indices, Parameters, and What They Represent", br(), br()),
                tags$ul(
                  tags$li(HTML("<b>Faith PD</b><br/><br/>Measures the total evolutionary history represented by the species in a cell. It sums the branch lengths of all lineages that occur in that area. Cells with longer total branch lengths contain more evolutionary diversity.<br/><i>Parameter:</i> <code>faith_method</code> (root, node, exclude; default = root). <br/><i>Sensitive to:</i> total branch length in the tree, richness (more species usually increase PD), long deep branches, tree scaling/units, and uneven sampling across clades.<br/><br/>")),
                  tags$li(HTML("<b>MPD (Mean Pairwise Distance)</b><br/><br/>Measures how evolutionarily distant species are from each other on average. Large MPD values indicate phylogenetic overdispersion; low values suggest close relatives dominate.<br/><i>Parameters:</i> <code>mpd_method</code> (root, node, exclude; default = root), <code>mpd_abundance</code> (default = FALSE).<br/><i>Sensitive to:</i> deep/basal structure of the tree (older splits), presence of distantly related lineages, abundance weighting (if used), singletons handling, and tree topology uncertainty.<br/><br/>")),
                  tags$li(HTML("<b>MNTD (Mean Nearest-Taxon Distance)</b><br/><br/>Focuses on the closest relatives of each species in a cell and highlights fine-scale clustering of closely related species.<br/><i>Parameters:</i> <code>mntd_method</code> (root, node, exclude; default = exclude), <code>mntd_abundance</code> (default = FALSE).<br/><i>Sensitive to:</i> tip-level structure (recent splits), local clustering/overdispersion, occurrence of sister taxa, minimum distances in unbalanced clades, and how singletons are treated.<br/><br/>")),
                  tags$li(HTML("<b>PE (Phylogenetic Endemism)</b><br/><br/>Identifies areas where unique evolutionary lineages have small geographic ranges. High PE indicates restricted and distinctive evolutionary history.<br/><i>Parameters:</i> <code>pe_RPE</code> (Relative PE), <code>pe_faith_method</code>.<br/><i>Sensitive to:</i> branch lengths weighted by geographic range size (short-ranged lineages), spatial grain/extent, range-map bias, and tree scaling.<br/><br/>")),
                  tags$li(HTML("<b>Rao’s Q (\u03B1)</b><br/><br/>Reflects the average phylogenetic difference between individuals or species in a community, combining species abundances and phylogenetic distances.<br/><i>Parameters:</i> <code>rao_abundance</code> (FALSE), <code>rao_scale01</code> (TRUE).<br/><i>Sensitive to:</i> abundance structure (if enabled), distribution of pairwise phylogenetic distances, dominance/evenness, and the choice to rescale to 0–1.<br/><br/>")),
                  tags$li(HTML("<b>NRI (Net Relatedness Index)</b><br/><br/>Shows whether a community’s species are more or less related than expected by chance. Negative values suggest overdispersion; positive values indicate clustering. See the null model section for parameter details.<br/><i>Parameters:</i> <code>nri_mpd_method</code>, <code>nri_abundance</code>, <code>nri_iterations</code>, <code>nri_sampling</code>, <code>nri_n_directions</code>, <code>nri_regional_weight</code>.<br/><i>Sensitive to:</i> MPD behavior (deep structure), the <em>null model</em> and constraints (sampling scheme, spatial constraints, pool definition), richness, and spatial grain.<br/><br/>")),
                  tags$li(HTML("<b>NNI (Nearest-Neighbor Index)</b><br/><br/>Similar to NRI but based on nearest relatives rather than all pairwise distances; useful for detecting fine-scale clustering of close relatives. See the null model section for parameter details.<br/><i>Parameters:</i> <code>nti_mntd_method</code>, <code>nti_abundance</code>, <code>nti_iterations</code>, <code>nti_sampling</code>, <code>nti_n_directions</code>, <code>nti_regional_weight</code>.<br/><i>Sensitive to:</i> MNTD behavior (tip structure), <em>null model</em> choices and spatial constraints, presence of sister species, richness, and cell size.<br/><br/>"))
                ),

                h4("Method Options (root, node, exclude)", br(), br()),
                tags$ul(
                  tags$li(HTML("<b>exclude</b> – skips any cell that has only one species.")),
                  tags$li(HTML("<b>root</b> – if a cell has a single species, its value is calculated using the entire evolutionary path from that species to the root of the tree, giving weight to its full evolutionary history.")),
                  tags$li(HTML("<b>node</b> – if a cell has a single species, its value is based only on the last branch that connects it to its closest ancestor, focusing on recent evolutionary differences.<br/><br/>"))
                ),

                p(strong("Output — "), code("vilma.pd")),
                tags$ul(
                  tags$li("Cell table and raster layer with index values"),
                  tags$li("Printable summary statistics")
                ),
                div(class="divider"),

                # ---------- 3) Null Models ----------
                h3(strong("3) Null Models (for \u03B1-Diversity)"), br(), br()),
                p(strong("Purpose:"), " to test whether your observed patterns are different from random expectations. ",
                  "The index is locked to the \u03B1-diversity you ran (Faith PD, MPD, MNTD, PE, RaoQ).", br()),

                h4("Global vs Cell", br(), br()),
                tags$ul(
                  tags$li(HTML("<b>Global method</b> tests whether the overall region differs from randomness. Useful for detecting large-scale phylogenetic structure. <i>Output:</i> one observed value, null distribution, SES, and p-value.<br/>")),
                  tags$li(HTML("<b>Cell method</b> tests each cell individually. Highlights hotspots (significant clustering or overdispersion). <i>Output:</i> SES and p-value rasters.<br/><br/>"))
                ),

                h4("Sampling Strategies — How and Why", br(), br()),
                tags$ul(
                  tags$li(HTML("<b>taxa.label</b><br/><br/>Randomizes species labels on the phylogenetic tree.<br/><i>Goal:</i> Removes the link between traits and evolutionary history to test whether the pattern depends on phylogenetic relationships.<br/><br/>")),
                  tags$li(HTML("<b>range</b><br/><br/>Swaps occurrences among species while keeping species’ range sizes and cell richness constant.<br/><i>Goal:</i> Tests whether the pattern is driven by how widespread species are, rather than by phylogeny.<br/><br/>")),
                  tags$li(HTML("<b>neighbor</b><br/><br/>Exchanges species only between neighboring cells (based on rook, bishop, or queen adjacency).<br/><i>Goal:</i> Keeps spatial structure realistic, testing whether local species replacement shapes the pattern.<br/><br/>")),
                  tags$li(HTML("<b>regional</b><br/><br/>Randomizes within a defined regional species pool using different weighting schemes.<br/><i>Goal:</i> Tests whether community structure is driven by regional biogeography instead of random species assembly.<br/><br/>"))
                ),

                # Parameter details table (alpha/null models)
                h5(strong("Parameter details"), br(), br()),
                tags$table(
                  tags$thead(
                    tags$tr(
                      tags$th("Parameter"),
                      tags$th("Applies when"),
                      tags$th("Options (default)"),
                      tags$th("Purpose / Notes")
                    )
                  ),
                  tags$tbody(
                    tags$tr(
                      tags$td(code("iterations")),
                      tags$td("All null model runs"),
                      tags$td("Positive integer (default = 999)"),
                      tags$td("How many randomizations to build the null distribution. Higher values give more stable SES and p-values but take longer.")
                    ),
                    tags$tr(
                      tags$td(code("method")),
                      tags$td("All null model runs"),
                      tags$td("cell, global (choose one)"),
                      tags$td(HTML("<b>cell</b>: test each cell independently (hotspots map). <b>global</b>: test a single regional value."))
                    ),
                    tags$tr(
                      tags$td(code("sampling")),
                      tags$td("All null model runs"),
                      tags$td("taxa.label, range, neighbor, regional"),
                      tags$td(HTML("<b>taxa.label</b>: shuffle species labels on the tree; breaks trait–phylogeny link. <br/><b>range</b>: swap occurrences while keeping species range size and cell richness the same. <br/><b>neighbor</b>: swap only among adjacent cells (spatially realistic). <br/><b>regional</b>: draw from a defined regional species pool with weights."))
                    ),
                    tags$tr(
                      tags$td(code("n_directions")),
                      tags$td(HTML("When <code>sampling</code> = <b>neighbor</b>")),
                      tags$td("rook, bishop, queen (default = rook)"),
                      tags$td(HTML("<b>rook</b>: 4-neighbors (N,S,E,W). <b>bishop</b>: diagonals only. <b>queen</b>: 8-neighbors (rook + bishop). Controls which cells can exchange species."))
                    ),
                    tags$tr(
                      tags$td(code("regional_weight")),
                      tags$td(HTML("Only when <code>sampling</code> = <b>regional</b>")),
                      tags$td("uniform, frequency, range (default = uniform)"),
                      tags$td(HTML("<b>uniform</b>: all species equally likely. <br/><b>frequency</b>: species drawn in proportion to how often they occur in the region (common species more likely). <br/><b>range</b>: species drawn in proportion to their range size (widespread species more likely)."))
                    ),
                    tags$tr(
                      tags$td(code("rao_abundance")),
                      tags$td("Index is RaoQ (alpha Rao only)"),
                      tags$td("TRUE / FALSE (default = FALSE)"),
                      tags$td("Include abundance in Rao’s Q calculations. Ignored for other indices.")
                    ),
                    tags$tr(
                      tags$td(code("rao_scale01")),
                      tags$td("Index is RaoQ (alpha Rao only)"),
                      tags$td("TRUE / FALSE (default = TRUE)"),
                      tags$td("Rescale Rao’s Q to 0–1 for easier comparison. Ignored for other indices.")
                    )
                  )
                ),
                p(br(),em("Notes:"), " NRI uses MPD nulls; NNI uses MNTD nulls. Choose the corresponding method and sampling accordingly.", br(), br()),

                p(strong("Output — "), code("vilma.null"), br()),
                tags$ul(
                  tags$li("Global: regional SES and p-value"),
                  tags$li("Cell: SES and p-value rasters with index map for reference")
                ),
                div(class="divider"),

                # ---------- 4) Beta Diversity ----------
                h3(strong("4) \u03B2-Diversity (between-cell)"), br(), br()),
                h4("Indices, Parameters, and What They Represent", br(), br()),
                tags$ul(
                  tags$li(HTML("<b>Phylo Beta</b><br/><br/>Compares how evolutionary history is shared or replaced between cells; helps identify turnover of lineages across the landscape.<br/><i>Parameters:</i> <code>pb_method</code>, <code>pb_normalize</code>.<br/><i>Sensitive to:</i> partitioning of shared vs. unique branch lengths, spatial turnover vs. nestedness, normalization choices, and variation in richness among cells.<br/><br/>")),
                  tags$li(HTML("<b>Phylosor</b><br/><br/>Measures phylogenetic similarity between communities, combining richness and shared branch lengths. High similarity means communities share much of their evolutionary history.<br/><i>Parameters:</i> <code>ps_method</code>, <code>ps_singleton_overlap</code>, <code>ps_abundance</code>, <code>ps_normalize</code>.<br/><i>Sensitive to:</i> shared deep branches, richness differences (can inflate similarity), overlap of rare/singleton lineages (controlled by options), and abundance weighting (if enabled).<br/><br/>")),
                  tags$li(HTML("<b>Rao \u03B2</b><br/><br/>Evaluates the average evolutionary difference between individuals or species from different cells; useful for understanding how distinct communities are in their evolutionary composition.<br/><i>Parameters:</i> <code>rb_abundance</code>, <code>rb_scale01</code>.<br/><i>Sensitive to:</i> abundance/evenness differences, the distribution of between-site phylogenetic distances, and rescaling choices (0–1).<br/><br/>")),
                  tags$li(HTML("<b>UniFrac</b><br/><br/>Measures how much unique evolutionary history is found in one community versus another; helps detect phylogenetic turnover across space.<br/><i>Parameters:</i> <code>uf_method</code> (unweighted/weighted).<br/><i>Sensitive to:</i> unique vs. shared branch lengths (deep vs. tip contributions), detection of rare lineages, and abundance weighting in the weighted version.<br/><br/>")),
                  tags$li(HTML("<b>\u03B2MPD</b><br/><br/>Captures the average evolutionary distance between species in different cells; shows whether nearby sites are composed of close relatives or distant ones.<br/><i>Parameters:</i> <code>bmpd_method</code>, <code>bmpd_abundance</code>, <code>bmpd_exclude_conspecific</code>, <code>bmpd_normalize</code>, <code>bmpd_scale01</code>.<br/><i>Sensitive to:</i> deep/basal structure across communities, richness imbalance (can bias distances), abundance weighting, and handling of conspecific matches.<br/><br/>")),
                  tags$li(HTML("<b>\u03B2MNTD</b><br/><br/>Focuses on the closest relatives between communities, emphasizing fine-scale turnover of species; useful for identifying gradual vs. abrupt lineage replacement.<br/><i>Parameters:</i> <code>bmntd_method</code>, <code>bmntd_abundance</code>, <code>bmntd_exclude_conspecific</code>, <code>bmntd_normalize</code>, <code>bmntd_scale01</code>.<br/><i>Sensitive to:</i> tip-level replacements near community boundaries, occurrence of sister taxa across sites, abundance weighting, and treatment of conspecifics.<br/><br/>"))
                ),

                # ---- NEW: Detailed Parameter Table for Beta ----
                h5(strong("Parameter details — Beta"), br(), br()),
                tags$table(
                  tags$thead(
                    tags$tr(
                      tags$th("Parameter"),
                      tags$th("Applies when"),
                      tags$th("Options (default)"),
                      tags$th("Purpose / Notes")
                    )
                  ),
                  tags$tbody(
                    # ---- Phylo Beta ----
                    tags$tr(
                      tags$td(code("pb_method")),
                      tags$td("Phylo Beta"),
                      tags$td("unweighted, weighted (default = unweighted)"),
                      tags$td(HTML("How shared vs. unique branch lengths are partitioned between cells. <b>weighted</b> can incorporate abundance if supported; <b>unweighted</b> uses presence/absence."))
                    ),
                    tags$tr(
                      tags$td(code("pb_normalize")),
                      tags$td("Phylo Beta"),
                      tags$td("TRUE / FALSE (default = TRUE)"),
                      tags$td("Rescale metric to a common range (e.g., 0–1) to improve comparability across regions.")
                    ),

                    # ---- PhyloSor ----
                    tags$tr(
                      tags$td(code("ps_method")),
                      tags$td("Phylosor"),
                      tags$td("root, node, exclude (default = root)"),
                      tags$td("How to treat singletons or shallow branches when computing shared branch lengths.")
                    ),
                    tags$tr(
                      tags$td(code("ps_singleton_overlap")),
                      tags$td("Phylosor"),
                      tags$td("TRUE / FALSE (default = FALSE)"),
                      tags$td("If TRUE, allows overlap via singleton handling when communities share very few lineages.")
                    ),
                    tags$tr(
                      tags$td(code("ps_abundance")),
                      tags$td("Phylosor"),
                      tags$td("TRUE / FALSE (default = FALSE)"),
                      tags$td("If TRUE, weights shared branch lengths by species abundances.")
                    ),
                    tags$tr(
                      tags$td(code("ps_normalize")),
                      tags$td("Phylosor"),
                      tags$td("TRUE / FALSE (default = TRUE)"),
                      tags$td("Normalize similarity (or distance) to 0–1 to compare across richness gradients.")
                    ),

                    # ---- UniFrac ----
                    tags$tr(
                      tags$td(code("uf_method")),
                      tags$td("UniFrac"),
                      tags$td("unweighted, weighted (default = unweighted)"),
                      tags$td(HTML("<b>unweighted</b>: presence/absence, sensitive to rare lineages; <b>weighted</b>: incorporates abundance, down-weights rare branches if they are low-abundance."))
                    ),

                    # ---- βMPD ----
                    tags$tr(
                      tags$td(code("bmpd_method")),
                      tags$td("\u03B2MPD"),
                      tags$td("exclude, root, node (default = exclude)"),
                      tags$td("Controls how within-cell distances are treated when normalizing pairwise distances across cells.")
                    ),
                    tags$tr(
                      tags$td(code("bmpd_abundance")),
                      tags$td("\u03B2MPD"),
                      tags$td("TRUE / FALSE (default = FALSE)"),
                      tags$td("If TRUE, uses abundance-weighted pairwise distances between communities.")
                    ),
                    tags$tr(
                      tags$td(code("bmpd_exclude_conspecific")),
                      tags$td("\u03B2MPD"),
                      tags$td("TRUE / FALSE (default = TRUE)"),
                      tags$td("If TRUE, ignores matches of the same species across sites so the metric reflects phylogenetic (not taxonomic) distance.")
                    ),
                    tags$tr(
                      tags$td(code("bmpd_normalize")),
                      tags$td("\u03B2MPD"),
                      tags$td("TRUE / FALSE (default = TRUE)"),
                      tags$td("Standardizes distances for comparability (e.g., divide by max or expected).")
                    ),
                    tags$tr(
                      tags$td(code("bmpd_scale01")),
                      tags$td("\u03B2MPD"),
                      tags$td("TRUE / FALSE (default = FALSE)"),
                      tags$td("Optional rescaling to 0–1 after normalization; useful for plotting.")
                    ),

                    # ---- βMNTD ----
                    tags$tr(
                      tags$td(code("bmntd_method")),
                      tags$td("\u03B2MNTD"),
                      tags$td("exclude, root, node (default = exclude)"),
                      tags$td("Same idea as \u03B2MPD but applied to nearest-taxon distances (tip-level focus).")
                    ),
                    tags$tr(
                      tags$td(code("bmntd_abundance")),
                      tags$td("\u03B2MNTD"),
                      tags$td("TRUE / FALSE (default = FALSE)"),
                      tags$td("If TRUE, abundance-weights nearest-taxon distances between communities.")
                    ),
                    tags$tr(
                      tags$td(code("bmntd_exclude_conspecific")),
                      tags$td("\u03B2MNTD"),
                      tags$td("TRUE / FALSE (default = TRUE)"),
                      tags$td("If TRUE, ignores conspecific matches across communities to emphasize phylogenetic replacement.")
                    ),
                    tags$tr(
                      tags$td(code("bmntd_normalize")),
                      tags$td("\u03B2MNTD"),
                      tags$td("TRUE / FALSE (default = TRUE)"),
                      tags$td("Standardizes nearest-taxon distances to limit effects of richness and tree scaling.")
                    ),
                    tags$tr(
                      tags$td(code("bmntd_scale01")),
                      tags$td("\u03B2MNTD"),
                      tags$td("TRUE / FALSE (default = FALSE)"),
                      tags$td("Optional 0–1 rescaling for visualization and comparability.")
                    ),

                    # ---- Rao β ----
                    tags$tr(
                      tags$td(code("rb_abundance")),
                      tags$td("Rao \u03B2"),
                      tags$td("TRUE / FALSE (default = FALSE)"),
                      tags$td("If TRUE, incorporates community evenness/dominance into between-site phylogenetic differences.")
                    ),
                    tags$tr(
                      tags$td(code("rb_scale01")),
                      tags$td("Rao \u03B2"),
                      tags$td("TRUE / FALSE (default = TRUE)"),
                      tags$td("Rescale to 0–1 for easier interpretation and mapping.")
                    )
                  )
                ),
                p(em("Tips:"), " For any metric with both ", code("normalize"), " and ", code("scale01"),
                  ", normalization standardizes the metric (e.g., by expected or max), while ",
                  code("scale01"), " linearly rescales the standardized values into 0–1 for display.", br(), br()),

                p(strong("Output — "), code("vilma.beta"), br(), br()),
                tags$ul(
                  tags$li("Pairwise dissimilarity matrix and optional raster summaries"),
                  tags$li("Downloadable .rda with metadata")
                ),
                div(class="divider"),

                # ---------- 5) Export ----------
                h3(strong("5) Export"), br(), br()),
                tags$table(
                  tags$thead(tags$tr(tags$th("Object"), tags$th("Contents"), tags$th("Files"))),
                  tags$tbody(
                    tags$tr(tags$td(code("vilma.dist")), tags$td("Presence matrix, richness/abundance rasters"), tags$td(".csv, .tif")),
                    tags$tr(tags$td(code("vilma.pd")),   tags$td("\u03B1-index per cell, raster"), tags$td(".csv, .tif")),
                    tags$tr(tags$td(code("vilma.null")), tags$td("SES/p (global or cell), rasters"), tags$td(".csv, .tif")),
                    tags$tr(tags$td(code("vilma.beta")), tags$td("\u03B2-distance matrix, optional raster summaries"), tags$td(".csv, .tif"))
                  )
                ),

                br(), br(),
                div(style="text-align:center;",
                    actionButton("continue_howto","Continue \u2192 Points to Raster",
                                 class="btn btn-lg",
                                 style="background-color:#ff9999; border:none; color:white; font-weight:bold; padding:10px 25px;"))
              )
            )
          )
        ),

        # ----------------- Existing Analysis Tabs ---------------
        tabPanel("Points to Raster",       mod_points_to_raster_ui("ptr")),
        tabPanel("Phylogenetic Diversity", mod_alpha_ui("alpha")),
        tabPanel("Null Models",            mod_nulls_ui("nulls")),
        tabPanel("Beta Diversity",         mod_beta_ui("beta")),
        tabPanel("Write / Export",         mod_write_export_ui("write"))
      )
  )
)

# ============================= SERVER ==============================
server <- function(input, output, session) {

  # Welcome -> About
  observeEvent(input$continue_button, {
    updateTabsetPanel(session, "tabs", selected = "About")
  }, ignoreInit = TRUE)

  # About -> How to Use
  observeEvent(input$continue_about, {
    updateTabsetPanel(session, "tabs", selected = "How to Use")
  }, ignoreInit = TRUE)

  # How to Use -> Points to Raster
  observeEvent(input$continue_howto, {
    updateTabsetPanel(session, "tabs", selected = "Points to Raster")
  }, ignoreInit = TRUE)

  # ---- Modules ----
  ptr <- mod_points_to_raster_server("ptr")
  alpha <- mod_alpha_server("alpha", vilma_dist = ptr$vilma_dist)
  nulls <- mod_nulls_server(
    "nulls",
    vilma_dist   = ptr$vilma_dist,
    alpha_result = alpha$alpha_result,
    alpha_tree   = alpha$alpha_tree,
    alpha_index  = alpha$alpha_index
  )
  beta <- mod_beta_server(
    "beta",
    vilma_dist = ptr$vilma_dist,
    beta_tree  = alpha$alpha_tree
  )

  # -------------------- Navigation --------------------
  observeEvent(ptr$go_next(), {
    if (!is.null(ptr$vilma_dist())) {
      updateTabsetPanel(session, "tabs", selected = "Phylogenetic Diversity")
    } else {
      showNotification("Build a vilma.dist first.", type = "warning")
    }
  }, ignoreInit = TRUE)

  observeEvent(alpha$go_back(), {
    updateTabsetPanel(session, "tabs", selected = "Points to Raster")
  }, ignoreInit = TRUE)

  # Only after a successful Run; decide on Next. If last run was NRI/NNI, show modal here.
  observeEvent(alpha$go_next(), {
    if (!isTRUE(alpha$alpha_has_run())) {
      showNotification("Please run the analysis first.", type = "warning")
      return()
    }
    ran <- alpha$alpha_last_run_idx()
    if (is.null(ran)) {
      showNotification("No completed alpha result found.", type = "warning")
      return()
    }

    if (ran %in% c("nri","nti")) {
      showModal(modalDialog(
        title = "NRI / NNI routing",
        easyClose = TRUE,
        size = "m",
        tagList(
          p(strong("NRI"), "and", strong("NNI"),
            "don’t have standalone null-model configuration in the Alpha page."),
          p("They’re derived from MPD (for NRI) and MNTD (for NNI)."),
          p("Do you want to continue to the Beta page or stay on the Alpha page?")
        ),
        footer = tagList(
          modalButton("Stay on Alpha"),
          actionButton("alpha_next_confirm_beta", "Continue to Beta", class = "btn btn-primary")
        )
      ))
    } else {
      updateTabsetPanel(session, "tabs", selected = "Null Models")
    }
  }, ignoreInit = TRUE)

  observeEvent(input$alpha_next_confirm_beta, {
    removeModal()
    updateTabsetPanel(session, "tabs", selected = "Beta Diversity")
  }, ignoreInit = TRUE)

  observeEvent(nulls$go_back(), {
    updateTabsetPanel(session, "tabs", selected = "Phylogenetic Diversity")
  }, ignoreInit = TRUE)

  observeEvent(nulls$go_next(), {
    updateTabsetPanel(session, "tabs", selected = "Beta Diversity")
  }, ignoreInit = TRUE)

  observeEvent(nulls$go_beta(), {
    updateTabsetPanel(session, "tabs", selected = "Beta Diversity")
  }, ignoreInit = TRUE)

  observeEvent(beta$go_back(), {
    updateTabsetPanel(session, "tabs", selected = "Null Models")
  }, ignoreInit = TRUE)

  observeEvent(beta$go_next(), {
    updateTabsetPanel(session, "tabs", selected = "Write / Export")
  }, ignoreInit = TRUE)

  # ----------------- Results bundle for Write -----------------
  safeCall <- function(f) {
    out <- NULL
    if (is.null(f)) return(NULL)
    try(out <- f(), silent = TRUE)
    out
  }

  write_results <- reactive({
    list(
      dist    = safeCall(ptr$vilma_dist),
      alpha   = safeCall(alpha$alpha_result),
      beta    = safeCall(beta$beta_result),
      null    = safeCall(nulls$null_result),
      log     = NULL,
      summary = NULL
    )
  })

  mod_write_export_server("write", results = write_results)

  # ---- Close app ----
  observeEvent(input$`write-close_app`, {
    session$sendCustomMessage("vilma-close", list())
    later::later(function() {
      try(shiny::stopApp(), silent = TRUE)
      invisible(gc())
    }, 0.12)
  }, ignoreInit = TRUE)

  session$onSessionEnded(function() {
    try(shiny::stopApp(), silent = TRUE)
    invisible(gc())
  })
}

# ============================ RUN APP =============================
shinyApp(ui, server)

