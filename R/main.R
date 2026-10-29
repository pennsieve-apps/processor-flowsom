# processor-flowsom: FlowSOM clustering and metaclustering of transformed FCS files
# (Van Gassen et al. 2015), following the Nature Protocols practice (Quintelier et
# al. 2021): build the self-organising map on an aggregate of equal-sized,
# seeded per-sample subsamples, then map every event of every file onto it.
#
#   INPUT_DIR/fcs/<sample_id>.fcs        (or *.fcs anywhere under INPUT_DIR)
#   INPUT_DIR/metadata/samples.csv       optional
#   INPUT_DIR/metadata/panel.csv         clustering / state markers
#
#   OUTPUT_DIR/clusters/  counts, medians, labels, model, figures, provenance
#   OUTPUT_DIR/cells/     cells.csv.gz plotting subsample + cells.json
#   OUTPUT_DIR/fcs/       FCS files with cluster channels (WRITE_FCS=true)
#
# Memory: files are read one at a time, twice (pass 1: subsample + ranges;
# pass 2: map + count). Between files only the training aggregate, per-sample
# tables, the plotting subsample and fixed-size histograms are kept.

suppressPackageStartupMessages(suppressWarnings({
  library(flowCore)
  library(FlowSOM)
  library(ggplot2)
}))

# Figures go to files only; without this, functions that open the default
# device would try to write Rplots.pdf into the (read-only) working directory.
grDevices::pdf(NULL)

APP_NAME <- "processor-flowsom"
OUTPUTS <- c("clusters", "cells", "fcs")
HIST_BINS <- 4096L

app_dir <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- Sys.getenv("FLOWSOM_APP_DIR", "")
  if (d != "") normalizePath(d) else if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else normalizePath(".")
})
source(file.path(app_dir, "R", "labels.R"))

log_line <- function(fmt, ...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", sprintf(fmt, ...), "\n")
  flush(stdout())
}

env_str <- function(name, default = "") {
  v <- Sys.getenv(name, unset = "")
  if (identical(v, "")) default else v
}
env_int <- function(name, default) {
  v <- env_str(name, "")
  if (v == "") return(as.integer(default))
  x <- suppressWarnings(as.integer(v))
  if (is.na(x)) stop(sprintf("%s must be an integer, got '%s'", name, v))
  x
}
env_bool <- function(name, default = FALSE) {
  v <- tolower(env_str(name, ""))
  if (v == "") return(default)
  v %in% c("1", "true", "yes", "y", "on")
}

read_params <- function() {
  p <- list(
    XDIM = env_int("XDIM", 10), YDIM = env_int("YDIM", 10),
    N_METACLUSTERS = env_int("N_METACLUSTERS", 20),
    CELLS_PER_SAMPLE = env_int("CELLS_PER_SAMPLE", 20000),
    CELLS_FOR_PLOT = env_int("CELLS_FOR_PLOT", 5000),
    SEED = env_int("SEED", 42),
    MARKERS = env_str("MARKERS", ""),
    CELL_TYPE_RULES = env_str("CELL_TYPE_RULES", ""),
    WRITE_FCS = env_bool("WRITE_FCS", FALSE)
  )
  if (p$XDIM < 2 || p$YDIM < 2) stop("XDIM and YDIM must be at least 2")
  if (p$N_METACLUSTERS < 2 || p$N_METACLUSTERS > p$XDIM * p$YDIM)
    stop("N_METACLUSTERS must be between 2 and XDIM*YDIM")
  if (p$CELLS_PER_SAMPLE < 100) stop("CELLS_PER_SAMPLE must be at least 100")
  if (p$CELLS_FOR_PLOT < 0) stop("CELLS_FOR_PLOT must be >= 0")
  p
}

# ---- inputs -----------------------------------------------------------------

# Files by name under dir, following symlinks (INPUT_DIR entries are symlinks).
find_files <- function(dir, pattern) {
  if (!dir.exists(dir)) return(character(0))
  sort(list.files(dir, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE))
}

# The port folder when it exists (linked inputs), else anywhere under INPUT_DIR.
find_input <- function(input_dir, port, pattern) {
  hits <- find_files(file.path(input_dir, port), pattern)
  if (length(hits)) hits else find_files(input_dir, pattern)
}

rel_path <- function(path, root) {
  root <- normalizePath(root, mustWork = FALSE)
  p <- normalizePath(path, mustWork = FALSE)
  if (startsWith(p, paste0(root, "/"))) substring(p, nchar(root) + 2) else path
}

as_flag <- function(x) tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")

read_panel <- function(path, first_fcs) {
  if (!is.null(path)) {
    panel <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
    need <- c("channel", "marker", "type")
    if (!all(need %in% names(panel))) stop(sprintf("panel.csv needs columns %s", paste(need, collapse = ", ")))
    blank <- is.na(panel$marker) | panel$marker == ""
    panel$marker[blank] <- panel$channel[blank]
    panel$clustering <- if ("clustering" %in% names(panel)) as_flag(panel$clustering) else panel$type == "marker"
    panel$state <- if ("state" %in% names(panel)) as_flag(panel$state) else FALSE
    panel <- panel[panel$type == "marker", c("channel", "marker", "clustering", "state")]
    if (!nrow(panel)) stop("panel.csv has no rows of type 'marker'")
    return(list(panel = panel, source = "panel.csv"))
  }
  # No panel: every channel with a $PnS, other than scatter and time, is a
  # marker, and all of them cluster.
  ff <- read.FCS(first_fcs, which.lines = 1, transformation = FALSE, truncate_max_range = FALSE)
  pd <- pData(parameters(ff))
  desc <- as.character(pd$desc)
  keep <- !is.na(desc) & desc != "" & !grepl("^(FSC|SSC|Time)", pd$name, ignore.case = TRUE)
  if (!any(keep)) stop("no panel.csv and no channel with a marker name ($PnS)")
  panel <- data.frame(channel = as.character(pd$name[keep]), marker = desc[keep],
                      clustering = TRUE, state = FALSE, stringsAsFactors = FALSE)
  list(panel = panel, source = "fcs-keywords")
}

select_markers <- function(panel, spec) {
  if (trimws(spec) == "") {
    m <- panel$marker[panel$clustering]
    if (!length(m)) stop("the panel marks no marker with clustering=true; set MARKERS")
    return(m)
  }
  wanted <- trimws(strsplit(spec, ",")[[1]])
  wanted <- wanted[wanted != ""]
  idx <- vapply(wanted, function(w) {
    i <- which(norm_marker(panel$marker) == norm_marker(w))
    if (!length(i)) i <- which(panel$channel == w)
    if (!length(i)) stop(sprintf("MARKERS: '%s' is not a marker in the panel (%s)", w, paste(panel$marker, collapse = ", ")))
    i[1]
  }, integer(1))
  unique(panel$marker[idx])
}

sample_id_of <- function(path) gsub("[^A-Za-z0-9_-]", "_", sub("\\.fcs$", "", basename(path), ignore.case = TRUE))

# Marker matrix of one file, columns named by marker.
read_markers <- function(path, panel) {
  ff <- read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
  miss <- setdiff(panel$channel, colnames(ff))
  if (length(miss)) stop(sprintf("channels missing: %s", paste(miss, collapse = ", ")))
  x <- exprs(ff)[, panel$channel, drop = FALSE]
  colnames(x) <- panel$marker
  list(ff = ff, x = x)
}

# ---- outputs ----------------------------------------------------------------

write_csv <- function(df, path) write.csv(df, path, row.names = FALSE, na = "")

write_json <- function(x, path) {
  writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA), path)
}

# Median from a histogram (linear interpolation within the bin).
hist_median <- function(counts, lo, width) {
  n <- sum(counts)
  if (n == 0) return(NA_real_)
  cs <- cumsum(counts)
  b <- which(cs >= n / 2)[1]
  before <- if (b > 1) cs[b - 1] else 0
  lo + width * (b - 1 + (n / 2 - before) / counts[b])
}

plot_heatmap <- function(medians, labels, frac, markers_order, clustering, path) {
  k <- nrow(medians)
  m <- medians[, markers_order, drop = FALSE]
  s <- scale_medians(m, 0)
  rows <- sprintf("%d  %s  (%.1f%%)", seq_len(k), labels, 100 * frac)
  df <- data.frame(
    metacluster = factor(rep(rows, times = ncol(s)), levels = rev(rows)),
    marker = factor(rep(colnames(s), each = k), levels = colnames(s)),
    value = as.vector(s), median = as.vector(m))
  g <- ggplot(df, aes(marker, metacluster, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f", median)), size = 2.2,
              colour = ifelse(is.na(df$value) | df$value > 0.6, "black", "white")) +
    scale_fill_viridis_c(name = "scaled\nmedian", limits = c(0, 1), na.value = "grey90") +
    labs(x = NULL, y = NULL,
         title = "FlowSOM metacluster medians (arcsinh), all events",
         subtitle = "Numbers: median. Colour: median scaled 0-1 across metaclusters per marker. Clustering markers left of the line.") +
    geom_vline(xintercept = length(clustering) + 0.5, colour = "black", linewidth = 0.6) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid = element_blank())
  ggsave(path, g, width = max(8, 0.32 * ncol(s) + 4), height = max(5, 0.28 * k + 2), dpi = 120, bg = "white")
}

plot_stars <- function(fsom, labels, path) {
  # Columns are already named by marker; drop FlowSOM's "<channel>" suffix.
  fsom$prettyColnames <- setNames(names(fsom$prettyColnames), names(fsom$prettyColnames))
  mc <- as.integer(fsom$metaclustering)
  lv <- sprintf("%d %s", seq_along(labels), labels)
  bg <- factor(lv[mc], levels = lv)
  g <- PlotStars(fsom, backgroundValues = bg,
                 title = "FlowSOM tree of the training aggregate: star = node medians of the clustering markers, background = metacluster")
  ggsave(path, g, width = 13, height = 10, dpi = 110, bg = "white")
}

tool_versions <- function() {
  pk <- c("FlowSOM", "flowCore", "ConsensusClusterPlus", "igraph", "ggplot2", "jsonlite")
  v <- lapply(pk, function(p) as.character(packageVersion(p)))
  names(v) <- pk
  c(list(R = paste(R.version$major, R.version$minor, sep = ".")), v)
}

secs_since <- function(t) round(as.numeric(difftime(Sys.time(), t, units = "secs")), 2)

# ---- main -------------------------------------------------------------------

main <- function() {
  t0 <- Sys.time()
  input_dir <- env_str("INPUT_DIR", "")
  output_dir <- env_str("OUTPUT_DIR", "")
  if (input_dir == "" || output_dir == "") stop("INPUT_DIR and OUTPUT_DIR must be set")
  p <- read_params()
  version <- trimws(readLines(file.path(app_dir, "VERSION"), warn = FALSE)[1])
  log_line("%s %s: FlowSOM %s, INPUT_DIR=%s OUTPUT_DIR=%s", APP_NAME, version,
           as.character(packageVersion("FlowSOM")), input_dir, output_dir)

  out <- setNames(file.path(output_dir, OUTPUTS), OUTPUTS)
  for (d in out) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  fcs_files <- find_input(input_dir, "fcs", "\\.fcs$")
  if (!length(fcs_files)) stop(sprintf("no .fcs files found under %s", input_dir))
  samples_path <- find_input(input_dir, "metadata", "^samples\\.csv$")
  panel_path <- find_input(input_dir, "metadata", "^panel\\.csv$")
  pn <- read_panel(if (length(panel_path)) panel_path[1] else NULL, fcs_files[1])
  panel <- pn$panel
  markers <- panel$marker
  clustering <- select_markers(panel, p$MARKERS)
  rules <- load_rules(p$CELL_TYPE_RULES, app_dir)
  log_line("%d FCS files; panel from %s: %d markers; clustering on %d: %s",
           length(fcs_files), pn$source, length(markers), length(clustering), paste(clustering, collapse = ", "))

  ids <- vapply(fcs_files, sample_id_of, character(1), USE.NAMES = FALSE)
  ord <- order(ids, fcs_files)
  fcs_files <- fcs_files[ord]
  ids <- ids[ord]
  dup <- duplicated(ids)
  if (any(dup)) log_line("WARNING: duplicate sample ids, skipping: %s", paste(fcs_files[dup], collapse = ", "))
  fcs_files <- fcs_files[!dup]
  ids <- ids[!dup]
  status <- setNames(lapply(ids, function(id) list(sample_id = id, status = "pending")), ids)

  # Per-sample seeds drawn once from SEED, in sample-id order.
  set.seed(p$SEED)
  seeds <- setNames(sample.int(.Machine$integer.max - 1L, length(ids)), ids)

  # ---- pass 1: per-sample training subsample + marker ranges ----
  t1 <- Sys.time()
  train <- list()
  lo <- setNames(rep(Inf, length(markers)), markers)
  hi <- setNames(rep(-Inf, length(markers)), markers)
  for (i in seq_along(ids)) {
    id <- ids[i]
    ts <- Sys.time()
    r <- tryCatch({
      x <- read_markers(fcs_files[i], panel)$x
      n <- nrow(x)
      if (n == 0) stop("no events")
      set.seed(seeds[[id]])
      take <- if (n > p$CELLS_PER_SAMPLE) sort(sample.int(n, p$CELLS_PER_SAMPLE)) else seq_len(n)
      rg <- apply(x, 2, function(v) { v <- v[is.finite(v)]; if (length(v)) range(v) else c(Inf, -Inf) })
      list(sub = x[take, , drop = FALSE], n = n, rg = rg)
    }, error = function(e) e)
    if (inherits(r, "error")) {
      status[[id]] <- list(sample_id = id, status = "failed", stage = "read", error = conditionMessage(r))
      log_line("%s: FAILED reading: %s", id, conditionMessage(r))
      next
    }
    lo <- pmin(lo, r$rg[1, ])
    hi <- pmax(hi, r$rg[2, ])
    train[[id]] <- r$sub
    status[[id]] <- list(sample_id = id, status = "read", n_events = r$n, n_training = nrow(r$sub))
    log_line("%s: %d events, %d into the training aggregate (%.1fs)", id, r$n, nrow(r$sub), secs_since(ts))
  }
  if (!length(train)) stop("every sample failed to read")
  agg <- do.call(rbind, train)
  n_training <- nrow(agg)
  rm(train)
  invisible(gc(verbose = FALSE))

  # ---- build the SOM + metaclusters on the aggregate ----
  t2 <- Sys.time()
  if (n_training < p$XDIM * p$YDIM) stop(sprintf("only %d training cells for a %dx%d grid", n_training, p$XDIM, p$YDIM))
  log_line("building a %dx%d SOM on %d cells from %d samples, %d metaclusters (seed %d)",
           p$XDIM, p$YDIM, n_training, sum(vapply(status, function(s) s$status == "read", logical(1))),
           p$N_METACLUSTERS, p$SEED)
  fsom <- FlowSOM(agg, compensate = FALSE, transform = FALSE, scale = FALSE,
                  colsToUse = clustering, xdim = p$XDIM, ydim = p$YDIM,
                  nClus = p$N_METACLUSTERS, seed = p$SEED, silent = TRUE)
  n_nodes <- p$XDIM * p$YDIM
  meta_of_node <- as.integer(fsom$metaclustering)
  k <- nlevels(fsom$metaclustering)
  rm(agg)
  t3 <- Sys.time()
  log_line("SOM and metaclusters built in %.1fs", as.numeric(difftime(t3, t2, units = "secs")))

  # ---- pass 2: map every event of every file ----
  width <- (hi - lo) / HIST_BINS
  bad <- !is.finite(width) | width <= 0
  width[bad] <- 1
  lo[!is.finite(lo)] <- 0
  hist <- array(0, dim = c(k, length(markers), HIST_BINS))
  cluster_counts <- list()
  meta_counts <- list()
  by_sample <- list()
  cells <- list()
  for (i in seq_along(ids)) {
    id <- ids[i]
    if (status[[id]]$status != "read") next
    ts <- Sys.time()
    r <- tryCatch({
      rd <- read_markers(fcs_files[i], panel)
      x <- rd$x
      mapped <- suppressWarnings(NewData(fsom, x, silent = TRUE))
      cl <- as.integer(GetClusters(mapped))
      far <- sum(mapped$outliers$perCluster$Number_of_outliers)
      rm(mapped)
      mc <- meta_of_node[cl]
      h <- array(0, dim = c(k, length(markers), HIST_BINS))
      for (j in seq_along(markers)) {
        b <- floor((x[, j] - lo[j]) / width[j]) + 1
        b[!is.finite(b)] <- 1
        b <- pmin(pmax(b, 1), HIST_BINS)
        h[, j, ] <- tabulate((b - 1) * k + mc, nbins = k * HIST_BINS)
      }
      med <- matrix(NA_real_, k, length(markers), dimnames = list(NULL, markers))
      grp <- split(seq_along(mc), mc)
      for (g in names(grp)) med[as.integer(g), ] <- apply(x[grp[[g]], , drop = FALSE], 2, median)
      set.seed(seeds[[id]] + 1L)
      take <- if (nrow(x) > p$CELLS_FOR_PLOT) sort(sample.int(nrow(x), p$CELLS_FOR_PLOT)) else seq_len(nrow(x))
      if (p$WRITE_FCS) {
        ff <- fr_append_cols(rd$ff, cbind(FlowSOM_cluster = as.numeric(cl), FlowSOM_metacluster = as.numeric(mc)))
        keyword(ff)[["PENNSIEVE_FLOWSOM"]] <- sprintf("xdim=%d;ydim=%d;metaclusters=%d;seed=%d", p$XDIM, p$YDIM, k, p$SEED)
        suppressWarnings(write.FCS(ff, file.path(out[["fcs"]], paste0(id, ".fcs"))))
        rm(ff)
      }
      list(cc = tabulate(cl, nbins = n_nodes), mcc = tabulate(mc, nbins = k), med = med, h = h,
           cells = list(cluster = cl[take], x = x[take, , drop = FALSE]), n = nrow(x), far = far)
    }, error = function(e) e)
    if (inherits(r, "error")) {
      status[[id]] <- modifyList(status[[id]], list(status = "failed", stage = "map", error = conditionMessage(r)))
      log_line("%s: FAILED mapping: %s", id, conditionMessage(r))
      next
    }
    hist <- hist + r$h
    cluster_counts[[id]] <- r$cc
    meta_counts[[id]] <- r$mcc
    by_sample[[id]] <- r$med
    cells[[id]] <- r$cells
    secs <- secs_since(ts)
    status[[id]] <- modifyList(status[[id]], list(status = "ok", n_mapped = r$n, n_far_from_node = r$far, seconds_mapping = secs))
    top <- order(r$mcc, decreasing = TRUE)[seq_len(min(3, k))]
    log_line("%s: mapped %d events (%.1fs); largest metaclusters %s", id, r$n, secs,
             paste(sprintf("%d:%.1f%%", top, 100 * r$mcc[top] / r$n), collapse = " "))
    rm(r)
    invisible(gc(verbose = FALSE))
  }
  ok_ids <- names(meta_counts)
  if (!length(ok_ids)) stop("every sample failed")
  t4 <- Sys.time()

  # ---- medians over all events + labels ----
  medians <- matrix(NA_real_, k, length(markers), dimnames = list(NULL, markers))
  for (j in seq_along(markers)) for (m in seq_len(k)) medians[m, j] <- hist_median(hist[m, j, ], lo[j], width[j])
  rm(hist)
  n_meta <- Reduce(`+`, meta_counts)
  used <- n_meta > 0
  lab <- assign_labels(medians[used, , drop = FALSE], rules)
  labels <- rep("empty", k)
  labels[used] <- lab$label
  rule_idx <- rep(NA_integer_, k)
  rule_idx[used] <- lab$rule
  scaled <- matrix(NA_real_, k, length(markers), dimnames = list(NULL, markers))
  scaled[used, ] <- lab$scaled
  if (length(lab$skipped)) log_line("rules skipped (marker not in the data): %s", paste(lab$skipped, collapse = "; "))
  profile <- vapply(seq_len(k), function(m) if (!used[m]) "" else marker_profile(scaled[m, ], clustering, rules), character(1))
  for (m in seq_len(k)) log_line("metacluster %2d: %-30s %9d events  %s", m, labels[m], n_meta[m], profile[m])

  # ---- clusters/ ----
  cd <- out[["clusters"]]
  write_csv(do.call(rbind, lapply(ok_ids, function(id) {
    n <- meta_counts[[id]]
    data.frame(sample_id = id, metacluster = seq_len(k), label = labels, n = n, frac = round(n / sum(n), 8))
  })), file.path(cd, "counts.csv"))
  write_csv(do.call(rbind, lapply(ok_ids, function(id) {
    n <- cluster_counts[[id]]
    data.frame(sample_id = id, cluster = seq_len(n_nodes), metacluster = meta_of_node, n = n, frac = round(n / sum(n), 8))
  })), file.path(cd, "cluster_counts.csv"))
  write_csv(data.frame(metacluster = seq_len(k), label = labels, n_cells = n_meta, round(medians, 4), check.names = FALSE),
            file.path(cd, "medians.csv"))
  write_csv(do.call(rbind, lapply(ok_ids, function(id) {
    data.frame(sample_id = id, metacluster = rep(seq_len(k), each = length(markers)),
               label = rep(labels, each = length(markers)), marker = rep(markers, times = k),
               median = round(as.vector(t(by_sample[[id]])), 4),
               n = rep(meta_counts[[id]], each = length(markers)))
  })), file.path(cd, "medians_by_sample.csv"))
  write_csv(data.frame(metacluster = seq_len(k), label = labels, n_cells = n_meta,
                       frac = round(n_meta / sum(n_meta), 8), rule = rule_idx, profile = profile),
            file.path(cd, "labels.csv"))
  # The model without the training events (fsom$data); the node assignments of
  # the training events stay, NewData compares new data's distances with them.
  model <- fsom
  model$data <- NULL
  saveRDS(model, file.path(cd, "model.rds"))
  rm(model)
  heat_order <- c(clustering, setdiff(markers, clustering))
  fig_errors <- list()
  e <- tryCatch({ plot_heatmap(medians, labels, n_meta / sum(n_meta), heat_order, clustering, file.path(cd, "heatmap.png")); NULL },
                error = function(e) conditionMessage(e))
  if (!is.null(e)) { fig_errors$heatmap <- e; log_line("WARNING: heatmap failed: %s", e) }
  e <- tryCatch({ plot_stars(fsom, labels, file.path(cd, "star_tree.png")); NULL }, error = function(e) conditionMessage(e))
  if (!is.null(e)) { fig_errors$star_tree <- e; log_line("WARNING: star tree failed: %s", e) }

  # ---- cells/ ----
  gz <- gzfile(file.path(out[["cells"]], "cells.csv.gz"), "w")
  n_cells <- 0
  for (j in seq_along(ok_ids)) {
    d <- cells[[ok_ids[j]]]
    df <- data.frame(sample_id = rep(ok_ids[j], length(d$cluster)), cluster = d$cluster,
                     metacluster = meta_of_node[d$cluster], label = labels[meta_of_node[d$cluster]],
                     round(d$x, 4), check.names = FALSE)
    write.table(df, gz, sep = ",", row.names = FALSE, col.names = (j == 1), na = "", qmethod = "double")
    n_cells <- n_cells + nrow(df)
  }
  close(gz)
  write_json(list(clustering_markers = clustering, state_markers = panel$marker[panel$state],
                  markers = markers, cells_per_sample = p$CELLS_FOR_PLOT, seed = p$SEED,
                  n_cells = n_cells, n_samples = length(ok_ids),
                  values = "marker values as in the input FCS files (arcsinh-transformed by fcs-prep)"),
             file.path(out[["cells"]], "cells.json"))

  # ---- provenance ----
  t5 <- Sys.time()
  secs <- function(a, b) round(as.numeric(difftime(b, a, units = "secs")), 2)
  inputs_seen <- c(vapply(fcs_files, rel_path, character(1), root = input_dir, USE.NAMES = FALSE),
                   if (length(samples_path)) rel_path(samples_path[1], input_dir),
                   if (length(panel_path)) rel_path(panel_path[1], input_dir))
  prov <- list(
    app = APP_NAME, version = version, tools = tool_versions(),
    parameters = p,
    method = list(
      training = sprintf("SOM built on an aggregate of up to CELLS_PER_SAMPLE=%d seeded random events per sample (%d events in all); every event of every file then mapped to its nearest node with FlowSOM::NewData",
                         p$CELLS_PER_SAMPLE, n_training),
      metaclustering = "FlowSOM default: consensus hierarchical clustering (ConsensusClusterPlus) of the node codes into N_METACLUSTERS",
      transform = "none: inputs are already arcsinh-transformed by fcs-prep; scale=FALSE",
      medians = sprintf("medians.csv: per metacluster over all mapped events, from a %d-bin histogram per marker spanning its full range (error below one bin width, (max - min) / %d); medians_by_sample.csv: exact, per file", HIST_BINS, HIST_BINS),
      labels = rules$description,
      label_thresholds = list(high = rules$high, low = rules$low, min_range = rules$min_range),
      rules = lapply(rules$rules, function(r) list(label = r$label, markers = as.list(r$markers))),
      rules_skipped = lab$skipped),
    clustering_markers = clustering,
    state_markers = panel$marker[panel$state],
    panel_source = pn$source,
    inputs = unname(inputs_seen),
    figure_errors = fig_errors,
    timings = list(pass1_read_seconds = secs(t1, t2), som_seconds = secs(t2, t3),
                   pass2_map_seconds = secs(t3, t4), outputs_seconds = secs(t4, t5), total_seconds = secs(t0, t5)),
    samples = unname(status),
    created = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  for (o in OUTPUTS) {
    pr <- c(list(output = o), prov)
    if (o == "fcs") pr$write_fcs <- p$WRITE_FCS
    write_json(pr, file.path(out[[o]], "provenance.json"))
  }
  n_failed <- sum(vapply(status, function(s) s$status == "failed", logical(1)))
  log_line("done: %d samples ok, %d failed, %d metaclusters, %.1fs", length(ok_ids), n_failed, k, secs(t0, t5))
  invisible(0)
}

if (sys.nframe() == 0L && !identical(Sys.getenv("FLOWSOM_NO_MAIN"), "1")) {
  rc <- tryCatch({ main(); 0L }, error = function(e) {
    log_line("ERROR: %s", conditionMessage(e))
    1L
  })
  quit(status = rc, save = "no")
}
