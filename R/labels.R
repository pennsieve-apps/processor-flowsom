# Metacluster labels from a marker table (CELL_TYPE_RULES).
#
# A transparent heuristic, not a classifier. For every marker, the metacluster
# medians are scaled to 0..1 across the metaclusters (0 = the lowest metacluster
# median, 1 = the highest). A marker whose metacluster medians span less than
# `min_range` (arcsinh units) separates nothing, so every metacluster counts as
# negative for it. A rule is a list of marker conditions: "+" needs a scaled
# value >= `high`, "-" needs <= `low`, "int" needs a value strictly between,
# "pos" needs > `low` ("+" or "int"). Rules are tried in order; the first rule
# whose every condition holds names the metacluster, else it is "unassigned".
# A rule naming a marker that is not in the data is skipped (and reported).

norm_marker <- function(x) gsub("[^A-Z0-9]", "", toupper(x))

default_rules_path <- function(app_dir) file.path(app_dir, "R", "cell_type_rules.json")

# spec: "" / "default" for the built-in table, else a JSON object
# {high, low, min_range, rules: [{label, markers: {CD3: "+", ...}}]} or a bare
# JSON array of rules. Missing high/low/min_range take the built-in values.
load_rules <- function(spec, app_dir) {
  base <- jsonlite::fromJSON(default_rules_path(app_dir), simplifyVector = FALSE)
  spec <- trimws(if (is.null(spec)) "" else spec)
  if (spec == "" || tolower(spec) == "default") {
    cfg <- base
  } else {
    user <- jsonlite::fromJSON(spec, simplifyVector = FALSE)
    if (is.null(names(user))) user <- list(rules = user)
    cfg <- base
    for (k in names(user)) cfg[[k]] <- user[[k]]
    if (is.null(user$rules)) stop("CELL_TYPE_RULES: no 'rules' given")
  }
  rules <- lapply(cfg$rules, function(r) {
    if (is.null(r$label) || is.null(r$markers) || length(r$markers) == 0)
      stop("CELL_TYPE_RULES: every rule needs a label and at least one marker condition")
    conds <- unlist(r$markers)
    if (!all(conds %in% c("+", "-", "int", "pos")))
      stop(sprintf("CELL_TYPE_RULES: rule '%s': conditions must be '+', '-', 'int' or 'pos'", r$label))
    list(label = as.character(r$label), markers = conds)
  })
  list(high = as.numeric(cfg$high), low = as.numeric(cfg$low),
       min_range = as.numeric(cfg$min_range), rules = rules,
       description = if (is.null(cfg$description)) "" else cfg$description)
}

# medians: metacluster x marker matrix (colnames = marker names).
scale_medians <- function(medians, min_range) {
  s <- medians
  for (j in seq_len(ncol(medians))) {
    x <- medians[, j]
    r <- range(x, na.rm = TRUE)
    s[, j] <- if (!all(is.finite(r)) || diff(r) < min_range) 0 else (x - r[1]) / diff(r)
  }
  s
}

assign_labels <- function(medians, cfg) {
  s <- scale_medians(medians, cfg$min_range)
  have <- norm_marker(colnames(medians))
  labels <- rep("unassigned", nrow(medians))
  rule_used <- rep(NA_integer_, nrow(medians))
  skipped <- character(0)
  usable <- vapply(seq_along(cfg$rules), function(i) {
    miss <- setdiff(norm_marker(names(cfg$rules[[i]]$markers)), have)
    if (length(miss)) skipped <<- c(skipped, sprintf("%s (no %s)", cfg$rules[[i]]$label, paste(miss, collapse = ", ")))
    length(miss) == 0
  }, logical(1))
  for (k in seq_len(nrow(medians))) {
    for (i in which(usable)) {
      conds <- cfg$rules[[i]]$markers
      v <- s[k, match(norm_marker(names(conds)), have)]
      ok <- ifelse(conds == "+", v >= cfg$high,
            ifelse(conds == "-", v <= cfg$low,
            ifelse(conds == "pos", v > cfg$low, v > cfg$low & v < cfg$high)))
      if (all(ok)) {
        labels[k] <- cfg$rules[[i]]$label
        rule_used[k] <- i
        break
      }
    }
  }
  list(label = labels, rule = rule_used, scaled = s, skipped = skipped)
}

# "CD3+ CD4+ CD8- CD45RAint ..." (int: between low and high) for the given markers, from scaled values.
marker_profile <- function(scaled_row, markers, cfg) {
  v <- scaled_row[markers]
  tag <- ifelse(v >= cfg$high, "+", ifelse(v <= cfg$low, "-", "int"))
  paste0(markers, tag, collapse = " ")
}
