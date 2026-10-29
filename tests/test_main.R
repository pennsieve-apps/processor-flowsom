# Unit tests for processor-flowsom on a tiny synthetic fixture (no Docker needed
# when R with FlowSOM is installed; `make test` otherwise runs this in the image).
#
#   Rscript tests/test_main.R
#
# The fixture: three samples of four populations on arcsinh-like scales
# (CD4 T, CD8 T, B cells, NK cells), FCS 3.1 float32 as fcs-prep writes them.

suppressPackageStartupMessages(suppressWarnings({
  library(flowCore)
  library(FlowSOM)
}))

repo <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  normalizePath(file.path(dirname(f[1]), ".."))
})
Sys.setenv(FLOWSOM_APP_DIR = repo, FLOWSOM_NO_MAIN = "1")
source(file.path(repo, "R", "main.R"))

failures <- 0L
check <- function(name, ok) {
  if (isTRUE(ok)) cat("ok   -", name, "\n") else { cat("FAIL -", name, "\n"); failures <<- failures + 1L }
}

MARKERS <- c(CD3 = "BUV563-A", CD4 = "BUV395-A", CD8 = "BUV805-A", CD19 = "BUV661-A",
             CD56 = "BUV496-A", CD16 = "BV480-A", HLADR = "BV750-A", CD69 = "BV605-A")
POPS <- list(
  "CD4 T" = c(CD3 = 2.5, CD4 = 2.5, CD8 = 0, CD19 = 0, CD56 = 0, CD16 = 0, HLADR = 0.2),
  "CD8 T" = c(CD3 = 2.5, CD4 = 0, CD8 = 3, CD19 = 0, CD56 = 0, CD16 = 0, HLADR = 0.2),
  "B"     = c(CD3 = 0, CD4 = 0, CD8 = 0, CD19 = 2.5, CD56 = 0, CD16 = 0, HLADR = 3),
  "NK"    = c(CD3 = 0, CD4 = 0, CD8 = 0.3, CD19 = 0, CD56 = 2.5, CD16 = 2, HLADR = 0.1))

make_sample <- function(path, n, seed, stimulated = FALSE) {
  set.seed(seed)
  sizes <- round(n * c(0.4, 0.25, 0.2, 0.15))
  x <- do.call(rbind, lapply(seq_along(POPS), function(i) {
    mu <- POPS[[i]]
    m <- sapply(names(mu), function(k) rnorm(sizes[i], mu[[k]], 0.25))
    cbind(m, CD69 = rnorm(sizes[i], if (stimulated && i <= 2) 2 else 0, 0.25))
  }))
  x <- x[sample(nrow(x)), ]
  ex <- cbind(`FSC-A` = runif(nrow(x), 5e4, 2e5), x[, names(MARKERS)], Time = seq_len(nrow(x)))
  colnames(ex) <- c("FSC-A", unname(MARKERS), "Time")
  storage.mode(ex) <- "double"
  ff <- flowFrame(ex)
  pd <- pData(parameters(ff))
  pd$desc <- c(NA, names(MARKERS), NA)
  pData(parameters(ff)) <- pd
  suppressWarnings(write.FCS(ff, path))
}

make_fixture <- function(root, layout = "linked", n = 1500) {
  fcs <- if (layout == "linked") file.path(root, "fcs") else root
  meta <- if (layout == "linked") file.path(root, "metadata") else root
  dir.create(fcs, recursive = TRUE, showWarnings = FALSE)
  dir.create(meta, recursive = TRUE, showWarnings = FALSE)
  ids <- c("s1_unstim", "s2_unstim", "s3_stim")
  for (i in seq_along(ids)) make_sample(file.path(fcs, paste0(ids[i], ".fcs")), n, 100 + i, grepl("stim$", ids[i]) && !grepl("unstim", ids[i]))
  write.csv(data.frame(sample_id = ids, file = paste0(ids, ".fcs"), batch = c("b1", "b1", "b2"),
                       condition = c("unstimulated", "unstimulated", "stimulated"), n_events = n),
            file.path(meta, "samples.csv"), row.names = FALSE)
  write.csv(data.frame(channel = c("FSC-A", unname(MARKERS), "Time"),
                       marker = c("FSC-A", names(MARKERS), "Time"),
                       type = c("scatter", rep("marker", length(MARKERS)), "time"),
                       cofactor = c("", rep("6000", length(MARKERS)), ""),
                       clustering = c("false", rep("true", 7), "false", "false"),
                       state = c("false", rep("false", 6), "true", "true", "false")),
            file.path(meta, "panel.csv"), row.names = FALSE)
  ids
}

run_app <- function(input, output, ...) {
  env <- c(INPUT_DIR = input, OUTPUT_DIR = output, XDIM = "4", YDIM = "4", N_METACLUSTERS = "6",
           CELLS_PER_SAMPLE = "600", CELLS_FOR_PLOT = "200", SEED = "42", MARKERS = "",
           CELL_TYPE_RULES = "", WRITE_FCS = "false")
  extra <- c(...)
  env[names(extra)] <- extra
  do.call(Sys.setenv, as.list(env))
  tryCatch({ main(); "ok" }, error = function(e) conditionMessage(e))
}

tmp <- tempfile("flowsom-test-")
dir.create(tmp)

# ---- rules ----
cfg <- load_rules("", repo)
check("built-in rules load", length(cfg$rules) >= 10 && cfg$high == 0.5)
med <- rbind(c(CD3 = 3, CD4 = 3, CD8 = 0), c(CD3 = 3, CD4 = 0, CD8 = 4), c(CD3 = 3, CD4 = 0, CD8 = 1.6), c(CD3 = 0, CD4 = 0.1, CD8 = 0))
lab <- assign_labels(med, load_rules('[{"label":"T4","markers":{"CD3":"+","CD4":"+"}},{"label":"T8","markers":{"CD3":"+","CD8":"+"}},{"label":"T8dim","markers":{"CD3":"+","CD8":"int"}},{"label":"X","markers":{"FOO":"+"}}]', repo))
check("rules: first match wins, int, unmatched is unassigned", identical(lab$label, c("T4", "T8", "T8dim", "unassigned")))
lab3 <- assign_labels(med, load_rules('[{"label":"T8pos","markers":{"CD3":"+","CD8":"pos"}}]', repo))
check("rules: pos is + or int", identical(lab3$label, c("unassigned", "T8pos", "T8pos", "unassigned")))
check("rules: a rule naming an absent marker is skipped", length(lab$skipped) == 1)
lab2 <- assign_labels(cbind(med, FOXP3 = c(0.1, 0.2, 0.15, 0.3)), load_rules('{"rules":[{"label":"Treg","markers":{"FOXP3":"+"}}]}', repo))
check("rules: a marker below min_range is negative everywhere", all(lab2$label == "unassigned"))
check("rules: HLA-DR matches HLADR", identical(norm_marker("HLA-DR"), norm_marker("HLADR")))
check("rules: bad condition rejected", inherits(tryCatch(load_rules('[{"label":"a","markers":{"CD3":"++"}}]', repo), error = function(e) e), "error"))

# ---- linked layout, default parameters ----
ids <- make_fixture(file.path(tmp, "in-linked"), "linked")
out <- file.path(tmp, "out-linked")
res <- run_app(file.path(tmp, "in-linked"), out)
check("linked: run succeeds", identical(res, "ok"))
cd <- file.path(out, "clusters")
for (f in c("counts.csv", "cluster_counts.csv", "medians.csv", "medians_by_sample.csv", "labels.csv", "model.rds", "heatmap.png", "star_tree.png", "provenance.json"))
  check(paste("linked: clusters/", f), file.exists(file.path(cd, f)))
check("linked: cells/cells.csv.gz", file.exists(file.path(out, "cells", "cells.csv.gz")))
check("linked: fcs/ has only provenance without WRITE_FCS", identical(list.files(file.path(out, "fcs")), "provenance.json"))
counts <- read.csv(file.path(cd, "counts.csv"))
check("counts: long format with zeros, sums to events", nrow(counts) == 3 * 6 && all(tapply(counts$n, counts$sample_id, sum) == 1500))
check("counts: frac sums to 1 per sample", all(abs(tapply(counts$frac, counts$sample_id, sum) - 1) < 1e-6))
cc <- read.csv(file.path(cd, "cluster_counts.csv"))
check("cluster_counts: 16 nodes per sample", nrow(cc) == 3 * 16 && all(tapply(cc$n, cc$sample_id, sum) == 1500))
labels <- read.csv(file.path(cd, "labels.csv"), stringsAsFactors = FALSE)
found <- unique(labels$label)
check(sprintf("labels: CD4 T, CD8 T, B and NK found (%s)", paste(found, collapse = " | ")),
      any(grepl("^CD4 T", found)) && any(grepl("^CD8 T", found)) && "B cells" %in% found && any(grepl("^NK", found)))
cd4 <- grepl("^CD4 T", labels$label)
check("labels: CD4 T is ~40% of events", abs(sum(labels$n_cells[cd4]) / sum(labels$n_cells) - 0.4) < 0.05)
med <- read.csv(file.path(cd, "medians.csv"), check.names = FALSE)
check("medians: metacluster x all markers", nrow(med) == 6 && all(names(MARKERS) %in% names(med)))
b <- med[med$label == "B cells", ]
check("medians: B cells CD19 ~2.5 (histogram median)", nrow(b) >= 1 && all(abs(b$CD19 - 2.5) < 0.2))
mbs <- read.csv(file.path(cd, "medians_by_sample.csv"))
check("medians_by_sample: long sample x metacluster x marker", nrow(mbs) == 3 * 6 * length(MARKERS) && all(c("sample_id", "metacluster", "marker", "median") %in% names(mbs)))
cd69 <- mbs[mbs$marker == "CD69" & grepl("^CD4 T", mbs$label), ]
check("medians_by_sample: stimulated CD4 T have higher CD69", mean(cd69$median[cd69$sample_id == "s3_stim"]) > mean(cd69$median[cd69$sample_id != "s3_stim"]) + 1)
cells <- read.csv(gzfile(file.path(out, "cells", "cells.csv.gz")), check.names = FALSE)
check("cells: 200 per sample, columns by marker name", nrow(cells) == 600 && identical(names(cells)[1:4], c("sample_id", "cluster", "metacluster", "label")) && all(names(MARKERS) %in% names(cells)))
cj <- jsonlite::fromJSON(file.path(out, "cells", "cells.json"))
check("cells.json names the clustering markers", identical(cj$clustering_markers, c("CD3", "CD4", "CD8", "CD19", "CD56", "CD16", "HLADR")))
prov <- jsonlite::fromJSON(file.path(cd, "provenance.json"), simplifyVector = FALSE)
check("provenance: app, version, tools, samples", prov$app == "processor-flowsom" && prov$version == "0.1.0" && !is.null(prov$tools$FlowSOM) && length(prov$samples) == 3)
check("provenance: inputs relative", all(!startsWith(unlist(prov$inputs), "/")))
model <- readRDS(file.path(cd, "model.rds"))
again <- NewData(model, read_markers(file.path(tmp, "in-linked", "fcs", "s1_unstim.fcs"), read_panel(file.path(tmp, "in-linked", "metadata", "panel.csv"))$panel)$x, silent = TRUE)
check("model.rds maps new data with NewData", length(GetClusters(again)) == 1500)

# ---- reproducible ----
out2 <- file.path(tmp, "out-linked-2")
invisible(run_app(file.path(tmp, "in-linked"), out2))
check("same SEED, same counts", identical(read.csv(file.path(out2, "clusters", "counts.csv")), counts))

# ---- flat layout (no links), symlinked entries, MARKERS, WRITE_FCS ----
ids <- make_fixture(file.path(tmp, "real"), "flat")
flat <- file.path(tmp, "in-flat")
dir.create(flat)
for (f in list.files(file.path(tmp, "real"), full.names = TRUE)) file.symlink(f, file.path(flat, basename(f)))
out3 <- file.path(tmp, "out-flat")
res <- run_app(flat, out3, MARKERS = "CD3, CD4,CD8,CD19,CD56,CD16,HLA-DR", WRITE_FCS = "true")
check("flat + symlinks: run succeeds", identical(res, "ok"))
check("flat: three samples counted", length(unique(read.csv(file.path(out3, "clusters", "counts.csv"))$sample_id)) == 3)
fo <- list.files(file.path(out3, "fcs"), pattern = "\\.fcs$")
check("WRITE_FCS writes one FCS per sample", length(fo) == 3)
ff <- read.FCS(file.path(out3, "fcs", fo[1]), transformation = FALSE, truncate_max_range = FALSE)
check("WRITE_FCS adds cluster channels and keeps events", all(c("FlowSOM_cluster", "FlowSOM_metacluster") %in% colnames(ff)) && nrow(ff) == 1500 &&
      !is.null(keyword(ff)[["PENNSIEVE_FLOWSOM"]]))
res <- run_app(flat, file.path(tmp, "out-bad"), MARKERS = "CD3,NOPE")
check("MARKERS with an unknown marker fails", grepl("NOPE", res))

# ---- a broken sample is skipped; all broken fails ----
writeLines("not an fcs file", file.path(tmp, "in-linked", "fcs", "zz_broken.fcs"))
out4 <- file.path(tmp, "out-broken")
res <- run_app(file.path(tmp, "in-linked"), out4)
p4 <- jsonlite::fromJSON(file.path(out4, "clusters", "provenance.json"), simplifyVector = FALSE)
st <- vapply(p4$samples, function(s) s$status, character(1))
check("one broken sample: run succeeds, sample marked failed", identical(res, "ok") && sum(st == "failed") == 1 && sum(st == "ok") == 3)
bad <- file.path(tmp, "in-allbad")
dir.create(bad)
writeLines("junk", file.path(bad, "a.fcs"))
res <- run_app(bad, file.path(tmp, "out-allbad"), MARKERS = "CD3")
check("every sample broken: run fails", !identical(res, "ok"))
res <- run_app(file.path(tmp, "empty-does-not-exist"), file.path(tmp, "out-none"))
check("no input: run fails", grepl("no .fcs files", res))

unlink(tmp, recursive = TRUE)
cat(sprintf("\n%d failure(s)\n", failures))
quit(status = if (failures) 1L else 0L, save = "no")
