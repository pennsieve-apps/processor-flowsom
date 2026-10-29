# processor-flowsom

Unsupervised clustering of flow / spectral / mass cytometry data with **FlowSOM**
(Van Gassen et al. 2015), the clustering step of the Immune Profiling & Differential
Abundance workflow (`fcs-prep → peacoqc → cleanet → flowsom → {cytometry-umap, diffcyt} →
cytometry-report`).

FlowSOM trains a self-organising map (by default a 10 × 10 grid, so 100 small clusters,
"nodes") on the cells' marker expression, then groups the nodes into a smaller number of
**metaclusters** (consensus hierarchical clustering) that usually correspond to cell
populations. The app follows the practice of the FlowSOM Nature Protocols paper
(Quintelier et al. 2021): the map is built on an **aggregate of equal-sized, seeded random
subsamples of every sample** (so a large sample does not dominate), and **every event of
every file** is then mapped to its nearest node, one file at a time. The outputs are what
the next steps need: per-sample cell counts per metacluster (differential abundance),
per-sample marker medians per metacluster (differential state), a metacluster × marker
heatmap, the FlowSOM star tree, and a per-sample cell subsample for UMAP.

Each metacluster also gets a **label** (`CD4 T CD45RA-`, `B cells`, `NK CD56bright`, …) from
a small, published-in-the-repo marker table. This is a transparent heuristic to make the
tables readable, not a cell-type classifier; check the heatmap.

## Ports

| Port | Kind | Files | Notes |
|---|---|---|---|
| `fcs` | input (subfolder) | `<sample_id>.fcs` | FCS 3.1, marker channels arcsinh-transformed by fcs-prep (after peacoqc / cleanet). Looked for in `INPUT_DIR/fcs/`, else anywhere under `INPUT_DIR` (symlinks followed). The sample id is the file stem. Values are used as they are: no compensation, transform or scaling. |
| `metadata` | input (subfolder, optional) | `panel.csv`, `samples.csv` | `panel.csv` (`channel, marker, type, clustering, state`) chooses the clustering markers (`clustering == true`) and lists the markers reported (`type == marker`). Without it every channel with a `$PnS` name other than scatter/time is a marker and clusters. `samples.csv` is not needed here. |
| `clusters` | output | see below | Tables, model, figures, `provenance.json`. Input of diffcyt and cytometry-report. |
| `cells` | output | `cells.csv.gz`, `cells.json`, `provenance.json` | Plotting subsample. Input of cytometry-umap. |
| `fcs` | output | `<sample_id>.fcs` (only with `WRITE_FCS=true`), `provenance.json` | Every event, all original channels and keywords, plus `FlowSOM_cluster` and `FlowSOM_metacluster` channels and keyword `PENNSIEVE_FLOWSOM=xdim=…;ydim=…;metaclusters=…;seed=…`. Declared so it can be linked to a later step (e.g. gating on metaclusters); with the default `WRITE_FCS=false` it holds only `provenance.json`, which keeps the default run small (the demo data would add 1.3 GB). |

## Parameters

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `XDIM`, `YDIM` | integer | `10`, `10` | SOM grid size: XDIM × YDIM clusters (nodes). |
| `N_METACLUSTERS` | integer | `20` | Metaclusters the nodes are grouped into (FlowSOM's `nClus`, ConsensusClusterPlus). |
| `CELLS_PER_SAMPLE` | integer | `20000` | Events drawn at random from each sample for the training aggregate (all of a smaller sample). 24 samples × 20,000 = 480,000 training events on the demo data. |
| `CELLS_FOR_PLOT` | integer | `5000` | Events per sample in `cells/cells.csv.gz`. |
| `SEED` | integer | `42` | Seeds every random step: per-sample seeds are drawn from it in sample-id order, and FlowSOM gets it for the SOM and metaclustering. Same inputs + same seed = same outputs. |
| `MARKERS` | string | empty | Comma list of markers (panel `marker`, case and punctuation ignored, so `HLA-DR` = `HLADR`; or a channel name) to cluster on, overriding the panel. |
| `CELL_TYPE_RULES` | string | empty | JSON marker table for the labels (format below). Empty: the built-in table. |
| `WRITE_FCS` | boolean | `false` | Also write each FCS with the cluster channels to the `fcs` output. |

## Outputs

`clusters/`

| File | Columns / content |
|---|---|
| `counts.csv` | `sample_id, metacluster, label, n, frac` — long format, one row per sample × metacluster, zeros included; `frac` = n / events of the sample. |
| `cluster_counts.csv` | `sample_id, cluster, metacluster, n, frac` — the same per SOM node (1 … XDIM·YDIM). |
| `medians.csv` | `metacluster, label, n_cells`, then one column per marker (all `type == marker` channels, by marker name): the median over **all mapped events** of the run. Computed from a 4096-bin histogram per marker over its full range, accumulated file by file (error below one bin width, (max − min) / 4096, about 0.002 arcsinh units). |
| `medians_by_sample.csv` | `sample_id, metacluster, label, marker, median, n` — exact per-sample medians (long), for diffcyt's differential-state tests; `median` is empty where a sample has no events in a metacluster. |
| `labels.csv` | `metacluster, label, n_cells, frac, rule, profile` — `rule` is the index of the matching rule (empty = unassigned); `profile` is the metacluster's clustering-marker profile in the rule terms (`CD3+ CD4+ CD8- CD45RAint …`). |
| `model.rds` | The FlowSOM object without the training events (`fsom$data` removed): `FlowSOM::NewData(readRDS("model.rds"), x)` maps new data (a matrix with columns named by marker). |
| `heatmap.png` | Metacluster × marker medians (numbers), coloured by the median scaled 0–1 across metaclusters per marker; clustering markers left of the line; rows named with label and % of events. |
| `star_tree.png` | FlowSOM's minimal spanning tree of the nodes (training aggregate): each star shows the node's clustering-marker medians, the background colour its metacluster. |
| `provenance.json` | App and tool versions, parameters, method notes, the rule table, inputs seen, timings, per-sample status (events, training events, mapped, events far from their node, seconds). |

`cells/cells.csv.gz`: `sample_id, cluster, metacluster, label`, then every marker by marker
name, values as in the input (arcsinh). `cells/cells.json`: `clustering_markers`,
`state_markers`, `markers`, `cells_per_sample`, `seed`, `n_cells`.

## Metacluster labels

For each marker, the metacluster medians (`medians.csv`) are rescaled to 0–1 across the
metaclusters: 0 is the lowest metacluster, 1 the highest. A marker whose metacluster medians
span less than `min_range` (0.5 arcsinh units) separates nothing and counts as negative
everywhere. A rule is a list of conditions:

- `+`: scaled value ≥ `high` (0.5)
- `-`: scaled value ≤ `low` (0.3)
- `int`: in between
- `pos`: > `low` (`+` or `int`)

Rules are tried in order and the **first rule whose every condition holds** names the
metacluster; otherwise it is `unassigned`. Several metaclusters can share a label. Rules
use every marker in the panel, not only the clustering ones (CD25 and HLA-DR are state
markers in the demo panel, but still help to name populations). A rule that names a marker
the data lacks is skipped (listed in the log and provenance). Marker names are matched
ignoring case and punctuation.

Built-in table (`R/cell_type_rules.json`), for human blood leukocytes pre-gated
CD45+CD14− (mostly T, B, NK cells, DCs, basophils):

| # | Label | Conditions |
|---|---|---|
| 1 | T/B doublet-like | CD3+ CD19+ |
| 2 | CD4 T CD25+FOXP3+ (Treg or activated) | CD3+ CD4+ CD8− CD25+ FOXP3+ |
| 3 | CD4 T CD45RA+ | CD3+ CD4+ CD8− CD45RA+ |
| 4 | CD4 T CD45RA− | CD3+ CD4+ CD8− CD45RA− |
| 5 | CD4 T | CD3+ CD4+ CD8− |
| 6 | CD8 T CD45RA+ | CD3+ CD8+ CD4− CD45RA+ |
| 7 | CD8 T CD45RA− | CD3+ CD8+ CD4− CD45RA− |
| 8 | CD8 T | CD3+ CD8+ CD4− |
| 9 | CD4+CD8+ T | CD3+ CD4+ CD8+ |
| 10 | CD8dim T | CD3+ CD4− CD8 int |
| 11 | DN T (gd-like) | CD3+ CD4− CD8− |
| 12 | B cells | CD3− CD19+ |
| 13 | NK CD56bright | CD3− CD19− CD56+ CD16− |
| 14 | NK CD56dim CD16+ | CD3− CD19− CD16+ HLA-DR− |
| 15 | NK CD56dim CD16int | CD3− CD19− CD56 pos, CD16 int, HLA-DR− |
| 16 | CD16+ HLA-DR+ (monocyte-like) | CD3− CD19− CD56− CD16+ HLA-DR+ |
| 17 | CD14+ HLA-DR+ (monocyte/DC-like) | CD3− CD19− CD14+ HLA-DR pos |
| 18 | DC (HLA-DR+ lin−) | CD3− CD19− CD56− CD16− HLA-DR+ |
| 19 | Basophil-like (HLA-DR− lin−) | CD3− CD19− CD56− CD16− HLA-DR− |

Things the table cannot know: FOXP3 and CD25 also rise on activated conventional CD4 T
cells (hence "Treg or activated"; the panel has no CD127); CD45RA+ T cells are naive *or*
TEMRA (CCR7 and CD27 tell them apart, see the heatmap); "DN T (gd-like)" has no TCRγδ
marker behind it; "basophil-like" has no CD123/FcεRI. Replace the table for another panel
or tissue:

```json
{"high": 0.5, "low": 0.3, "min_range": 0.5,
 "rules": [{"label": "CD4 T", "markers": {"CD3": "+", "CD4": "+", "CD8": "-"}},
           {"label": "B cells", "markers": {"CD19": "+"}}]}
```

(or just the list of rules; missing thresholds take the built-in values).

## Results on the demo data

fcs-prep output of the CytoBatchFlagR dataset (24 files, 6 batches × 2 unstimulated + 2
anti-CD3/CD28 + IL-2 stimulated, 10,212,112 events, 25 markers, 14 clustering markers:
HELIOS, CCR7, CD4, CD56, CD3, CD19, CD8, CD103, CD16, CD14, CD27, CCR6, CD45RA, FOXP3),
default parameters. All 20 metaclusters got a label; mean % of events per sample:

| Label | metaclusters | unstimulated % | stimulated % |
|---|---|---|---|
| B cells | 1 | 20.9 | 15.8 |
| CD4 T CD45RA− | 3 | 18.2 | 25.8 |
| CD4 T CD45RA+ | 1 | 16.7 | 12.3 |
| CD8 T CD45RA+ | 2 | 11.7 | 6.9 |
| NK CD56dim CD16int | 1 | 10.9 | 3.5 |
| DN T (gd-like) | 2 | 7.8 | 1.1 |
| CD8 T CD45RA− | 2 | 4.6 | 6.6 |
| CD14+ HLA-DR+ (monocyte/DC-like) | 2 | 2.4 | 1.2 |
| CD16+ HLA-DR+ (monocyte-like) | 1 | 2.0 | 0.0 |
| NK CD56dim CD16+ | 1 | 1.8 | 19.3 |
| CD4 T CD25+FOXP3+ (Treg or activated) | 1 | 1.2 | 5.1 |
| NK CD56bright | 1 | 1.2 | 2.2 |
| DC (HLA-DR+ lin−) | 1 | 0.4 | 0.1 |
| T/B doublet-like | 1 | < 0.05 | < 0.05 |

The stimulated samples look activated: the CD25+FOXP3+ CD4 metacluster grows four-fold,
T cells shift from CD45RA+ to CD45RA−, granzyme B rises about 2.5-fold in T-cell
metaclusters (median 0.7–0.8 → 1.8–2.1 arcsinh) and HLA-DR, PD-1 and TIGIT rise modestly;
CD16+ monocyte-like cells disappear from the cultures. The HLA-DR+ CD4+ CD45RA+ lin−
metacluster labelled DC is pDC-like.

## Runtime

`computeTypes: [standard]`, 2 vCPU / 8 GiB, 1 h timeout. FlowSOM is single-threaded.

Measured on the full demo data as fcs-prep wrote it (10.2 M events, 24 files of 25 k–603 k events, 1.3 GB of
float32 FCS) in Docker, **amd64 emulated on an Apple M-series Mac**, 4 GB memory limit:
**313 s** wall time (reading + subsampling 13 s, SOM + metaclustering 26 s, mapping every
event 267 s, outputs 7 s), **peak memory 3.07 GB** (cgroup `memory.peak`). Memory is
bounded by the largest file (each file is read and mapped on its own; FlowSOM's `NewData`
keeps a few copies of it) plus the training aggregate, not by the total number of events.
Native amd64 should be several times faster. The 20,000-events-per-file fixture runs in
about 70 s.

The image is 2.17 GB unpacked (558 MB compressed): R 4.6 and FlowSOM's dependency tree
(ggpubr, igraph, flowCore/cytolib).

## Running it locally

```sh
make fixture   # data/fixture/demo20k/: 20,000 events per demo file, arcsinh, panel + samples (Python: numpy, flowio)
make run       # builds the image (amd64) and runs it over the fixture into data/output/
make test      # R tests on a tiny synthetic fixture (local Rscript with FlowSOM, else inside the image)
```

By hand:

```sh
docker build --platform linux/amd64 -t processor-flowsom .
docker run --rm --platform linux/amd64 -e INPUT_DIR=/data/input -e OUTPUT_DIR=/data/output \
  -v "$PWD/data/fixture/demo20k:/data/input:ro" -v "$PWD/data/output:/data/output" processor-flowsom
```

## Limitations

- The labels are a heuristic on relative expression; a population missing from the data
  (no B cells, say) shifts the scale of its marker. Check `heatmap.png`.
- `N_METACLUSTERS` is fixed, not estimated; choosing somewhat more metaclusters than the
  populations you expect, and merging afterwards, is safer than too few.
- Inputs must already be transformed (fcs-prep); the app does not compensate, transform
  or scale.
- `medians.csv` is a histogram median (see above); `medians_by_sample.csv` is exact.
- The star tree shows the training aggregate, not all events.

## Licences

Our code is Apache-2.0. The image contains GPL and AGPL R packages loaded in-process
(FlowSOM GPL (>= 2), cytolib AGPL-3.0-only), so the image is distributed under those terms;
see `THIRD_PARTY.md`.

## Citations

> Van Gassen S, Callebaut B, Van Helden MJ, Lambrecht BN, Demeester P, Dhaene T, Saeys Y.
> **FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry data.**
> Cytometry A 87(7):636-645 (2015). https://doi.org/10.1002/cyto.a.22625, PMID 25573116

> Quintelier K, Couckuyt A, Emmaneel A, Aerts J, Saeys Y, Van Gassen S.
> **Analyzing high-dimensional cytometry data using FlowSOM.**
> Nat Protoc 16(8):3775-3801 (2021). https://doi.org/10.1038/s41596-021-00550-0, PMID 34172973

> Wilkerson MD, Hayes DN. **ConsensusClusterPlus: a class discovery tool with confidence
> assessments and item tracking.** Bioinformatics 26(12):1572-1573 (2010).
> https://doi.org/10.1093/bioinformatics/btq170, PMID 20427518 (metaclustering)

> Hahne F, LeMeur N, Brinkman RR, Ellis B, Haaland P, Sarkar D, Spidlen J, Strain E,
> Gentleman R. **flowCore: a Bioconductor package for high throughput flow cytometry.**
> BMC Bioinformatics 10:106 (2009). https://doi.org/10.1186/1471-2105-10-106, PMID 19358741

Demo data: Koenig Z, Andorf S, Tilburgs T, Eswar S. **CytoBatchFlagR: Spectral flow
cytometry dataset.** Zenodo (2025). https://doi.org/10.5281/zenodo.15388817 (CC BY 4.0)
