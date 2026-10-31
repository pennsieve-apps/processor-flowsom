"""
Build an input fixture for processor-flowsom from the raw demo FCS files, the way
the upstream steps (fcs-prep -> peacoqc -> cleanet) would hand them over, until
their real outputs exist:

  <out>/fcs/<sample_id>.fcs        FCS 3.1, float32, marker channels arcsinh(x / 6000),
                                   scatter and Time as acquired, CD28 (BV785-A) dropped
                                   (missing from batches 05-06 upstream)
  <out>/metadata/samples.csv       sample_id, file, batch, condition, n_events, control
  <out>/metadata/panel.csv         channel, marker, type, cofactor, clustering, state

Events are a seeded random subsample of each file (--events, 0 = all events).
The clustering/state flags below are a stand-in for fcs-prep's defaults: lineage
markers cluster, functional markers are state markers; the gating markers CD45
and CD14 (every event is CD45+CD14-) and CD163 do not cluster.

  python -m tests.make_fixture <raw-dir> <out-dir> [--events 20000] [--seed 42] [--sheet tests/demo-samples.csv]

Needs numpy and flowio (requirements-dev.txt).
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import flowio
import numpy as np

COFACTOR = 6000.0
DROP_CHANNELS = {"BV785-A"}  # CD28: absent from batches 05-06
CLUSTERING = {"CD3", "CD4", "CD8", "CD19", "CD56", "CD16", "HLADR", "CD45RA", "CCR7",
              "CD27", "CD25", "FOXP3", "CD95", "CCR6", "CD103"}
STATE = {"CD69", "PD1", "TIGIT", "CD39", "GZMB", "GNLY", "PFN", "HELIOS", "GITR",
         "HLADR", "CD25", "CD95", "CD103"}


def channel_type(name: str, marker: str) -> str:
    if name == "Time":
        return "time"
    if re.match(r"^(FSC|SSC)", name):
        return "scatter"
    if name == "AF-A" or not marker:
        return "other"
    return "marker"


def read_fcs(path: Path):
    fd = flowio.FlowData(str(path))
    n = int(fd.text["par"])
    names = [fd.text[f"p{i}n"] for i in range(1, n + 1)]
    markers = [fd.text.get(f"p{i}s", "") or "" for i in range(1, n + 1)]
    events = np.asarray(fd.as_array(preprocess=False), dtype=np.float32).reshape(-1, n)
    return names, markers, events, fd.text


def write_fcs(path: Path, names, markers, events: np.ndarray, keywords: dict[str, str]):
    with open(path, "wb") as fh:
        flowio.create_fcs(fh, events.astype(np.float32).ravel(), names,
                          opt_channel_names=markers, metadata_dict=keywords)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("raw_dir", type=Path)
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--events", type=int, default=20000, help="events per file (0 = all)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--sheet", type=Path, help="the sample sheet (default: <raw_dir>/samples.csv), e.g. tests/demo-samples.csv")
    args = ap.parse_args()

    sheet = {}
    with open(args.sheet or args.raw_dir / "samples.csv", newline="") as fh:
        for row in csv.DictReader(fh):
            sheet[row["file"]] = row
    fcs_dir = args.out_dir / "fcs"
    meta_dir = args.out_dir / "metadata"
    fcs_dir.mkdir(parents=True, exist_ok=True)
    meta_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    panel = None
    for i, (fname, row) in enumerate(sorted(sheet.items())):
        names, markers, ev, text = read_fcs(args.raw_dir / fname)
        keep = [j for j, n in enumerate(names) if n not in DROP_CHANNELS]
        names = [names[j] for j in keep]
        markers = [markers[j] for j in keep]
        ev = ev[:, keep]
        rng = np.random.default_rng(args.seed + i)
        if args.events and ev.shape[0] > args.events:
            ev = ev[np.sort(rng.choice(ev.shape[0], args.events, replace=False))]
        types = [channel_type(n, m) for n, m in zip(names, markers)]
        for j, t in enumerate(types):
            if t == "marker":
                ev[:, j] = np.arcsinh(ev[:, j] / COFACTOR)
        if panel is None:
            panel = []
            for n, m, t in zip(names, markers, types):
                mk = m or n
                panel.append({
                    "channel": n, "marker": mk, "type": t,
                    "cofactor": f"{COFACTOR:g}" if t == "marker" else "",
                    "clustering": str(t == "marker" and mk in CLUSTERING).lower(),
                    "state": str(t == "marker" and mk in STATE).lower(),
                })
        sid = row["sample_id"]
        write_fcs(fcs_dir / f"{sid}.fcs", names, markers, ev, {
            "$DATE": text.get("date", ""),
            "PENNSIEVE_PREP": f"transform=arcsinh;cofactor={COFACTOR:g}",
            "PENNSIEVE_FIXTURE": f"subsample={args.events};seed={args.seed + i}",
        })
        rows.append({"sample_id": sid, "file": fname, "batch": row["batch"],
                     "condition": row["condition"], "n_events": ev.shape[0],
                     "control": row.get("control", "")})
        print(f"{sid}: {ev.shape[0]} events, {len(names)} channels", flush=True)

    with open(meta_dir / "samples.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(meta_dir / "panel.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(panel[0]))
        w.writeheader()
        w.writerows(panel)


if __name__ == "__main__":
    main()
