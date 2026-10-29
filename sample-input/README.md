# sample-input

processor-flowsom takes the transformed FCS files and metadata that fcs-prep writes (after
QC by peacoqc and cleanet in the workflow), so the real input is an upstream step's output.

## Demo data

The CytoBatchFlagR spectral flow dataset (Koenig Z, Andorf S, Tilburgs T, Eswar S, Zenodo
2025, https://doi.org/10.5281/zenodo.15388817, CC BY 4.0): 24 FCS files, healthy donor
blood leukocytes pre-gated CD45+CD14−, 6 batches × {2 unstimulated, 2 anti-CD3/CD28 + IL-2}.
Download the files listed in the record (or `manifest.txt` of the workflow repo) into one
folder together with a `samples.csv` sample sheet (`file, sample_id, batch, condition`).

Then either run fcs-prep on it and point `INPUT_DIR` at its output (`fcs/`, `metadata/`),
or make a stand-in with the fixture script, which does what fcs-prep would (arcsinh with
cofactor 6000 on marker channels, CD28 dropped because batches 05–06 lack it, `panel.csv`
and `samples.csv`) on a seeded subsample:

```sh
python3 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
.venv/bin/python -m tests.make_fixture <folder-with-the-fcs-files> data/fixture/demo20k --events 20000
# --events 0 keeps every event (1.3 GB)
make run
```

`tests/test_main.R` builds its own tiny synthetic FCS files and needs no download.
