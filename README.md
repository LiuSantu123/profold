# ProFold

ProFold is a structure-prediction and structure-screening pipeline for designed
proteins and protein complexes. It separates fast first-pass folding from
multi-model confirmation and preserves reproducible manifests, model metrics,
logs, and resumable Slurm runs.

This repository does not perform solubility, expression, or thermal-stability
prediction. The former Level 3 sequence-screening module has been removed.

## Pipeline

| Stage | Purpose | Models |
| --- | --- | --- |
| Level 1 | High-throughput foldback and first-pass complex confidence | ESMFold + AlphaFold 3 |
| Level 2 | Multi-model confirmation for a selected subset | Protenix 2 + Boltz-2 + OpenDDE |

Each stage accepts a TSV manifest and writes an independent run directory. Use
`common/select_manifest.py` to choose designs between stages; no hidden state is
required from an earlier run.

## Install and configure

```bash
git clone https://github.com/LiuSantu123/profold.git
cd profold
python -m pip install -r requirements.txt
python scripts/configure.py --help
python scripts/check_environment.py --config config.local.json --level level1
```

Copy `config.example.json` or use the environment template in
`environments/README.md`. Model weights, databases, GPU runtimes, and machine-
specific paths are intentionally kept outside the repository.

## Run a stage

Prepare and run Level 1:

```bash
python level1/prepare_level1.py \
  --manifest input.tsv --outdir runs/level1 --shard-size 32 \
  --config config.local.json
python level1/run_level1.py --run-dir runs/level1 --shard-index 0
python level1/aggregate_level1.py --run-dir runs/level1
```

Select successful designs for Level 2:

```bash
python common/select_manifest.py \
  --result-manifest runs/level1/results/result_manifest.tsv \
  --top-n 100 --sort-field af3_ranking_score --require-success \
  --output runs/level1_selected.tsv
```

Then prepare and run Level 2 with `level2/prepare_level2.py`,
`level2/run_level2.py`, and `level2/aggregate_level2.py`, or submit the
corresponding `submit_level*.sh` wrapper to Slurm. Repeating a worker enables
resume for complete results; use `--no-resume` when a fresh prediction is
required.

## Input and outputs

The required manifest key is `design_id`; `sequence` or `fasta_path` provides
the protein sequence. Level 1 may additionally use AF3 templates or CID
inputs. Level 2 accepts model-specific Protenix JSON, Boltz YAML, and OpenDDE
JSON inputs for complex systems.

Every run contains an input snapshot, shard tasks, model artifacts, per-model
metrics, a status table, a result manifest, and logs. The result status describes
whether the run and parsing completed; it is not a biological pass/fail label.
Inspect pLDDT, pTM, ipTM, PAE/PDE, ranking, and interface metrics separately.

For structure-prediction outputs, complete metric extraction before cleanup.
Keep the top-ranked structure and its full confidence data per design, together
with run-level `summary.json`/`ranking.csv`; archive intermediate samples and
model caches rather than treating them as deliverables.

## Validation and release

Run the regression suite with:

```bash
python -m unittest discover -s tests -v
```

The current release is `0.0.37` (`v0.0.37`). To build the release archive:

```bash
bash scripts/publish_release.sh LiuSantu123/profold v0.0.37
```

The command expects authenticated `gh` and a clean worktree. It creates a
tagged source archive and SHA256 checksums; it does not include weights,
databases, local configuration, or prediction results.
