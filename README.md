# ProFold

ProFold is a structure-prediction and structure-screening pipeline for designed
proteins and protein complexes. It separates fast first-pass folding from
multi-model confirmation and preserves reproducible manifests, model metrics,
logs, and resumable Slurm runs.

This repository does not perform solubility, expression, or thermal-stability
prediction. The former Level 3 sequence-screening module has been removed.

**详细中文说明：** [ProFold 使用与输出说明](docs/PROFOLD_GUIDE.md)

> **v38.1 验证状态：** Level 1 ESMFold2/AF3 和 Level 2 的主要批处理流程已完成局部真实调试；跨 Level 1→Level 2→Rosetta→最终交付的完整端到端测试尚未完成。

## Pipeline

| Stage | Purpose | Models |
| --- | --- | --- |
| Level 1 | High-throughput foldback and first-pass complex confidence | ESMFold + AlphaFold 3 |
| Level 2 | Multi-model confirmation for a selected subset | Protenix 2 + Boltz-2 + OpenDDE + ESMFold2 |

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
specific paths are intentionally kept outside the repository. Copy the template
to `config.local.json` and set every executable, wrapper, model, database,
USalign, and Rosetta path for the current machine before running.

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

For a user-facing sequence CSV, use the campaign entry point:

```bash
python scripts/submit_campaign.py \
  --sequence-csv designs.csv \
  --config config.local.json \
  --level level2 \
  --outdir runs/campaign \
  --dry-run
```

Remove `--dry-run` after checking the generated command. The CSV must contain
`design_id` and either `sequence` or `fasta_path`; an optional
`reference_structure_path` is carried into USalign. The JSON configuration can
record `conda_env`, the Slurm `run_script` (normally `submit_level2.sh`), model
wrapper paths, and `usalign_path`.

Level 2 submission now schedules a dependent closeout job after aggregation.
Its default delivery directory is `<run>/review`; override it with
`--delivery-dir`. The closeout copies structures, confidence/full-data files,
and metric tables, then moves raw artifacts, inputs, and logs into a reversible
`_archive/` directory.

## Input and outputs

The required manifest key is `design_id`; `sequence` or `fasta_path` provides
the protein sequence. Level 1 may additionally use AF3 templates or CID
inputs. Level 2 accepts model-specific Protenix JSON, Boltz YAML, and OpenDDE
JSON inputs for complex systems, plus one per-design FASTA for ESMFold2.
Within each shard, each predictor receives a directory batch and is started
once; results are mapped back to design IDs after prediction.
Each shard also records predictor wall time in `artifacts/shard_*/stage_timings.json`;
the aggregate stage collects these into `metrics/stage_timings.tsv` so shard size,
kernel settings, and model-specific throughput can be tuned from real runs.
Set `level2_mode` to `"cid"` to run apo and holo states for every Level 2
model; provide state-specific inputs such as `apo_protenix_input_path` and
`holo_protenix_input_path` when the states differ in ligand or complex content.

Every run contains an input snapshot, shard tasks, model artifacts, per-model
metrics, a status table, a result manifest, and logs. The result status describes
whether the run and parsing completed; it is not a biological pass/fail label.
Inspect pLDDT, pTM, ipTM, PAE/PDE, ranking, and interface metrics separately.

Level 2 also writes `metrics/canonical_metrics.csv` with stable columns for
complex ipTM/pTM/pLDDT, mean PAE, cross-chain interface PAE, IPSAE, per-chain
and per-chain-pair JSON metrics, USalign values, and optional Rosetta filter
values. Model-specific fields remain in `metrics/level2.csv`.

When Level 2 aggregation runs, `geometry_filters.json` is also produced from
the structures on CPU. It records radius of gyration, chain contacts, contact
density, and optional Shrake-Rupley SASA/ΔSASA. These geometry values are
separate from official Rosetta values; configure `rosetta_filter_script` when
Rosetta shape-complementarity or site-specific filters are required.
The Rosetta batch script receives `--input result_manifest.tsv --output
rosetta_filters.json` and must return a JSON list (or `{"rows": [...]}`) with
`design_id` and numeric filter fields. The output is schema-checked before it
is merged into canonical metrics.

For the final handoff, run
`level2/compact_level2.py --run-dir runs/level2 --out-dir runs/level2/review --delivery-only`;
the result contains only structures, confidence/full-data files, `summary.tsv`,
and metric tables. The default mode also copies status and result manifests
for audit review.

For structure-prediction outputs, complete metric extraction before cleanup.
Keep the top-ranked structure and its full confidence data per design, together
with run-level `summary.json`/`ranking.csv`; archive intermediate samples and
model caches rather than treating them as deliverables.

## Validation and release

Run the regression suite with:

```bash
python -m unittest discover -s tests -v
```

The current development version is `0.0.38.1` (`v0.0.38.1`). To build the release archive:

```bash
bash scripts/publish_release.sh LiuSantu123/profold v0.0.38.1
```

The command expects authenticated `gh` and a clean worktree. It creates a
tagged source archive and SHA256 checksums; it does not include weights,
databases, local configuration, or prediction results.
