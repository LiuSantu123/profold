# 0.0.38.2

ProFold release for designed-protein structure prediction and screening.

- Added the compact v38.2 facade with one normalized design manifest, machine
  config, resumable execution, and a stable `summary.tsv` output.
- Added Level 2 aggregate PDF reporting for large campaigns, including metric
  distributions, quantiles, model summaries, and Top-N structure figures.
- Added VMD/Tachyon multi-view rendering with chain colors, transparent
  surface, nucleic-acid/metal representations, ligand sticks, and PPI views.
- Added f101 Protenix wrapper and regression coverage for the facade/report.

- Level 1: ESMFold and AlphaFold 3 first-pass structure prediction.
- Level 2: Protenix 2, Boltz-2, and OpenDDE multi-model confirmation.
- Removed the former Level 3 solubility and thermal-stability pipeline,
  including its runners, wrappers, parsers, fixtures, and environment entries.
- Updated documentation and release tooling for the `profold` project name.
- Preserved manifest-driven execution, Slurm sharding, resume behavior, and
  structure-prediction metric aggregation.

Validation covers the repository regression suite and configuration/shell
checks. Model weights, databases, machine-specific paths, and prediction
outputs are not included in the source release.
