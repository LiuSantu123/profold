# Bundled adapters and model dependencies

`tools/af3/` snapshots the local AF3 CID wrapper, monitor, multichain metrics
and archive modules (2026-09-07). `tools/monitors/` contains the existing
multimodel confidence scanner. These are workflow tools, not AlphaFold3
weights or an AlphaFold3 source distribution.

`wrappers/protenix_100.py` retains its ByteDance copyright and Apache-2.0
header; the license is in `wrappers/LICENSE.protenix.txt`. This is the existing
local inference adapter derived from Protenix, not an assertion of
compatibility with every newer Protenix release.

NetSolP extraction, TemBERTure scoring and OpenDDE conversion/launcher scripts
come from existing local workflow adapters. Model packages, datasets,
pretrained parameters and third-party executables are not bundled. Install
them separately under their upstream terms. No new blanket license is assigned
to third-party software by this release.

The legacy multimodel scanner and AF3 multichain_v2 have different historical
metric conventions; their ipSAE thresholds are not interchangeable.
