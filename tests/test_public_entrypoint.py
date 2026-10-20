import csv
import json
import tempfile
import unittest
from pathlib import Path

from common.metric_schema import CANONICAL_FIELDS, canonical_metric_row
from common.sequence_csv import read_sequence_csv, write_manifest_from_csv
from common.confidence_metrics import extract_confidence_metrics
from common.model_parsers import parse_summary_model
from common.predictor_registry import PREDICTORS, get_predictor
from level2.structure_filters import structure_metrics
from level2.compact_level2 import compact


class PublicEntrypointTests(unittest.TestCase):
    def test_predictor_registry_has_batch_contracts(self):
        self.assertEqual(tuple(PREDICTORS), ("protenix", "boltz2", "opendde", "esmfold2"))
        self.assertEqual(get_predictor("boltz2").input_suffix, ".yaml")

    def test_geometry_filters_compute_radius_and_contacts(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "complex.pdb"
            path.write_text(
                "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n"
                "ATOM      2  CA  ALA B   1       3.000   0.000   0.000  1.00 20.00           C  \n"
                "END\n", encoding="utf-8"
            )
            metrics = structure_metrics(path)
            self.assertEqual(metrics["chain_count"], 2)
            self.assertGreater(metrics["radius_of_gyration"], 0)
            self.assertIn("A:B", metrics["interface_pairs"])
    def test_sequence_csv_normalizes_paths_and_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "designs.csv"
            source.write_text(
                "design_id,sequence,reference_structure_path\n"
                "d1,AAA,ref/model.pdb\n", encoding="utf-8"
            )
            rows = read_sequence_csv(source)
            self.assertEqual(rows[0]["sequence"], "AAA")
            self.assertEqual(rows[0]["reference_structure_path"], str((root / "ref/model.pdb").resolve()))
            manifest = write_manifest_from_csv(source, root / "run/input.tsv")
            self.assertTrue(manifest.is_file())
            self.assertIn("reference_structure_path", manifest.read_text())

    def test_canonical_metric_row_keeps_predictor_values(self):
        row = canonical_metric_row(
            {
                "design_id": "d1",
                "protenix_status": "success",
                "protenix_iptm": 0.8,
                "protenix_ptm": 0.7,
                "protenix_mean_plddt": 82,
                "protenix_structure_path": "/tmp/d1.cif",
            },
            model="protenix",
        )
        self.assertEqual(row["complex_iptm"], 0.8)
        self.assertEqual(row["complex_mean_plddt"], 82)
        self.assertEqual(row["structure_path"], "/tmp/d1.cif")

    def test_confidence_metrics_reduce_complex_and_chain_values(self):
        metrics = extract_confidence_metrics({
            "plddt": [80, 90, 70, 60],
            "token_pair_pae": [[1, 2, 8, 9], [2, 1, 7, 8], [8, 7, 1, 2], [9, 8, 2, 1]],
            "token_asym_id": ["A", "A", "B", "B"],
        })
        self.assertEqual(metrics["mean_plddt"], 75.0)
        self.assertEqual(metrics["mean_interface_pae"], 8.0)
        self.assertEqual(metrics["chain_metrics"]["A"]["length"], 2)
        self.assertIn("A:B", metrics["chain_pair_metrics"])
        self.assertGreater(metrics["chain_pair_metrics"]["A:B"]["ipsae"], 0)

    def test_summary_parser_attaches_full_data_and_chain_ptm(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "d1_summary_confidence_sample_0.json").write_text(json.dumps({
                "ptm": 0.7, "iptm": 0.8, "chain_ptm": [0.6, 0.5],
            }), encoding="utf-8")
            (root / "d1_full_data_sample_0.json").write_text(json.dumps({
                "token_pair_pae": [[1, 2, 8, 9], [2, 1, 7, 8], [8, 7, 1, 2], [9, 8, 2, 1]],
                "token_asym_id": ["A", "A", "B", "B"], "plddt": [80, 90, 70, 60],
            }), encoding="utf-8")
            (root / "d1_sample_0.cif").write_text("data_d1\n", encoding="utf-8")
            parsed = parse_summary_model(root, "d1", "protenix")
            self.assertEqual(parsed["protenix_iptm"], 0.8)
            self.assertEqual(parsed["protenix_mean_interface_pae"], 8.0)
            self.assertEqual(parsed["protenix_chain_metrics"]["A"]["ptm"], 0.6)

    def test_delivery_rewrites_paths_and_keeps_only_selected_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run = root / "run"
            (run / "results").mkdir(parents=True)
            (run / "metrics").mkdir(parents=True)
            raw = run / "raw"
            raw.mkdir()
            structure = raw / "d1_model_0.cif"
            full = raw / "d1_full_data.json"
            structure.write_text("data_d1\n", encoding="utf-8")
            full.write_text("{}\n", encoding="utf-8")
            with (run / "results/result_manifest.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["design_id", "protenix_structure_path", "protenix_full_data_json"], delimiter="\t")
                writer.writeheader()
                writer.writerow({"design_id": "d1", "protenix_structure_path": str(structure), "protenix_full_data_json": str(full)})
            row = {field: "" for field in CANONICAL_FIELDS}
            row.update({"design_id": "d1", "model": "protenix", "status": "success", "structure_path": str(structure), "full_data_path": str(full)})
            with (run / "metrics/canonical_metrics.csv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(row))
                writer.writeheader()
                writer.writerow(row)
            (run / "summary.json").write_text(json.dumps({"rows": [row]}), encoding="utf-8")
            destination = root / "review"
            compact(run, destination, delivery_only=True)
            with (destination / "summary.tsv").open(encoding="utf-8", newline="") as handle:
                summary = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(summary["structure_path"], "structures/d1_protenix.cif")
            self.assertEqual(summary["full_data_path"], "confidences/d1_protenix_full_data.json")
            self.assertFalse(Path(summary["structure_path"]).is_absolute())
            self.assertTrue((destination / "summary.json").is_file())
