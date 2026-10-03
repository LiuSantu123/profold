import csv
import json
import tempfile
import unittest
from pathlib import Path

from scripts.profold import _canonical_rows, _fingerprint, _summary_rows, _write_summary


class FacadeTests(unittest.TestCase):
    def test_compact_alias_and_fixed_summary_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            template = root / "template.json"
            template.write_text(json.dumps({"sequences": []}), encoding="utf-8")
            source = root / "designs.tsv"
            source.write_text(
                "design_id\tsequence\tmode\ttemplate\n"
                "d1\tMKT\tprotein\ttemplate.json\n",
                encoding="utf-8",
            )
            row = _canonical_rows(source)[0]
            self.assertEqual(row["mode"], "protein")
            self.assertEqual(row["target_spec_path"], str(template.resolve()))
            output = root / "summary.tsv"
            _write_summary(output, [{"design_id": "d1", "stage": "level1", "status": "success"}])
            with output.open(encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                self.assertEqual(reader.fieldnames, [
                    "design_id", "stage", "status", "structure_path", "plddt", "ptm",
                    "iptm", "interface_pae", "model_consensus", "selection", "error",
                ])


if __name__ == "__main__":
    unittest.main()
