from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path


V37 = Path(__file__).resolve().parents[1]
if str(V37) not in sys.path:
    sys.path.insert(0, str(V37))

from common.model_parsers import parse_af3  # noqa: E402


class Af3ParserTests(unittest.TestCase):
    def test_for_wj_cid_keeps_both_variants_and_legacy_aliases(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            archive = root / "archive"
            confidence = root / "af3_confidence_json"
            archive.mkdir()
            confidence.mkdir()

            for variant, iptm in (("with_lig", 0.94), ("without_lig", 0.91)):
                stem = f"d1_af3_{variant}"
                (archive / f"{stem}.cif").write_text("data_test\n", encoding="utf-8")
                (confidence / f"{stem}_summary_confidences.json").write_text(
                    json.dumps({"iptm": iptm, "ptm": 0.92, "ranking_score": 0.93}),
                    encoding="utf-8",
                )
            (confidence / "d1_af3_with_lig_confidences.json").write_text(
                json.dumps({"plddt": [95.0, 97.0]}),
                encoding="utf-8",
            )

            result = parse_af3(root, "d1", mode="for_wj_cid")

            self.assertEqual(result["af3_status"], "success")
            self.assertEqual(result["af3_with_lig_status"], "success")
            self.assertEqual(result["af3_without_lig_status"], "success")
            self.assertEqual(result["af3_with_lig_iptm"], 0.94)
            self.assertEqual(result["af3_without_lig_iptm"], 0.91)
            self.assertEqual(result["af3_iptm"], 0.94)
            self.assertEqual(result["af3_mean_plddt"], 96)
            self.assertTrue(str(result["af3_with_lig_structure_path"]).endswith("d1_af3_with_lig.cif"))

    def test_for_wj_cid_is_partial_when_one_variant_is_missing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "d1_af3_with_lig.cif").write_text("data_test\n", encoding="utf-8")
            (root / "d1_af3_with_lig_summary_confidences.json").write_text(
                json.dumps({"iptm": 0.8}),
                encoding="utf-8",
            )

            result = parse_af3(root, "d1", mode="for_wj_cid")

            self.assertEqual(result["af3_status"], "partial")
            self.assertEqual(result["af3_with_lig_status"], "success")
            self.assertEqual(result["af3_without_lig_status"], "failed")


if __name__ == "__main__":
    unittest.main()
