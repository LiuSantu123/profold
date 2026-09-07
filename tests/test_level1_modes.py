from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


V37 = Path(__file__).resolve().parents[1]
if str(V37) not in sys.path:
    sys.path.insert(0, str(V37))

from level1.run_level1 import _default_command, _write_af3_fasta  # noqa: E402
from common.model_parsers import parse_esmfold  # noqa: E402


class Level1ModeTests(unittest.TestCase):
    def test_native_esmfold_cli(self) -> None:
        command = _default_command({}, "esmfold", {})
        self.assertEqual(command, ["/xcfhome/yzmeng/miniconda3/envs/zb/bin/esm-fold",
                                   "-i", "{esm_input}", "-o", "{esm_outdir}"])
        command = _default_command({"esmfold_bin": "/custom/esm-fold", "esmfold_num_recycles": 0,
                                    "esmfold_chunk_size": 64, "esmfold_cpu_offload": True}, "esmfold", {})
        self.assertEqual(command, ["/custom/esm-fold", "-i", "{esm_input}", "-o", "{esm_outdir}",
                                   "--num-recycles", "0", "--chunk-size", "64", "--cpu-offload"])

    def test_native_pdb_confidence_and_resume_contract(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "d1.pdb").write_text(
                "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 70.00           N  \n"
                "ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00 90.00           C  \nEND\n"
            )
            result = parse_esmfold(root, "d1", native=True)
            self.assertEqual(result["esmfold_status"], "success")
            self.assertEqual(result["esmfold_mean_plddt"], 80)
            self.assertEqual(result["esmfold_summary_path"], "")
            self.assertNotEqual(parse_esmfold(root, "d2", native=True)["esmfold_status"], "success")
            (root / "d2.cif").write_text("data_legacy_esmfold2\n")
            (root / "esmfold2_summary.csv").write_text("design_id,mean_plddt\nd2,99\n")
            self.assertNotEqual(parse_esmfold(root, "d2", native=True)["esmfold_status"], "success")
            (root / "d1.pdb").write_text("END\n")
            self.assertEqual(parse_esmfold(root, "d1", native=True)["esmfold_status"], "failed")

    def test_for_wj_cid_fasta_is_normalized_to_two_colon_separated_chains(self) -> None:
        rows = [{"design_id": "d1", "sequence": "AAA|BBB", "chain_ids": "X,Y"}]
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "af3_inputs.fasta"
            count = _write_af3_fasta(rows, output, "for_wj_cid")
            self.assertEqual(count, 1)
            self.assertEqual(output.read_text(encoding="utf-8"), ">d1\nAAA:BBB\n")

    def test_for_wj_cid_requires_an_explicit_wrapper_command(self) -> None:
        context = {"af3_mode": "for_wj_cid"}
        self.assertIsNone(_default_command({}, "af3", context))
        command = _default_command(
            {"commands": {"af3": ["wrapper", "{af3_fasta}"]}},
            "af3",
            {"af3_mode": "for_wj_cid", "af3_fasta": "/tmp/in.fasta"},
        )
        self.assertEqual(command, ["wrapper", "/tmp/in.fasta"])


if __name__ == "__main__":
    unittest.main()
