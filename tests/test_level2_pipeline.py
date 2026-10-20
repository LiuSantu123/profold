from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path


V37 = Path(__file__).resolve().parents[1]
if str(V37) not in sys.path:
    sys.path.insert(0, str(V37))

from common.pipeline import prepare_run  # noqa: E402
from level2.run_level2 import LEVEL2_MODELS, _default_command, run_shard  # noqa: E402


class Level2PipelineTests(unittest.TestCase):
    def test_level2_declares_esmfold2_and_has_a_default_adapter(self) -> None:
        self.assertEqual(LEVEL2_MODELS, ("protenix", "boltz2", "opendde", "esmfold2"))
        command = _default_command(
            {}, "esmfold2", {"esmfold2_input_dir": "/tmp/in", "esmfold2_outdir": "/tmp/out"}
        )
        self.assertEqual(command[0:2], ["bash", str(V37 / "wrappers/run_esmfold2.sh")])

    def test_fake_run_executes_all_four_level2_models(self) -> None:
        fake_wrapper = V37 / "tests/fixtures/fake_wrappers.py"
        self.assertTrue(fake_wrapper.is_file())
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "input.tsv"
            manifest.write_text("design_id\tsequence\nd1\tMKTAA\n", encoding="utf-8")
            command_prefix = [sys.executable, str(fake_wrapper), "--kind"]
            config = {
                "commands": {
                    "esmfold": command_prefix + [
                        "esmfold",
                        "--input",
                        "{esm_input}",
                        "--outdir",
                        "{esm_outdir}",
                    ],
                    "protenix": command_prefix + ["protenix", "{protenix_input}", "{protenix_outdir}"],
                    "boltz2": command_prefix
                    + ["boltz2", "predict", "{boltz_input}", "--out_dir", "{boltz_outdir}"],
                    "opendde": command_prefix
                    + ["opendde", "-i", "{opendde_input}", "-o", "{opendde_outdir}"],
                    "esmfold2": command_prefix
                    + ["esmfold", "--input", "{esmfold2_input}", "--outdir", "{esmfold2_outdir}"],
                }
            }
            run_dir, _ = prepare_run(manifest, root / "run", level="level2", shard_size=1, config=config)
            rows = run_shard(run_dir, 0)

            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["status"], "success")
            self.assertEqual(rows[0]["esmfold2_status"], "success")
            model_dirs = {path.name for path in (run_dir / "artifacts/shard_00000/models").iterdir()}
            self.assertEqual(model_dirs, set(LEVEL2_MODELS))
            self.assertEqual(
                sorted(path.name for path in (run_dir / "logs").glob("level2_*.log")),
                ["level2_boltz2_batch.log", "level2_esmfold2_batch.log", "level2_opendde_batch.log", "level2_protenix_batch.log"],
            )
            for log in (run_dir / "logs").glob("level2_*_batch.log"):
                self.assertEqual(log.read_text(encoding="utf-8").count("command="), 1)
            monitor = run_dir / "metrics/esmfold2_monitor_shard_00000.csv"
            self.assertTrue(monitor.is_file())
            self.assertIn("d1", monitor.read_text(encoding="utf-8"))

    def test_cid_run_executes_apo_and_holo_for_all_models(self) -> None:
        fake_wrapper = V37 / "tests/fixtures/fake_wrappers.py"
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "input.tsv"
            manifest.write_text("design_id\tsequence\nd1\tMKTAA\n", encoding="utf-8")
            command_prefix = [sys.executable, str(fake_wrapper), "--kind"]
            config = {
                "level2_mode": "cid",
                "commands": {
                    "esmfold2": command_prefix + ["esmfold", "--input", "{esmfold2_input}", "--outdir", "{esmfold2_outdir}"],
                    "protenix": command_prefix + ["protenix", "{protenix_input}", "{protenix_outdir}"],
                    "boltz2": command_prefix + ["boltz2", "predict", "{boltz_input}", "--out_dir", "{boltz_outdir}"],
                    "opendde": command_prefix + ["opendde", "-i", "{opendde_input}", "-o", "{opendde_outdir}"],
                },
            }
            run_dir, _ = prepare_run(manifest, root / "run", level="level2", shard_size=1, config=config)
            rows = run_shard(run_dir, 0)
            self.assertEqual(rows[0]["status"], "success")
            for model in LEVEL2_MODELS:
                for state in ("apo", "holo"):
                    self.assertEqual(rows[0][f"{model}_{state}_status"], "success")
                    self.assertTrue(
                        (run_dir / "artifacts/shard_00000/models" / model / state).exists()
                    )
            self.assertTrue((run_dir / "metrics/esmfold2_monitor_apo_shard_00000.csv").is_file())
            self.assertTrue((run_dir / "metrics/esmfold2_monitor_holo_shard_00000.csv").is_file())


if __name__ == "__main__":
    unittest.main()
