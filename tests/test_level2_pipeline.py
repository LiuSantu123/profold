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
    def test_level2_has_only_the_three_declared_models(self) -> None:
        self.assertEqual(LEVEL2_MODELS, ("protenix", "boltz2", "opendde"))
        with self.assertRaisesRegex(ValueError, "unsupported level2 model"):
            _default_command({}, "esmfold2", {})

    def test_fake_run_does_not_execute_or_emit_esmfold(self) -> None:
        fake_wrapper = V37 / ".tmp_fake_wrappers.py"
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
                }
            }
            run_dir, _ = prepare_run(manifest, root / "run", level="level2", shard_size=1, config=config)
            rows = run_shard(run_dir, 0)

            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["status"], "success")
            self.assertNotIn("esmfold2_status", rows[0])
            model_dirs = {path.name for path in (run_dir / "artifacts/shard_00000/models").iterdir()}
            self.assertEqual(model_dirs, set(LEVEL2_MODELS))
            self.assertFalse((run_dir / "artifacts/shard_00000/models/esmfold2").exists())
            self.assertEqual(
                sorted(path.name for path in (run_dir / "logs").glob("level2_*.log")),
                ["level2_boltz2_d1.log", "level2_opendde_d1.log", "level2_protenix_d1.log"],
            )


if __name__ == "__main__":
    unittest.main()
