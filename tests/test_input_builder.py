from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path


V37 = Path(__file__).resolve().parents[1]
if str(V37) not in sys.path:
    sys.path.insert(0, str(V37))

from common.input_builder import (  # noqa: E402
    build_af3_json,
    build_boltz_yaml,
    build_protenix_json,
)


class InputBuilderTests(unittest.TestCase):
    def _template(self, path: Path) -> None:
        path.write_text(
            json.dumps(
                {
                    "name": "template",
                    "sequences": [
                        {
                            "protein": {
                                "id": ["A"],
                                "sequence": "OLD_A",
                                "templates": [{"old": True}],
                                "unpairedMsa": "old.a3m",
                                "pairedMsa": "old.paired.a3m",
                            }
                        },
                        {
                            "protein": {
                                "id": ["B"],
                                "sequence": "OLD_B",
                                "templates": [{"old": True}],
                                "unpairedMsa": "old.a3m",
                            }
                        },
                        {"ligand": {"id": "C", "ccdCodes": ["ATP"]}},
                    ],
                }
            ),
            encoding="utf-8",
        )

    def test_af3_multichain_replacement_clears_msa_and_keeps_ligand(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "template.json"
            output = root / "out.json"
            self._template(source)
            row = {
                "design_id": "d1",
                "sequence": "NEW_A:NEW_B",
                "chain_ids": "A,B",
                "af3_json_path": str(source),
                "metadata_json": json.dumps({"clear_msa": True}),
            }
            build_af3_json(row, output, seed=7)
            payload = json.loads(output.read_text(encoding="utf-8"))
            proteins = [item["protein"] for item in payload["sequences"] if "protein" in item]
            self.assertEqual([protein["sequence"] for protein in proteins], ["NEW_A", "NEW_B"])
            self.assertEqual(proteins[0]["templates"], [])
            self.assertEqual(proteins[0]["unpairedMsa"], [])
            self.assertEqual(proteins[0]["pairedMsa"], [])
            self.assertEqual(payload["sequences"][2]["ligand"]["ccdCodes"], ["ATP"])

    def test_af3_missing_chain_is_rejected_instead_of_leaving_template_sequence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "template.json"
            self._template(source)
            row = {
                "design_id": "d1",
                "sequence": "NEW_A:NEW_B",
                "chain_ids": "A,B",
                "af3_json_path": str(source),
                "metadata_json": json.dumps({"af3_sequences": {"A": "ONLY_A"}}),
            }
            with self.assertRaisesRegex(ValueError, "missing chain id"):
                build_af3_json(row, root / "out.json", seed=7)

    def test_af3_sequence_list_count_is_strict(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "template.json"
            self._template(source)
            row = {
                "design_id": "d1",
                "sequence": "NEW_A:NEW_B",
                "chain_ids": "A,B",
                "af3_json_path": str(source),
                "metadata_json": json.dumps({"af3_sequences": ["ONLY_A"]}),
            }
            with self.assertRaisesRegex(ValueError, "count does not match"):
                build_af3_json(row, root / "out.json", seed=7)

    def test_explicit_protenix_msa_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "input.json"
            source.write_text(
                json.dumps(
                    [
                        {
                            "name": "d1",
                            "sequences": [
                                {"proteinChain": {"id": ["A"], "sequence": "AAA", "msa": "x.a3m"}}
                            ],
                        }
                    ]
                ),
                encoding="utf-8",
            )
            row = {"design_id": "d1", "protenix_input_path": str(source), "sequence": "AAA"}
            with self.assertRaisesRegex(ValueError, "non-empty MSA"):
                build_protenix_json(row, root / "copy.json")

    def test_explicit_boltz_msa_is_rejected_but_empty_is_allowed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "input.yaml"
            source.write_text(
                "version: 1\nsequences:\n  - protein:\n      id: [A]\n      sequence: AAA\n      msa: x.a3m\n",
                encoding="utf-8",
            )
            row = {"design_id": "d1", "boltz_input_path": str(source), "sequence": "AAA"}
            with self.assertRaisesRegex(ValueError, "non-empty MSA"):
                build_boltz_yaml(row, root / "copy.yaml")

            source.write_text(
                "version: 1\nsequences:\n  - protein:\n      id: [A]\n      sequence: AAA\n      msa: empty\n",
                encoding="utf-8",
            )
            output = Path(build_boltz_yaml(row, root / "copy_empty.yaml"))
            self.assertTrue(output.is_file())

    def test_generated_level2_inputs_are_no_msa(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            row = {"design_id": "d1", "sequence": "AAA:BBB", "chain_ids": "A,B"}
            protenix = Path(build_protenix_json(row, root / "d1.json"))
            boltz = Path(build_boltz_yaml(row, root / "d1.yaml"))
            self.assertNotIn("msa", protenix.read_text(encoding="utf-8").lower())
            self.assertIn("msa: empty", boltz.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
