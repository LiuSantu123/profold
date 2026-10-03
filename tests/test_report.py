import csv
import tempfile
import unittest
from pathlib import Path

from report.generate_report import generate_report


class ReportTests(unittest.TestCase):
    def test_report_writes_plots_ranking_and_pdf_without_renderer(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp)
            (run / "metrics").mkdir()
            rows = [
                {"design_id": "d1", "model": "protenix", "status": "success", "plddt": "90", "ptm": ".7", "iptm": ".6", "ranking_score": ".65", "confidence_score": ".8", "complex_plddt": "90", "structure_path": ""},
                {"design_id": "d1", "model": "boltz2", "status": "success", "plddt": "85", "ptm": ".65", "iptm": ".55", "ranking_score": ".6", "confidence_score": ".75", "complex_plddt": "85", "structure_path": ""},
                {"design_id": "d2", "model": "protenix", "status": "success", "plddt": "70", "ptm": ".5", "iptm": ".2", "ranking_score": ".4", "confidence_score": ".5", "complex_plddt": "70", "structure_path": ""},
            ]
            with (run / "metrics" / "level2.csv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
                writer.writeheader()
                writer.writerows(rows)
            pdf = generate_report(run, top_n=1, render_structures=False)
            self.assertTrue(pdf.is_file())
            report = run / "report"
            for name in ("metric_distributions.png", "metric_scatter.png", "model_heatmap.png", "model_bars.png", "ranking.tsv", "report_manifest.json"):
                self.assertTrue((report / name).is_file(), name)
            self.assertIn("d1", (report / "ranking.tsv").read_text())


if __name__ == "__main__":
    unittest.main()
