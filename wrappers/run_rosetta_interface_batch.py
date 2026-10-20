#!/usr/bin/env python3
"""Run Rosetta InterfaceAnalyzer on Level2 structures in one manifest.

The wrapper deliberately consumes the aggregate manifest, so AF3 structures
can be supplied as ``af3_structure_path`` and are scored alongside predictor
structures.  Rosetta output is converted to the stable ProFold JSON contract.
"""
from __future__ import annotations

import argparse, csv, json, os, re, subprocess, tempfile
from pathlib import Path

ROOT = Path("/nvme_af3/softwares/rosetta3.9/rosetta_2018.09/main")
DEFAULT_BIN = ROOT / "source/build/src/release/linux/3.10/64/x86/gcc/4.8/static/InterfaceAnalyzer.static.linuxgccrelease"
DEFAULT_DB = ROOT / "database"

def _chains(path: Path) -> list[str]:
    try:
        import gemmi
        s = gemmi.read_structure(str(path)); return [str(c.name) for c in s[0]]
    except Exception:
        return []

def _pdb(src: Path, dst: Path) -> None:
    import gemmi
    s = gemmi.read_structure(str(src)); s.write_pdb(str(dst))

def _score(path: Path, outdir: Path, binary: Path, database: Path) -> dict[str, object]:
    chains = _chains(path)
    if len(chains) < 2:
        return {"status": "partial", "error": "fewer than two chains"}
    interface = f"{chains[0]}_{'_'.join(chains[1:])}"
    pdb = outdir / (path.stem + ".pdb"); score = outdir / (path.stem + ".sc")
    _pdb(path, pdb)
    cmd = [str(binary), "-database", str(database), "-s", str(pdb), "-interface", interface,
           "-pack_input", "false", "-out:file:scorefile", str(score), "-overwrite"]
    run = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    if run.returncode:
        return {"status": "failed", "error": (run.stderr or run.stdout)[-1200:]}
    lines = score.read_text(errors="replace").splitlines() if score.is_file() else []
    header = next((line.split() for line in lines if line.startswith("SEQUENCE:") is False and line.startswith("SCORE:") and "description" in line), [])
    value = next((line.split() for line in lines if line.startswith("SCORE:") and "description" not in line), [])
    if not header or not value: return {"status":"failed", "error":"Rosetta scorefile has no SCORE row"}
    data = dict(zip(header[1:], value[1:]))
    def num(*names):
        for name in names:
            try: return float(data[name])
            except (KeyError, TypeError, ValueError): pass
        return ""
    return {"status":"success", "interface":interface,
            "rosetta_interface_complementarity":num("sc_value", "shape_complementarity"),
            "rosetta_sasa":num("complexed_sasa", "sasa_total"),
            "rosetta_delta_sasa":num("dSASA_int", "dSASA_int total"),
            "rosetta_radius_of_gyration":num("rg"),
            "rosetta_total_score":num("total_score"),
            "rosetta_dG_separated":num("dG_separated")}

def main() -> int:
    ap=argparse.ArgumentParser(); ap.add_argument("--input", required=True); ap.add_argument("--output", required=True)
    ap.add_argument("--binary", default=os.environ.get("ROSETTA_INTERFACE_ANALYZER", str(DEFAULT_BIN)))
    ap.add_argument("--database", default=os.environ.get("ROSETTA_DATABASE", str(DEFAULT_DB)))
    a=ap.parse_args(); manifest=Path(a.input); out=Path(a.output); out.parent.mkdir(parents=True, exist_ok=True)
    rows=[]; work=Path(tempfile.mkdtemp(prefix="profold_rosetta_", dir=str(out.parent)))
    with manifest.open() as h:
        for row in csv.DictReader(h, delimiter="\t"):
            candidates=[(k[0:-len("_structure_path")], Path(v)) for k,v in row.items() if k.endswith("_structure_path") and v]
            for model,path in candidates:
                if not path.is_file(): continue
                item=_score(path,work,Path(a.binary),Path(a.database)); item.update({"design_id":row.get("design_id",""),"model":model,"structure_path":str(path)})
                rows.append(item)
    out.write_text(json.dumps(rows, indent=2, ensure_ascii=False)+"\n")
    return 0
if __name__ == "__main__": raise SystemExit(main())
