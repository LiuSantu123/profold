#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
multimodel_monitor_v1.py
Protenix / Boltz-2 / OpenDDE 统一置信度监控脚本（对齐 multichain_af3_metrics_v1.py 列口径，
便于后续与 AF3 做多模型一致性验证）。

每个设计输出一行（与 AF3 脚本一致的列名）：
  全局      : iptm / ptm / ranking_score / mean_plddt(统一 0-100) / has_clash / disorder
  软件特有  : boltz2 -> confidence_score / ligand_iptm / protein_iptm / complex_plddt /
              complex_iplddt / complex_pde / complex_ipde / affinity_pred_value /
              affinity_probability_binary（预测小分子时）
  每链      : {c}_mean_plddt / {c}_iptm / {c}_ptm / {c}_len_tokens
  每两两链  : {X}-{Y}_iptm / {X}-{Y}_ipae(min) / {X}-{Y}_ipae_mean / {X}-{Y}_ipsae_max/min/mean /
              {X}-{Y}_avg_if_pae
  整体复合物 : overall_ipae / overall_ipsae_max/min/mean / overall_avg_if_pae

各软件输出约定：
  Protenix : <name>_summary_confidence_sample_<N>.json（默认有）
             <name>_full_data_sample_<N>.json（需配置 need_atom_confidence=True 才有，
             configs/configs_inference.py 默认 False）—— 含 token_pair_pae + token_asym_id，
             无此文件时 ipae/ipsae 标 N/A（仅 summary 级 iptm/ptm/plddt）。
  Boltz-2  : predictions/<id>/confidence_<id>_model_<rank>.json + pae_<id>_model_<rank>.npz +
             plddt_<id>_model_<rank>.npz + <id>_model_<rank>.cif + affinity_<id>.json（可选）。
             链划分从 CIF 解析（pae/plddt 为 token 级，与 CIF 链序对齐）。
  OpenDDE  : <name>_summary_confidence_sample_<N>.json + <name>_full_data_sample_<N>.json（默认有）。

用法：
  python3 multimodel_monitor_v1.py --software protenix --scan <root> -o protenix.csv
  python3 multimodel_monitor_v1.py --software boltz2   --scan <root> -o boltz2.csv
  python3 multimodel_monitor_v1.py --software opendde  --scan <root> -o opendde.csv
"""
import os, sys, json, glob, csv, argparse, itertools
import numpy as np

try:
    from Bio.PDB import MMCIFParser
    _HAS_BIO = True
except Exception:
    _HAS_BIO = False


# ---------------- 通用工具 ----------------
def _num(x, nd=4):
    try:
        if x is None:
            return "N/A"
        return round(float(x), nd)
    except (TypeError, ValueError):
        return "N/A"


def _plddt100(x):
    """统一 plddt 到 0-100（Boltz/OpenDDE atom 级为 0-1，AF3/Protenix 为 0-100）。"""
    try:
        if x is None:
            return "N/A"
        v = float(x)
        return round(v * 100.0 if v <= 1.1 else v, 3)
    except (TypeError, ValueError):
        return "N/A"


def _avg_pair(mat, i, j):
    try:
        n = len(mat)
        if i >= n or j >= n:
            return "N/A"
        a, b = mat[i][j], mat[j][i]
        vals = [v for v in (a, b) if isinstance(v, (int, float))]
        if not vals:
            return "N/A"
        return round(float(np.mean(vals)), 4)
    except (TypeError, IndexError):
        return "N/A"


def calc_ipsae_from_indices(pae, idx_a, idx_b, cutoff=15.0):
    """界面 PAE 指标（与 AF3 脚本同公式）。返回 (max_sae, min_sae, mean_sae, avg_if_pae)。"""
    if pae is None or not idx_a or not idx_b:
        return "N/A", "N/A", "N/A", "N/A"
    max_dim = pae.shape[0]
    va = [i for i in idx_a if i < max_dim]
    vb = [i for i in idx_b if i < max_dim]
    if not va or not vb:
        return "N/A", "N/A", "N/A", "N/A"
    sub = pae[np.ix_(va, vb)]
    valid = sub < cutoff
    if not np.any(valid):
        return 0.0, 99.0, 0.0, 99.0
    L = np.sum(valid)
    d0 = 1.24 * (L ** (1/3)) - 1.8 if L >= 27 else 1.0
    pae_if = sub[valid]
    scores = 1.0 / (1.0 + (pae_if / d0) ** 2)
    return (round(float(np.max(scores)), 4), round(float(np.min(scores)), 4),
            round(float(np.mean(scores)), 4), round(float(np.mean(pae_if)), 3))


def pair_metrics_from_pae(pae, chain_labels, tok_idx, cutoff):
    """由 PAE 矩阵 + 链 token 索引计算两两/整体 ipae、ipsae。写入 dict m。"""
    m = {}
    pair_ipae, pair_ipae_mean = [], []
    for id1, id2 in itertools.combinations(chain_labels, 2):
        key = f"{id1}-{id2}"
        i1, i2 = tok_idx.get(id1, []), tok_idx.get(id2, [])
        if pae is not None and i1 and i2:
            sub = pae[np.ix_(i1, i2)]
            m[f"{key}_ipae"] = _num(float(np.min(sub)))
            m[f"{key}_ipae_mean"] = _num(float(np.mean(sub)))
            mx, mn, mean_s, avg = calc_ipsae_from_indices(pae, i1, i2, cutoff)
            m[f"{key}_ipsae_max"], m[f"{key}_ipsae_min"] = mx, mn
            m[f"{key}_ipsae_mean"], m[f"{key}_avg_if_pae"] = mean_s, avg
            pair_ipae.append(float(np.min(sub)))
            pair_ipae_mean.append(float(np.mean(sub)))
        else:
            for k in (f"{key}_ipae", f"{key}_ipae_mean", f"{key}_ipsae_max",
                      f"{key}_ipsae_min", f"{key}_ipsae_mean", f"{key}_avg_if_pae"):
                m[k] = "N/A"
    # 整体
    if pae is not None and len(chain_labels) >= 2:
        n = pae.shape[0]
        mask = np.ones((n, n), dtype=bool)
        for c in chain_labels:
            if tok_idx[c]:
                mask[np.ix_(tok_idx[c], tok_idx[c])] = False
        np.fill_diagonal(mask, False)
        cross = pae[mask]
        valid = cross < cutoff
        if np.any(valid):
            L = int(np.sum(valid))
            d0 = 1.24 * (L ** (1/3)) - 1.8 if L >= 27 else 1.0
            pae_if = cross[valid]
            sc = 1.0 / (1.0 + (pae_if / d0) ** 2)
            m["overall_ipsae_max"] = round(float(np.max(sc)), 4)
            m["overall_ipsae_min"] = round(float(np.min(sc)), 4)
            m["overall_ipsae_mean"] = round(float(np.mean(sc)), 4)
            m["overall_avg_if_pae"] = round(float(np.mean(pae_if)), 3)
        else:
            m.update({"overall_ipsae_max": 0.0, "overall_ipsae_min": 99.0,
                      "overall_ipsae_mean": 0.0, "overall_avg_if_pae": 99.0})
    else:
        for k in ("overall_ipsae_max", "overall_ipsae_min", "overall_ipsae_mean", "overall_avg_if_pae"):
            m[k] = "N/A"
    m["overall_ipae"] = _num(np.mean(pair_ipae)) if pair_ipae else "N/A"
    return m


# ---------------- Protenix ----------------
def parse_protenix(summary_path):
    p = {}
    p["base"] = os.path.basename(summary_path).replace("_summary_confidence_sample_", "__SPLIT__")
    name, rank = p["base"].split("__SPLIT__")
    p["name"] = name
    p["rank"] = rank.replace(".json", "")
    p["model_name"] = f"{os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(summary_path))))}/{name}_sample_{p['rank']}"
    d = json.load(open(summary_path, encoding="utf-8"))
    p["iptm"], p["ptm"] = d.get("iptm"), d.get("ptm")
    p["ranking_score"] = d.get("ranking_score")
    p["plddt"] = d.get("plddt")          # 0-100
    p["has_clash"], p["disorder"] = d.get("has_clash"), d.get("disorder")
    p["chain_ptm"] = d.get("chain_ptm") or []
    p["chain_iptm"] = d.get("chain_iptm") or []
    p["chain_plddt"] = d.get("chain_plddt") or []   # 0-100
    p["chain_pair_iptm"] = d.get("chain_pair_iptm") or []
    # full data（need_atom_confidence=True 才有）
    full = summary_path.replace(f"_summary_confidence_sample_{p['rank']}.json",
                                f"_full_data_sample_{p['rank']}.json")
    p["full_path"] = full if os.path.exists(full) else None
    return p


def load_full_data_protenix(p, chains):
    if not p["full_path"]:
        return None, {}, None, None
    d = json.load(open(p["full_path"], encoding="utf-8"))
    pae = np.asarray(d.get("token_pair_pae"), dtype=float) if d.get("token_pair_pae") else None
    asym = d.get("token_asym_id") or []
    chain_labels = []
    for a in asym:
        c = chr(ord('A') + int(a)) if isinstance(a, (int, float)) else str(a)
        if c not in chain_labels:
            chain_labels.append(c)
    if chains:
        keep = [c for c in chains if c in chain_labels]
        if keep:
            chain_labels = keep
    tok_idx = {c: [] for c in chain_labels}
    for a, i in zip(asym, range(len(asym))):
        c = chr(ord('A') + int(a)) if isinstance(a, (int, float)) else str(a)
        if c in tok_idx:
            tok_idx[c].append(i)
    atom_plddt = d.get("atom_plddt") or None  # 0-1 (原子级)
    # atom -> 链：通过 atom_to_token_idx 映射到 token 级链
    atom_chain_ids = None
    if atom_plddt and d.get("atom_to_token_idx"):
        a2t = d["atom_to_token_idx"]
        tok_chain = {}
        for t_i, a in enumerate(asym):
            c = chr(ord('A') + int(a)) if isinstance(a, (int, float)) else str(a)
            tok_chain[t_i] = c
        atom_chain_ids = [tok_chain.get(int(t), "?") for t in a2t]
    return pae, tok_idx, atom_plddt, atom_chain_ids


# ---------------- Boltz-2 ----------------
def parse_boltz2(conf_path):
    p = {}
    base = os.path.basename(conf_path).replace("confidence_", "").replace("_model_", "__MODEL__")
    name, rank = base.split("__MODEL__")
    p["name"] = name
    p["rank"] = rank.replace(".json", "")
    p["model_name"] = f"{os.path.basename(os.path.dirname(conf_path))}/{name}_model_{p['rank']}"
    d = json.load(open(conf_path, encoding="utf-8"))
    for k in ("confidence_score", "iptm", "ptm", "ligand_iptm", "protein_iptm",
              "complex_plddt", "complex_iplddt", "complex_pde", "complex_ipde"):
        p[k] = d.get(k)
    p["ranking_score"] = None   # boltz 无 ranking_score
    p["has_clash"], p["disorder"] = None, None
    chains_ptm = d.get("chains_ptm") or {}
    p["chains_ptm"] = {int(k): v for k, v in chains_ptm.items()}
    pair = d.get("pair_chains_iptm") or {}
    p["pair_chains_iptm"] = {int(k): {int(k2): val for k2, val in v.items()} for k, v in pair.items()}
    # 链标签（cif 顺序）
    cif = os.path.join(os.path.dirname(conf_path), f"{name}_model_{p['rank']}.cif")
    p["cif_path"] = cif if os.path.exists(cif) else None
    p["pae_path"] = os.path.join(os.path.dirname(conf_path), f"pae_{name}_model_{p['rank']}.npz")
    if not os.path.exists(p["pae_path"]):
        p["pae_path"] = None
    p["plddt_path"] = os.path.join(os.path.dirname(conf_path), f"plddt_{name}_model_{p['rank']}.npz")
    if not os.path.exists(p["plddt_path"]):
        p["plddt_path"] = None
    # affinity（小分子亲和力，可缺）
    aff = os.path.join(os.path.dirname(conf_path), f"affinity_{name}.json")
    p["affinity"] = None
    if os.path.exists(aff):
        a = json.load(open(aff, encoding="utf-8"))
        p["affinity"] = a
    return p


def load_full_data_boltz2(p, chains):
    """从 CIF 构建链 token 索引；加载 pae/plddt npz。"""
    pae = None
    if p["pae_path"]:
        try:
            pae = np.load(p["pae_path"])["pae"].astype(float)
        except Exception:
            pae = None
    plddt = None
    if p["plddt_path"]:
        try:
            plddt = np.load(p["plddt_path"])["plddt"].astype(float)
        except Exception:
            plddt = None
    # 链划分：CIF 链序（token 级，与 pae/plddt 对齐）
    chain_labels, tok_idx = [], {}
    if _HAS_BIO and p["cif_path"] and os.path.exists(p["cif_path"]):
        try:
            s = MMCIFParser(QUIET=True).get_structure("x", p["cif_path"])
            model = s[0]
            labels_in_cif = [c.id for c in model]
            if chains:
                keep = [c for c in chains if c in labels_in_cif]
                if keep:
                    labels_in_cif = keep
            off = 0
            for c in labels_in_cif:
                chain = model[c]
                nres = len([r for r in chain])
                chain_labels.append(c)
                tok_idx[c] = list(range(off, off + nres))
                off += nres
        except Exception:
            pass
    return pae, tok_idx, plddt, None   # boltz plddt 为 token 级，无 atom_chain_ids


# ---------------- OpenDDE ----------------
def parse_opendde(summary_path):
    p = {}
    base = os.path.basename(summary_path).replace("_summary_confidence_sample_", "__SPLIT__")
    name, rank = base.split("__SPLIT__")
    p["name"] = name
    p["rank"] = rank.replace(".json", "")
    p["model_name"] = f"{os.path.basename(os.path.dirname(os.path.dirname(summary_path)))}/{name}_sample_{p['rank']}"
    d = json.load(open(summary_path, encoding="utf-8"))
    p["iptm"], p["ptm"] = d.get("iptm"), d.get("ptm")
    p["ranking_score"] = d.get("ranking_score")
    p["plddt"] = d.get("plddt")           # 0-100（opendde 全局已转 0-100）
    p["has_clash"], p["disorder"] = d.get("has_clash"), d.get("disorder")
    p["chain_ptm"] = d.get("chain_ptm") or []
    p["chain_iptm"] = d.get("chain_iptm") or []
    p["chain_plddt"] = d.get("chain_plddt") or []  # 0-1（opendde chain 级是 0-1）
    p["chain_pair_iptm"] = d.get("chain_pair_iptm") or []
    p["chain_pair_gpde"] = d.get("chain_pair_gpde") or []  # 界面 PAE 类
    full = summary_path.replace(f"_summary_confidence_sample_{p['rank']}.json",
                                f"_full_data_sample_{p['rank']}.json")
    p["full_path"] = full if os.path.exists(full) else None
    return p


def load_full_data_opendde(p, chains):
    if not p["full_path"]:
        return None, {}, None, None
    d = json.load(open(p["full_path"], encoding="utf-8"))
    pae = np.asarray(d.get("token_pair_pae"), dtype=float) if d.get("token_pair_pae") else None
    asym = d.get("token_asym_id") or []
    chain_labels = []
    for a in asym:
        c = chr(ord('A') + int(a)) if isinstance(a, (int, float)) else str(a)
        if c not in chain_labels:
            chain_labels.append(c)
    if chains:
        keep = [c for c in chains if c in chain_labels]
        if keep:
            chain_labels = keep
    tok_idx = {c: [] for c in chain_labels}
    for a, i in zip(asym, range(len(asym))):
        c = chr(ord('A') + int(a)) if isinstance(a, (int, float)) else str(a)
        if c in tok_idx:
            tok_idx[c].append(i)
    atom_plddt = d.get("atom_plddt") or None  # 0-1 (原子级)
    atom_chain_ids = None
    if atom_plddt and d.get("atom_to_token_idx"):
        a2t = d["atom_to_token_idx"]
        tok_chain = {}
        for t_i, a in enumerate(asym):
            c = chr(ord('A') + int(a)) if isinstance(a, (int, float)) else str(a)
            tok_chain[t_i] = c
        atom_chain_ids = [tok_chain.get(int(t), "?") for t in a2t]
    return pae, tok_idx, atom_plddt, atom_chain_ids


# ---------------- 统一指标计算 ----------------
def compute_metrics(p, software, chains, pae_cutoff):
    r = {"model_name": p["model_name"], "software": software,
         "status": "success", "error": ""}
    # 全局
    r["iptm"] = _num(p.get("iptm"))
    r["ptm"] = _num(p.get("ptm"))
    r["ranking_score"] = _num(p.get("ranking_score"))
    r["has_clash"] = p.get("has_clash")
    r["disorder"] = _num(p.get("disorder"))
    # 软件特有
    if software == "boltz2":
        for k in ("confidence_score", "ligand_iptm", "protein_iptm",
                  "complex_iplddt", "complex_pde", "complex_ipde"):
            r[k] = _num(p.get(k))
        r["complex_plddt"] = _plddt100(p.get("complex_plddt"))
        if p.get("affinity"):
            for k, v in p["affinity"].items():
                r[f"affinity_{k}"] = _num(v)
    # 加载 full（pae / 链 token 索引 / 原子或 token 级 plddt）
    if software == "protenix":
        pae, tok_idx, atom_plddt, atom_chain_ids = load_full_data_protenix(p, chains)
    elif software == "boltz2":
        pae, tok_idx, atom_plddt, atom_chain_ids = load_full_data_boltz2(p, chains)
    else:
        pae, tok_idx, atom_plddt, atom_chain_ids = load_full_data_opendde(p, chains)

    # 链顺序：优先 full 的 token_asym_id 顺序；否则 summary 数组长度
    if tok_idx:
        chain_labels = list(tok_idx.keys())
    else:
        n = len(p.get("chain_ptm") or [])
        chain_labels = [chr(ord('A') + i) for i in range(n)]
    if not chain_labels:
        chain_labels = [chr(ord('A') + i) for i in range(len(p.get("chain_iptm") or []))]
    r["chain_order"] = ",".join(chain_labels)

    # 每链指标
    cp = p.get("chain_ptm") or []
    ci = p.get("chain_iptm") or []
    cpl = p.get("chain_plddt") or []
    bz_cp = p.get("chains_ptm") or {}      # boltz: {chain_idx: ptm}
    bz_ci = p.get("pair_chains_iptm") or {} # boltz: 对角线 = 链自身 iptm
    # 每链 plddt 索引：原子级用 atom_chain_ids；token 级(boltz)仅当长度与 token 总数一致时用
    chain_atom_idx = None
    n_tok_total = sum(len(v) for v in tok_idx.values()) if tok_idx else 0
    if atom_plddt is not None:
        if atom_chain_ids is not None:
            chain_atom_idx = {c: [i for i, a in enumerate(atom_chain_ids) if a == c]
                              for c in chain_labels}
        elif len(atom_plddt) == n_tok_total:
            chain_atom_idx = {c: list(tok_idx.get(c, [])) for c in chain_labels}
    for pos, c in enumerate(chain_labels):
        r[f"{c}_len_tokens"] = len(tok_idx.get(c, [])) if tok_idx else "N/A"
        if chain_atom_idx is not None:
            vals = [atom_plddt[i] for i in chain_atom_idx.get(c, [])]
            if vals:
                r[f"{c}_mean_plddt"] = _plddt100(float(np.mean(vals)))
            else:
                r[f"{c}_mean_plddt"] = "N/A"
        else:
            r[f"{c}_mean_plddt"] = _plddt100(cpl[pos]) if pos < len(cpl) else "N/A"
        if software == "boltz2":
            r[f"{c}_ptm"] = _num(bz_cp.get(pos))
            r[f"{c}_iptm"] = _num(bz_ci.get(pos, {}).get(pos))
        else:
            r[f"{c}_iptm"] = _num(ci[pos]) if pos < len(ci) else "N/A"
            r[f"{c}_ptm"] = _num(cp[pos]) if pos < len(cp) else "N/A"

    # 全局 plddt
    if atom_plddt is not None and len(atom_plddt) > 0:
        r["mean_plddt"] = _plddt100(float(np.mean(atom_plddt)))
    else:
        r["mean_plddt"] = _plddt100(p.get("plddt"))

    # 两两/整体（先从 summary 写 iptm，再补 pae 派生指标）
    pair_mat = None
    if software == "boltz2":
        pair_mat = [[None] * len(chain_labels) for _ in range(len(chain_labels))]
        pc = p.get("pair_chains_iptm") or {}
        for i in range(len(chain_labels)):
            for j in range(len(chain_labels)):
                pair_mat[i][j] = pc.get(i, {}).get(j)
    else:
        pair_mat = p.get("chain_pair_iptm") or []
    for id1, id2 in itertools.combinations(chain_labels, 2):
        key = f"{id1}-{id2}"
        i1, i2 = chain_labels.index(id1), chain_labels.index(id2)
        r[f"{key}_iptm"] = _avg_pair(pair_mat, i1, i2)
    r.update(pair_metrics_from_pae(pae, chain_labels, tok_idx, pae_cutoff))
    return r


# ---------------- 扫描 / 主入口 ----------------
def scan(root, software):
    if software == "boltz2":
        pat = os.path.join(root, "**", "confidence_*_model_*.json")
    else:
        pat = os.path.join(root, "**", "*_summary_confidence_sample_*.json")
    return sorted(glob.glob(pat, recursive=True))


def main():
    ap = argparse.ArgumentParser(description="Protenix/Boltz2/OpenDDE 统一置信度监控")
    ap.add_argument("--software", required=True, choices=["protenix", "boltz2", "opendde"])
    ap.add_argument("--scan", required=True, help="递归扫描根目录")
    ap.add_argument("--chains", default=None, help="逗号分隔链集合，如 A,B,C")
    ap.add_argument("--pae-cutoff", type=float, default=15.0)
    ap.add_argument("-o", "--out", default="multimodel_metrics.csv")
    args = ap.parse_args()

    files = scan(args.scan, args.software)
    if not files:
        print(f"在 {args.scan} 下未找到 {args.software} 的置信度输出")
        sys.exit(1)
    chains = [c.strip() for c in args.chains.split(",") if c.strip()] if args.chains else None

    parser = {"protenix": parse_protenix, "boltz2": parse_boltz2, "opendde": parse_opendde}[args.software]
    rows = []
    for f in files:
        try:
            p = parser(f)
            rows.append(compute_metrics(p, args.software, chains, args.pae_cutoff))
        except Exception as e:
            rows.append({"model_name": os.path.basename(f), "software": args.software,
                         "status": "failed", "error": str(e)[:300]})

    order = ["model_name", "software", "status", "error", "chain_order",
             "iptm", "ptm", "ranking_score", "mean_plddt", "has_clash", "disorder",
             "overall_ipae", "overall_ipsae_max", "overall_ipsae_min",
             "overall_ipsae_mean", "overall_avg_if_pae"]
    rest = [k for r in rows for k in r if k not in order]
    cols = order + list(dict.fromkeys(rest))
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})
    ok = sum(1 for r in rows if r["status"] == "success")
    print(f"完成: {ok}/{len(rows)} 成功 ({args.software}), 输出 -> {args.out}")
    for r in rows:
        if r["status"] != "success":
            print(f"  ❌ {r['model_name']}: {r['error']}")


if __name__ == "__main__":
    main()
