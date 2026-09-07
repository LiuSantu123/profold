#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AF3 监控脚本 v16 (完整改进版)
核心改进：
1. 严格检查：只有四个结构都有，才在CSV中输出该设计
2. completion_status字段：AF3双变体为2/2，启用ESMFold比对时为(N+2)/(N+2)
3. 最终扫描：AF3结束后，完整扫描所有设计并生成缺失报告
4. 对比FA和CSV：确保设计个数完全对应
5. 增量JSON保存：每个设计完成时立即保存JSON到指定目录（不等待批处理）
6. 自动清理：设计完成后立即删除af3_output下对应的原始文件夹，节省磁盘空间
"""

import os
import sys
import json
import time
import logging
import argparse
import multiprocessing as mp
import csv
import numpy as np
import shutil
import re
import subprocess
import tempfile
from pathlib import Path
from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa
from af3_multichain_archive import process_task, load_bundle, cleanup_bundle

# ======================== 配置 ========================
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('monitor_runtime.log', encoding='utf-8'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# ======================== 核心检查函数（新增） ========================
def get_fa_design_list(fa_path):
    """从FA文件中提取所有设计名，按顺序"""
    design_names = []
    try:
        with open(fa_path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('>'):
                    name = line[1:].strip().split()[0]
                    design_names.append(name)
        logger.info(f"✅ 从FA解析出 {len(design_names)} 个设计")
    except Exception as e:
        logger.error(f"❌ 读取FA文件失败: {e}")
    return design_names

def check_design_completion(design_name, esm_results, af3_tasks_dict, require_esm=False):
    tasks = af3_tasks_dict.get(design_name, {})
    missing = [v for v in ("with_lig", "without_lig") if tasks.get(v, {}).get("status") != "success"]
    expected = 2
    if require_esm:
        order = tasks.get("with_lig", {}).get("protein_chain_order", "")
        chains = order.split(",") if order else []
        if not chains:
            missing.append("protein_chain_order")
            expected += 1
        for cid in chains:
            expected += 1
            item = esm_results.get(f"{design_name}_{cid.lower()}", {})
            if not item.get("pdb_path") or not Path(item["pdb_path"]).is_file():
                missing.append(f"esm_{cid}")
    return not missing, f"{expected - len(missing)}/{expected}", missing


def complete_status(value):
    parts = str(value).split("/")
    return len(parts) == 2 and parts[0].isdigit() and parts[0] == parts[1]


def generate_final_report(fa_path, csv_path, esm_output_dir, af3_output_dir):
    """
    最终检查：对比FA和CSV，生成详细的缺失报告
    这个函数在AF3完全结束时调用，进行终极质量检查
    """
    logger.info("\n" + "="*70)
    logger.info("🔍 最终质量检查：验证所有设计的完整性")
    logger.info("="*70)
    
    # 1. 获取FA中的所有设计
    fa_designs = get_fa_design_list(fa_path)
    fa_designs_set = set(fa_designs)
    logger.info(f"📌 FA文件设计总数: {len(fa_designs)}")
    
    # 2. 获取CSV中已完成的设计
    csv_designs = {}
    if os.path.exists(csv_path):
        try:
            with open(csv_path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    name = row.get('design_name')
                    status = row.get('completion_status', 'unknown')
                    if name:
                        csv_designs[name] = status
        except Exception as e:
            logger.error(f"❌ 读取CSV失败: {e}")
    
    csv_designs_set = set(csv_designs.keys())
    logger.info(f"✅ CSV中已完成的设计数: {len(csv_designs_set)}")
    
    # 3. 找出缺失的设计
    missing_designs = fa_designs_set - csv_designs_set
    incomplete_designs = {k: v for k, v in csv_designs.items() if not complete_status(v)}
    
    # 4. 详细报告
    logger.info("\n--- 结果汇总 ---")
    if not missing_designs and not incomplete_designs:
        logger.info(f"🎉 完美！全部 {len(fa_designs)} 个设计都已完成！")
        print(f"\n✨ 成功完成所有设计: {len(fa_designs)}/{len(fa_designs)}")
        return {"status": "COMPLETE", "total": len(fa_designs), "completed": len(fa_designs), "missing": 0}
    
    # 5. 列出缺失和不完整的
    logger.warning(f"\n⚠️  发现问题:")
    
    if missing_designs:
        logger.warning(f"  ❌ 完全缺失: {len(missing_designs)} 个设计")
        for name in sorted(missing_designs):
            logger.warning(f"     - {name}")
    
    if incomplete_designs:
        logger.warning(f"  ⚠️  不完整: {len(incomplete_designs)} 个设计")
        for name, status in sorted(incomplete_designs.items()):
            logger.warning(f"     - {name} ({status})")
    
    # 6. 生成详细报告文件
    report_path = csv_path.replace(".csv", "_completion_report.txt")
    try:
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("="*70 + "\n")
            f.write("CID设计完整性检查最终报告\n")
            f.write("="*70 + "\n\n")
            
            f.write(f"FA文件设计总数: {len(fa_designs)}\n")
            f.write(f"CSV已完成数: {len(csv_designs_set)}\n")
            f.write(f"完全缺失: {len(missing_designs)}\n")
            f.write(f"不完整: {len(incomplete_designs)}\n\n")
            
            f.write("详情:\n")
            f.write("-"*70 + "\n")
            
            for design_name in fa_designs:
                if design_name in csv_designs_set:
                    status = csv_designs.get(design_name, "unknown")
                    if complete_status(status):
                        f.write(f"✅ {design_name:<50} {status}\n")
                    else:
                        f.write(f"⚠️  {design_name:<50} {status}\n")
                else:
                    f.write(f"❌ {design_name:<50} 0/4 (缺失)\n")
        
        logger.info(f"\n📄 详细报告已保存: {report_path}")
    except Exception as e:
        logger.error(f"❌ 保存报告失败: {e}")
    
    # 7. 返回结果
    return {
        "status": "INCOMPLETE",
        "total": len(fa_designs),
        "completed": len(csv_designs_set),
        "incomplete": len(incomplete_designs),
        "missing": len(missing_designs),
        "missing_list": sorted(missing_designs),
        "incomplete_list": sorted(incomplete_designs.keys())
    }

# ======================== 原有辅助函数 ========================
class ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id
    def accept_chain(self, chain):
        return chain.id == self.chain_id

def extract_chain_to_pdb(cif_path, chain_id, temp_dir):
    try:
        parser = MMCIFParser(QUIET=True)
        struct = parser.get_structure('x', str(cif_path))
        io = PDBIO()
        io.set_structure(struct)
        tmp_file = tempfile.mktemp(suffix=f"_{chain_id}.pdb", dir=temp_dir)
        io.save(tmp_file, ChainSelect(chain_id))
        return tmp_file
    except Exception as e:
        logger.debug(f"提取链失败: {e}")
        return None

def run_usalign_and_parse(esm_pdb_path, af3_chain_pdb_path, usalign_path):
    try:
        cmd = [usalign_path, esm_pdb_path, af3_chain_pdb_path]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        output = result.stdout
        
        rmsd = "N/A"
        seq_id = "N/A"
        tm_score = "N/A"
        
        m_rmsd = re.search(r'RMSD=\s*([\d.]+)', output)
        if m_rmsd:
            rmsd = round(float(m_rmsd.group(1)), 3)
        
        m_seqid = re.search(r'Seq_ID=n_identical/n_aligned=\s*([\d.]+)', output)
        if m_seqid:
            seq_id = round(float(m_seqid.group(1)), 3)
        
        m2 = re.search(r'TM-score=\s*([\d\.]+)\s*\(normalized by length of Structure_2', output)
        if m2:
            tm_score = round(float(m2.group(1)), 4)
        
        return rmsd, seq_id, tm_score
        
    except Exception as e:
        logger.debug(f"USalign 运行失败: {e}")
        return "N/A", "N/A", "N/A"

def get_entities_from_json(model_name, json_input_dir):
    if not json_input_dir:
        return {}, {}, []
    try:
        json_file = os.path.join(json_input_dir, f"{model_name}.json")
        if not os.path.exists(json_file):
            candidates = list(Path(json_input_dir).glob(f"*{model_name}*.json"))
            if candidates:
                json_file = str(candidates[0])
            else:
                return {}, {}, []
        
        with open(json_file, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        seqs = {}
        types = {}
        ids_order = []
        
        for s in data.get('sequences', []):
            if 'protein' in s:
                chain_id = s['protein']['id']
                seqs[f"seq_{chain_id}"] = s['protein']['sequence']
                types[chain_id] = 'protein'
                ids_order.append(chain_id)
            elif 'ligand' in s:
                chain_id = s['ligand']['id']
                ligand_name = s['ligand'].get('smiles', s['ligand'].get('ccdCode', 'LIGAND'))
                seqs[f"seq_{chain_id}"] = ligand_name
                types[chain_id] = 'ligand'
                ids_order.append(chain_id)
        
        return seqs, types, ids_order
    except Exception as e:
        logger.debug(f"读取JSON实体失败: {e}")
        return {}, {}, []

def get_chain_indices_and_lengths(structure, json_chain_ids, json_entity_types):
    try:
        model = structure[0]
        indices = {}
        lengths = {}
        
        pae_idx = 0
        for chain in model:
            c_id = chain.id
            indices[c_id] = []
            elements = []
            for res in chain:
                if is_aa(res):
                    elements.append(pae_idx)
                    pae_idx += 1
                else:
                    elements.append(pae_idx)
                    pae_idx += 1
            
            indices[c_id] = elements
            lengths[f"len_{c_id}"] = len(elements)
            
        return indices, lengths, model
    except Exception as e:
        logger.debug(f"解析结构索引失败: {e}")
        return {}, {}, None

def calculate_ip_sae_metrics(pae_matrix, a_idx, b_idx, cutoff=15.0):
    if pae_matrix is None or not a_idx or not b_idx:
        return 0.0, 99.0, 0.0
    
    max_dim = pae_matrix.shape[0]
    valid_a = [i for i in a_idx if i < max_dim]
    valid_b = [i for i in b_idx if i < max_dim]
    
    if not valid_a or not valid_b:
        return 0.0, 99.0, 0.0

    sub = pae_matrix[np.ix_(valid_a, valid_b)]
    valid = sub < cutoff
    if not np.any(valid):
        return 0.0, 99.0, 0.0
    
    L = np.sum(valid)
    d0 = 1.24 * (L ** (1/3)) - 1.8 if L >=27 else 1.0
    
    scores = []
    pae_values_at_interface = []
    for i,j in np.argwhere(valid):
        scores.append(1.0/(1.0 + (sub[i,j]/d0)**2))
        pae_values_at_interface.append(sub[i,j])
    
    if not scores:
        return 0.0, 99.0, 0.0
        
    max_score = round(max(scores), 4)
    min_score = round(min(scores), 4)
    mean_pae_at_interface = round(np.mean(pae_values_at_interface), 3)
    
    return max_score, min_score, mean_pae_at_interface

# ======================== 核心分析 ========================
def process_single_af3_task(folder_path, chain_ids, pae_cutoff, archive_dir, keep_source, json_input_dir, noseqs, json_archive_dir=None):
    return process_task(folder_path, chain_ids, pae_cutoff, archive_dir, keep_source,
                        json_input_dir, noseqs, json_archive_dir)


def process_design(design_name, with_lig_result, without_lig_result, esm_results, usalign_path):
    dr = {"design_name": design_name, "completion_status": "2/2", "metrics_schema": "multichain_v2"}
    for variant, result in (("with_lig", with_lig_result), ("without_lig", without_lig_result)):
        for key, value in result.items():
            if key not in ("model_name", "status", "error", "af3_cif_path"):
                dr[f"{variant}_{key}"] = value
        dr[f"{variant}_cif"] = result.get("af3_cif_path")
    chains = with_lig_result.get("protein_chain_order", "").split(",")
    for cid in filter(None, chains):
        label = cid.lower()
        esm = esm_results.get(f"{design_name}_{label}")
        if not esm:
            continue
        dr.update({f"{label}_esm_plddt": esm.get("esm_mean_plddt"),
                   f"{label}_esm_ptm": esm.get("esm_ptm"), f"{label}_esm_pdb": esm.get("pdb_path")})
        if usalign_path:
            for variant in ("with_lig", "without_lig"):
                with tempfile.TemporaryDirectory(prefix="af3_chain_") as temp_dir:
                    chain_pdb = extract_chain_to_pdb(dr[f"{variant}_cif"], cid, temp_dir)
                    if chain_pdb:
                        rmsd, seqid, tm = run_usalign_and_parse(esm["pdb_path"], chain_pdb, usalign_path)
                        dr.update({f"{variant}_{label}_esm_rmsd": rmsd,
                                   f"{variant}_{label}_esm_tmscore": tm})
    return dr

# ======================== 加载已有结果 ========================
def load_existing_csv(csv_path):
    """加载已有的CSV，把之前的结果都读进来，绝对不丢！"""
    existing = {}
    if os.path.exists(csv_path):
        try:
            with open(csv_path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    name = row.get('design_name')
                    if name:
                        existing[name] = row
            logger.info(f"📂 已加载已有CSV中的 {len(existing)} 个历史结果，不会覆盖它们！")
        except Exception as e:
            logger.warning(f"⚠️  加载已有CSV失败: {e}，将从头开始")
    return existing

# ======================== 结果保存 ========================
def save_results(all_designs, csv_path, noseqs):
    """保存所有的设计结果，永远不覆盖，每次都把所有的写进去"""
    
    # 排序，按带小分子的iptm排序（处理非数字值如"N/A"）
    def get_iptm_score(design):
        iptm = design.get("with_lig_iptm", 0)
        if isinstance(iptm, (int, float)):
            return iptm
        try:
            return float(iptm)
        except (ValueError, TypeError):
            return 0  # 如果是"N/A"或其他非数字，返回0
    
    final_results = sorted(all_designs.values(), key=get_iptm_score, reverse=True)
    
    # 列顺序
    def column_sort_key(col):
        if col == 'design_name': return 0
        if col == 'completion_status': return 1  # 新增
        if col.startswith('with_lig_A-B_'): return 2
        if col.startswith('without_lig_A-B_'): return 3
        if col.startswith('with_lig_') and not col.startswith('with_lig_A-B_') and not col.startswith('with_lig_a_') and not col.startswith('with_lig_b_'): return 4
        if col.startswith('without_lig_') and not col.startswith('without_lig_A-B_') and not col.startswith('without_lig_a_') and not col.startswith('without_lig_b_'): return 5
        if col.startswith('a_esm_') or col.startswith('b_esm_'): return 6
        if col.startswith('with_lig_a_esm_') or col.startswith('with_lig_b_esm_'): return 7
        if col.startswith('without_lig_a_esm_') or col.startswith('without_lig_b_esm_'): return 8
        if col.startswith('len_') or col.startswith('seq_'): return 9
        return 10
    
    all_keys = {'design_name', 'completion_status'}
    for r in final_results:
        all_keys.update(r.keys())
    sorted_fields = sorted(list(all_keys), key=lambda col: (column_sort_key(col), col))
    if noseqs:
        sorted_fields = [f for f in sorted_fields if not f.startswith('seq_') and '_seq_' not in f]

    # 永远是写模式，但是我们把所有的设计都写进去，所以不会丢
    Path(csv_path).parent.mkdir(parents=True, exist_ok=True)
    tmp_path = str(csv_path) + f'.{os.getpid()}.tmp'
    with open(tmp_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=sorted_fields, restval='N/A', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(final_results)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp_path, csv_path)
    
    logger.info(f"📄 汇总CSV完成: {csv_path}, 共 {len(final_results)} 个完整设计")

# ======================== 主循环 ========================
def monitor_loop(args):
    noseqs = args.noseqs.lower() in ('true', '1', 'yes')
    require_esm = args.run_esm_compare.lower() in ('true', '1', 'yes')
    args.archive_dir = args.archive_dir or str(Path(args.af3_output_dir) / "archive")
    processed = load_existing_csv(args.csv_output)
    esm = {}
    if require_esm and args.esm_results_json:
        with open(args.esm_results_json) as handle:
            esm = json.load(handle)
    with mp.Pool(processes=args.processes) as pool:
        while True:
            # Archive-first recovery also works after raw task directories are deleted.
            paths = {}
            for root in (Path(args.archive_dir), Path(args.af3_output_dir)):
                if root.exists():
                    for folder in root.iterdir():
                        if folder.is_dir() and folder.name.endswith(("_af3_with_lig", "_af3_without_lig")):
                            if root == Path(args.archive_dir) and not (folder / 'result.json').is_file():
                                continue
                            paths.setdefault(folder.name, str(folder))
            tasks = [(p, args.chain_ids, args.pae_cutoff, args.archive_dir, args.keep_source,
                      args.json_input_dir, noseqs, args.json_archive_dir) for p in paths.values()]
            results = pool.starmap(process_single_af3_task, tasks)
            grouped = {}
            for result in results:
                name = result["model_name"]
                variant = "without_lig" if name.endswith("_af3_without_lig") else "with_lig"
                design = name[:-len("_af3_" + variant)]
                grouped.setdefault(design, {})[variant] = result
                if result["status"] != "success":
                    logger.warning("%s: %s", name, result.get("error"))
                    processed.pop(design, None)
            completed = []
            for name, variants in grouped.items():
                ready, status, missing = check_design_completion(name, esm, grouped, require_esm)
                if not ready:
                    processed.pop(name, None)
                    continue
                row = process_design(name, variants["with_lig"], variants["without_lig"], esm, args.usalign_path)
                row["completion_status"] = status
                processed[name] = row
                completed.append(name)
            save_results(processed, args.csv_output, noseqs)
            # Re-read the committed CSV before any manifest-listed source removal.
            committed = load_existing_csv(args.csv_output)
            if not args.keep_source:
                for name in completed:
                    if name not in committed:
                        continue
                    for variant in ("with_lig", "without_lig"):
                        try:
                            cleanup_bundle(Path(args.archive_dir) / f"{name}_af3_{variant}", args.af3_output_dir)
                        except Exception as exc:
                            logger.warning("cleanup retained %s: %s", name, exc)
            if args.once:
                return
            time.sleep(args.interval)

def find_unprocessed_tasks(af3_dir, processed_names, archive_dir=None):
    """找到所有未处理的AF3任务
    
    方案B改进：同时支持af3_dir和archive_dir目录
    这样即使AF3设置了--archive-dir参数，Monitor也能找到所有已完成任务
    """
    unprocessed = []
    
    # 首先扫描af3_output目录
    if os.path.exists(af3_dir):
        for item in os.listdir(af3_dir):
            item_path = os.path.join(af3_dir, item)
            if os.path.isdir(item_path) and item not in processed_names:
                has_core = any(Path(item_path).glob("*_summary_confidences.json")) or any(Path(item_path).glob("*model.cif"))
                if has_core:
                    unprocessed.append(item_path)
                    processed_names.add(item)
    
    # 方案B：如果指定了归档目录，也扫描它
    if archive_dir and os.path.exists(archive_dir):
        for item in os.listdir(archive_dir):
            item_path = os.path.join(archive_dir, item)
            # 归档目录中通常只有CIF文件，但我们也检查匹配的task目录
            if os.path.isdir(item_path) and item not in processed_names:
                has_core = any(Path(item_path).glob("*.cif"))
                if has_core:
                    unprocessed.append(item_path)
                    processed_names.add(item)
    
    return unprocessed

def backup_json_files(af3_output_dir, backup_dir):
    """
    自动备份AF3输出的JSON文件到af3_confidence_json目录
    保护JSON文件在CIF文件被归档后仍可访问
    """
    os.makedirs(backup_dir, exist_ok=True)
    backed_up = 0
    
    if not os.path.exists(af3_output_dir):
        return backed_up
    
    # 扫描af3_output中的所有JSON文件
    for item in os.listdir(af3_output_dir):
        item_path = os.path.join(af3_output_dir, item)
        if os.path.isdir(item_path):
            # 查找该目录中的JSON文件
            for json_file in Path(item_path).glob("*.json"):
                # 复制JSON文件到备份目录，保留原名
                backup_json_path = os.path.join(backup_dir, f"{item}_{json_file.name}")
                try:
                    shutil.copy2(json_file, backup_json_path)
                    backed_up += 1
                except Exception as e:
                    logger.debug(f"备份JSON失败 {json_file}: {e}")
    
    if backed_up > 0:
        logger.info(f"💾 已备份 {backed_up} 个JSON文件到 {backup_dir}")

    return backed_up

def main():
    parser = argparse.ArgumentParser(description="AF3监控脚本 v15 改进版")
    parser.add_argument('--af3-output-dir', required=True)
    parser.add_argument('--csv-output', required=True)
    parser.add_argument('--chain-ids', nargs='+')
    parser.add_argument('--json-input-dir', help='JSON目录')
    parser.add_argument('--fasta-path', help='FA文件路径，用于最终检查')
    parser.add_argument('--interval', type=int, default=30)
    parser.add_argument('--processes', type=int, default=4)
    parser.add_argument('--pae-cutoff', type=float, default=15.0)
    parser.add_argument('--archive-dir')
    parser.add_argument('--json-archive-dir', help='JSON备份目录')
    parser.add_argument('--keep-source', action='store_true')
    parser.add_argument('--once', action='store_true')
    parser.add_argument('--noseqs', default='false')
    parser.add_argument('--run-esm-compare', default='false')
    parser.add_argument('--esm-results-json')
    parser.add_argument('--usalign-path')
    args = parser.parse_args()
    
    if not args.chain_ids:
        args.chain_ids = []

    # 设置默认JSON备份目录
    if not args.json_archive_dir:
        args.json_archive_dir = os.path.join(os.path.dirname(args.af3_output_dir), "af3_confidence_json")

    monitor_loop(args)
    
    # AF3完全结束后，进行最终完整性检查
    if args.once and args.fasta_path and os.path.exists(args.fasta_path):
        esm_output_dir = os.path.join(os.path.dirname(args.af3_output_dir), "esm_output")
        report = generate_final_report(
            args.fasta_path,
            args.csv_output,
            esm_output_dir,
            args.af3_output_dir
        )
        
        # 打印总结
        logger.info("\n" + "="*70)
        logger.info("📊 最终统计")
        logger.info("="*70)
        logger.info(f"总设计数: {report.get('total', 0)}")
        logger.info(f"已完成: {report.get('completed', 0)}")
        logger.info(f"缺失: {report.get('missing', 0)}")
        logger.info(f"不完整: {report.get('incomplete', 0)}")
        
        if report['status'] != 'COMPLETE':
            logger.warning(f"\n⚠️  仍有问题需要解决，详见完整报告")
            raise SystemExit(1)

if __name__ == "__main__":
    main()
