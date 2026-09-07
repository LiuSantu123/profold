#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AlphaFold3 CID工作流一体化脚本 v18 (统一输出、增量处理版本)
优化了目录结构和处理流程，所有输出都在指定的路径下

核心改进：
- 所有输出(json、esm_output、af3_output、cifs、json备份)都在一个目录树下
- 用户只需指定一个输出路径即可
- 与v16监控脚本配合，支持增量JSON保存和af3_output自动清理
"""

import json
import argparse
import random
import os
import sys
import logging
import subprocess
import time
import signal
import re
import shutil
import hashlib
import tempfile
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser, Polypeptide
from af3_multichain_core import read_designs, protein_ids, fill_template
from af3_multichain_archive import load_bundle

# ======================== 日志配置 ========================
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

MONITOR_PID = None

# ======================== pTM提取工具 ========================
def find_all_ptm_values(log_text):
    """从ESMFold日志中提取所有pTM值"""
    pattern = r'Predicted structure for (\w+) with length \d+, pLDDT \d+\.\d+, pTM (\d+\.\d+)'
    matches = re.findall(pattern, log_text)
    return {protein_id: float(ptm) for protein_id, ptm in matches}

# ======================== 核心工具函数 ========================
def parse_fasta_cid(fasta_path):
    return read_designs(fasta_path)


def batch_generate_json_cid(seqs, template_with_lig, template_without_lig, json_out_dir, chain_ids=None):
    templates = {}
    for variant, path in (("with_lig", template_with_lig), ("without_lig", template_without_lig)):
        with open(path) as handle:
            templates[variant] = json.load(handle)
    order = chain_ids or protein_ids(templates["with_lig"])
    if set(order) != set(protein_ids(templates["without_lig"])):
        raise ValueError("with/without ligand templates must have identical protein chain IDs")
    payloads = []
    for name, *sequences in seqs:
        seed = int(hashlib.sha256(name.encode()).hexdigest()[:8], 16) % 10000 + 1
        for variant, template in templates.items():
            payload = fill_template(template, f"{name}_af3_{variant}", sequences, order, seed)
            payloads.append(payload)
    Path(json_out_dir).mkdir(parents=True, exist_ok=True)
    for payload in payloads:
        path = Path(json_out_dir) / (payload["name"] + ".json")
        tmp = path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(payload, indent=2) + "\n")
        tmp.replace(path)
    return json_out_dir

def run_esmfold_cid_clean(args, seqs, esm_output_dir):
    """
    运行ESMFold预测各蛋白链单体，保留结果和指标
    """
    logger.info(f"\n===== 启动 ESMFold 单体预测 =====")
    logger.info(f"📌 ESM输出目录: {esm_output_dir}")

    esm_fasta = os.path.join(esm_output_dir, "cid_monomers.fasta")
    with open(esm_fasta, 'w') as f:
        for name, *sequences in seqs:
            for cid, sequence in zip(args.sequence_chain_ids, sequences):
                f.write(f">{name}_{cid.lower()}\n{sequence}\n")

    log_file = os.path.join(esm_output_dir, "esmfold.log")
    esm_cmd = [
        args.esmfold_path,
        "-i", esm_fasta,
        "-o", esm_output_dir
    ]
    logger.info(f"📌 运行: {' '.join(esm_cmd)}")

    try:
        with open(log_file, 'w', encoding='utf-8') as log_f:
            subprocess.run(esm_cmd, check=True, stdout=log_f, stderr=subprocess.STDOUT)
        logger.info(f"✅ ESMFold预测完成")
    except Exception as e:
        logger.error(f"❌ ESMFold运行失败: {e}")
        return {}

    # 解析结果
    esm_results = {}
    pdb_files = list(Path(esm_output_dir).glob("*.pdb"))
    parser = PDBParser(QUIET=True)

    logger.info(f"📊 处理 {len(pdb_files)} 个ESMFold结果...")

    for pdb_file in pdb_files:
        model_name = pdb_file.stem
        try:
            struct = parser.get_structure('x', str(pdb_file))
            model = struct[0]

            b_factors = []
            for chain in model:
                for res in chain:
                    if Polypeptide.is_aa(res):
                        for atom in res:
                            b_factors.append(atom.get_bfactor())
                            break

            if b_factors:
                mean_plddt = round(np.mean(b_factors), 3)
                esm_results[model_name] = {
                    "pdb_path": str(pdb_file),
                    "esm_mean_plddt": mean_plddt,
                    "esm_ptm": None
                }
        except Exception as e:
            logger.debug(f"处理PDB失败 {model_name}: {e}")

    # 从日志提取pTM
    try:
        with open(log_file, 'r', encoding='utf-8') as f:
            log_content = f.read()
        ptm_dict = find_all_ptm_values(log_content)

        for model_name, ptm in ptm_dict.items():
            if model_name in esm_results:
                esm_results[model_name]["esm_ptm"] = round(ptm, 3)
                logger.info(f"  ✓ {model_name}: pLDDT={esm_results[model_name]['esm_mean_plddt']}, pTM={ptm:.3f}")
    except Exception as e:
        logger.error(f"❌ 提取pTM失败: {e}")

    # 清理中间文件：只保留FASTA和LOG，删除其他临时文件
    logger.info(f"🧹 清理ESMFold中间文件...")
    for item in Path(esm_output_dir).iterdir():
        if item.is_file() and item.suffix in ['.pkl', '.pth', '.bak']:
            try:
                item.unlink()
            except:
                pass

    return esm_results

def start_monitor_process(args, json_input_dir, esm_results, fasta_path=None):
    """启动监控进程"""
    global MONITOR_PID
    monitor_script = args.monitor_script
    if not os.path.exists(monitor_script):
        raise FileNotFoundError(f"监控脚本不存在: {monitor_script}")

    current_af3_out = str(args.af3_output_dir)
    csv_path = str(args.csv_output)

    monitor_cmd = [
        str(args.af3_env),
        str(monitor_script),
        "--af3-output-dir", current_af3_out,
        "--csv-output", csv_path,
        "--interval", str(args.check_interval),
        "--processes", str(args.monitor_processes),
        "--pae-cutoff", str(args.pae_cutoff),
        "--json-input-dir", str(json_input_dir)
    ]

    if fasta_path:
        monitor_cmd.extend(["--fasta-path", str(fasta_path)])

    if args.run_esmfold:
        monitor_cmd.extend(["--run-esm-compare", "true"])
        json_esm_results = os.path.join(args.temp_dir, "esm_results.json")
        monitor_cmd.extend(["--esm-results-json", str(json_esm_results)])
        monitor_cmd.extend(["--usalign-path", str(args.usalign_path)])

    if args.keep_source:
        monitor_cmd.append('--keep-source')

    if args.analysis_chains:
        monitor_cmd.extend(["--chain-ids", *[str(c) for c in args.analysis_chains]])

    if args.noseqs:
        monitor_cmd.extend(["--noseqs", "true"])

    if args.archive_dir:
        monitor_cmd.extend(["--archive-dir", args.archive_dir])
        logger.info(f"📦 CIF归档目录: {args.archive_dir}")

    # JSON备份目录
    if args.json_archive_dir:
        monitor_cmd.extend(["--json-archive-dir", str(args.json_archive_dir)])

    monitor_log = os.path.join(current_af3_out, "monitor.log")
    Path(os.path.dirname(monitor_log)).mkdir(parents=True, exist_ok=True)

    logger.info(f"\n===== 启动监控进程 =====")
    if args.json_archive_dir:
        logger.info(f"📌 JSON备份目录: {args.json_archive_dir}")

    with open(monitor_log, 'w', encoding='utf-8') as log_f:
        monitor_process = subprocess.Popen(
            monitor_cmd + ([] if '--keep-source' in monitor_cmd else ['--keep-source']),
            stdout=log_f, stderr=log_f, preexec_fn=os.setsid
        )

    MONITOR_PID = monitor_process.pid
    logger.info(f"✅ 监控已启动 (PID: {MONITOR_PID})")
    return monitor_process, monitor_cmd

def stop_monitor_process():
    """停止监控进程"""
    global MONITOR_PID
    if MONITOR_PID is not None:
        try:
            os.killpg(os.getpgid(MONITOR_PID), signal.SIGTERM)
            logger.info(f"📌 已终止监控 (PID: {MONITOR_PID})")
            MONITOR_PID = None
        except Exception as e:
            logger.warning(f"⚠️  终止监控失败: {e}")

def run_af3_prediction(args, json_input_dir, af3_output_dir):
    archive = Path(args.archive_dir or (Path(af3_output_dir) / 'archive'))
    with tempfile.TemporaryDirectory(prefix='af3_pending_', dir=args.temp_dir) as pending:
        count = 0
        for source in Path(json_input_dir).glob('*.json'):
            bundle = archive / source.stem
            if (bundle / 'result.json').exists():
                load_bundle(bundle)
                if not (bundle / 'input.json').is_file():
                    raise ValueError(f'archive lacks verifiable input for {source.stem}; use a new output directory')
                saved = json.loads((bundle / 'input.json').read_text())
                if saved == json.loads(source.read_text()):
                    continue
                raise ValueError(f'archived input differs for {source.stem}; use a new output directory')
            shutil.copy2(source, Path(pending) / source.name)
            count += 1
        if not count:
            logger.info('All AF3 inputs have validated archived outputs')
            return 0
        return _run_af3_prediction(args, pending, af3_output_dir)


def _run_af3_prediction(args, json_input_dir, af3_output_dir):
    """运行AF3预测"""
    af3_cmd = [
        args.af3_env, args.af3_script,
        "--input_dir", json_input_dir,
        "--output_dir", af3_output_dir,
        "--model_dir", args.model_dir,
        "--db_dir", args.db_dir
    ]
    logger.info(f"\n===== 启动AF3预测 =====")
    logger.info(f"📌 输入: {json_input_dir}")
    logger.info(f"📌 输出: {af3_output_dir}")

    try:
        result = subprocess.run(af3_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        logger.info(f"✅ AF3预测正常结束")
        return 0
    except subprocess.CalledProcessError as e:
        logger.error(f"❌ AF3预测失败 (返回码: {e.returncode})")
        if e.stderr:
            logger.error(f"错误信息: {e.stderr[:500]}...")
        return e.returncode

# ======================== 主程序 ========================
def main():
    parser = argparse.ArgumentParser(description="AF3 CID工作流一体化脚本 v17-unified")

    # 必需参数
    parser.add_argument('--fasta', required=True, help='输入FASTA (格式: >name\\nseq_a:seq_b)')
    parser.add_argument('--af3-script', required=True, help='run_alphafold.py路径')
    parser.add_argument('--model-dir', required=True, help='AF3模型目录')
    parser.add_argument('--db-dir', required=True, help='AF3数据库目录')
    parser.add_argument('--template-with-lig', required=True, help='含小分子的模板JSON')
    parser.add_argument('--template-without-lig', required=True, help='不含小分子的模板JSON')

    # 输出相关
    parser.add_argument('--json-output-dir', required=True, help='JSON输出目录')
    parser.add_argument('--af3-output-dir', required=True, help='AF3输出目录')
    parser.add_argument('--csv-output', required=True, help='CSV汇总输出路径')
    parser.add_argument('--archive-dir', help='CIF归档目录')
    parser.add_argument('--json-archive-dir', help='JSON备份目录')
    parser.add_argument('--esm-output-dir', help='ESMFold输出目录')
    parser.add_argument('--temp-dir', help='临时文件目录')

    # 监控相关
    parser.add_argument('--monitor-script', required=True, help='监控脚本路径')
    parser.add_argument('--af3-env', default='/xcfhome/zpzeng/originrepo/alphafold3/alphafold3_env/bin/python', help='AF3 Python环境')
    parser.add_argument('--check-interval', type=int, default=300, help='监控间隔(秒)')
    parser.add_argument('--monitor-processes', type=int, default=4, help='监控进程数')
    parser.add_argument('--pae-cutoff', type=float, default=15.0, help='PAE截断值')

    # ESMFold相关
    parser.add_argument('--run-esmfold', action='store_true', help='是否运行ESMFold')
    parser.add_argument('--esmfold-path', default='/xcfhome/yzmeng/miniconda3/envs/zb/bin/esm-fold', help='ESMFold命令路径')
    parser.add_argument('--usalign-path', default='/xcfhome/yhliu/002_software/009_TMalign/USalign', help='USalign路径')

    # 链分析
    parser.add_argument('--analysis-chains', nargs='+', default=None, help='分析链ID；默认所有输出链')

    parser.add_argument('--sequence-chain-ids', nargs='+', help='FASTA链顺序；默认按with-lig模板蛋白链顺序')

    # 其他选项
    parser.add_argument('--keep-source', action='store_true', help='保留原始AF3文件')
    parser.add_argument('--noseqs', action='store_true', help='CSV中不输出序列')

    args = parser.parse_args()

    # 设置默认临时目录
    if not args.temp_dir:
        args.temp_dir = os.path.join(os.path.dirname(args.af3_output_dir), ".temp")

    Path(args.temp_dir).mkdir(parents=True, exist_ok=True)

    def signal_handler(signum, frame):
        logger.info("\n⚠️  收到终止信号，正在清理...")
        stop_monitor_process()
        sys.exit(1)

    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    try:
        # 第1步: 解析输入
        logger.info("\n" + "="*70)
        logger.info("第1步: 解析输入FASTA")
        logger.info("="*70)
        seqs = parse_fasta_cid(args.fasta)
        with open(args.template_with_lig) as handle:
            args.sequence_chain_ids = args.sequence_chain_ids or protein_ids(json.load(handle))
        if len({c.lower() for c in args.sequence_chain_ids}) != len(args.sequence_chain_ids):
            raise ValueError('protein chain IDs must remain unique when lowercased')
        # Validate every design before launching expensive prediction.
        batch_generate_json_cid(seqs, args.template_with_lig, args.template_without_lig,
                                args.json_output_dir, args.sequence_chain_ids)

        # 第2步: ESMFold预测
        esm_results = {}
        if args.run_esmfold:
            logger.info("\n" + "="*70)
            logger.info("第2步: ESMFold单体预测")
            logger.info("="*70)

            if not args.esm_output_dir:
                args.esm_output_dir = os.path.join(os.path.dirname(args.af3_output_dir), "esm_output")

            Path(args.esm_output_dir).mkdir(parents=True, exist_ok=True)
            esm_results = run_esmfold_cid_clean(args, seqs, args.esm_output_dir)

            esm_results_file = os.path.join(args.temp_dir, "esm_results.json")
            with open(esm_results_file, 'w') as f:
                json.dump(esm_results, f)

            logger.info(f"✅ ESMFold结果保存: {len(esm_results)} 个")

        # 第3步: 生成AF3 JSON
        logger.info("\n" + "="*70)
        logger.info("第3步: 生成AF3输入JSON")
        logger.info("="*70)
        json_dir = args.json_output_dir

        # 第4步: 启动监控和AF3
        logger.info("\n" + "="*70)
        logger.info("第4步: 启动AF3预测和监控")
        logger.info("="*70)

        monitor_process, monitor_base_cmd = start_monitor_process(
            args, json_input_dir=json_dir, esm_results=esm_results, fasta_path=args.fasta
        )

        code = run_af3_prediction(args, json_dir, args.af3_output_dir)

        # 第5步: 最终收尾
        logger.info("\n" + "="*70)
        logger.info("第5步: 最终处理和汇总")
        logger.info("="*70)

        try:
            if monitor_process.poll() is None:
                stop_monitor_process()
                monitor_process.wait()

            logger.info(f"📌 运行最终一次性扫描...")
            once_cmd = monitor_base_cmd + ['--once']
            if code and '--keep-source' not in once_cmd:
                once_cmd.append('--keep-source')
            subprocess.run(once_cmd, check=True, capture_output=False)
            logger.info(f"✅ 所有任务处理完成")
        except Exception as e:
            logger.error(f"❌ 最终处理出错: {e}", exc_info=True)
            code = code or 1

        logger.info("\n" + "="*70)
        logger.info(f"🎉 全部完成！")
        logger.info(f"📄 结果汇总: {args.csv_output}")
        if args.json_archive_dir:
            logger.info(f"📊 JSON备份: {args.json_archive_dir}")
        if args.esm_output_dir:
            logger.info(f"🔬 ESM输出: {args.esm_output_dir}")
        logger.info("="*70)

        sys.exit(code)

    except Exception as e:
        stop_monitor_process()
        logger.error(f"\n❌ 程序异常: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
