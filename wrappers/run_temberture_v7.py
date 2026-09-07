#!/usr/bin/env python3
"""
TemBERTure蛋白质热稳定性预测脚本
修复分类模型返回结果处理问题
"""

from temBERTure import TemBERTure
import pandas as pd
from Bio import SeqIO
import argparse
import sys
import time
from pathlib import Path
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def parse_fasta(fasta_path):
    """解析FASTA文件，返回序列列表"""
    sequences = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences.append({
            "id": record.id,
            "seq": str(record.seq).replace("\n", "")
        })
    return sequences

def predict_single_sequence(seq_id, sequence, tm_models, cls_model):
    """
    预测单个序列的Tm值和热稳定性分数
    
    参数:
        seq_id: 序列ID
        sequence: 蛋白质序列
        tm_models: Tm预测模型列表
        cls_model: 热稳定性分类模型
    
    返回:
        dict: 包含预测结果的字典
    """
    # 1. Tm值预测
    tm_predictions = []
    for i, model in enumerate(tm_models, 1):
        try:
            result = model.predict(sequence)
            
            # 从不同格式的返回结果中提取数值
            if isinstance(result, (list, np.ndarray)):
                tm_value = float(result[0])
            elif isinstance(result, (int, float)):
                tm_value = float(result)
            else:
                # 尝试直接转换为浮点数
                tm_value = float(result)
            
            tm_predictions.append(tm_value)
            
        except Exception as e:
            print(f"警告: 序列 {seq_id} 的Tm预测失败 (模型{i}): {e}")
            tm_predictions.append(None)
    
    # 计算平均Tm值（忽略None值）
    valid_tm = [tm for tm in tm_predictions if tm is not None]
    if valid_tm:
        avg_tm = sum(valid_tm) / len(valid_tm)
    else:
        avg_tm = None
        print(f"警告: 序列 {seq_id} 所有Tm预测都失败")
    
    # 2. 热稳定性分类预测
    thermophilicity_score = None
    try:
        # 分类模型返回的是列表：[分类标签列表, 分数数组]
        result = cls_model.predict(sequence)
        
        # 调试输出
        print(f"序列 {seq_id} 的分类模型返回结果: {result}")
        print(f"返回结果类型: {type(result)}")
        
        # 提取分数
        thermophilicity_score = extract_score_from_result(result)
        
        if thermophilicity_score is None:
            print(f"警告: 无法从分类结果中提取分数: {result}")
        else:
            print(f"提取的热稳定性分数: {thermophilicity_score}")
            
    except Exception as e:
        print(f"警告: 序列 {seq_id} 的热稳定性预测失败: {e}")
    
    return {
        "sid": seq_id,
        "fasta": sequence,
        "tem_TM": round(avg_tm, 4) if avg_tm is not None else None,
        "tem_SC": round(thermophilicity_score, 4) if thermophilicity_score is not None else None
    }

def extract_score_from_result(result):
    """从分类模型返回结果中提取分数"""
    if result is None:
        return None
    
    try:
        # 处理列表格式：[分类标签列表, 分数数组]
        if isinstance(result, list) and len(result) >= 2:
            # 第二个元素应该是分数数组
            score_element = result[1]
            
            # 如果是numpy数组
            if isinstance(score_element, np.ndarray):
                if score_element.size > 0:
                    return float(score_element[0])
            
            # 如果是列表
            elif isinstance(score_element, list):
                if len(score_element) > 0:
                    return float(score_element[0])
            
            # 如果是标量
            elif isinstance(score_element, (int, float, np.float32, np.float64)):
                return float(score_element)
            
            # 尝试直接转换
            else:
                return float(score_element)
        
        # 如果是元组格式
        elif isinstance(result, tuple) and len(result) >= 2:
            score_element = result[1]
            if isinstance(score_element, (np.ndarray, list)):
                if isinstance(score_element, np.ndarray) and score_element.size > 0:
                    return float(score_element[0])
                elif isinstance(score_element, list) and len(score_element) > 0:
                    return float(score_element[0])
            else:
                return float(score_element)
        
        # 如果是直接的分数值
        elif isinstance(result, (int, float, np.float32, np.float64)):
            return float(result)
        
        # 如果是单个数组
        elif isinstance(result, np.ndarray):
            if result.size > 0:
                return float(result[0])
        
        # 如果是单个列表
        elif isinstance(result, list):
            if len(result) > 0:
                # 列表中的第一个元素可能是分数
                return float(result[0])
        
        # 其他情况
        else:
            print(f"警告: 未知的结果格式: {type(result)}, 值: {result}")
            return None
            
    except Exception as e:
        print(f"警告: 提取分数时出错: {e}")
        print(f"出错的结果类型: {type(result)}, 值: {result}")
        return None

def predict_tm_ensemble(fasta_path, output_csv, tm_model_paths, cls_model_path, device='cuda'):
    """
    批量预测FASTA文件中的所有序列
    
    参数:
        fasta_path: 输入FASTA文件路径
        output_csv: 输出CSV文件路径
        tm_model_paths: Tm预测模型路径列表
        cls_model_path: 热稳定性分类模型路径
        device: 计算设备
    """
    # 检查输入文件
    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        print(f"错误: 输入文件不存在: {fasta_path}")
        sys.exit(1)
    
    print(f"读取FASTA文件: {fasta_path}")
    sequences = parse_fasta(fasta_path)
    
    if not sequences:
        print("错误: FASTA文件中没有找到序列")
        sys.exit(1)
    
    print(f"找到 {len(sequences)} 个序列")
    
    # 初始化模型
    print("初始化Tm预测模型...")
    tm_models = []
    for i, model_path in enumerate(tm_model_paths, 1):
        print(f"  加载Tm模型 {i}/{len(tm_model_paths)}: {model_path}")
        try:
            model = TemBERTure(
                adapter_path=model_path,
                device=device,
                batch_size=16,
                task='regression'
            )
            tm_models.append(model)
        except Exception as e:
            print(f"错误: 无法加载Tm模型 {model_path}: {e}")
            sys.exit(1)
    
    print("初始化热稳定性分类模型...")
    print(f"  加载分类模型: {cls_model_path}")
    try:
        cls_model = TemBERTure(
            adapter_path=cls_model_path,
            device=device,
            batch_size=1,
            task='classification'
        )
    except Exception as e:
        print(f"错误: 无法加载分类模型 {cls_model_path}: {e}")
        sys.exit(1)
    
    print(f"所有模型加载完成，使用设备: {device}")
    print("-" * 60)
    
    # 批量预测
    print("开始预测...")
    start_time = time.time()
    results = []
    
    for i, seq_data in enumerate(sequences, 1):
        seq_id = seq_data["id"]
        sequence = seq_data["seq"]
        
        print(f"预测序列 {i}/{len(sequences)}: {seq_id[:30]}...")
        
        result = predict_single_sequence(seq_id, sequence, tm_models, cls_model)
        results.append(result)
        
        # 显示当前进度
        if i % 5 == 0 or i == len(sequences):
            elapsed = time.time() - start_time
            avg_time = elapsed / i
            print(f"进度: {i}/{len(sequences)} | 已用时: {elapsed:.1f}s | 平均: {avg_time:.1f}s/序列")
    
    # 创建DataFrame
    df = pd.DataFrame(results)
    
    # 检查预测结果
    tm_failed = df[df['tem_TM'].isna()]
    sc_failed = df[df['tem_SC'].isna()]
    
    if not tm_failed.empty:
        print(f"警告: {len(tm_failed)} 个序列的Tm值预测失败")
        for _, row in tm_failed.iterrows():
            print(f"  - {row['sid']}")
    
    if not sc_failed.empty:
        print(f"警告: {len(sc_failed)} 个序列的热稳定性预测失败")
        for _, row in sc_failed.iterrows():
            print(f"  - {row['sid']}")
    
    # 保存结果
    output_path = Path(output_csv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"保存结果到: {output_path}")
    df.to_csv(output_path, index=False)
    
    # 显示统计信息
    end_time = time.time()
    total_time = end_time - start_time
    
    print("\n" + "="*60)
    print("预测完成!")
    print("="*60)
    print(f"总序列数: {len(df)}")
    print(f"成功预测Tm值: {len(df) - len(tm_failed)}")
    print(f"成功预测热稳定性: {len(df) - len(sc_failed)}")
    print(f"总耗时: {total_time:.1f}秒")
    print(f"平均每个序列: {total_time/len(df):.1f}秒")
    
    # 显示结果预览
    if not df.empty:
        print("\n结果预览 (前5行):")
        print(df.head().to_string(index=False))
    
    return df

def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='TemBERTure蛋白质热稳定性预测脚本',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  python predict_tm.py input.fasta -o predictions.csv
  python predict_tm.py sequences.fasta -o results.csv
  python predict_tm.py data.fasta -o /path/to/output.csv --device cpu

默认模型路径:
  Tm模型: /xcfhome/yhliu/002_software/075_TemBERTure/temBERTure/temBERTure_TM/replica1/ 等
  分类模型: /xcfhome/yhliu/002_software/075_TemBERTure/temBERTure/temBERTure_CLS/
        """
    )
    
    # 必需参数
    parser.add_argument('fasta_file', type=str, 
                       help='输入FASTA文件路径')
    
    # 可选参数
    parser.add_argument('-o', '--output', type=str, required=True,
                       help='输出CSV文件路径 (必需)')
    parser.add_argument('--tm_models', type=str, nargs='+',
                       default=[
                           '/xcfhome/yhliu/002_software/075_TemBERTure/temBERTure/temBERTure_TM/replica1/',
                           '/xcfhome/yhliu/002_software/075_TemBERTure/temBERTure/temBERTure_TM/replica2/',
                           '/xcfhome/yhliu/002_software/075_TemBERTure/temBERTure/temBERTure_TM/replica3/'
                       ],
                       help='Tm预测模型路径 (默认: 三个replica模型)')
    parser.add_argument('--cls_model', type=str,
                       default='/xcfhome/yhliu/002_software/075_TemBERTure/temBERTure/temBERTure_CLS/',
                       help='热稳定性分类模型路径')
    parser.add_argument('--device', type=str, default='cuda',
                       choices=['cuda', 'cpu'],
                       help='计算设备: cuda 或 cpu (默认: cuda)')
    
    args = parser.parse_args()
    
    # 检查输出文件路径
    output_path = Path(args.output)
    if output_path.suffix.lower() != '.csv':
        print("注意: 输出文件扩展名不是.csv，已自动添加.csv后缀")
        output_path = output_path.with_suffix('.csv')
    
    # 执行预测
    predict_tm_ensemble(
        fasta_path=args.fasta_file,
        output_csv=str(output_path),
        tm_model_paths=args.tm_models,
        cls_model_path=args.cls_model,
        device=args.device
    )

if __name__ == "__main__":
    main()