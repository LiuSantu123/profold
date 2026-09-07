#!/usr/bin/env python3
"""
NetSolP结果处理脚本
用于提取和处理NetSolP预测结果，保留指定字段和格式
作者：David Baker实验室
"""

import pandas as pd
import numpy as np
import argparse
import sys
import os
from pathlib import Path

def process_netsolp_results(input_file, output_file, delimiter=','):
    """
    处理NetSolP预测结果CSV文件
    
    参数:
    ----------
    input_file : str
        输入CSV文件路径
    output_file : str
        输出CSV文件路径
    delimiter : str
        CSV文件分隔符，默认为逗号
    """
    
    print(f"正在处理文件: {input_file}")
    
    try:
        # 1. 读取CSV文件
        # 使用pandas的read_csv函数读取文件
        # 参数说明:
        # - sep: 分隔符
        # - header=0: 第一行作为列名
        # - low_memory=False: 避免类型推断时的内存警告
        df = pd.read_csv(input_file, sep=delimiter, header=0, low_memory=False)
        
        print(f"成功读取文件，共 {len(df)} 行数据")
        print(f"原始列名: {list(df.columns)}")
        
        # 2. 检查必要的列是否存在
        # 定义我们需要的列名
        required_columns = ['sid', 'fasta', 'predicted_solubility', 'predicted_usability']
        
        # 但NetSolP输出可能有不同的列名，我们需要尝试不同的可能性
        # 常见列名变体
        column_variants = {
            'sid': ['sid', 'id', 'ID', 'protein_id', 'seq_id'],
            'fasta': ['fasta', 'sequence', 'seq', 'protein_sequence'],
            'predicted_solubility': ['predicted_solubility', 'solubility', 'avg_solubility'],
            'predicted_usability': ['predicted_usability', 'usability', 'avg_usability']
        }
        
        # 查找实际列名
        actual_columns = []
        for required_col, variants in column_variants.items():
            found = False
            for variant in variants:
                if variant in df.columns:
                    actual_columns.append(variant)
                    found = True
                    print(f"找到列 {required_col}: 实际列名为 '{variant}'")
                    break
            if not found:
                print(f"警告: 未找到 {required_col} 列，尝试的列名: {variants}")
                print(f"可用列名: {list(df.columns)}")
                sys.exit(1)
        
        # 3. 提取需要的列
        # 使用实际找到的列名
        extracted_df = df[actual_columns].copy()
        
        # 4. 重命名列为标准名称
        # 创建重命名映射
        rename_dict = {}
        for i, std_name in enumerate(required_columns):
            rename_dict[actual_columns[i]] = std_name
        
        extracted_df = extracted_df.rename(columns=rename_dict)
        
        # 5. 处理数值列，保留6位小数
        # 检查哪些列是数值类型
        numeric_columns = ['predicted_solubility', 'predicted_usability']
        
        for col in numeric_columns:
            if col in extracted_df.columns:
                # 转换为数值类型，错误值转换为NaN
                extracted_df[col] = pd.to_numeric(extracted_df[col], errors='coerce')
                # 四舍五入到6位小数
                extracted_df[col] = extracted_df[col].round(6)
                print(f"已处理 {col} 列，最小值: {extracted_df[col].min():.6f}, 最大值: {extracted_df[col].max():.6f}")
        
        # 6. 检查并处理缺失值
        missing_count = extracted_df.isnull().sum().sum()
        if missing_count > 0:
            print(f"警告: 发现 {missing_count} 个缺失值")
            # 可以选择填充或删除，这里显示信息
            for col in extracted_df.columns:
                col_missing = extracted_df[col].isnull().sum()
                if col_missing > 0:
                    print(f"  - {col}: {col_missing} 个缺失值")
        
        # 7. 保存到新的CSV文件
        # 设置pandas显示格式，避免科学计数法
        pd.set_option('display.float_format', '{:.6f}'.format)
        
        # 保存为CSV，不保存索引
        extracted_df.to_csv(output_file, index=False, float_format='%.6f')
        
        print(f"处理完成!")
        print(f"原始数据形状: {df.shape}")
        print(f"提取后数据形状: {extracted_df.shape}")
        print(f"输出文件已保存: {output_file}")
        print(f"输出列: {list(extracted_df.columns)}")
        
        # 8. 显示前几行作为验证
        print("\n输出文件前5行预览:")
        print(extracted_df.head().to_string())
        
        return extracted_df
        
    except FileNotFoundError:
        print(f"错误: 输入文件 '{input_file}' 未找到")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"错误: 输入文件 '{input_file}' 为空")
        sys.exit(1)
    except Exception as e:
        print(f"处理文件时出错: {str(e)}")
        sys.exit(1)

def main():
    """主函数，解析命令行参数并调用处理函数"""
    
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(
        description='NetSolP结果处理器 - 提取和处理NetSolP预测结果',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  %(prog)s netsolp_results.csv processed_results.csv
  %(prog)s input.csv output.csv --delimiter ","
  %(prog)s --input data.csv --output processed.csv

David Baker实验室提示:
  1. 确保输入文件是逗号分隔的CSV格式
  2. 脚本会自动检测列名变体
  3. 数值结果会四舍五入到6位小数
  4. 输出文件可直接用于下游分析
        """
    )
    
    # 添加参数
    # 位置参数
    parser.add_argument('input_file', nargs='?', 
                       help='输入CSV文件路径')
    parser.add_argument('output_file', nargs='?', 
                       help='输出CSV文件路径')
    
    # 可选参数
    parser.add_argument('--input', '-i', 
                       help='输入CSV文件路径（替代位置参数）')
    parser.add_argument('--output', '-o', 
                       help='输出CSV文件路径（替代位置参数）')
    parser.add_argument('--delimiter', '-d', default=',',
                       help='CSV文件分隔符（默认: ,）')
    parser.add_argument('--version', '-v', action='version',
                       version='NetSolP Processor v1.0 (David Baker Lab)')
    
    # 解析参数
    args = parser.parse_args()
    
    # 确定输入输出文件
    input_file = args.input if args.input else args.input_file
    output_file = args.output if args.output else args.output_file
    
    # 验证参数
    if not input_file:
        print("错误: 必须指定输入文件")
        parser.print_help()
        sys.exit(1)
    
    if not output_file:
        # 如果没有指定输出文件，则基于输入文件名生成
        input_path = Path(input_file)
        output_file = input_path.stem + "_processed.csv"
        print(f"注意: 未指定输出文件，将使用: {output_file}")
    
    # 检查输入文件是否存在
    if not os.path.exists(input_file):
        print(f"错误: 输入文件 '{input_file}' 不存在")
        sys.exit(1)
    
    # 调用处理函数
    process_netsolp_results(input_file, output_file, args.delimiter)

if __name__ == "__main__":
    main()