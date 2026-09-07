# Level1 输入规范

版本：2026-09-07。TSV 使用 UTF-8，首行为字段名，真实 TAB 分隔，每行一个完整设计。
禁止重复列名、重复 design_id 和列数错位；空字段保留 TAB，不使用 `NA` 代替空路径。
输入文件与模型配置分开：TSV 描述设计，config JSON 描述模型环境和运行模式。

## 模式

| af3_mode | AF3 行为 | 模板输入 |
| --- | --- | --- |
| direct | 每个设计预测一次 | af3_json_path 或 target_spec_path |
| cid | 同一组蛋白链分别预测 apo/holo，分别统计 | cid_apo_json_path + cid_holo_json_path |
| for_wj_cid | 兼容旧 wrapper 配置 | commands.af3 内的 with/without-lig 参数 |

CID 的 apo 指同一组蛋白链去除配体后的复合物，不是拆成单体。
ESMFold stage 只运行一次；apo/holo 双预测属于 AF3 stage。
不自动用 holo 结构作为 apo 初始构象，也不自动做 apo/holo 结构差异比对。

## CID 必填字段

| 字段 | 规范 |
| --- | --- |
| design_id | 唯一，首位字母或数字，其余限字母、数字、点、下划线、连字符 |
| sequence | 所有蛋白链，标准写法 `SEQ_A:SEQ_B:SEQ_C`；单体只有一段，同聚体也逐链列出 |
| chain_ids | 与序列一一对应、逗号分隔，如 `A,B,C`；唯一的大写字母 ID，支持 `AA` |
| cid_apo_json_path | AF3 JSON 模板，包含全部蛋白链，不含 ligand 实体 |
| cid_holo_json_path | AF3 JSON 模板，包含相同蛋白链 ID，至少一个 ligand 实体 |

序列允许 20 种标准氨基酸及 X；小写转大写，`|` 转为 `:`，不接受空链、空格、终止符或 gap。
配体用 JSON 的 CCD/SMILES 等实体定义，不能写进 sequence。
两份模板可以使用不同实体顺序，序列始终按 chain_ids 匹配，不按位置替换。
一批 run 只接受一对模板及一种链 ID 顺序；不同聚体数、配体或模板体系分批运行。
模板蛋白序列会被 TSV 替换，MSA 清空；apo/holo 的同一设计使用相同 seed。
模板其他实体、修饰和键定义由用户保证有效，预检查不替代 AF3 完整输入校验。

路径相对输入 TSV 所在目录解析，随后保存为绝对路径。推荐直接填写 sequence；
兼容字段 fasta_path 不应指向“多个独立设计共用的 FASTA”，需先展开为一行一设计。

## 可选字段

`parent_design_id`、`target_id`、`ligand_id`：来源与分组标签，不控制预测。
`structure_path`、`chain_map`、`metadata_json`：参考结构与后处理信息，不是预测必需参数。
额外自定义字段会保留到结果 manifest。普通 direct 模式不要求 CID 两个模板字段。

## 示例与运行

[input.tsv](../examples/cid/input.tsv) 配有 apo/holo JSON，内容为三链短序列加 ATP 的格式示例，
仅用于验证输入格式，不是有科学意义的 CID 候选。
[config.cid.example.json](../config.cid.example.json) 可使用当前机器的模型入口，
`cid` 模式会自动构造 wrapper 命令，无需填写 commands.af3。

```bash
V37=/xcfhome/yhliu/14_magpcr/z_zcodex/pipeline_scripts/00_latest/pipeline/v37
python "$V37/level1/prepare_level1.py" --manifest "$V37/examples/cid/input.tsv" \
  --config "$V37/config.cid.example.json" --outdir runs/cid_l1 --shard-size 32
# 正式预测时将 manifest 和模板替换为实际候选；下列命令需要 GPU 环境。
python "$V37/level1/run_level1.py" --run-dir runs/cid_l1 --shard-index 0
python "$V37/level1/aggregate_level1.py" --run-dir runs/cid_l1
```

多个 shard 必须全部运行后再 aggregate，或通过 submit_level1.sh 提交 array。

## CID 输出

- `metrics/cid_apo.csv`、`metrics/cid_holo.csv`：两张独立指标表，每设计一行，包含各自状态。
- `results/result_manifest.tsv`：`af3_apo_*` 和 `af3_holo_*` 两组指标，均含整体、每链、逐对指标。
- `metrics/level1.csv`：保留 ESMFold/AF3 的全层汇总。
- `status/status.tsv`：分别记录 af3_apo、af3_holo，两个状态都成功才算 AF3 成功。
- 原有 `af3_with_lig_*` / `af3_without_lig_*` 和底层归档名称保留兼容。
  新 cid 模式的无状态前缀 `af3_*` 主指标仅对应 holo，不用 apo 替代缺失 holo。

ipTM 从 AF3 原生链对矩阵提取，ipAE 为跨链 PAE 双向均值，ipSAE 为蛋白链对 d0res。
无链对或非蛋白 ipSAE 留空；两种状态独立统计，不平均、不混用、不自动设筛选阈值。
归档校验、CSV 原子提交、按清单清理和 resume 沿用多链 CID wrapper 契约。
