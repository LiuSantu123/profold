# ProFold 使用与输出说明

ProFold 是面向设计蛋白和蛋白复合物的结构预测筛选流水线。它把结构预测拆成两个相互独立的阶段：Level 1 用于高通量初筛，Level 2 用于对人工选择的子集进行多模型复核。

ProFold 当前只负责结构预测、结构质量评估和结果整理，不负责可溶性、表达量、热稳定性或其他湿实验前序列性质预测。原项目的 Level 3 序列筛选模块已经移除。

## 目录

1. [流水线概览](#流水线概览)
2. [安装和环境](#安装和环境)
3. [配置文件](#配置文件)
4. [输入 manifest](#输入-manifest)
5. [Level 1](#level-1)
6. [Level 2](#level-2)
7. [Slurm 提交](#slurm-提交)
8. [输出目录](#输出目录)
9. [监控与状态](#监控与状态)
10. [指标解释](#指标解释)
11. [resume 行为](#resume-行为)
12. [文件整理和清理](#文件整理和清理)
13. [常见问题](#常见问题)
14. [运行建议](#运行建议)

## 流水线概览

| 阶段 | 主要用途 | 默认模型 | 默认分片大小 |
| --- | --- | --- | ---: |
| Level 1 | 快速 foldback、初步结构和复合物置信度筛选 | ESMFold + AlphaFold 3 | 32 designs/shard |
| Level 2 | 对人工选择的子集进行多模型结构复核 | Protenix、Boltz-2、OpenDDE | 8 designs/shard |

两个阶段通过 TSV manifest 传递 design ID、序列、输入文件路径和结果路径。Level 2 不会隐式读取某个 Level 1 目录；通常由用户使用 `common/select_manifest.py` 从 Level 1 结果中选择设计，然后生成新的 Level 2 manifest。

推荐工作流：

```text
设计序列/复合物输入
        |
        v
Level 1: ESMFold + AF3
        |
        v
metrics/level1.csv + result_manifest.tsv
        |
人工选择或按指标选择
        |
        v
Level 2: Protenix + Boltz-2 + OpenDDE
        |
        v
metrics/level2.csv + result_manifest.tsv
        |
统一检查结构、界面指标和模型共识
```

预测运行状态不等于科学筛选通过。`success` 只表示必需输出存在并且能够解析；是否进入下一阶段仍需要根据 pLDDT、pTM、ipTM、PAE/PDE、界面指标、结构一致性和人工规则判断。

## 安装和环境

### 获取代码

```bash
git clone https://github.com/LiuSantu123/profold.git
cd profold
```

### 安装控制器依赖

```bash
python -m pip install -r requirements.txt
```

控制器依赖负责：

- 读取和校验 TSV manifest
- 生成模型输入
- 启动外部模型命令
- 保存日志和分片结果
- 解析结构和 confidence summary
- 聚合 metrics、status 和 downstream manifest

模型本身需要独立环境、模型权重、数据库和 GPU。ProFold 不把这些大文件放入仓库，也不假设所有机器拥有相同的绝对路径。

### 推荐环境布局

典型布局如下，实际路径可以不同：

```text
conda/envs/
├── profold/       # 控制器、配置、测试
├── esmfold/       # ESMFold CLI
├── af3/           # AlphaFold 3 + JAX
├── protenix/      # Protenix
├── boltz/         # Boltz-2
└── opendde/       # OpenDDE

software/
├── alphafold3/
└── Protenix/

data/
├── alphafold3/models/
└── opendde/
```

### 生成本机配置

仓库提供 `environments/site.example.json` 模板。使用 `scripts/configure.py` 将 `${...}` 路径变量展开：

```bash
python scripts/configure.py \
  --conda-root /opt/miniconda3 \
  --software-root /opt/software \
  --data-root /data/models \
  --output config.local.json
```

脚本不会覆盖已有配置。如果目标文件已经存在，请指定新的输出路径，确认无误后再替换旧配置。

检查 Level 1：

```bash
python scripts/check_environment.py \
  --config config.local.json \
  --level level1
```

检查 Level 2：

```bash
python scripts/check_environment.py \
  --config config.local.json \
  --level level2
```

环境检查只检查 Python 模块、可执行文件和配置路径，不会验证 GPU 显存、模型权重完整性、实际推理兼容性或预测质量。正式批量运行前应使用少量输入进行 GPU smoke test。

## 配置文件

配置文件是 JSON 对象。设计输入不应写进配置文件，而应写入 manifest。

### Level 1 配置示例

```json
{
  "esmfold_bin": "/opt/conda/envs/esmfold/bin/esm-fold",
  "esmfold_num_recycles": 4,
  "esmfold_max_tokens_per_batch": 1024,
  "af3_mode": "direct",
  "af3_python": "/opt/conda/envs/af3/bin/python",
  "af3_script": "/opt/software/alphafold3/run_alphafold.py",
  "af3_model_dir": "/data/models/alphafold3/models",
  "af3_db_dir": "/data/models/alphafold3"
}
```

`af3_mode` 有三个取值：

- `direct`：每个 design 使用一个 AF3 JSON 输入。
- `cid`：同一组蛋白链分别运行 apo 和 holo 两个 AF3 状态。
- `for_wj_cid`：兼容历史 wrapper，需要在 `commands.af3` 中显式配置完整命令。

### Level 2 配置示例

```json
{
  "boltz_bin": "/opt/conda/envs/boltz/bin/boltz",
  "boltz_no_kernels": false,
  "commands": {
    "protenix": [
      "env",
      "PROTENIX_PYTHON=/opt/conda/envs/protenix/bin/python",
      "PROTENIX_SOURCE=/opt/software/Protenix",
      "bash",
      "/opt/software/profold/wrappers/run_protenix.sh",
      "{protenix_input}",
      "{protenix_outdir}"
    ],
    "opendde": [
      "env",
      "CONDA_BASE=/opt/miniconda3",
      "OPENDDE_ENV=opendde",
      "OPENDDE_ROOT_DIR=/data/models/opendde",
      "bash",
      "/opt/software/profold/wrappers/opendde_predict.sh",
      "-i",
      "{opendde_input}",
      "-o",
      "{opendde_outdir}"
    ]
  }
}
```

命令模板中的占位符由 worker 在运行时展开。常用占位符包括：

| 占位符 | 含义 |
| --- | --- |
| `{run_dir}` | 当前 run 根目录 |
| `{shard_dir}` | 当前 shard 的 artifacts 目录 |
| `{design_id}` | 当前设计 ID |
| `{input}` / `{model_input}` | 当前模型输入文件 |
| `{outdir}` / `{model_outdir}` | 当前模型输出目录 |
| `{protenix_input}` | Protenix JSON |
| `{protenix_outdir}` | Protenix 输出目录 |
| `{boltz_input}` | Boltz YAML |
| `{boltz_outdir}` | Boltz 输出目录 |
| `{opendde_input}` | OpenDDE JSON 或转换输入 |
| `{opendde_outdir}` | OpenDDE 输出目录 |

配置中的命令使用参数数组时，不需要自行处理 shell quoting。字符串形式会经过 `shlex.split`，路径包含空格时优先使用数组形式。

### Protenix 的 PAE 数据

Protenix 默认 summary 文件通常只包含 ipTM、pTM 和 pLDDT 等汇总指标。要计算 PAE、ipAE 或 ipSAE，必须启用 atom-level confidence。当前仓库的 `wrappers/run_protenix.sh` 已使用：

```text
--need_atom_confidence true
```

因此应检查输出中是否同时存在：

```text
*_summary_confidence_sample_*.json
*_full_data_sample_*.json
```

如果只有 summary 文件，不能假设 PAE 或 ipSAE 存在。

## 输入 manifest

manifest 使用 UTF-8、首行字段名和真实 TAB 分隔：

```tsv
design_id<TAB>sequence<TAB>chain_ids<TAB>metadata_json
d001<TAB>MKTAA...<TAB>A<TAB>{"target":"protein_x"}
```

命令行中可以用普通制表符文件创建，不建议用空格代替 TAB。

### 必填字段

| 字段 | 说明 |
| --- | --- |
| `design_id` | 稳定主键，首字符必须是字母或数字；允许字母、数字、`.`、`_`、`-` |
| `sequence` 或 `fasta_path` | 一个蛋白序列或多链序列 |

`design_id` 不能重复。一个 design 的序列可以写成多链格式：

```text
MSEQUENCEA:MSEQUENCEB:MSEQUENCEC
```

配合：

```text
chain_ids=A,B,C
```

如果省略 `chain_ids`，系统会按顺序自动使用 A、B、C 等链 ID。

### 常用可选字段

| 字段 | 用途 |
| --- | --- |
| `parent_design_id` | 记录设计来源或母设计 |
| `target_id` | 靶标标识 |
| `ligand_id` | 配体标识 |
| `structure_path` | 参考结构路径 |
| `chain_map` | 结构比较时的链映射信息 |
| `metadata_json` | 额外 JSON 元数据 |
| `af3_json_path` | direct AF3 输入或模板 |
| `cid_apo_json_path` | CID apo 模板 |
| `cid_holo_json_path` | CID holo 模板 |
| `protenix_input_path` | Protenix 专用输入 |
| `boltz_input_path` | Boltz YAML 输入 |
| `opendde_input_path` | OpenDDE 专用输入 |
| `target_spec_path` | 其他目标/模板输入 |

路径相对 manifest 文件所在目录解析，而不是相对当前 shell 目录。准备 run 后，输入路径会被规范化为绝对路径并写入 `input/input_manifest.tsv`。

### Level 1 direct 示例

```tsv
design_id	sequence	af3_json_path
d001	MKTAA...	templates/d001.json
d002	MKTTT...	templates/d002.json
```

### Level 1 CID 示例

```tsv
design_id	sequence	chain_ids	cid_apo_json_path	cid_holo_json_path
cid_001	MKT...:GAV...	A,B	templates/apo.json	templates/holo.json
```

CID 的 apo/holo 模板必须包含相同的蛋白链 ID。apo 通常不含 ligand，holo 至少包含目标 ligand。模板中的蛋白序列会被 manifest 序列替换，MSA 会按无 MSA 路线处理；模板中的其他实体、修饰和键定义由用户负责校验。

### Level 2 复杂体系

对于蛋白-小分子、蛋白-核酸或多实体体系，推荐显式提供每个模型的输入文件：

```tsv
design_id	sequence	protenix_input_path	boltz_input_path	opendde_input_path
ligand_001	MKT...	inputs/ligand_001.json	inputs/ligand_001.yaml	inputs/ligand_001.json
```

不要把 SMILES、CCD ligand 或核酸实体直接拼进 `sequence` 字段。它们应该写在对应模型的输入 JSON/YAML 中。

## Level 1

### 准备 run

```bash
python level1/prepare_level1.py \
  --manifest input.tsv \
  --outdir runs/level1 \
  --shard-size 32 \
  --config config.local.json
```

准备阶段会：

1. 读取并校验 manifest。
2. 检查 design ID 是否重复或非法。
3. 检查序列是否存在。
4. 为每行写入 `task_index`、`shard_index` 和输入 fingerprint。
5. 保存原始输入快照。
6. 按 shard size 生成任务文件。
7. 保存本次 run 的 `config.json`。

### 运行单个 shard

```bash
python level1/run_level1.py \
  --run-dir runs/level1 \
  --shard-index 0
```

worker 会：

1. 读取当前 shard。
2. 生成 ESMFold FASTA。
3. 生成或复制 AF3 输入。
4. 调用 ESMFold。
5. 调用 AF3 或 CID wrapper。
6. 将标准输出和错误输出写入 `logs/`。
7. 解析结构、summary 和 confidence 文件。
8. 写入 `results/parts/part_XXXXX.tsv`。

强制重新运行：

```bash
python level1/run_level1.py \
  --run-dir runs/level1 \
  --shard-index 0 \
  --no-resume
```

### 聚合 Level 1

所有 shard 完成后运行：

```bash
python level1/aggregate_level1.py --run-dir runs/level1
```

聚合器不会重新运行模型。它会读取所有 `results/parts/part_*.tsv`，按照任务 manifest 恢复设计顺序，然后生成：

```text
runs/level1/results/result_manifest.tsv
runs/level1/metrics/level1.csv
runs/level1/status/status.tsv
```

CID 模式还会生成：

```text
runs/level1/metrics/cid_apo.csv
runs/level1/metrics/cid_holo.csv
```

### Level 1 CID 监控

当 `af3_mode=cid` 时，Level 1 的 AF3 命令由 `level1/cid.py` 生成，并调用 bundled AF3 wrapper。wrapper 会启动：

```text
tools/af3/af3_monitor_v16.py
```

该 monitor 负责轮询 AF3 输出、提取指标、归档结构和 confidence，并在校验完成后清理任务级中间产物。CID 的 archive 逻辑位于：

```text
tools/af3/af3_multichain_archive.py
```

如果 monitor 或 archive 校验失败，源文件应保留，不能把任务标记为完整成功。

## Level 2

### 准备 run

```bash
python level2/prepare_level2.py \
  --manifest level1_selected.tsv \
  --outdir runs/level2 \
  --shard-size 8 \
  --config config.local.json
```

### 运行 shard

```bash
python level2/run_level2.py \
  --run-dir runs/level2 \
  --shard-index 0
```

Level 2 对每个 design 依次处理：

```text
Protenix -> Boltz-2 -> OpenDDE
```

每个模型有独立的输入目录、输出目录、日志和解析状态。某个模型失败不会阻止其他模型运行；worker 会把该模型记为 `failed` 或 `missing`，同时保留其他模型结果。

### 聚合 Level 2

```bash
python level2/aggregate_level2.py --run-dir runs/level2
```

输出：

```text
runs/level2/results/result_manifest.tsv
runs/level2/metrics/level2.csv
runs/level2/status/status.tsv
```

### Level 2 当前监控边界

仓库中有独立的：

```text
tools/monitors/multimodel_monitor_v1.py
```

它可以扫描 Protenix、Boltz-2 和 OpenDDE 的输出并计算更完整的 ipTM、ipAE、ipSAE、链级指标和 Boltz affinity。但当前 `level2/run_level2.py` 默认不自动启动该脚本，而是由 `common/model_parsers.py` 做轻量解析：

- 检查 summary 文件是否存在。
- 检查结构文件是否存在。
- 读取 summary 中已有指标。
- 在 Protenix full data 存在时记录其路径。

因此，当前 Level 2 的 monitor 是“结果解析器”，不是持续轮询的多模型 monitor。若需要完整 interface 指标，应单独运行 monitor：

```bash
python tools/monitors/multimodel_monitor_v1.py \
  --software protenix \
  --scan runs/level2/artifacts \
  -o runs/level2/metrics/protenix_monitor.csv
```

不同软件的参数和输出布局以脚本的 `--help` 及实际模型版本为准。正式批处理前应先用一个 design 验证 scan 根目录、文件名和链顺序。

## Slurm 提交

### Level 1

```bash
bash level1/submit_level1.sh \
  --manifest input.tsv \
  --outdir runs/level1 \
  --config config.local.json \
  --partition gpu \
  --time 04:00:00
```

脚本会：

1. 调用 `prepare_level1.py`。
2. 读取 `tasks/n_shards.txt`。
3. 提交 array worker。
4. 提交等待 array 完成的 aggregate job。

输出会打印：

```text
level1_array_job=...
level1_aggregate_job=...
```

### Level 2

```bash
bash level2/submit_level2.sh \
  --manifest level1_selected.tsv \
  --outdir runs/level2 \
  --config config.local.json \
  --partition gpu \
  --time 08:00:00 \
  --dependency LEVEL1_AGGREGATE_JOB_ID
```

`--dependency` 使用 `afterok`。它只约束下游 array；下游 aggregate 会在自己的 array 完成后运行。

排查 Slurm 任务时同时查看：

```bash
squeue -j JOB_ID
sacct -j JOB_ID --format=JobID,State,ExitCode,Elapsed,MaxRSS
```

以及：

```text
runs/levelN/status/status.tsv
runs/levelN/logs/
runs/levelN/results/parts/
```

## 输出目录

标准 run 目录如下：

```text
runs/levelN/
├── config.json
├── input/
│   └── input_manifest.tsv
├── tasks/
│   ├── task_manifest.tsv
│   ├── shard_00000.tsv
│   ├── shard_00001.tsv
│   └── n_shards.txt
├── artifacts/
│   └── shard_XXXXX/
│       ├── inputs/
│       ├── models/
│       ├── af3_inputs/
│       ├── esmfold/
│       └── af3/
├── results/
│   ├── parts/
│   │   └── part_XXXXX.tsv
│   └── result_manifest.tsv
├── metrics/
│   ├── level1.csv 或 level2.csv
│   └── cid_apo.csv / cid_holo.csv
├── status/
│   └── status.tsv
└── logs/
    └── *.log
```

### `input/input_manifest.tsv`

这是输入快照。它用于确认本次 run 实际使用的 design、序列和路径，不应在预测中途手工修改。

### `tasks/task_manifest.tsv`

在输入字段基础上增加：

- `task_index`
- `shard_index`
- `input_fingerprint`

### `results/parts/part_XXXXX.tsv`

每个 worker 生成一个 part 文件。它包含设计状态、模型状态、日志路径、结构路径和可解析指标。只有所有 shard 完成后，aggregate 才能得到完整结果。

### `metrics/levelN.csv`

每个 design/model 一行。缺失指标保持为空，不会用 0 填充。这样可以区分：

- 指标确实为 0。
- 指标缺失。
- 模型没有运行。
- 模型运行失败。

### `status/status.tsv`

既记录 design 总状态，也记录每个模型的状态：

```text
design_id  level   model    status   error   log_path
```

### `results/result_manifest.tsv`

这是给下一阶段或人工筛选使用的结果 manifest。它保留原始设计字段，并增加当前阶段的模型结果字段。

## 结构和 confidence 文件

不同模型的原始命名可能不同，但 ProFold 解析器通常识别以下类型：

### AF3

```text
*.cif
*_summary_confidences.json
*_confidences.json
```

### Protenix

```text
*_summary_confidence_sample_*.json
*_full_data_sample_*.json
*.cif
```

### Boltz-2

常见文件包括：

```text
confidence_*.json
pae_*.npz
plddt_*.npz
*_model_*.cif
affinity_*.json
```

### OpenDDE

```text
*_summary_confidence_sample_*.json
*_full_data_sample_*.json
*.cif
```

这些是模型原始输出，不等于最终交付格式。当前 Level 2 不会自动把所有模型结果扁平化到统一的 `structures/` 和 `confidences/` 目录，也不会自动删除所有非 best sample。

## 监控、状态和错误

### 状态含义

| 状态 | 含义 |
| --- | --- |
| `success` | 必需结构和 summary 存在且成功解析 |
| `partial` | 有部分结果，例如只有结构或只有 summary；或命令失败但留下可解析结果 |
| `failed` | 输出目录存在，但没有有效可识别结果，或 design ID 无法匹配 |
| `missing` | 模型未配置、未运行或输出目录不存在 |

### 典型错误定位顺序

1. 查看 `status/status.tsv`。
2. 找到对应模型的 `log_path`。
3. 查看 `artifacts/shard_XXXXX/` 中是否有输入文件。
4. 检查模型输出是否使用预期文件名。
5. 检查 summary 与 structure 是否属于同一个 design。
6. 检查 GPU 显存、CUDA、权重和数据库路径。
7. 重新运行单个 shard，必要时加 `--no-resume`。

### 不要只看进程退出码

某些模型命令可能返回非零，但仍留下结构或 summary。ProFold 会将这种情况标记为 `partial`，而不是直接当作成功。相反，命令返回 0 但没有结构和 summary 时，也不会被当作成功。

## resume 行为

Level 1 和 Level 2 worker 默认启用 resume：

```bash
python level2/run_level2.py --run-dir runs/level2 --shard-index 0
```

当已有结果被解析为 `success` 时，worker 会复用它，不重新运行对应模型。某个模型失败时，下一次运行通常只重新处理不完整模型。

强制重跑：

```bash
python level1/run_level1.py --run-dir runs/level1 --shard-index 0 --no-resume
```

如果更换了模型版本、输入模板、权重或关键参数，建议使用新的 run 目录，而不是覆盖旧 run。旧 run 应作为可追溯记录保留。

## 文件整理和清理

### 当前已经自动完成的整理

- 输入 manifest 快照。
- 任务分片。
- 运行配置快照。
- 模型日志集中保存。
- 每个 shard 的结果 part 原子写入。
- Level 1 CID wrapper 的 AF3 归档和清理。
- 统一的 metrics/status/result manifest 聚合。

### 当前不会自动完成的整理

Level 2 当前不会自动完成以下交付整理：

- 每个 design 选择 top-ranked best structure。
- 将 best 结构复制到统一 `structures/`。
- 将 confidence/full data 复制到统一 `confidences/`。
- 生成统一的 `summary.json`、`ranking.csv` 和 `index.json`。
- 归档后删除非 best diffusion sample。
- 删除 model cache、processed、MSA 或其他中间目录。

因此不要直接把 `artifacts/` 目录当作最终交付目录。对于正式结果，应先完成指标提取和汇总，再按照以下原则整理：

```text
run/
├── structures/
│   └── <design>_model.cif
├── confidences/
│   └── <design>_confidences.json 或 .npz
├── summary.json
├── ranking.csv
├── index.json
└── _archive/
```

清理必须遵循：

1. 先确认所有 design 的指标已提取。
2. 先生成 summary、ranking 和 checksum manifest。
3. 将非 best 中间文件移动到 `_archive/`。
4. 校验交付目录能够独立读取。
5. 最后再删除 archive，且重要结果应另存到 review/交付目录。

## 指标解释

### 蛋白-蛋白

建议共同查看：

- `ipTM`：界面/相互作用相关置信度。
- `pTM`：整体拓扑置信度。
- interface PAE 或 ipAE：界面相对定位不确定性。
- ipSAE：界面结构质量相关指标。
- interface pLDDT：界面区域局部置信度。
- 跨模型结构 RMSD：不同预测器是否收敛到相似界面。
- clash 或几何异常：是否存在明显不合理接触。

单独的高 pLDDT 不能证明蛋白-蛋白结合界面正确。应优先检查界面区域，而不是只看全链平均值。

### 蛋白-小分子

建议共同查看：

- ligand pose 是否位于合理结合口袋。
- ligand confidence 或等价配体置信度。
- interface PAE/PDE。
- Boltz affinity 预测及其概率字段。
- 不同模型的 ligand pose 是否一致。
- clash、键长、价态和配体几何。

affinity 预测不能替代实验亲和力测定。尤其不要把一个模型的 affinity 分数与另一个模型的 confidence 分数直接放在同一尺度上平均。

### 模型共识

推荐先按模型分别排序，再计算共识，而不是直接混合原始数值：

1. 每个模型内部进行 percentile 或 rank normalization。
2. 对缺失模型保持缺失，不填 0。
3. 记录参与共识的模型数量。
4. 对结构 pose 做 RMSD/接触图一致性分析。
5. 保存每个模型的原始分数和归一化分数。

## 常见问题

### 只有 summary，没有 PAE

检查模型 wrapper 是否启用了 full data 输出。Protenix 需要 `--need_atom_confidence true`。如果输出已经完成但没有 full data，通常需要重新运行，不能从 summary 可靠恢复 PAE。

### 结构存在但状态是 partial

ProFold 要求结构和 summary 同时存在。检查：

- 结构文件是否是 `.cif`、`.mmcif` 或 `.pdb`。
- summary 文件名是否符合 parser pattern。
- design ID 是否出现在文件名或父目录中。
- 结构和 summary 是否属于同一 design。

### Level 2 没有生成结果

检查：

```bash
cat runs/level2/status/status.tsv
find runs/level2/artifacts -maxdepth 5 -type f | sort | head -100
find runs/level2/logs -type f -maxdepth 1 -print
```

常见原因是 wrapper 路径、模型环境变量、输入格式或输出目录契约不匹配。

### 显存不足

不要让 Protenix、Boltz、OpenDDE 与其他占用显存的任务共用同一 GPU。为 Slurm job 显式申请 GPU，并将不同大模型任务分配到空闲设备。

### 任务重复运行

如果修改了输入模板、模型权重、模型版本或关键推理参数，使用新的 `--outdir`。不要依赖旧目录的 resume 结果混合不同版本的预测。

### 路径在另一台机器失效

运行目录会保存绝对路径。迁移机器时重新生成配置并重新 prepare；不要直接复用包含旧机器路径的 `config.json` 和 task manifest。

## 建议的正式交付检查

提交结构给下游分析或湿实验前，建议执行以下检查：

```text
[ ] 输入 manifest 和序列 fingerprint 已保存
[ ] config.json 已保存，且记录模型路径和关键参数
[ ] 所有 shard 都已 aggregate
[ ] status.tsv 中没有未解释的 partial/failed/missing
[ ] 每个 design 的结构和 confidence 文件可以匹配
[ ] Protenix 需要的 full data 已生成
[ ] PAE/ipAE/ipSAE 的计算口径一致
[ ] 已提取 run-level summary 和 ranking
[ ] 每个 design 已选择 best structure
[ ] 非 best sample 已归档，不会混入交付目录
[ ] 交付结构、confidence、metrics 已生成 checksum
[ ] 结果目录中没有未说明的旧版本输出
```

## 测试和版本

运行回归测试：

```bash
python -m unittest discover -s tests -v
```

当前 release：

```text
ProFold 0.0.37
Git tag: v0.0.37
```

模型权重、数据库、机器配置和预测结果不属于源码 release。正式报告中应额外记录模型源码 commit、权重版本、CUDA/驱动、输入 fingerprint 和实际命令行。
