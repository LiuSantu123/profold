# ProFold Code Wiki

> 适用版本：`0.0.38.2`（见 [VERSION](../VERSION)）。本文档面向开发者，系统描述 ProFold 仓库的整体架构、模块职责、关键类与函数、依赖关系与运行方式。用户视角的使用说明见 [PROFOLD_GUIDE.md](PROFOLD_GUIDE.md)。

---

## 目录

1. [项目概述](#1-项目概述)
2. [整体架构](#2-整体架构)
3. [目录结构](#3-目录结构)
4. [核心模块详解](#4-核心模块详解)
5. [关键类与函数索引](#5-关键类与函数索引)
6. [数据流与运行目录布局](#6-数据流与运行目录布局)
7. [配置系统](#7-配置系统)
8. [依赖关系](#8-依赖关系)
9. [项目运行方式](#9-项目运行方式)
10. [测试与发布](#10-测试与发布)
11. [设计原则与注意事项](#11-设计原则与注意事项)

---

## 1. 项目概述

**ProFold** 是一个面向设计蛋白质/蛋白复合物的**结构预测与结构筛选流水线**。它将"快速首轮折叠"与"多模型确认"分离为两个独立阶段，并以**清单（manifest）驱动、分片（shard）执行、可断点续跑（resume）**的方式组织所有运行：

| 阶段 | 目的 | 模型 |
| --- | --- | --- |
| **Level 1** | 高通量回折（foldback）与首轮复合物置信度评估 | ESMFold + AlphaFold 3 |
| **Level 2** | 对筛选子集做多模型交叉确认 | Protenix + Boltz-2 + OpenDDE |

核心特征：

- **manifest 是唯一阶段间接口**：`design_id` 是稳定主键；每个阶段接受一个 TSV 清单，写出独立 run 目录；阶段之间通过 [common/select_manifest.py](../common/select_manifest.py) 显式选样，无隐式状态传递。
- **控制器/模型分离**：本仓库只做编排、解析、归档与报告；所有 GPU 模型（权重、数据库、CUDA 环境）都在仓库外以独立 conda 环境安装，通过配置文件中的可执行路径接入。
- **可审计输出**：每次运行包含输入快照、分片任务、模型产物、逐模型指标、状态表、结果清单和日志；状态（success/partial/failed/missing）只描述"运行与解析是否完成"，**不是生物学意义上的通过/失败标签**。
- **报告与交付**：Level 2 聚合后可生成 BoltzGen 风格 PDF 报告（指标分布、模型散点、design×model 热图、Top-N 结构面板）。

> 本仓库不做溶解度、表达量或热稳定性预测（原 Level 3 序列筛选模块已移除）。

---

## 2. 整体架构

### 2.1 分层视图

```
┌─────────────────────────────────────────────────────────────────────┐
│ 用户入口层                                                           │
│   scripts/profold.py  (v38.2 紧凑 facade: prepare/run/aggregate/report)│
│   或旧版各 level 独立脚本 (prepare_*/run_*/aggregate_*)                │
├─────────────────────────────────────────────────────────────────────┤
│ 共享编排层  common/                                                   │
│   pipeline.py  runner.py  manifest.py  input_builder.py              │
│   model_parsers.py  select_manifest.py                               │
├──────────────────────────┬──────────────────────────────────────────┤
│ Level 1 (level1/)        │ Level 2 (level2/)                        │
│  ESMFold → AF3 批式执行   │  Protenix/Boltz-2/OpenDDE 逐设计执行      │
│  direct / cid / for_wj_cid│  单态 / CID apo+holo 双态                │
├──────────────────────────┴──────────────────────────────────────────┤
│ 模型接入层  tools/ + wrappers/                                        │
│   AF3 CID 一体化 wrapper + 监控器 (tools/af3/)                        │
│   多模型置信度监控器 (tools/monitors/)                                │
│   Protenix/OpenDDE 包装脚本 (wrappers/)                               │
├─────────────────────────────────────────────────────────────────────┤
│ 后处理层                                                             │
│   report/generate_report.py  (统计+绘图+PDF)                          │
│   structure_align/align_structures.py (USalign 骨架比对)              │
├─────────────────────────────────────────────────────────────────────┤
│ 外部模型环境 (仓库外, config.local.json 指定路径)                      │
│   esm-fold CLI │ AlphaFold3 (run_alphafold.py) │ Protenix │ Boltz-2  │
│   OpenDDE │ USalign │ VMD/Tachyon 渲染器                              │
└─────────────────────────────────────────────────────────────────────┘
```

### 2.2 执行模型

- **prepare**：读取用户 TSV 清单 → 规范化（校验 `design_id`、解析相对路径、补齐 sequence）→ 按 `shard_size` 切片 → 写出 run 目录（`input/`、`tasks/`、`config.json` 等）。
- **run**：每个 shard 独立可执行（Slurm array 一个任务对应一个 shard）。worker 负责构造模型输入、以批式（Level 1，每个预测器只启动一次）或逐设计（Level 2）方式调用外部模型，然后解析输出目录为统一的前缀化字段（如 `af3_iptm`、`protenix_mean_plddt`），写出 `results/parts/part_XXXXX.tsv`。
- **resume**：worker 重跑时先解析已有输出，已 success 的模型标记 `*_reused=1` 并跳过预测；`--no-resume` 强制重新预测。facade 还会用 sha256 指纹校验 input/config/VERSION 是否变化，变化则拒绝续跑。
- **aggregate**：合并所有 part → 与 task manifest 按 `design_id` 连接 → 写出 `results/result_manifest.tsv`、`metrics/levelN.csv`、`status/status.tsv`。
- **select**：`common/select_manifest.py` 显式选择进入下一阶段的设计（按 ID 或按字段 Top-N），并写 `selection.json` 审计记录。

### 2.3 AF3 的三种模式（Level 1）

`config.af3_mode` 决定 AF3 stage 行为（见 [level1/run_level1.py](../level1/run_level1.py) 与 [level1/INPUT_SCHEMA.md](../level1/INPUT_SCHEMA.md)）：

| 模式 | 行为 | 模板输入 |
| --- | --- | --- |
| `direct`（默认） | 每个设计预测一次，每个设计生成一份 AF3 JSON（`build_af3_json`） | `af3_json_path` 或 `target_spec_path` |
| `cid` | 同一组蛋白链分别预测 apo（无配体）与 holo（含配体）两个状态 | `cid_apo_json_path` + `cid_holo_json_path`（整批共用一对模板） |
| `for_wj_cid` | 兼容旧 wrapper：由 `commands.af3` 配置的外部 wrapper 自己生成 JSON | wrapper 参数 |

CID 模式下，`level1/cid.py:cid_command` 组装对 `tools/af3/af3_single_node_v18.py`（CID 一体化脚本）的调用；该脚本内部生成双态 JSON、运行 ESMFold 单体预测、调用 AF3，并配合 `af3_monitor_v16.py` 做增量监控与归档。

---

## 3. 目录结构

```
profold/
├── scripts/                 # 入口与机器配置
│   ├── profold.py           # ★ v38.2 紧凑 facade（推荐入口）
│   ├── configure.py         # site.example.json → config.local.json 渲染
│   ├── check_environment.py # 只读依赖/路径检查
│   └── publish_release.sh   # 发布打包（tag + SHA256）
├── common/                  # ★ 三个 level 共享的编排原语
│   ├── manifest.py          # TSV 清单读写、design_id 校验、指纹
│   ├── pipeline.py          # prepare/分片/聚合共享逻辑
│   ├── runner.py            # 命令渲染、日志执行、原子 JSON
│   ├── input_builder.py     # 各模型输入文件构造（无 MSA）
│   ├── model_parsers.py     # 各模型输出解析（防御式）
│   └── select_manifest.py   # 阶段间显式选样 CLI
├── level1/                  # ESMFold + AF3 高通量阶段
│   ├── prepare_level1.py / run_level1.py / aggregate_level1.py
│   ├── cid.py               # CID 清单契约 + AF3 命令组装
│   ├── str20.py             # 独立 STR 三态验证 runner（私有数据）
│   ├── str20_closeout.py    # STR 运行的紧凑归档/校验
│   ├── INPUT_SCHEMA.md      # Level 1 输入规范（含 CID）
│   └── submit_level1.sh / run_level1_array.sh / aggregate_level1.sh  # Slurm
├── level2/                  # Protenix + Boltz-2 + OpenDDE 确认阶段
│   ├── prepare_level2.py / run_level2.py / aggregate_level2.py
│   └── submit_level2.sh / run_level2_array.sh / aggregate_level2.sh
├── report/
│   └── generate_report.py   # Level 2 统计、绘图、结构面板、PDF
├── tools/
│   ├── af3/                 # AF3 CID 工作流工具（在 AF3 环境内运行）
│   │   ├── af3_multichain_core.py      # 链映射 + token 对齐置信度指标
│   │   ├── af3_multichain_archive.py   # 校验归档 bundle + 清理
│   │   ├── af3_single_node_v18.py      # CID 一体化 wrapper（v18）
│   │   └── af3_monitor_v16.py          # AF3 输出监控器（v16）
│   └── monitors/
│       └── multimodel_monitor_v1.py    # Protenix/Boltz-2/OpenDDE 统一监控
├── wrappers/                # 外部模型包装（多为集群内路径，需按本机适配）
│   ├── protenix_100.py      # Protenix 推理 runner（ByteDance, Apache-2.0）
│   ├── run_protenix.sh / run_protenix_f101.sh
│   ├── run_opendde.sh / opendde_predict.sh / fasta_to_opendde_json.py
│   └── LICENSE.protenix.txt
├── structure_align/
│   ├── align_structures.py  # 通用结构-骨架对齐（USalign + Kabsch）
│   └── README.md
├── tests/                   # unittest 回归套件 + fake wrappers fixtures
├── environments/            # controller.yml、site.example.json、说明
├── examples/cid/            # CID 五列格式示例（apo.json/holo.json/input.tsv）
├── docs/                    # PROFOLD_GUIDE.md（用户手册）、任务清单
├── .github/workflows/test.yml  # CI：unittest + shell 语法检查
├── config.example.json      # Level 1/2 通用配置模板（direct 模式）
├── config.cid.example.json  # CID 模式配置模板
├── requirements.txt         # 控制器依赖（不含任何 GPU 依赖）
├── README.md / RELEASE_NOTES.md / THIRD_PARTY_NOTICES.md / VERSION
```

---

## 4. 核心模块详解

### 4.1 `common/` — 共享编排层

#### [manifest.py](../common/manifest.py) — 清单契约

模块 docstring：*"The manifest is the only interface between levels."* 所有跨阶段数据交换都通过 TSV 清单完成，`design_id` 是稳定主键。

- 常量：`DESIGN_ID_RE`（合法 ID 正则）、`PATH_FIELDS`（需解析为绝对路径的字段集合）、`CORE_FIELDS`（列排序时的核心字段顺序）。
- `read_fasta(path)`：不依赖生物信息库的轻量 FASTA 读取。
- `sequence_from_row(row)`：优先取 `sequence` 列，否则从 `fasta_path` 读第一条记录。
- `normalize_row(raw, manifest_path)`：单行规范化——校验 `design_id` 非空且合法、相对路径相对 TSV 所在目录解析为绝对路径、补齐 `sequence`。
- `read_manifest(path, require_sequence=False)`：读入整个清单；强校验列名唯一/非空、列数对齐、`design_id` 唯一；`require_sequence=True` 时要求每行有序列。
- `write_tsv / write_csv / write_fasta`：**原子写**（先写 `.tmp` 再 rename），列序 = preferred + CORE_FIELDS + 出现序。
- `row_fingerprint(row)`：对排序后的行内容做 sha256，用于任务级输入指纹。
- `metadata(row)`：解析 `metadata_json` 列（必须是 JSON 对象）。

#### [pipeline.py](../common/pipeline.py) — prepare/shard/aggregate 共享逻辑

- `RUN_DIRS`：标准 run 目录骨架 `("input", "tasks", "results/parts", "metrics", "status/parts", "artifacts", "logs")`。
- `prepare_run(input_manifest, outdir, *, level, shard_size, config)`：读清单（level1+cid 模式会先调用 `level1.cid.validate_cid_rows` 预检）→ 建目录 → 每行加 `task_index`/`shard_index`/`input_fingerprint` → 写 `input/input_manifest.tsv`、`tasks/task_manifest.tsv`、每个分片 `tasks/shard_XXXXX.tsv` 与 `tasks/n_shards.txt` → 把 config 与运行元数据快照进 `config.json`（原子写）。
- `load_config(run_dir)` / `load_shard(run_dir, shard_index)`：读取快照配置与分片任务。
- `write_part(run_dir, shard_index, rows)`：写该分片结果 part；`read_part_files`：枚举所有 part。
- `join_results(run_dir, *, model_names)`：以 task manifest 为基准连接所有 part；`design_id` 重复直接报错；缺失 part 的行标记 `status=missing`；**part 中仅合并结果字段**（空值不覆盖基准行，避免 loader 默认空列冲掉原始输入）。
- `status_rows(rows, level, model_names)`：生成逐设计+逐模型的状态审计表。
- `metric_rows(rows, model_names)`：把 `f"{model}_*"` 前缀字段拆成"每设计×每模型一行"的长表（剔除 `_status/_error/_log_path/_output_dir`），**不填补缺失分数**。

#### [runner.py](../common/runner.py) — 命令与执行原语

- `atomic_json / read_json / tail_text`：原子 JSON 写、读、日志尾部提取。
- `command_template(value, context)`：把配置中的命令字符串/列表按 `{placeholder}` 渲染成 argv 列表；缺失占位符即抛 `KeyError`。
- `configured_command(config, name, context)`：从 `config["commands"][name]` 读取用户自定义命令（优先级最高）。
- `run_logged(command, log_path, *, cwd, env)`：以子进程执行命令，stdout/stderr 合流写入日志文件（首行记录完整命令，尾行记录 `returncode`），返回 `(rc, tail)`；`OSError` 记为 rc=127。
- `stage_result(model, status, *, error, log_path, **extra)`：构造 `{model}_status/{model}_error/{model}_log_path` 前缀化结果字典——这是所有 worker 结果行的统一约定。

#### [input_builder.py](../common/input_builder.py) — 模型输入构造

原则：**只构造 run 本地输入，绝不修改用户源文件**。

- `protein_records(row)`：核心解析函数——把"单序列（`A:B`/`A|B` 多链）或 FASTA"解析为有序 `(chain_id, sequence)` 列表；`chain_ids` 列可选（缺省按 A、B、C… 命名），数量/重复强校验；FASTA 支持 `__chain__` 命名约定。
- `build_af3_json(row, destination, *, seed)`：复制现成 AF3 JSON 或基于 `target_spec_path` 模板填充；`metadata_json.af3_sequences` 支持按链 ID 的 dict 或列表两种映射；写 `name`/`modelSeeds`，默认清空 MSA/templates（`clear_msa`）；模板链 ID 与序列链严格双向匹配。
- `build_protenix_json / build_boltz_yaml / build_opendde_fasta`：为三个 Level 2 模型构造**无 MSA** 最小输入；`ligand` 列支持 `CCD_XXX` 与 `smiles:...` 两种写法（逗号/分号分隔多个）。若用户直接给了 `*_input_path`，则复制该文件并先经 `_validate_no_msa_input` 校验（递归扫描 JSON/YAML，禁止非空 MSA 字段）。
- `status_from_models(statuses)`：多模型状态归并规则——全 success→`success`；任一 success/partial→`partial`；全 missing→`missing`；否则 `failed`。

#### [model_parsers.py](../common/model_parsers.py) — 输出解析

原则：*"defensive parsers for model output contracts"*——缺失细节文件记为 `partial`，**绝不伪造置信度数值**。

- 通用工具：`_files`（按模式递归收集）、`_choose`（把输出文件关联到 design——优先精确名匹配；目录里只有一个文件才允许模糊关联，否则宁缺毋滥）、`_status`（success/partial/failed/missing 判定）、`_number`（清洗 N/A、非有限数）、`_scaled_plddt`（0–1 的 pLDDT 统一放大到 0–100）。
- `parse_esmfold(root, design_id, *, prefix, native=False)`：`native=True` 时只认 `{design_id}.pdb`，并用 gemmi 校验坐标与 B-factor 形式的 pLDDT；否则解析自定义 summary CSV + 结构文件。
- `parse_af3(root, design_id, *, mode)`：`direct` 走 `_parse_af3_variant`（summary confidences JSON + 结构 CIF/PDB + 细节 confidences JSON，提取 iptm/ptm/ranking_score/fraction_disordered/has_clash/mean_plddt）；`for_wj_cid`/`cid` 模式解析 `with_lig`/`without_lig` 两个变体、校验归档 `result.json` 的 sha256、再合成 `af3_*`（主结果，cid 下强制取 holo/with_lig）与 `af3_apo_*`/`af3_holo_*` 双组字段。
- `parse_summary_model(root, design_id, model)`：Protenix/OpenDDE（`*_summary_confidence_sample_*.json`）与 Boltz-2（`confidence_*.json` + `affinity_*.json`）的统一解析；提取约 25 种置信度字段（iptm、ptm、ranking_score、complex_plddt、chain_*、chain_pair_*、affinity_* 等）。

#### [select_manifest.py](../common/select_manifest.py) — 阶段间选样 CLI

- 两种互斥选法：`--ids/--ids-file`（显式 ID，找不到即报错）或 `--top-n --sort-field`（数值降序优先，文本次之，缺失值永远排最后）。
- `--require-success` 只保留 `status=success` 行；输出写 TSV + `selection.json` 审计记录（来源、数量、ID、时间戳）。

### 4.2 `level1/` — 高通量阶段

#### [run_level1.py](../level1/run_level1.py) — 单分片执行（批式）

`run_shard(run_dir, shard_index, *, resume)` 流程：

1. 读取 `config.json` 与分片任务；cid 模式先 `validate_cid_rows`。
2. 逐行构造 AF3 输入（direct 模式 `build_af3_json`，seed 由 design_id sha256 确定性导出 `_seed()`）；失败行记 `af3_status=missing` 并继续。
3. 写 ESMFold 批 FASTA（多链用 `:` 连接）；cid/for_wj 模式另写 AF3 FASTA（`_write_af3_fasta`，要求整分片共享同一链顺序）。
4. `_context()` 组装占位符上下文（`esm_input`、`af3_outdir`、`af3_csv` 等）。
5. **先 ESMFold 后 AF3**，各自：先解析已有输出判断是否需要运行 → `_run_stage`（`_default_command` 提供默认命令，`commands.*` 配置可覆盖）→ `run_logged` 执行 → 解析输出 → rc≠0 时把 success 降级为 `partial`。
6. 汇总状态（`status_from_models`）写 part TSV。

`_default_command(config, name, context)` 的默认命令逻辑：

- `esmfold`：`{esmfold_bin} -i {esm_input} -o {esm_outdir}` + 可选 `-m/--num-recycles/--max-tokens-per-batch/--chunk-size/--cpu-only/--cpu-offload`。
- `af3` + cid：`level1.cid.cid_command`（调用 `af3_single_node_v18.py`）。
- `af3` + for_wj_cid：返回 None（必须显式配置 `commands.af3`）。
- `af3` + direct：需要 `af3_script`/`af3_model_dir`/`af3_db_dir` 三项齐全才生成 `run_alphafold.py` 命令。

#### [cid.py](../level1/cid.py) — CID 契约与命令组装

- `validate_cid_rows(rows)`：GPU 作业启动前的完整预检——CID 五列必填（`design_id/sequence/chain_ids/cid_apo_json_path/cid_holo_json_path`）、序列仅允许 20 标准氨基酸+X、链 ID 唯一大写且与序列链数匹配、**整批只能有一对 apo/holo 模板与一种链顺序**、apo 模板必须无 ligand 实体而 holo 必须至少一个、模板蛋白 ID 必须与 `chain_ids` 一致。
- `cid_command(config)`：组装 `af3_single_node_v18.py` 调用（`--fasta/--template-with-lig/--template-without-lig/--af3-script/--model-dir/--db-dir/--monitor-script/--json-output-dir/--af3-output-dir/--csv-output/--archive-dir/...`）。

#### [aggregate_level1.py](../level1/aggregate_level1.py)

`join_results` 合并 → 写 `results/result_manifest.tsv`、`metrics/level1.csv`；cid/for_wj 模式额外写 `metrics/cid_apo.csv`、`metrics/cid_holo.csv` 与 `status/status.tsv`（状态行包含 `af3_apo`/`af3_holo`）。

#### [str20.py](../level1/str20.py) / [str20_closeout.py](../level1/str20_closeout.py) — STR 三态验证（独立工具）

- `str20.py`：从私有历史 CSV（`STR20_SOURCE_DIR`）选出 20 对 STR 蛋白，构造 apo/holo/quat 三态共 60 个 AF3 任务；单 GPU 串行队列执行（`--wait-idle` 等待 GPU 空闲 3 次、preflight 检查 CUDA/CCD、连续 3 次失败停队列、文件锁防并发），任务级校验归档。
- `str20_closeout.py`：`plan/verify/closeout` 三段式紧凑归档——以 sha256 清单为检查点，保留每任务 CIF+完整置信度与合并日志 `log/run.log`，删除 inputs/JAX cache；全程防 symlink、防路径逃逸、防文件变更竞态。

### 4.3 `level2/` — 多模型确认阶段

#### [run_level2.py](../level2/run_level2.py) — 单分片执行（逐设计）

`run_shard` 对每个设计、每个状态（普通=`(None,)`；`mode=cid` 时=`(apo, holo)`）、每个模型（protenix → boltz2 → opendde）调用 `_run_one`：

- `_build_input`：按模型构造输入（protenix→JSON、boltz2→YAML、opendde→无配体 FASTA/含配体 JSON）。
- resume：已有 success 解析结果 → 标 `*_reused=1` 跳过。
- `_default_command`：protenix=`bash {protenix_wrapper} IN OUT`；opendde=`bash {opendde_wrapper} -i IN -o OUT`；boltz2=`{boltz_bin} predict IN --out_dir OUT --write_full_pae --write_full_pde`（`boltz_no_kernels` 可追加）。`commands.*` 配置可完全覆盖（boltz2 兼容 `commands.boltz` 别名）。
- `_state_row`：CID 双态的轻量物化——apo 清空 `ligand`，holo 保留 `ligand`（缺失时从 `cid_holo_json_path` 提取第一个 CCD code）；模型侧 `design_id` 变为 `{design_id}_{state}`；结果字段前缀为 `{model}_{state}_*`，并记 `states=apo,holo`。

#### [aggregate_level2.py](../level2/aggregate_level2.py)

合并 part → `results/result_manifest.tsv`、`metrics/level2.csv`（CID 双态时模型名展开为 `{model}_apo`/`{model}_holo`）与 `status/status.tsv`。

### 4.4 `report/` — 报告生成

[generate_report.py](../report/generate_report.py)：`generate_report(run_dir, top_n, render_structures, renderer, resolution, aa_samples, pdf_name)` 主入口。**只读规范化产物** `metrics/level2.csv`（永不解析模型私有文件，因此归档后仍可重出报告）。

- `_prepare_rows`：仅保留 success 行；计算复合 `score`（pLDDT/100、pTM、ipTM 均值；单体的结构性 ipTM=0 不计入）；为 pae/ipsae/apae 做多别名匹配。
- `_rank_designs`：按 design 分组取分数均值排序，写 `report/ranking.tsv`（含 rank/selected 标记）。
- `_summary_stats`：`summary_stats.tsv`（逐指标分位数统计）与 `model_stats.tsv`（逐模型均值/中位数）。
- `_plot_stats`：matplotlib 以 **SVG 后端渲染 + ImageMagick `convert` 栅格化**（规避特定环境 Agg/Pillow 兼容问题，convert 不可用时回退保留 SVG）；产出 `metric_distributions.png`、`metric_scatter.png`、`model_heatmap.png`（design×model 复合分热图，热图用 Rectangle 逐格绘制以规避 QuadMesh 兼容问题，设计数封顶 40）、`model_bars.png`。
- `_render_top`：对 Top-N 调外部渲染脚本（`--renderer`，VMD/Tachyon 多视角面板），清理 `.tga/.dat/.tcl` 中间产物，保留 4 张视角图。
- `_build_pdf`：reportlab 组装 A4 PDF——标题、模型统计表、指标汇总表、Top-N 排名表、图集、4×6 结构图网格。

### 4.5 `scripts/` — 入口与机器配置

#### [profold.py](../scripts/profold.py) — v38.2 紧凑 facade（推荐入口)

四个子命令，`SUMMARY_FIELDS` 固定 11 列（`design_id, stage, status, structure_path, plddt, ptm, iptm, interface_pae, model_consensus, selection, error`）：

- `prepare`：`_canonical_rows` 把用户友好别名（`mode`、`template`→`target_spec_path`、cid 路径绝对化）在边界处一次性翻译 → 写规范化 `input.tsv` → 调 `common.pipeline.prepare_run` → 计算 sha256 指纹（input+config+level+VERSION）写 `fingerprint.json`；任一行 `mode=cid` 会自动把 config 提升为 `af3_mode=cid`。
- `run`：校验指纹（不符则拒绝 resume）→ 按分片数顺序子进程调用对应 `run_levelN.py`；`--no-resume` 透传。
- `aggregate`：子进程调用 `aggregate_levelN.py` → `_summary_rows` 从 result manifest 提取（level1 优先 AF3 后 ESMFold；level2 CID 双态时用 `{model}_holo` 指标）→ 写 `summary.tsv`。
- `report`：直接委托 `report.generate_report.generate_report`。

#### [configure.py](../scripts/configure.py) / [check_environment.py](../scripts/check_environment.py)

- `configure.py`：以 `--conda-root/--software-root/--data-root` 展开 `environments/site.example.json` 中的 `${V37_ROOT}/${CONDA_ROOT}/...` 模板 → 生成 `config.local.json`（拒绝覆盖已有输出）。
- `check_environment.py`：**只读**检查——控制器 Python 模块（numpy/Bio/gemmi/pandas/yaml）、内置工具文件存在性、配置中各模型可执行/路径存在性、CID 模式的 wrapper/monitor 路径；明确声明不验证 CUDA、权重与模型版本。

### 4.6 `tools/` — 模型接入工具

#### `tools/af3/`（运行于 AF3 conda 环境）

- [af3_multichain_core.py](../tools/af3/af3_multichain_core.py)：链映射与指标核心。`read_designs`（CID FASTA 解析，ID 安全校验）、`protein_ids`（模板蛋白链 ID）、`fill_template`（按 chain_ids 把序列填入模板，清空 MSA/templates）、`confidence_metrics(summary, detail, chains, cutoff, protein_chains)`——基于 `token_chain_ids`/`atom_chain_ids` 做 **token 对齐**的全套指标：整体 iptm/ptm/ranking_score、每链 plddt/PAE/iptm、每链对 iptm/ipAE（双向均值）/ipSAE（Dunbrack d0res 公式 `directional_ipsae`，仅蛋白-蛋白对）、overall 汇总。
- [af3_multichain_archive.py](../tools/af3/af3_multichain_archive.py)：`digest`（流式 sha256）、`atomic_bytes/atomic_json`（fsync+rename）、`process_task`（把 AF3 输出压缩为校验 bundle：CIF+summary+detail+input.json，各自记录 sha256；bundle 必须在 raw 目录外）、`load_bundle`（读 bundle 并重验所有 checksum）、`cleanup_bundle`（按清单逐文件校验后删除 raw，保留被改动的文件并记录 `retained_changed_files`）。
- [af3_single_node_v18.py](../tools/af3/af3_single_node_v18.py)：CID 一体化脚本——读 CID FASTA、由 apo/holo 模板批量生成双态 JSON（`batch_generate_json_cid`，name 形如 `{design}_af3_{variant}`）、运行 ESMFold 单体预测、调 AF3，并拉起 `af3_monitor_v16.py` 做增量 JSON 保存/输出归档/CSV 汇总。
- `af3_monitor_v16.py`：AF3 任务监控器（配合 v18 的增量处理与归档；由 `cid_command --monitor-script` 指定）。

#### `tools/monitors/multimodel_monitor_v1.py`

Protenix/Boltz-2/OpenDDE 的**统一置信度监控器**，列口径与 AF3 脚本对齐（便于多模型一致性验证）：全局 iptm/ptm/ranking_score/mean_plddt(0–100)/has_clash/disorder；Boltz-2 特有 confidence_score/complex_*/affinity_*；每链与每链对指标。用法：`--software {protenix|boltz2|opendde} --scan <root> -o xxx.csv`。

### 4.7 `wrappers/` — 外部模型包装

| 文件 | 职责 |
| --- | --- |
| [protenix_100.py](../wrappers/protenix_100.py) | Protenix 推理 runner（ByteDance 上游拷贝，`InferenceRunner`）；默认模型 `protenix_base_20250630_v1.0.0`，含 PyTorch 2.6 `weights_only` 兼容补丁 |
| [run_protenix.sh](../wrappers/run_protenix.sh) | 环境变量驱动（`PROTENIX_PYTHON/PROTENIX_SOURCE/PROTENIX_SCRIPT/PROTENIX_MODEL/PROTENIX_SEED/...`）；固定 `--use_msa false --need_atom_confidence true`、5 diffusion samples |
| [run_opendde.sh](../wrappers/run_opendde.sh) | 激活 conda 环境 → FASTA 自动转 JSON → `opendde pred`（默认无 MSA 快速模式；`--msa/--template` 需 search_database） |
| [fasta_to_opendde_json.py](../wrappers/fasta_to_opendde_json.py) | FASTA→OpenDDE JSON：`:` 分隔多链、`--auto-type` 识别 DNA/RNA、`--ligand CCD`/`--smiles` 添加配体 |
| [opendde_predict.sh](../wrappers/opendde_predict.sh) | 独立 OpenDDE 推理包装（可调 seeds/cycles/steps/dtype/model） |

> 这些脚本内含集群绝对路径默认值（如 `/xcfhome/...`、`/pubhome/...`），跨机器部署必须通过配置项或环境变量覆盖。

### 4.8 `structure_align/` — 通用结构比对

[align_structures.py](../structure_align/align_structures.py)：独立于两个 level 的后处理工具，解决"一个骨架对应多条设计序列"的结构归属与对齐问题。输入是普通 CSV（`design_id, sequence_id, backbone_structure_path, predicted_structure_path` + 可选 `chain_map`/motif 字段），输出三张审计表：

- `alignment_summary.csv`：每 design 一行，复合物级 USalign 指标（RMSD、以骨架长度归一的 TM-score）+ 配体重原子 RMSD；
- `chain_pairs.csv`：每个骨架链×预测链对一行（蛋白/核酸/配体）；
- `motifs.csv`：指定残基/原子的 paired-atom Kabsch RMSD。

关键函数：`auto_chain_mapping`（按实体类型+序列相似度+长度自动配链）、`run_usalign`（多线程调用 USalign 并解析旋转矩阵）、`_kabsch`（局部刚体拟合）、`align_one_row`/`align_rows`/`main`。依赖 gemmi + numpy + 外部 USalign。

### 4.9 `tests/` — 回归套件

| 测试文件 | 覆盖范围 |
| --- | --- |
| `test_facade.py` | facade 别名翻译、指纹、summary 固定列 |
| `test_input_builder.py` | 多链解析、AF3 模板填充、无 MSA 校验 |
| `test_model_parsers.py` | 解析器对缺失/畸形输出的防御行为 |
| `test_level1_modes.py` | direct/cid/for_wj 命令与状态 |
| `test_level2_pipeline.py` | prepare→run→aggregate 全链路（fake wrapper） |
| `test_cid_mode.py` / `test_cid_multichain.py` | CID 契约校验与双态结果字段 |
| `test_str20.py` | STR20 closeout 校验逻辑 |
| `test_report.py` | 报告指标读取与排名 |
| `test_portability.py` | 无集群路径环境下的可移植性 |
| `fixtures/fake_wrappers.py` | 假模型 wrapper（测试中替代真实 GPU 模型） |

CI（[.github/workflows/test.yml](../.github/workflows/test.yml)）：Python 3.11 安装 requirements → `unittest discover` → 所有 `.sh` 做 `bash -n` 语法检查。

---

## 5. 关键类与函数索引

### 5.1 common 层（最重要的公共 API）

| 函数 | 位置 | 说明 |
| --- | --- | --- |
| `read_manifest(path, require_sequence)` | [manifest.py:131](../common/manifest.py#L131) | 读清单并全面校验（唯一 ID、列对齐、序列） |
| `write_tsv/write_csv(path, rows, preferred)` | [manifest.py:174](../common/manifest.py#L174) | 原子写 TSV/CSV，稳定列序 |
| `row_fingerprint(row)` | [manifest.py:218](../common/manifest.py#L218) | 行内容 sha256 |
| `prepare_run(input, outdir, level, shard_size, config)` | [pipeline.py:17](../common/pipeline.py#L17) | 创建 run 目录、切片、快照 config.json |
| `join_results(run_dir, model_names)` | [pipeline.py:95](../common/pipeline.py#L95) | part ↔ task manifest 连接（duplicate ID 即错） |
| `metric_rows(rows, model_names)` | [pipeline.py:157](../common/pipeline.py#L157) | 前缀字段 → 长 表（每设计×模型一行） |
| `command_template(value, context)` | [runner.py:39](../common/runner.py#L39) | `{placeholder}` 渲染为 argv |
| `run_logged(command, log_path)` | [runner.py:61](../common/runner.py#L61) | 子进程执行 + 全量日志 + returncode |
| `stage_result(model, status, ...)` | [runner.py:89](../common/runner.py#L89) | `{model}_*` 前缀化结果字典构造器 |
| `protein_records(row)` | [input_builder.py:27](../common/input_builder.py#L27) | 序列/FASTA → 有序 (chain_id, seq) 列表 |
| `build_af3_json / build_protenix_json / build_boltz_yaml / build_opendde_fasta` | [input_builder.py](../common/input_builder.py) | 四个模型输入构造器（无 MSA 约束） |
| `status_from_models(statuses)` | [input_builder.py:397](../common/input_builder.py#L397) | 多模型状态归并 |
| `parse_esmfold / parse_af3 / parse_summary_model` | [model_parsers.py](../common/model_parsers.py) | 三类模型输出的防御式解析 |
| `unavailable(model, reason)` | [model_parsers.py:483](../common/model_parsers.py#L483) | 统一 `missing` 结果构造 |

### 5.2 worker 层

| 函数 | 位置 | 说明 |
| --- | --- | --- |
| `run_shard(run_dir, shard_index, resume)` | [run_level1.py:179](../level1/run_level1.py#L179) | Level 1 分片：批式 ESMFold→AF3 |
| `_default_command(config, name, context)` | [run_level1.py:73](../level1/run_level1.py#L73) | Level 1 默认命令（esmfold/af3 三模式） |
| `validate_cid_rows(rows)` | [cid.py:12](../level1/cid.py#L12) | CID 五列契约 + apo/holo 模板校验 |
| `cid_command(config)` | [cid.py:61](../level1/cid.py#L61) | AF3 CID wrapper 命令组装 |
| `run_shard(...)` | [run_level2.py:183](../level2/run_level2.py#L183) | Level 2 分片：逐设计×状态×模型 |
| `_run_one(row, model, ..., state)` | [run_level2.py:113](../level2/run_level2.py#L113) | 单模型单态执行（构造输入→运行→解析→resume） |
| `_state_row(row, state)` | [run_level2.py:83](../level2/run_level2.py#L83) | CID apo/holo 行物化（ligand 差异 + ID 后缀） |

### 5.3 入口与报告

| 函数 | 位置 | 说明 |
| --- | --- | --- |
| `prepare/run/aggregate/report (args)` | [profold.py](../scripts/profold.py#L81) | facade 四个子命令实现 |
| `_canonical_rows(path)` | [profold.py:49](../scripts/profold.py#L49) | 用户别名→内部字段的一次性翻译 |
| `_fingerprint(input, config, level)` | [profold.py:72](../scripts/profold.py#L72) | resume 指纹（含 VERSION） |
| `_summary_rows(run_dir, level)` | [profold.py:126](../scripts/profold.py#L126) | result manifest → 固定 11 列 summary |
| `generate_report(run_dir, ...)` | [generate_report.py:400](../report/generate_report.py#L400) | 报告主入口（统计→排名→绘图→渲染→PDF） |
| `_score(row)` | [generate_report.py:69](../report/generate_report.py#L69) | 复合分数（容忍单体 ipTM=0） |
| `align_one_row(row, ...)` | [align_structures.py:732](../structure_align/align_structures.py#L732) | 单 design 的骨架对齐全流程 |

### 5.4 AF3 工具层

| 函数 | 位置 | 说明 |
| --- | --- | --- |
| `fill_template(template, name, sequences, chain_ids, seed)` | [af3_multichain_core.py:51](../tools/af3/af3_multichain_core.py#L51) | 按链 ID 填充 AF3 JSON 模板 |
| `confidence_metrics(summary, detail, chains, cutoff, protein_chains)` | [af3_multichain_core.py:98](../tools/af3/af3_multichain_core.py#L98) | token 对齐全套置信度指标 |
| `directional_ipsae(block, cutoff)` | [af3_multichain_core.py:85](../tools/af3/af3_multichain_core.py#L85) | Dunbrack ipSAE d0res 方向性计算 |
| `process_task(folder, ..., archive_dir, ...)` | [af3_multichain_archive.py:48](../tools/af3/af3_multichain_archive.py#L48) | AF3 输出 → 校验归档 bundle |
| `cleanup_bundle(bundle, raw_root)` | [af3_multichain_archive.py:122](../tools/af3/af3_multichain_archive.py#L122) | 按清单清理 raw（防误删变更文件） |
| `plan/verify/closeout(root)` | [str20_closeout.py](../level1/str20_closeout.py#L30) | STR 运行三段式紧凑归档 |

---

## 6. 数据流与运行目录布局

### 6.1 端到端数据流

```
designs.tsv ──prepare──▶ run 目录 (tasks/*.tsv, config.json)
                │
                ├─ run (shard 0..N) ──▶ 构造模型输入 → 调外部模型
                │                        → 解析输出 → results/parts/part_XXXXX.tsv
                │
                ├─ aggregate ──▶ results/result_manifest.tsv
                │                metrics/levelN.csv (+ cid_*.csv)
                │                status/status.tsv
                │
                ├─ [facade] aggregate ──▶ summary.tsv (固定 11 列)
                │
                └─ select_manifest ──▶ 下一阶段的输入 TSV (+ selection.json)

Level 2 聚合后:
   metrics/level2.csv ──report──▶ report/{PDF, ranking.tsv, *.png, structures/}
```

### 6.2 标准 run 目录

```
<run>/
├── config.json              # 输入/config/level/shard_size/n_tasks/n_shards 快照
├── fingerprint.json         # facade: resume 校验指纹
├── input.tsv                # facade: 规范化后的边界 manifest
├── input/input_manifest.tsv # 原始输入快照
├── tasks/
│   ├── task_manifest.tsv    # 全部任务（task_index/shard_index/input_fingerprint）
│   ├── shard_00000.tsv ...  # 每分片任务
│   └── n_shards.txt
├── results/
│   ├── parts/part_00000.tsv # 每分片 worker 结果
│   └── result_manifest.tsv  # 聚合后的完整结果（下游选择输入）
├── metrics/
│   ├── level1.csv / level2.csv      # 长表指标
│   ├── cid_apo.csv / cid_holo.csv   # CID 双态（level1）
│   └── stage_timings.tsv            # 每分片预测器墙钟时间
├── status/status.tsv        # 逐设计×逐模型状态审计
├── artifacts/shard_XXXXX/   # 模型输入/输出（esmfold/、af3/、models/、inputs/）
│   └── stage_timings.json
├── logs/                    # 逐分片或逐设计的模型日志
├── summary.tsv              # facade: 用户主输出
└── report/                  # report 命令产物
```

### 6.3 状态语义

| 状态 | 含义 |
| --- | --- |
| `success` | summary 与结构文件均解析成功 |
| `partial` | 部分产物缺失，或命令 rc≠0 但已有可解析输出 |
| `failed` | 输出目录存在但无可识别产物 |
| `missing` | 无输出目录/命令未配置/输入构造失败 |
| `{model}_reused=1` | resume 时复用已有成功结果，未重新预测 |

---

## 7. 配置系统

### 7.1 配置文件层次

```
environments/site.example.json   # 站点模板（${CONDA_ROOT} 等占位符）
        │ scripts/configure.py 展开
        ▼
config.local.json                # 本机实际配置（不入库）
config.example.json              # direct 模式示例 / config.cid.example.json  # CID 示例
        │ prepare 时快照
        ▼
<run>/config.json                # 每次运行的不可变配置快照
```

### 7.2 关键配置键

**Level 1（[config.example.json](../config.example.json) / [config.cid.example.json](../config.cid.example.json)）**

| 键 | 说明 |
| --- | --- |
| `esmfold_bin` | esm-fold CLI 路径 |
| `esmfold_num_recycles` / `esmfold_max_tokens_per_batch` / `esmfold_chunk_size` / `esmfold_cpu_only` / `esmfold_cpu_offload` | ESMFold 调参 |
| `af3_mode` | `direct` \| `cid` \| `for_wj_cid` |
| `af3_python` / `af3_script` / `af3_model_dir` / `af3_db_dir` | AF3 环境与权重/数据库（direct 模式必需三项齐全） |
| `cid_wrapper` / `cid_monitor` / `cid_pae_cutoff` / `cid_keep_source` | CID wrapper 覆盖项（默认用 tools/af3 内置脚本） |
| `commands.{esmfold,af3,...}` | 完全自定义命令（占位符见 worker `_context`），优先级最高 |

**Level 2**

| 键 | 说明 |
| --- | --- |
| `protenix_wrapper` | Protenix 包装脚本（默认指向 wrappers 逻辑/集群路径） |
| `boltz_bin` / `boltz_no_kernels` | Boltz-2 可执行与 kernel 开关 |
| `opendde_wrapper` | OpenDDE 包装脚本 |
| `level2_mode` | `single` \| `cid`（cid 时每个模型跑 apo+holo） |
| `level2_models` | 模型清单声明（当前 runner 固定执行 protenix/boltz2/opendde） |
| `usalign_path` / `usalign_workers` / `rosetta_filter_script` | 结构比对与 Rosetta 过滤（报告/交付用） |
| `conda_env` / `run_script` | Slurm 提交所用环境与脚本 |

### 7.3 占位符约定

- `${...}`：**configure 阶段**展开的站点路径（site 模板专用）。
- `{...}`：**worker 运行时**展开的占位符（如 `{esm_input}`、`{af3_outdir}`、`{protenix_input}`、`{boltz_outdir}`、`{opendde_input}` 等，定义在各 worker 的 `_context`）。不要把未展开的 site 模板直接传给 prepare。

---

## 8. 依赖关系

### 8.1 模块依赖图（仓库内）

```
scripts/profold.py ──▶ common.pipeline / common.manifest / report.generate_report
level1/run_level1.py ──▶ common.{input_builder,manifest,model_parsers,pipeline,runner} + level1.cid
level1/cid.py ──▶ tools/af3 (脚本路径引用)
level1/str20.py ──▶ tools.af3.{multichain_core,multichain_archive} + level1.str20_closeout
level2/run_level2.py ──▶ common.{input_builder,model_parsers,pipeline,runner}
report/generate_report.py ──▶ (独立; 读 metrics/level2.csv; 可选 matplotlib/reportlab/PIL)
tools/af3/af3_single_node_v18.py ──▶ tools.af3.{multichain_core,multichain_archive}
level1/str20_closeout.py ──▶ tools.af3.af3_multichain_archive
tools/af3/af3_multichain_archive.py ──▶ tools.af3.af3_multichain_core
structure_align/align_structures.py ──▶ (独立; gemmi+numpy+USalign)
```

依赖方向原则：`level* → common → (无仓库内依赖)`；`tools/af3` 自成一体（在 AF3 环境内运行）；`report`、`structure_align` 与 level 无反向耦合。

### 8.2 Python 依赖（[requirements.txt](../requirements.txt)）

**控制器（必需）**：`numpy>=1.24`、`biopython>=1.81`、`gemmi>=0.6`、`pandas>=2`、`PyYAML>=6`。
**报告阶段（可选）**：`matplotlib>=3.7`、`reportlab>=4`、`Pillow>=10`；另需系统 ImageMagick `convert`（SVG→PNG，缺失时保留 SVG 回退）。
**语言**：Python ≥ 3.10（代码使用 `str | Path` 原生联合类型语法）；CI 用 3.11。
**结构比对**：gemmi + numpy + 外部 [USalign](https://github.com/Terazus/USalign) 二进制。

### 8.3 外部模型环境（仓库外，按 environments/README.md 布局）

| 模型 | 调用契约 | 接入方式 |
| --- | --- | --- |
| ESMFold | `esm-fold -i FASTA -o PDB_DIR` | `esmfold_bin` 或 `commands.esmfold` |
| AlphaFold 3 | `run_alphafold.py --json_path ... --model_dir ... --db_dir ...`（v3.0.4 验证过） | `af3_python/af3_script/af3_model_dir/af3_db_dir` |
| Protenix | `protenix_100.py --input_json_path ... --dump_dir ...`（PyTorch 环境） | `protenix_wrapper` + `PROTENIX_*` 环境变量 |
| Boltz-2 | `boltz predict IN --out_dir OUT --write_full_pae --write_full_pde` | `boltz_bin` |
| OpenDDE | `opendde pred -i JSON -o OUT`（conda 环境 `opendde`） | `opendde_wrapper` |
| VMD/Tachyon 渲染器 | 外部 `render_structure_panel.py`（报告 Top-N 面板） | `--renderer` 参数 |

注意：CID wrapper/monitor 运行于 **AF3 Python** 环境（额外需要 numpy/biopython/gemmi）；控制器逻辑运行于控制器 Python。两条 Python 链路分离。

---

## 9. 项目运行方式

### 9.1 安装与配置

```bash
git clone https://github.com/LiuSantu123/profold.git
cd profold
python -m pip install -r requirements.txt

# 方式 A：从 site 模板生成本机配置
python scripts/configure.py \
  --conda-root /opt/miniconda3 --software-root /opt/software \
  --data-root /data/models --output config.local.json

# 方式 B：复制示例配置后手改路径
cp config.example.json config.local.json

# 只读环境检查（不执行模型、不下载权重）
python scripts/check_environment.py --config config.local.json --level level1
python scripts/check_environment.py --config config.local.json --level level2
```

模型环境（esmfold/af3/protenix/boltz/opendde）按 [environments/README.md](../environments/README.md) 单独安装；每个模型先跑一条真实 GPU smoke 再放量。

### 9.2 推荐路径：v38.2 紧凑 facade

```bash
# 1) 准备：designs.tsv + config.local.json → run 目录
python scripts/profold.py prepare \
  --input designs.tsv --config config.local.json \
  --level level1 --outdir runs/l1 --shard-size 32

# 2) 运行：顺序执行全部分片（可断点续跑）
python scripts/profold.py run --run-dir runs/l1 --level level1
#    强制重算: --no-resume

# 3) 聚合：结果清单 + 指标 + 状态 + summary.tsv
python scripts/profold.py aggregate --run-dir runs/l1 --level level1

# 4) 选样进入 Level 2
python common/select_manifest.py \
  --result-manifest runs/l1/results/result_manifest.tsv \
  --top-n 100 --sort-field af3_ranking_score --require-success \
  --output runs/l1_selected.tsv

# 5) Level 2 同样四步（prepare/run/aggregate 换 --level level2）

# 6) 报告（Level 2）：统计图 + Top-N 结构面板 + PDF
python scripts/profold.py report --run-dir runs/level2 --top-n 10
#    只要统计不要结构渲染: --no-structures
```

### 9.3 旧版直接调用（各 level 独立脚本）

```bash
# Level 1
python level1/prepare_level1.py --manifest input.tsv --outdir runs/level1 \
  --shard-size 32 --config config.local.json
python level1/run_level1.py --run-dir runs/level1 --shard-index 0
python level1/aggregate_level1.py --run-dir runs/level1

# Level 2（CID 双态：manifest 行加 mode=cid 与状态特异输入路径）
python level2/prepare_level2.py --manifest selected.tsv --outdir runs/level2 \
  --shard-size 8 --config config.local.json
python level2/run_level2.py --run-dir runs/level2 --shard-index 0
python level2/aggregate_level2.py --run-dir runs/level2

# 报告
python report/generate_report.py --run-dir runs/level2 --top-n 10
```

### 9.4 Slurm 提交

```bash
# prepare + array(0..N-1) + 依赖式聚合 一次提交
PYTHON_BIN=/opt/miniconda3/envs/profold/bin/python \
bash level1/submit_level1.sh --manifest input.tsv --outdir runs/level1 \
  --shard-size 32 --config config.local.json --partition gpu --time 24:00:00
bash level2/submit_level2.sh --manifest selected.tsv --outdir runs/level2 ...

# 也可以手动提交 array
sbatch --array=0-N level1/run_level1_array.sh runs/level1
```

- array 任务 ID ↔ shard index；`run_levelN_array.sh` 每任务申请 1 GPU/8 CPU。
- 重复提交 worker 即可恢复（complete 结果自动跳过）。

### 9.5 独立工具

```bash
# 结构-骨架对齐（独立于 level）
python structure_align/align_structures.py --input-csv align_input.csv \
  --outdir runs/align --workers 4 --usalign /path/to/USalign

# 多模型置信度监控（模型环境内）
python tools/monitors/multimodel_monitor_v1.py --software protenix --scan <root> -o protenix.csv

# OpenDDE FASTA→JSON 转换
python wrappers/fasta_to_opendde_json.py input.fasta -o out.json --ligand ATP

# STR20 三态验证（需要私有历史数据 STR20_SOURCE_DIR）
python level1/str20.py prepare --root runs/str20
python level1/str20.py run --root runs/str20 --gpu 1 --wait-idle
python level1/str20.py closeout --root runs/str20
```

### 9.6 测试与发布

```bash
python -m unittest discover -s tests -v          # 回归套件（无 GPU 依赖）
bash -n level1/*.sh level2/*.sh wrappers/*.sh    # shell 语法（CI 同款）
bash scripts/publish_release.sh LiuSantu123/profold v0.0.38.2   # 需已认证 gh + 干净 worktree
```

---

## 10. 测试与发布

- **测试策略**：所有测试不依赖 GPU/真实模型——`tests/fixtures/fake_wrappers.py` 提供假模型 wrapper，用真实流水线代码路径跑通 prepare→run→aggregate；解析器测试聚焦缺失/畸形输入的防御行为；`test_portability.py` 保证无集群路径环境可运行。
- **CI**：GitHub Actions（`Controller Tests`）：Python 3.11 + requirements + unittest + 全部 shell 脚本 `bash -n`。
- **发布**：`scripts/publish_release.sh <repo> <tag>` 创建带 tag 的源码归档与 SHA256 校验；发布包**不含**权重、数据库、本地配置、预测结果与 STR20 私有数据。

---

## 11. 设计原则与注意事项

1. **manifest 是唯一接口**：跨阶段只传 TSV；`design_id` 全局唯一且正则受限；额外自定义列会原样保留到结果 manifest。
2. **显式优于隐式**：阶段间选择必须走 `select_manifest.py`；模型命令优先级为 `commands.* > 内置默认`；缺配置宁可 `missing` 也不猜路径。
3. **不伪造数据**：解析器对缺失细节文件记 `partial`；缺分数留空；小分子不硬套 USalign TM-score。
4. **原子性与可恢复**：所有清单/JSON 写入都是 tmp+rename 原子操作；resume 基于输出解析而非日志；facade 指纹校验防止"输入变了还用旧结果"。
5. **归档即证据**：AF3 bundle/STR closeout 全程 sha256 校验、防 symlink/路径逃逸；closeout 清单不可手改。
6. **跨机器移植注意**：
   - 代码中存在集群默认绝对路径（如 [run_level1.py 的 `ESMFOLD_BIN`](../level1/run_level1.py#L23)、wrappers 内 `/xcfhome/...`、报告默认 `--renderer` 路径），部署时必须通过 `config.local.json` / CLI 参数覆盖。
   - README 提及的 `level2/compact_level2.py` 与 `scripts/submit_campaign.py` 属于部署侧可选脚本，**当前仓库快照中不存在**；`level2_models` 配置中的 `esmfold2` 也是预留项（未接入 runner）。
7. **状态≠生物学结论**：`status` 只描述执行与解析；生物学判断需单独检查 pLDDT/pTM/ipTM/PAE/ipSAE/USalign 等指标。
