# ProFold 简洁输入输出优化任务清单

更新时间：2026-10-02

## 目标

把 ProFold 的用户接口收敛为“一个设计输入表、一个运行配置、一个主结果表”。
内部模型输入、原始结果、日志和审计指标继续保留，但不再要求用户理解模型专属字段。

## 目标接口

### 输入

主输入 `designs.tsv`：

```text
design_id  sequence  mode  template  ligand  reference_structure
```

- `design_id`、`sequence` 必填。
- `mode` 取 `protein` 或 `cid`。
- `template`、`ligand`、`reference_structure` 按模式可选；CID 的 apo/holo 模板由一个清晰的配置对象表达，不暴露四个以上路径字段。
- 保留旧 manifest 字段的兼容读取，但新入口统一转换为内部标准 manifest。

运行配置 `config.local.json` 只描述机器和模型环境：模型命令、权重、数据库、Slurm 资源和输出保留策略，不写设计内容。

### 输出

用户主结果 `summary.tsv` 每个设计一行，固定字段：

```text
design_id  stage  status  structure_path  plddt  ptm  iptm  interface_pae  model_consensus  selection  error
```

内部结果另存为 `metrics/`，原始文件另存为 `raw/`，日志另存为 `logs/`。主表不暴露模型专属路径和长字段名。

推荐目录：

```text
run/
├── input.tsv
├── config.snapshot.json
├── summary.tsv
├── selected.tsv
├── structures/
├── metrics/
├── raw/
└── logs/
```

## 实施顺序

- [x] P0：定义 `designs.tsv` v1 schema、字段校验和错误信息。
- [x] P0：新增 `scripts/profold.py` facade，负责 prepare、run、aggregate 的统一调用。
- [x] P0：生成固定字段的 `summary.tsv`，把详细模型字段下沉到 `metrics/`。
- [x] P0：为每个输入、配置、版本生成 `fingerprint.json`；后续补强 resume 拒绝策略。
- [x] P1：统一 `protein`、`cid` 两种模式的结构、状态和选择字段（facade 已覆盖）。
- [x] P1：把 Protenix、Boltz-2、OpenDDE 的结构与 confidence 按 design/state/sample 成套绑定（单体、小分子真实 STR、PPI、CID apo/holo 均已在 f101 RTX 3090 实跑）。
- [ ] P1：增加 `selected.tsv` 和 `selection.json`，记录筛选规则、时间、版本和来源 summary。
- [ ] P1：实现统一 cleanup/closeout：先验证 summary 和结构完整，再移动 raw，不删除审计证据。
- [ ] P2：接入实时 progress/status，区分 queued、running、success、partial、failed、blocked。
- [ ] P2：统一多模型指标口径，生成 `model_consensus` 和可解释的筛选原因。
- [x] P2：加入 Level 2 统计报告：全量均值/分位数汇总 TSV、分布图、散点图、限幅设计×模型热图、模型均值柱形图、Top-N 等比例多视角结构面板和 PDF。
- [x] P2：更新 README、使用指南、示例和版本号，补充 f101 Level 2 配置与验收记录。
- [x] P2：补充 f101 RTX 3090 四类测试输入目录和本机配置模板。

## f101 RTX 3090 验收记录（2026-10-02）

- `test/monomer`：ESMFold + AF3，1/1 success。
- `test/small_molecule`：含 STR 配体的 AF3 模板，ESMFold + AF3，1/1 success。
- `test/ppi`：A/B 双链 AF3 模板，ESMFold + AF3，1/1 success。
- `test/cid`：AF3 apo/without_lig 与 holo/with_lig 两状态均 success，ESMFold 亦 success。
- 四类均生成固定 11 列 `run/summary.tsv`；3090 实测显存峰值约 18.4 GiB。
- CID 运行发现 AF3 环境缺少 Biopython；已将归档/监控中的 Biopython 解析改为可选，AF3 核心运行不再依赖该包。
- 该记录覆盖 Level 1；Level 2 三模型 GPU 验收随后已完成四类场景，详见下方记录。

## Level 2 f101 RTX 3090 验收记录（2026-10-03）

- `test/level2_monomer`：Protenix、Boltz-2、OpenDDE，3/3 success。
- `test/level2_small_molecule`：三模型，3/3 success；三套输入均包含真实 `CCD_STR`。
- `test/level2_ppi`：三模型，3/3 success。
- `test/level2_cid`：apo/holo 两状态各运行三模型，6/6 success。
- Protenix 使用仓库内 `wrappers/run_protenix_f101.sh`，固定 f101 的 Protenix Python 环境，解决旧 wrapper 的 conda 初始化问题。
- USalign `/xcfhome/yhliu/002_software/009_TMalign/USalign` 已实测可执行；Rosetta 脚本根目录固定为 `/xcfhome/yhliu/003_scripts/08_rosetta`。该目录尚无符合 ProFold JSON 过滤器契约的统一 `rosetta_filter_script`，因此暂不伪装成已接入。

## 验收标准

- 新用户只需准备一个 TSV 和一个本机配置即可启动完整 Level 1/Level 2。
- `summary.tsv` 不包含动态模型专属列，列顺序固定，可直接用 Excel/Pandas 查看。
- 从 `summary.tsv` 生成的 `selected.tsv` 可直接作为下一阶段输入。
- 改变序列、模板、模型版本或关键参数后，resume 自动拒绝旧结果。
- 中断、部分失败和重复运行不会覆盖已有可审计结果。
- 文档中的命令、文件名和实际代码全部一致。
