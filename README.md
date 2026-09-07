# v37 三层蛋白筛选 Pipeline

跨机器安装优先阅读 [环境配置](environments/README.md)，使用 `scripts/configure.py`
生成 `config.local.json`，再运行 `scripts/check_environment.py` 检查对应 level。
AF3 多链模块与 CID wrapper 已随仓库放在 `tools/af3/`，不再依赖外部 for_wj 目录。
[发布说明](RELEASE_NOTES.md)列出验证范围；模型权重、数据库和运行结果不在发布包内。
原 `config*.example.json` 和下文 `/xcfhome` 路径仅是历史机器参考，新机器使用环境模板。

GitHub 发布：认证 `gh auth login --hostname github.com`，创建或指定自己的仓库，提交代码后
运行 `bash scripts/publish_release.sh OWNER/REPO v37.0.0`。脚本校验测试和 tag，推送 main/tag，
创建包含源码包与 SHA256 校验文件的 release；不会覆盖已有 tag。不同机器 clone 后
使用各自的 `config.local.json`，代码修改通过正常分支、commit、push/pull 同步。

v37 将筛选拆成三个相互独立的 pipeline。每层都接受 TSV manifest，生成自己的 run directory、分片任务、结果 manifest、指标表和状态表；用户可以自行决定多少设计进入下一层，也可以用 Slurm `afterok` 将多个 level 串联起来。结构预测按 no-MSA 路线配置，具体由各模型 wrapper 和输入文件保证。

## 三个 level

| Level | 做什么 | 模型/分析 | 默认分片粒度 |
| --- | --- | --- | ---: |
| Level1 | 高通量第一轮结构筛选；快速 foldback，并预测复合物结构和置信度 | ESMFold + AlphaFold3 | 32 designs/shard |
| Level2 | 对用户选定的子集补充多模型结构预测，比较不同模型的一致性 | Protenix + Boltz-2 + OpenDDE | 8 designs/shard |
| Level3 | 对用户选定的子集做序列层面的可溶性和热稳定性分析 | NetSolP + TemBERTure | 64 designs/shard |

### Level1：ESMFold + AF3

ESMFold 用单序列快速预测设计结构，适合大规模初筛；AF3 用于复合物结构和置信度评估。AF3 可使用 `direct` 单次预测模式，或 `cid` 专属 apo/holo 双预测模式；`for_wj_cid` 保留为旧配置兼容入口。CID 支持任意蛋白链数，两种状态分别计算整体、逐链和逐对指标。

ESMFold 默认复用原 `for_wj` 的 `/xcfhome/yzmeng/miniconda3/envs/zb/bin/esm-fold -i FASTA -o PDB_DIR`，入口 shebang 已指定 Python 环境，无需先激活 conda。原脚本中的 `/xcfhome/bozhang/anaconda3/envs/protein_design` 用于 OmegaFold，该环境目前没有 `bin/esm-fold`。配置可设置 `esmfold_bin`、`esmfold_model_dir`、`esmfold_num_recycles`、`esmfold_max_tokens_per_batch`、`esmfold_chunk_size`、`esmfold_cpu_only`、`esmfold_cpu_offload`。默认模式读取精确的 `<design_id>.pdb`，校验坐标并提取原子平均 pLDDT；不要求 CSV，不复用旧 ESMFold2 CIF。多链 FASTA 统一为冒号分隔。已有 run 若显式配置 `commands.esmfold`，仍优先执行该命令，需自行核对配置；旧 `esmfold_model/esmfold_script` 等 ESMFold2 参数不再控制默认入口。

Level1 worker 的 Python 需要 `gemmi` 解析 PDB，可使用现有 base 环境（与 `structure_align` 相同依赖）。2026-09-07：base 环境 `python -m unittest discover -s tests -v` 通过 35 项测试；STR20 专用入口已完成 60/60 AF3 真实预测及收尾验证，ESMFold GPU 推理与完整 Level1 E2E 尚未验证。

### STR20 三态验证与自动收尾

`level1/str20.py` 是历史 STR CID 的专用验证入口：20 对蛋白分别预测 apo（A+B）、
holo（A+B+STR）、quat（A+B+STR+A8W），共60任务。该入口只运行AF3，不运行ESMFold；
AF3逐项使用 `--json_path`，不是一次 `--input_dir`。运行期间保留逐任务日志和JAX缓存。

全任务成功后自动调用 `level1/str20_closeout.py`，完成后的目录为：

```text
runs/str20/
  archive/<task>/<task>.cif
  archive/<task>/<task>_confidences.json
  log/run.log
  log/metrics/{all,apo,holo,quat}.csv
  log/closeout.json
  log/                     # 清单、设置、模板、状态等追溯记录
```

所有任务输入JSON和本run的JAX缓存自动删除；每任务归档仅留上述两个文件。
原任务、预检和队列日志按来源标记合并为 `log/run.log`。先原子写入集中校验清单，
校验结构、置信度、指标和日志后才清理；失败或未完成任务保留输入、缓存及诊断。
收尾中断可续接；精简归档的resume校验SHA256，不初始化CUDA，也不会重跑预测。

从v37目录执行：

```bash
# 运行或恢复；全成功后自动收尾，已经精简的结果只做校验。
python level1/str20.py run --root runs/str20 --gpu 0 --wait-idle
# 只对已完成的旧运行补做收尾；未完成会拒绝。
python level1/str20.py closeout --root runs/str20
```

已有STR20于2026-09-07 17:15:21在f101完成60/60，0失败；62份日志合并后，
总文件数由1283降至139。这套全run精简规则目前接入STR20入口，其他Level和通用CID
wrapper仍沿用各自归档格式，不能直接删除其summary/result/input文件。

### Level2：Protenix + Boltz-2 + OpenDDE

Level2 不要求所有 Level1 设计都进入，只运行 manifest 中实际传入的设计。三个模型按“一个 design、一个模型”独立执行，因此某个模型失败时可以单独恢复；输出包括结构、summary confidence，以及各模型可提供的 pLDDT、pTM、ipTM、PAE/PDE 和 Boltz-2 affinity 等指标。复杂配体、核酸或多实体体系应提供模型专用输入文件。

### Level3：可溶性 + 热稳定性

NetSolP 预测可溶性/可用性相关指标，TemBERTure 预测熔解温度或热稳定性分数。两者都是序列模型，不会从 Level1/Level2 结构中自动推导新序列；Level3 评分的是 manifest 中的 `level3_sequence` 或明确指定的序列。

三个 level 不依赖隐藏的上游状态，只通过 manifest 传递设计 ID、序列、输入路径和已产生的结果。因此可以单独运行，也可以按 `Level1 -> 人工选择 -> Level2 -> 人工选择 -> Level3` 串联。

## 成功、失败与筛选通过

这里的状态表示“运行和结果解析是否完整”，不表示 pLDDT、ipTM 或可溶性已经达到筛选阈值。科学筛选仍通过指标表和 `select_manifest.py` 单独完成。

单模型状态的基本口径：

- `success`：必需输出存在且能解析。Level1/2 要有结构和 summary；Level3 要找到对应设计行且至少有一个有效指标。
- `partial`：有部分可用结果，例如只生成了 summary 或结构，或命令返回非零但仍留下可解析结果。
- `failed`：输出目录存在，但没有找到可识别的有效结果，或设计 ID 不在结果表中。
- `missing`：输出目录/命令未提供，任务尚未产生模型结果。

一个 design 的 level 状态按该层所有模型合并：全部 `success` 才是 `success`；存在 `success`/`partial` 但未全部完成时是 `partial`；没有任何可用结果且包含失败时是 `failed`；全部模型都是 `missing` 时是 `missing`。`for_wj_cid` 的 AF3 只有 with-lig 和 without-lig 两套都成功时，AF3 总状态才是 `success`。

## 结构比对

`structure_align/align_structures.py` 是独立于三个 level 的通用后处理脚本，不修改结构、不运行预测、不产生 MSA。它读取一个 CSV，将设计结构与设计前 backbone 配对：

1. 先使用 `chain_map`，没有时按实体类型、序列相似度和长度自动映射链。
2. 蛋白/核酸做整体和逐链 USalign，记录 RMSD、序列一致性和 TM-score。
3. 小分子在聚合物 USalign 变换后按原子匹配计算 RMSD；没有聚合物时使用 ligand-only Kabsch。
4. 指定 `motif_json`/`motif_spec` 时，对 motif 做 paired-atom Kabsch 局部 RMSD。

输出为 `alignment_summary.csv`、`chain_pairs.csv` 和 `motifs.csv`。总表 `status=success` 表示结构可读且整体有有效的聚合物 USalign 或小分子坐标比较；USalign 超时、不可执行或没有可匹配实体时通常为 `partial`，输入结构解析异常或输入字段缺失时为 `failed`。链对和 motif 另有自己的状态，局部失败不会覆盖已经成功的整体比较。

## 运行监控模式

v37 的运行监控分为三层：

- **Slurm 层**：`submit_level*.sh` 提交 array worker，再提交等待 array 的 aggregate job。用 `squeue`/`sacct` 看排队、运行、退出和依赖状态；aggregate 只负责合并已完成的分片。
- **v37 worker 层**：每个 worker 调用模型 wrapper，将 stdout/stderr 写入 `logs/`，随后检查结构和 summary/CSV 并写入 `results/parts/part_*.tsv`。aggregate 再生成 `metrics/levelN.csv`、`status/status.tsv` 和下游 `results/result_manifest.tsv`。重复运行默认启用 resume，只重跑没有完整成功结果的模型。
- **模型后处理层**：`for_wj_cid` 模式由 `af3_single_node_v18.py` 调用 `af3_monitor_v16.py`，负责 AF3 两种变体的监控、汇总、归档和清理；v37 最后解析归档结果。`multimodel_monitor_v1.py` 则是 Protenix/Boltz-2/OpenDDE 的独立扫描器，读取模型输出目录并生成统一指标 CSV，目前不会自动控制 v37 的 Slurm 任务，也不是 level 成功状态的唯一判据。

因此，“任务成功”要同时看 Slurm 状态、`status/status.tsv` 和对应日志；“科学上值得进入下一层”要再看 `metrics/levelN.csv` 中的实际指标。

## 输入 manifest

字段规范见 [Level1 INPUT_SCHEMA.md](level1/INPUT_SCHEMA.md)。推荐 CID 输入采用五列：
`design_id`、`sequence`、`chain_ids`、`cid_apo_json_path`、`cid_holo_json_path`；
使用 [config.cid.example.json](config.cid.example.json) 的 `af3_mode=cid` 即可自动调用 wrapper。
[可直接检查的 TSV 示例](examples/cid/input.tsv) 附有两份模板。
CID 的 apo/holo 分别写入 `metrics/cid_apo.csv`、`metrics/cid_holo.csv`，
结果 manifest 分别保留 `af3_apo_*`、`af3_holo_*`，状态也分别记录。

普通 direct 模式 TSV 示例：

```text
design_id  sequence  af3_json_path  metadata_json
d001        MKT...    templates/d001.json  {"clear_msa":true}
```

常用字段：

| 字段 | 用途 |
| --- | --- |
| `design_id` | 稳定主键；只能含字母、数字、`.`、`_`、`-` |
| `sequence` | 单链序列；多链可写 `CHAIN_A:CHAIN_B` 或 `CHAIN_A\|CHAIN_B` |
| `chain_ids` | 多链时显式指定链 ID，例如 `A,B` |
| `fasta_path` | 没有 `sequence` 时使用的 FASTA；路径相对 manifest 所在目录解析 |
| `af3_json_path` / `target_spec_path` | Level1 的 AF3 JSON 或模板；蛋白-配体复杂输入建议显式提供 |
| `protenix_input_path` / `boltz_input_path` / `opendde_input_path` | Level2 的模型专用输入；不提供时仅自动生成简单蛋白输入 |
| `level3_sequence` | Level3 序列模型实际评分的序列；多链时建议显式提供 |
| `metadata_json` | `af3_sequences`、`clear_msa` 等额外设置 |

所有路径字段会在读取 manifest 时转为绝对路径。带小分子、核酸或特殊链拓扑的输入不要依赖自动 builder，应在对应模型字段中提供现有 JSON/YAML/FASTA。

## 配置

[`config.example.json`](config.example.json) 是当前机器路径的示例。复制后按实际环境修改：

```bash
cp config.example.json my_config.json
```

不提供 `--config` 时，Level2 使用代码中的 Protenix/Boltz/OpenDDE 默认路径，Level3 使用 v37 自带的 NetSolP wrapper 和历史 TemBERTure 脚本；Level1 默认 `af3_mode=direct`，需要在配置中提供 `af3_script`、`af3_model_dir`、`af3_db_dir`。命令可通过 `commands` 覆盖，支持这些占位符：

```text
Level1: {esm_input} {esm_outdir} {af3_input_dir} {af3_outdir} {af3_fasta}
Level2: {protenix_input} {protenix_outdir} {boltz_input} {boltz_outdir}
        {opendde_input} {opendde_outdir} {design_id}
Level3: {level3_fasta} {netsolp_raw} {netsolp_csv} {temberture_csv}
通用:   {run_dir} {shard_dir} {input} {outdir}
```

直接 AF3 模式由 v37 为每个 design 生成一个 JSON。若复用 `for_wj/3_filter_v2/af3_single_node_v18.py`，应设置 `af3_mode=for_wj_cid`；该模式支持每个 design 含 1 条或更多蛋白链，wrapper 根据规范化 FASTA 自己生成 with-lig/without-lig 两套 JSON，v37 不再重复生成 AF3 输入。同一批次的蛋白链数和链 ID 必须匹配两份模板；不同聚体数分批运行。v37 自动传入 manifest 的链顺序，默认分析全部输出链：

```json
{
  "af3_mode": "for_wj_cid",
  "commands": {
    "af3": [
      "/xcfhome/zpzeng/originrepo/alphafold3/alphafold3_env/bin/python",
      "/xcfhome/yhliu/14_magpcr/for_wj/3_filter_v2/af3_single_node_v18.py",
      "--fasta", "{af3_fasta}",
      "--af3-script", "/xcfhome/zpzeng/originrepo/alphafold3/run_alphafold.py",
      "--af3-env", "/xcfhome/zpzeng/originrepo/alphafold3/alphafold3_env/bin/python",
      "--model-dir", "/xcfhome/pubdata/folding/alphafold3/models",
      "--db-dir", "/xcfhome/pubdata/folding/alphafold3",
      "--monitor-script", "/xcfhome/yhliu/14_magpcr/for_wj/3_filter_v2/af3_monitor_v16.py",
      "--template-with-lig", "/path/to/template_with_lig.json",
      "--template-without-lig", "/path/to/template_without_lig.json",
      "--json-output-dir", "{af3_input_dir}",
      "--af3-output-dir", "{af3_outdir}",
      "--csv-output", "{af3_csv}",
      "--archive-dir", "{af3_archive_dir}"
    ]
  }
}
```

这个 wrapper 自己负责按 FASTA 生成 AF3 JSON 和监控；v37 仍会单独运行自己的 ESMFold stage，因此不要在该 command 中加 `--run-esmfold`，除非你明确要重复计算。完成后，结果会同时记录 `af3_with_lig_*` 与 `af3_without_lig_*`；原有 `af3_*` 字段默认映射 with-lig，便于继续使用旧的 selector。

## 直接运行一个 level

下面以 Level1 为例，Level2/3 只需替换目录和脚本名：

```bash
V37=/xcfhome/yhliu/14_magpcr/z_zcodex/pipeline_scripts/00_latest/pipeline/v37

python "$V37/level1/prepare_level1.py" \
  --manifest input.tsv \
  --outdir runs/l1 \
  --shard-size 32 \
  --config my_config.json

python "$V37/level1/run_level1.py" --run-dir runs/l1 --shard-index 0
python "$V37/level1/aggregate_level1.py" --run-dir runs/l1
```

恢复时重复执行同一个 `run_level*.py` 即可；默认 `resume` 会复用已存在且解析完整的结果，若要强制重跑使用 `--no-resume`。

## 手动选择下游数量

选择指定 ID：

```bash
python "$V37/common/select_manifest.py" \
  --result-manifest runs/l1/results/result_manifest.tsv \
  --ids-file l1_selected_ids.txt \
  --require-success \
  --output runs/l1_selected.tsv
```

按指标取前 N 个（默认降序；低 PAE/PDE 指标加 `--ascending`）：

```bash
python "$V37/common/select_manifest.py" \
  --result-manifest runs/l1/results/result_manifest.tsv \
  --top-n 100 \
  --sort-field af3_ranking_score \
  --require-success \
  --output runs/l1_top100.tsv
```

选择结果会同时写入 `runs/selection.json`，记录输入、输出、时间和设计 ID。

## Slurm 提交与串联

每个 `submit_level*.sh` 都提交一个 array 和一个聚合 job，并打印两个 job ID：

```bash
bash "$V37/level1/submit_level1.sh" \
  --manifest input.tsv --outdir runs/l1 --config my_config.json \
  --partition nitrogen --time 04:00:00
# 记下 level1_aggregate_job=...

# 聚合完成后人工选择，再启动下一层
python "$V37/common/select_manifest.py" \
  --result-manifest runs/l1/results/result_manifest.tsv \
  --top-n 100 --sort-field af3_ranking_score --require-success \
  --output runs/l1_top100.tsv

bash "$V37/level2/submit_level2.sh" \
  --manifest runs/l1_top100.tsv --outdir runs/l2 \
  --config my_config.json --partition nitrogen --time 04:00:00 \
  --dependency LEVEL1_AGGREGATE_JOBID
```

Level3 同理：从 `runs/l2/results/result_manifest.tsv` 选择后，将输出 manifest 传给 `submit_level3.sh`。`--dependency` 只约束 array；每个 submit 脚本内部还会让 aggregate 等待自己的 array。若人工选择发生在上游完成之后，依赖 job 已完成时 `afterok` 仍可直接满足。

## 输出目录

每个 run directory 都包含：

```text
input/input_manifest.tsv       原始输入快照
tasks/                         总任务表、分片表、n_shards.txt
artifacts/                     模型输入、模型原始输出
results/parts/                 每个 shard 的原子结果
results/result_manifest.tsv    聚合后的下游 manifest
metrics/levelN.csv             每个 design/model 一行的指标表
status/status.tsv              总状态和逐模型状态
logs/                          wrapper 日志
```

状态值为 `success`、`partial`、`failed`、`missing`。解析器不会用默认值填充缺失置信度；多链 Level3 若没有 `level3_sequence`，默认会明确报错，只有配置 `level3_multichain_mode` 为 `first` 或 `concat` 才会采用隐式策略。

## 已知限制

- Level1 默认是本地原版 ESMFold + 直接 AF3；`for_wj` CID wrapper 需要显式设置 `af3_mode=for_wj_cid`，并提供二聚体 FASTA 语义对应的 with-lig/without-lig 模板。
- Level2 自动 builder 只覆盖简单蛋白输入；复杂配体、核酸、多实体输入应提供模型专用输入文件。
- Level3 的 NetSolP 与 TemBERTure 是序列模型，不会自动从结构结果推导新序列。
- fake wrapper 的三层链路 smoke 已通过；真实 CID Level1 smoke 已提交，完整真实 Level1→Level2→Level3 E2E 尚未验收。正式提交前仍需确认各模型环境、GPU/分区和输出命名契约。
- `multimodel_monitor_v1.py` 当前作为独立的 Level2 指标扫描器运行，尚未自动接入 v37 的 aggregate/cleanup 生命周期。
- v37 自带的 Level3 wrapper 位于 `wrappers/`：NetSolP 固定使用 `screening` 环境，OpenDDE 固定使用 `/xcfhome` 下的 `opendde` 环境；如环境路径不同，请在配置中覆盖对应 wrapper。
