# 跨机器环境配置

控制器与模型分开安装。`controller.yml` 仅安装 manifest、指标解析、归档和结构比对所需的
轻量依赖，不包含 PyTorch/JAX/CUDA 或模型权重。模型按各自上游说明建立独立环境；
本 release 不提供所有 GPU 模型的通用锁文件。

推荐布局（路径可自由调整）：

```text
miniconda3/envs/
  profold/       # 控制器、监控、测试
  esmfold/       # 原版 esm-fold CLI
  af3/           # AlphaFold3 + JAX
  protenix/      # Protenix + PyTorch
  boltz/         # Boltz-2
  opendde/       # OpenDDE
software/
  alphafold3/ Protenix/
data/
  alphafold3/models/ opendde/
profold/                 # 本仓库
```

| 环境 | 上游安装来源 | ProFold 调用契约与限制 |
| --- | --- | --- |
| esmfold | [ESM](https://github.com/facebookresearch/esm) | 需要 `esm-fold -i FASTA -o PDB_DIR`；仅安装 fair-esm 不保证提供此 CLI，需按现有安装验证或通过 `commands.esmfold` 适配 |
| af3 | [AlphaFold3](https://github.com/google-deepmind/alphafold3) | 配置 Python、run_alphafold.py、数据库和权重；本地验证过 v3.0.4 + 自定义 patch，权重按上游许可自行获取 |
| protenix | [Protenix](https://github.com/bytedance/Protenix) | 保留原 `protenix_100.py` 及 `protenix_base_20250630_v1.0.0` 默认模型；不是已验证的 Protenix-v2 adapter。升级需匹配 `PROTENIX_SCRIPT`/模型名或 `commands.protenix` 并实测输出 |
| boltz | [Boltz](https://github.com/jwohlwend/boltz) | `boltz predict --write_full_pae --write_full_pde`；GPU kernel 由本机安装决定 |
| opendde | [OpenDDE](https://github.com/DeepPrinciple/OpenDDE) | 自带 `opendde_predict.sh` 与 FASTA 转换器，环境需提供 `opendde pred`；配置数据目录 |
| esmfold2 | 预留 | 当前 Level2 未接入，不提供假定可执行的环境 |

## 配置和检查

从仓库根目录执行：

```bash
conda env create -f environments/controller.yml
conda activate profold
python scripts/configure.py \
  --conda-root /opt/miniconda3 \
  --software-root /opt/software \
  --data-root /data/models \
  --output config.local.json
python scripts/check_environment.py --config config.local.json --level level1
python -m unittest discover -s tests -v
```

在生成的 `config.local.json` 中修改实际路径、环境名称和版本。
默认 `af3_mode=cid`；普通 AF3 单次预测改为 `direct`。
`site.example.json` 的 `${...}` 由 configure 展开；`{...}` 是 worker 的运行占位符，
不要把未生成的 site 模板直接传给 prepare。生成器不覆盖已有配置。
检查器仅检查依赖和路径，不执行模型、不下载权重、不证明 CUDA 可用。
每个模型先做一条真实 GPU smoke，再扩大任务。

prepare/run/aggregate 使用控制器 Python；模型使用配置里的独立 executable。
CID wrapper/monitor 使用 AF3 Python，AF3 环境还需提供 numpy、biopython、gemmi。
Slurm 可设置 `PYTHON_BIN=/opt/miniconda3/envs/v37/bin/python`，队列与 GPU 参数按本机配置。
结构比对另需 USalign，调用时传 `--usalign /path/to/USalign`。

## 复现与升级

记录模型源码 commit、权重版本、CUDA/驱动、环境导出和实际命令；不同硬件不要直接
复用带绝对 prefix 的 conda export。新环境验证后再替换路径。
已有 run/config.json 是快照；新机器重新 prepare 新运行目录。搬迁已有归档的路径迁移
不在本 release 范围内，STR20 closeout 中的校验清单不能随意编辑。

STR20 prepare 使用私有历史选样数据，发布包不含这些数据。
通过 `STR20_SOURCE_DIR` 指定脚本要求的 CSV/template 目录；AF3 路径通过
`AF3_PYTHON`、`AF3_SCRIPT`、`AF3_MODEL_DIR`、`AF3_DB_DIR` 覆盖。
一般 CID 新任务使用 examples/cid 五列格式，无需 STR20 历史数据。
