# v37.0.0

首次独立发布 v37 三层蛋白筛选脚本，加入跨机器环境配置及仓库内 AF3 多链监控依赖。

- Level1: 原版 ESMFold + AF3；CID 支持多蛋白链、apo/holo 分别预测和统计。
- Level2: Protenix、Boltz-2、OpenDDE 独立预测与结果解析。
- Level3: NetSolP、TemBERTure；支持人工子集、Slurm 分片与 resume。
- AF3 整体/逐链/逐对指标、校验归档、失败诊断保留、STR20 自动收尾。
- 独立 controller 环境、逐模型安装说明、配置生成器、只读环境检查器。
- 仅包含代码、格式示例和测试；不含运行结果、实验序列、缓存、权重或数据库。

验证：源码包在独立目录解压、干净 Python 3.11 环境中通过 38 项 unittest；配置生成、
CID 示例 prepare 与 shell 语法检查通过。已补齐正式测试 fixture 和 PyYAML 依赖。
历史 STR20 AF3 60/60 成功，本 release 未新增 GPU 推理验证。
ESMFold 单体与完整三层真实 E2E 尚未验收；ESMFold2 未接入；Level2 统一监控
尚未集成自动归档/删除；Protenix-v2 不是本次验证的默认模型。
仅 STR20 runner 实现全 run 精简收尾，通用 CID 保留原归档契约。

clone 后使用本地分支修改，通过 git push/pull 在机器间同步；
固定复现可 checkout v37.0.0。模型、权重与 config.local.json 单独管理。
