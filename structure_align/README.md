# v37 通用结构对齐

`align_structures.py` 独立于三个 screening level。它解决“一个骨架对应多条设计序列”时的结构归属和对齐问题，不修改输入结构、不运行预测，也不会产生 MSA。

## 输入 CSV

最小字段：

```csv
design_id,sequence_id,backbone_structure_path,predicted_structure_path
cid01_seq01,cid01,../backbones/cid01.cif,../level1/design01_af3.cif
```

路径相对 CSV 所在目录解析。可选字段：

- `chain_map`：骨架链到预测链，例如 `A:A,B:B,C:C` 或 `{"A":"B","B":"A"}`；不提供时按实体类型、序列相似度和长度自动映射。
- `motif_json` / `motif_spec`：内嵌 JSON、JSON 文件路径，或紧凑写法 `site:A:10-14,20;other:B:5-8`。
- `motif_name,motif_chain,motif_residues,motif_atom`：单个 motif 的简写；残基编号为结构中的 residue number，`motif_atom` 默认蛋白 `CA`、核酸 `C4'`。

完整 motif JSON 示例：

```json
[
  {"name": "active_site", "backbone_chain": "A", "predicted_chain": "A",
   "backbone_residues": [10, 12, 13], "predicted_residues": [10, 12, 13], "atom": "CA"},
  {"name": "interface", "chain": "B", "residues": "4-12", "atom": "CA"}
]
```

也可以用 `backbone_positions` / `predicted_positions` 指定各自聚合物链的 1-based 序列位置。这样能处理设计后 residue number 改变的结构。

## 输出

```text
alignment_summary.csv   每个 design 一行；overall_* 是预测结构→骨架的复合物结果
chain_pairs.csv          每个骨架链×每个预测链一行，包括蛋白、核酸和 ligand
motifs.csv              指定局部残基/原子的 paired-atom Kabsch RMSD/TM-like score
raw/                     仅在 --keep-raw 时保存 USalign 原始文本
```

蛋白/核酸链和复合物的 `rmsd`、`tm_score_norm_backbone` 来自 USalign；`overall_tmscore_norm_backbone` 的归一化长度是骨架。小分子不强行伪造 USalign TM-score：`overall_ligand_rmsd` 是在聚合物 USalign 变换后按 residue/atom name 匹配的重原子 RMSD；没有聚合物时使用 ligand-only Kabsch。

## 运行

```bash
V37=/xcfhome/yhliu/14_magpcr/z_zcodex/pipeline_scripts/00_latest/pipeline/v37
python "$V37/structure_align/align_structures.py" \
  --input-csv align_input.csv \
  --outdir runs/align \
  --workers 4 \
  --usalign /xcfhome/yhliu/002_software/009_TMalign/USalign
```

这个脚本可以直接消费 Level1/Level2 生成的结构路径，也可以单独用于已有结构与设计前骨架的通用比较。
