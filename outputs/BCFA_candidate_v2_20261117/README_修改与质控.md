# 羊膻味脂肪酸候选宿主合成通路 v2

日期：2026-11-17。依据用户上传的 BCFA_pathway.rar 重绘。原始 19 个文件逐字节保留在仓库 originals/，没有执行会覆盖旧图的原始程序。附件说明作为资料读取，没有据其设置额外审批或放行流程。

## 文件与依赖检查

已读取全部 5 个 Python 文件：draw_pathway.py、BCFA_mechanism_source.py、BCFA_shared_mechanism_source.py、BCFA_rumen_host_mechanism_source.py、create_host_mechanism.py；阅读两份 Markdown 说明；检查全部 4 套 SVG/PDF/PNG 参考图。其中 create_host_mechanism.py 用字符串替换生成并执行另一脚本，原始脚本均将图片写入自身目录，存在覆盖同名正式图的风险。原包没有 requirements.txt 或 pyproject.toml。

原代码依赖 pathlib、Matplotlib、RDKit；使用 RDKit 坐标后手工画化学键、氧原子及 C4 标记，且用固定原子索引校验。环境已有 Matplotlib 3.10.8，但起初没有 RDKit。本次在独立 venv 安装依赖；新脚本不依赖 Matplotlib，使用 RDKit MolDraw2DSVG 和 CairoSVG，Pillow 用于文字尺寸核查。RAR 解压另用 libarchive-c 5.3，不是绘图依赖。需要系统 Cairo 和 DejaVu Sans 字体；脚本内字体路径为 Debian/Linux 路径，移植时按系统修改。

## 修改内容

- 重排为瘤胃发酵与上皮吸收、脂肪细胞胞质 ACC 前体网络、三条并列 FASN 候选组装路线。保留与原图一致的英文图面标签。
- 分开乙酸外周供应、丙酸主要肝脏摄取和丁酸上皮转化为 β-羟丁酸；不将上皮代谢产物画成已证实直接供应胞质支链 CoA。
- 用实线、虚线、点线区分已知生化/生理步骤、候选反应、间接或不确定供应。线粒体 PCC–MCEE–MMUT 与胞质分开；保留 ECHDC1 修复支路。
- 逐条列出具体 ACP 中间体和总碳数：MOA 2→4→7→9，MNA 3→5→8→10，EOA 2→4→8→10。甲基丙二酰基供体脱羧保留 3 碳，乙基丙二酰基保留 4 碳；最后 +2 碳使支链位置从 C2 变为 C4。
- 用 RDKit 原生 SVG 绘制产物的单键、羧基和双键，高亮支链并标注 C4。通过羧基 SMARTS 找起点、最长碳链和支链分析代替硬编码原子编号。三种产物立体构型均未指定。
- 为延长供体标签添加背景，解决箭头穿字；扩大产品结构和分区留白。新代码及所有结果位于 outputs/BCFA_candidate_v2_20261117/。

## 校验结果

| 产物 | 分子式 | 母链碳数 | 支链 | 总碳数 |
|---|---|---:|---|---:|
| MOA：4-甲基辛酸 | C9H18O2 | 8 | C4 甲基 | 9 |
| MNA：4-甲基壬酸 | C10H20O2 | 9 | C4 甲基 | 10 |
| EOA：4-乙基辛酸 | C10H20O2 | 8 | C4 乙基 | 10 |

程序校验分子式、总碳数、母链长度、支链位置与长度；91 个作者添加文字标签无包围框重叠或越界。RDKit 结构内文字不包含在该检测中，已额外目视检查结构和 PDF 栅格预览；未见明显标签、碳链或排版问题。SVG 解析通过，PDF 已用 Poppler 渲染检查，PNG 为 5400×7020 像素。SVG 保留图面文字；化学结构字形由 RDKit 绘成矢量路径。validation.json 记录结构、版本及结果文件 SHA256。

## 科学边界

这是候选机制图，不是已证实的羊体内完整通路。胞质 ACC 旁反应、支链供体来源、FASN 对供体的接受及短链提前释放仍需羊组织与通量实验证据；未确定终止酶身份。原始资料中的文献线索及证据边界保存在 originals/ 两份说明中，本次没有进行新的系统文献检索或逐节点原文复核。图面实线表示一般生化/生理基础，不表示已经证明羊体内目标产物通量。

## 重新生成

```bash
python -m venv .venv
.venv/bin/pip install -r outputs/BCFA_candidate_v2_20261117/requirements.txt
.venv/bin/python outputs/BCFA_candidate_v2_20261117/draw_bcfa_v2.py
```

脚本只写同目录的 v2 结果。如需保留已有 v2，请复制整个版本目录后运行。原始正式图不会被修改。
