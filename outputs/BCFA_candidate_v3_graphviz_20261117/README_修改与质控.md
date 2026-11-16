# BCFA 候选通路 v3：RDKit + PyGraphviz + SVG

2026-11-17。按用户选定的三层绘图方案新增版本；originals/ 与 v2 均未修改。

## 实现与变化

- RDKit：原生 SVG 结构绘制，保留 C4 标记及支链颜色；由羧基 SMARTS 和最长碳路径核查结构，不依赖固定原子编号。
- PyGraphviz 2.0.4 / Graphviz dot：自动布置 28 个节点、24 条连线及瘤胃、上皮、脂肪细胞、FASN 区室。前体网络和 FASN 路线分别布局，避免将跨区 CoA 供应画成已证实反应。两个布局面板通过文字说明关联，未绘制前体面板至 FASN 面板的直接连线。
- SVG 精排：统一两个自动布局面板的宽度，添加生物学标题、图例、科学边界和背景遮罩，避免区室标题或边标签被箭头穿过。RDKit 分子结构按 dot 输出的产物节点坐标嵌入；不手写反应节点坐标。
- 新增 precursor_layout.dot、fasn_layout.dot（图定义），以及 *_positioned.dot（布局后的坐标与曲线）。可复用或直接用 Graphviz 查看。
- 导出 BCFA_candidate_v3.svg、.pdf、.png；PNG 5400×7485。矢量文件适合后续编辑。

## 验证

- MOA：C9H18O2，母链 8 碳、C4 甲基。
- MNA：C10H20O2，母链 9 碳、C4 甲基。
- EOA：C10H20O2，母链 8 碳、C4 乙基。
- 独立碳数算术校验：2+2+3+2=9；3+2+3+2=10；2+2+4+2=10。
- dot 实际节点矩形无相互重叠。该检测不是全部文字和箭头的自动无碰撞证明。
- 已目视检查 PDF 的 Poppler 渲染预览，修正跨区箭头经过区室标题及边标签的问题；未见明显文字、结构或裁切问题。
- validation.json 记录结构、布局检测、包版本、分辨率和图片 SHA256。

## 运行

```bash
python -m venv .venv
.venv/bin/pip install -r outputs/BCFA_candidate_v3_graphviz_20261117/requirements.txt
.venv/bin/python outputs/BCFA_candidate_v3_graphviz_20261117/draw_bcfa_v3.py
```

本次 Linux 环境使用 PyPI 的 PyGraphviz 二进制 wheel（包含 Graphviz 运行组件）。其他平台若没有相应 wheel，需按 PyGraphviz 官方说明安装 Graphviz；CairoSVG 需要系统 Cairo。标签测量使用 DejaVu Sans 的 Linux 字体路径，跨平台时需调整。

## 科学边界

继续保留候选模型：ACC 旁反应、脂肪细胞支链供体来源、FASN 供体接受性、短链提前水解未据此证明为羊体内完整通量。三种产物均不指定立体构型。线粒体分解与胞质 CoA 池分开，ECHDC1 修复背景保留。本次仅升级绘图架构，没有新增文献证据审查。原始证据线索见 originals/，原 v2 的检查说明见相邻版本目录。
