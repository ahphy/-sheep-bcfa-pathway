# v5：单一共享主线，R/X差异对应三种产物

2026-11-17。根据用户要求，第3部分从三条反应路线合并为一条；去掉中英文标题中的“候选”，删除图内小字、缩写说明、线型图例、证据边界和长解释，另存到 figure_caption_zh.txt 与 figure_caption_en.txt。标题按用户要求简化；证据等级没有因标题改变而升级，科学边界在图注中保留。

## 新版内容

- 保留前两部分：瘤胃发酵/上皮吸收，以及脂肪细胞ACC供体生成。丙酸主要肝摄取、丁酸上皮β-羟丁酸代谢、供应区室等信息移入图注。
- 第3部分仅有5个节点、4条边：引物→首次延长→C2支链引入→C4支链形成→游离酸。
- 用R表示引物烃基，X表示支链。下方三张卡只列不同的引物、支链供体及产物结构，不再重复路线。
- 仅展示 MOA：R=CH3/X=CH3；MNA：R=CH2CH3/X=CH3；EOA：R=CH3/X=CH2CH3 三个配对。
- 图面不放正文式小字，保留反应标签、化学载体标识与R/X选择作为主信息。
- RDKit原生绘制5个通式及3种产物，统一真实产物的羧基方向和支链朝向。PyGraphviz布局唯一主线，SVG完成图面排版，CairoSVG导出PDF/PNG。

## 共享通式核对

R–CO–SCoA
→ R–(CH2)2–CO–SACP
→ R–(CH2)3–CH(X)–CO–SACP
→ R–(CH2)3–CH(X)–(CH2)2–CO–SACP
→ R–(CH2)3–CH(X)–(CH2)2–COOH

支链之前有3个CH2，支链与羧基之间有2个CH2。R/X表示完整取代基，不作为单个碳原子计数。

程序检查通式中不计R/X的固定碳数为1→3→5→7→7；支链相对羰基的位置为C2→C4→C4。三种R/X展开结构与各自目标酸的规范SMILES完全一致；路线碳数分别2→4→7→9、3→5→8→10、2→4→8→10。真实产物的分子式、母链长度、总碳数、C4支链及甲/乙基长度均通过校验。结构不指定立体构型。

编辑层标签包围框未检测到相互重叠或越界；该检查不涵盖所有RDKit字形或路径。已目视检查中英文PNG及Poppler渲染PDF，包括通式与产物，并进行独立结构/版式复查。SVG/PDF为矢量，PNG 6600×3630；validation_*.json记录结构、R/X展开核对、版本和SHA256。

## 文件与运行

- BCFA_shared_v5_zh/en.svg、.pdf、.png：中英文图。
- figure_caption_zh/en.txt：对应独立图注，可直接用于论文或综述草稿。
- preview_zh/en.png：较小预览。
- shared_*.dot 与 *_positioned.dot：共享主线的图定义及自动布局坐标。
- draw_bcfa_shared_v5.py、requirements.txt：代码及依赖。

```bash
python -m venv .venv
.venv/bin/pip install -r outputs/BCFA_shared_v5_20261117/requirements.txt
.venv/bin/python outputs/BCFA_shared_v5_20261117/draw_bcfa_shared_v5.py --lang both
```

环境需Cairo、Graphviz及Noto Sans CJK字体。本环境PyGraphviz使用自带Graphviz的二进制wheel。字体路径在脚本FP/FB中，跨系统使用时按实际位置调整。脚本只改写本v5目录里的生成结果，originals及v2/v3/v4均不变。没有新增文献检索或羊体内通量证据。
