# 98头羊三种脂肪酸与加权气味指数可视化

使用 MOA（4-甲基辛酸）、EOA（4-乙基辛酸）、MNA（4-甲基壬酸）及 OAV 数据，按品种绘制浓度范围、四分位区间和逐头动物分布。另提供控制性别、月龄后的品种间总 OAV 比较图。

> **隐私提示**：这是公开 GitHub 仓库，只存储程序与说明。原始 Excel、逐头样本 CSV、模型统计 JSON、交互 HTML 及含样本数据的结果图 **不提交到 GitHub**。请只在本地处理实验数据。

## 1. 环境与数据文件

Python 3.10+；安装依赖：

```bash
python -m pip install -r requirements.txt
```

将以下两个私有输入文件放在 `sheep_visualization.py` 所在目录（或者使用命令行参数指定位置）：

- `sheep_plot_data.csv`：至少包含 `动物编号`、`品种`、`性别`、`月龄`、`批次`、`MOA μg/g`、`EOA μg/g`、`MNA μg/g`、`MOA状态`、`EOA状态`、`MNA状态`、`MOA指数`、`EOA指数`、`MNA指数`、`总指数`。在当前项目中是98头动物的整理数据。
- `sheep_model_summary.json`：包括 `raw_breed_comparison` 与 `adjusted_breed_prediction`，是原Excel统计结果的整理版，**并非由此脚本重新拟合模型**。

数据持有人可使用之前私下生成的完整项目包中的文件；仓库**不提供**原始实验数据。

## 2. 一键生成总览图 + 5套独立高清图

```bash
python sheep_visualization.py
```

首次运行后自动生成：

- 根目录：`sheep_fatty_acid_dashboard.html`（四酸/OAV综合交互面板、校正对比和统计表）；`sheep_fatty_acid_overview.png`（2×2总览，300 dpi）。
- `exports/` 目录：`panel_MOA`、`panel_EOA`、`panel_MNA`、`panel_OAV`、`panel_Adjusted_OAV` 五套独立图。

**每套独立图默认包含：**

| 文件名后缀 | 类型 | 用途 |
|---|---|---|
| `_300dpi.png` | 300 dpi PNG | PPT、报告、一般论文图片 |
| `_600dpi.png` | 600 dpi PNG | 打印、高清排版 |
| `_600dpi.tif` | 600 dpi，LZW 压缩 TIFF | 期刊可能要求的格式 |
| `.pdf` | 矢量 PDF | 论文排版，放大不失真 |
| `.svg` | 矢量 SVG | 排版与图形后期处理 |
| `_interactive.html` | 独立交互 HTML | 各面板分别悬停查看个体信息 |

总览另导出 `exports/sheep_fatty_acid_overview_300dpi.png`、`exports/sheep_fatty_acid_overview.pdf` 和 `exports/sheep_fatty_acid_overview.svg`。

独立的高分辨率图是**重新绘制的面板**，不是从总图截屏/裁切；不会出现悬停提示框遮挡图形的问题。

## 3. 自定义高清导出设置

```bash
# 指定PNG dpi（TIFF使用其中最大的dpi）
python sheep_visualization.py --dpi 300 600

# 只生成300 dpi图片、跳过TIFF、跳过独立交互HTML（减少磁盘空间）
python sheep_visualization.py --dpi 300 --no-tiff --no-single-html

# 指定文件位置及导出目录
python sheep_visualization.py --data ./private/sheep_plot_data.csv --summary ./private/sheep_model_summary.json --output-dir ./exports
```

程序会尝试自动选择中文字体：Windows 优先使用微软雅黑/黑体，macOS 可使用苹方，Linux 可安装 `Noto Sans CJK`。如果输出中中文显示为方框，请先安装可用的中文系统字体。

## 4. 图形含义与数据限制

- **灰色细线**：该品种数据的最小—最大值；**彩色粗线**：Q1–Q3；**菱形**：中位数；**实心圆**：单头动物的值；**空心三角形**：涉及低值替代。
- 三酸浓度图的橙色虚线为 **1.6667 μg/g 最低校准点**（不是经验证的检出限）；部分低值在源数据中以 0.83335 μg/g 代替，不能解释为真实测得的浓度。
- OAV 是根据原表选定的理论阈值加权指标，公式：`MOA/0.02 + EOA/0.006 + MNA/0.65`。其物理意义、阈值对脂肪基质的适用性需要独立验证。
- 浓度图与 OAV 图默认使用对数横轴；含极端样本时比线性坐标更能看出分布。单个品种 n=1 不能推断其分布。
- 校正后比较图的数据直接取自模型汇总 JSON；它把原始平均值和控制性别/月龄后的模型预测平均值清楚区分。**没有模型置信区间数据，不绘制虚构的置信区间**。
- 不同品种的样本量差距较大，且山羊/绵羊物种差异与品种分组存在潜在混杂。需要谨慎解读组间差异。

## 5. 使用的开源库

- [Plotly Python](https://github.com/plotly/plotly.py)——交互式图表。
- [Matplotlib](https://github.com/matplotlib/matplotlib)——PNG/TIFF/PDF/SVG 高清导出。

本程序 **不** 上传文件、不修改原始数据，也不自动向 GitHub 推送任何本地生成结果。