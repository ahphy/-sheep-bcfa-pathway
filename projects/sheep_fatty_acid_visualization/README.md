# 98头羊三种脂肪酸与OAV分布可视化

这是 `ahphy/-sheep-bcfa-pathway` 仓库中的独立可视化子项目。使用每头动物的MOA、EOA、MNA浓度及OAV绘制交互分布图和静态PNG。

> **仓库为公开仓库，因此只提交代码和使用说明。** 不公开98头羊逐个样本的CSV、模型统计JSON、原始Excel、PNG或含个体信息的HTML。数据持有人需在本目录本地放置对应输入数据后运行。

## 快速开始

```bash
python -m pip install -r requirements.txt
python sheep_visualization.py
```

运行前在同一目录准备：

- `sheep_plot_data.csv`：先前从原始Excel整理出的个体数据，需保留中文列名（MOA μg/g、EOA μg/g、MNA μg/g、MOA指数、EOA指数、MNA指数、总指数、动物编号、品种、性别、月龄、批次、MOA状态、EOA状态、MNA状态）。
- `sheep_model_summary.json`：先前生成的分组和调整模型统计汇总；包含 `raw_breed_comparison` 和 `adjusted_breed_prediction`。

**上述两个输入文件不在公开GitHub仓库中；请从本次对话提供的可视化完整项目包取出。**


基于源文件 `98头羊_更新A_H规则_OAV加权控制品种性别月龄.xlsx` 中8张工作表中的数据制作。

## 交互成果

- `sheep_fatty_acid_dashboard.html`：不需要联网即可打开的独立HTML；4个分布子图 + 控制性别和月龄后的品种比较 + 浓度范围汇总表。鼠标悬停可查看动物编号、品种、性别、月龄、批次及数值；右上角切换线性/对数坐标。
- `sheep_fatty_acid_overview.png`：适合幻灯片/报告的高分辨率静态概览。

## 运行源码

- `sheep_visualization.py`：绘图源码（Python + Plotly + Matplotlib），其余两个数据文件作为输入。
- `sheep_plot_data.csv`：98头动物的绘图明细，原始列和状态保留。
- `sheep_model_summary.json`：原文件的分组统计、模型结果、两两比较、敏感性分析、三酸贡献与说明。

运行：

```bash
pip install plotly matplotlib
python sheep_visualization.py
```

## 图形编码

灰色细横线：该品种用于分析的全部数值的最小至最大值；粗横线：25%–75%分位数；◆：中位数；实心圆：单动物值；空心三角：低值替代涉及的动物。三酸浓度图竖虚线为1.6667 μg/g最低校准点（**不是经过验证的LOD**）。

## 注意

1. 来源表以0.83335 μg/g代替“<1.6667 μg/g”；这不是精确实测值，也不是可靠的检出限。此类观测使用三角标记。EOA 98例中的54例是该类替代值。
2. `OAV = MOA/0.02 + EOA/0.006 + MNA/0.65`，按照原表计算逻辑；这是基于用户选定阈值的加权理论指标，阈值适用于实际基质的前提仍需核实。
3. 桑塔因斯绵羊n=1，仅显示单个观测，不用于分布推断；主模型n=96，剔除了此单例及月龄未知的萨福克个体。
4. 原始分布与校正后模型输出分开展示。没有凭空生成95%置信区间；原始均值与调整后的平均预测值不是同一个统计量。
5. 来源表格提醒部分数据需进一步复核峰质量、最低校准点、批次及样本来源等方法问题。

图形采用 Plotly（https://github.com/plotly/plotly.py）及 Matplotlib。源数据文件与含个体数据的生成结果未公开提交到GitHub。