#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Recreate the 98-animal fatty-acid distribution dashboard.
Input: sheep_plot_data.csv + sheep_model_summary.json (exported from the supplied Excel).
Install: pip install plotly matplotlib
Run: python sheep_visualization.py
"""
from __future__ import annotations
import argparse
import csv
import html
import json
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio

BASE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description="98头羊三种脂肪酸及OAV的交互图与高清单图导出")
parser.add_argument("--data", type=Path, default=BASE / "sheep_plot_data.csv",
                    help="个体数据CSV，默认与脚本位于同一目录")
parser.add_argument("--summary", type=Path, default=BASE / "sheep_model_summary.json",
                    help="模型结果JSON，默认与脚本位于同一目录")
parser.add_argument("--output-dir", type=Path, default=BASE / "exports",
                    help="独立图导出目录（默认：exports）")
parser.add_argument("--dpi", type=int, nargs="+", default=[300, 600],
                    help="单图PNG的输出分辨率，例如 --dpi 300 600")
parser.add_argument("--no-tiff", action="store_true", help="跳过单图TIFF输出")
parser.add_argument("--no-single-html", action="store_true", help="跳过4个单独的Plotly交互HTML")
args = parser.parse_args()
DPI_VALUES = sorted(set(args.dpi))
if not DPI_VALUES or min(DPI_VALUES) <= 0:
    parser.error("--dpi 必须是正整数")
OUT = args.output_dir
OUT.mkdir(parents=True, exist_ok=True)
for p in (args.data, args.summary):
    if not p.is_file():
        parser.error(f"没有找到输入文件：{p}")

records = []
with args.data.open(encoding="utf-8-sig", newline="") as fd:
    for row in csv.DictReader(fd):
        for c in ("MOA μg/g", "EOA μg/g", "MNA μg/g",
                  "MOA指数", "EOA指数", "MNA指数", "总指数"):
            row[c] = float(row[c])
        records.append(row)
meta = json.loads(args.summary.read_text(encoding="utf-8"))
order = [r["品种"] for r in meta["raw_breed_comparison"]]
by = defaultdict(list)
for r in records:
    by[r["品种"]].append(r)
COLORS = {"MOA": "#008C82", "EOA": "#E59D2D", "MNA": "#7C66B5",
          "OAV": "#315D80"}
facets = [
    ("MOA μg/g", "MOA：4-甲基辛酸 (μg/g)", "MOA", "MOA状态"),
    ("EOA μg/g", "EOA：4-乙基辛酸 (μg/g)", "EOA", "EOA状态"),
    ("MNA μg/g", "MNA：4-甲基壬酸 (μg/g)", "MNA", "MNA状态"),
    ("总指数", "加权气味指数：总 OAV (无量纲)", "OAV", None),
]
fig = make_subplots(rows=2, cols=2, horizontal_spacing=.11, vertical_spacing=.17,
                    subplot_titles=[f[1] for f in facets])
# Legend once, via common HTML note; each panel uses consistent glyphs.
for j, (column, _, key, statuscol) in enumerate(facets):
    rr, cc = divmod(j, 2)
    rr += 1
    cc += 1
    color = COLORS[key]
    for yi, breed in enumerate(order):
        sample = by[breed]
        vals = sorted(float(r[column]) for r in sample)
        n = len(vals)
        low, high = vals[0], vals[-1]
        med = statistics.median(vals)
        def quantile(frac):
            if n == 1:
                return vals[0]
            t = (n-1)*frac
            lo = math.floor(t)
            hi = math.ceil(t)
            return vals[lo] + (vals[hi]-vals[lo])*(t-lo)
        q1, q3 = quantile(.25), quantile(.75)
        summary = (f"{html.escape(breed)}（n={n}）<br>"
                   f"用于分析的范围：{low:,.3f} ~ {high:,.3f}<br>"
                   f"中位数：{med:,.3f}<br>IQR：{q1:,.3f} ~ {q3:,.3f}")
        if n > 1:
            fig.add_trace(go.Scatter(x=[low, high], y=[yi, yi],
                          mode="lines", line=dict(color="#A9B6C3", width=2),
                          hoverinfo="skip", showlegend=False), row=rr, col=cc)
            fig.add_trace(go.Scatter(x=[q1, q3], y=[yi, yi],
                          mode="lines", line=dict(color=color, width=10),
                          opacity=.8, hoverinfo="skip", showlegend=False),
                          row=rr, col=cc)
        fig.add_trace(go.Scatter(x=[med], y=[yi], mode="markers",
                      marker=dict(symbol="diamond", size=11, color=color,
                                  line=dict(color="white", width=1.2)),
                      hovertemplate=summary+"<extra></extra>", showlegend=False),
                      row=rr, col=cc)
        for kind, symbol in [(False, "circle"), (True, "triangle-up-open")]:
            selected = []
            for k, item in enumerate(sample):
                is_substitute = (
                    (statuscol is not None and
                     str(item[statuscol]).startswith("<")) or
                    (statuscol is None and any(
                        str(item[v]).startswith("<") for v in
                        ("MOA状态","EOA状态","MNA状态")))
                )
                if bool(is_substitute) == kind:
                    selected.append((k,item))
            if not selected:
                continue
            rand = random.Random(98231 + yi * 761 + j * 37)
            xx, yy, ht = [], [], []
            for k, item in selected:
                xx.append(item[column])
                yy.append(yi + rand.uniform(-.25, .25))
                details = [
                    f"<b>{html.escape(str(item['动物编号']))} · {html.escape(breed)}</b>",
                    f"{column}：{item[column]:,.4f}",
                    f"性别：{html.escape(str(item['性别']))}；月龄：{html.escape(str(item['月龄']))}",
                    f"批次：{html.escape(str(item['批次']))}",
                ]
                if statuscol:
                    details.append("浓度状态：" + html.escape(str(item[statuscol])))
                elif is_substitute:
                    details.append("至少一种脂肪酸使用低值替代")
                ht.append("<br>".join(details))
            fig.add_trace(go.Scatter(x=xx, y=yy, mode="markers",
                          marker=dict(symbol=symbol,
                                      size=7 if not kind else 9,
                                      opacity=.72,
                                      color=color if not kind else "white",
                                      line=dict(color=color, width=1.45)),
                          text=ht, hovertemplate="%{text}<extra></extra>",
                          showlegend=False), row=rr, col=cc)
    fig.update_yaxes(
        tickvals=list(range(len(order))),
        ticktext=[f"{b}  n={len(by[b])}" for b in order],
        range=[len(order)-.40, -.57], tickfont=dict(size=12),
        zeroline=False, showgrid=False, row=rr, col=cc)
    fig.update_xaxes(type="log", title_text="浓度 (μg/g)" if j < 3 else "总气味指数 (OAV)",
                     showgrid=True, gridcolor="#E6EBF0", zeroline=False,
                     ticks="outside", tickfont=dict(size=11),
                     row=rr, col=cc)
    if j < 3:
        fig.add_vline(x=1.6667, line_color="#CF8D48", line_dash="dot",
                      line_width=1.3, row=rr, col=cc)
fig.update_layout(template="plotly_white", width=1550, height=900,
                  paper_bgcolor="#FFFFFF", plot_bgcolor="#FFFFFF",
                  margin=dict(l=160,r=25,t=92,b=80),
                  font=dict(family="Arial,'Noto Sans CJK SC','Microsoft Yahei',sans-serif",
                            color="#283343",size=12),
                  hoverlabel=dict(bgcolor="white",font_size=12),
                  title=dict(text="98头动物｜三酸浓度范围、四分位分布与加权气味指数",x=.04,
                             font=dict(size=22,color="#1C334C")))
# A log scale remains the best default; the user can toggle linear across panels.
fig.update_layout(updatemenus=[dict(type="buttons", direction="right",
    x=.99,y=1.12,xanchor="right",yanchor="top",showactive=True,
    bgcolor="#EEF2F6",
    buttons=[
      dict(label="对数坐标",method="relayout",
           args=[{"xaxis.type":"log","xaxis2.type":"log","xaxis3.type":"log","xaxis4.type":"log"}]),
      dict(label="线性坐标",method="relayout",
           args=[{"xaxis.type":"linear","xaxis2.type":"linear","xaxis3.type":"linear","xaxis4.type":"linear"}])
    ])])

# Adjusted OAV: raw means and sex/age standardized predicted means are conceptually
# distinct, so they are labeled clearly and displayed separately from raw distribution.
adj = meta["adjusted_breed_prediction"]
raw = {r["品种"]:r for r in meta["raw_breed_comparison"]}
fig2=go.Figure()
for i, r in enumerate(adj):
    breed=r["品种"]
    mean=float(raw[breed]["总指数均值"])
    pred=float(r["标准化预测均值"])
    fig2.add_trace(go.Scatter(
        x=[min(mean,pred), max(mean,pred)], y=[i,i], mode="lines",
        line=dict(color="#B5BFCA",width=2),hoverinfo="skip",showlegend=False))
    fig2.add_trace(go.Scatter(x=[mean],y=[i],mode="markers",
        marker=dict(color="#ADB6C0",size=12,symbol="circle"),
        hovertemplate=f"{breed}<br>原始平均总指数：{mean:,.1f}<extra></extra>",
        name="原始均值",showlegend=i==0))
    fig2.add_trace(go.Scatter(x=[pred],y=[i],mode="markers",
        marker=dict(color="#266C95",size=14,symbol="diamond"),
        hovertemplate=f"{breed}<br>控制性别和月龄后的预测均值：{pred:,.1f}<extra></extra>",
        name="校正后预测均值",showlegend=i==0))
fig2.update_layout(
    template="plotly_white",width=1350,height=470,
    title=dict(text="校正性别、月龄后的品种平均总气味指数（模型有效 n=96）",
               x=.04,font=dict(size=21,color="#1C334C")),
    font=dict(family="Arial,'Noto Sans CJK SC','Microsoft Yahei',sans-serif",
              color="#283343",size=12),
    yaxis=dict(tickvals=list(range(len(adj))),ticktext=[f"{r['品种']}  n={r['模型有效n']}" for r in adj],
               range=[len(adj)-.6,-.6],showgrid=False),
    xaxis=dict(type="log",title="平均加权气味指数（对数坐标）",gridcolor="#E6EBF0"),
    legend=dict(orientation="h",y=-.22,x=.45),
    paper_bgcolor="white",plot_bgcolor="white",
    margin=dict(l=180,r=45,t=75,b=100)
)

def render_stats():
    headings=["品种","n","MOA 范围 μg/g","EOA 范围 μg/g","MNA 范围 μg/g",
              "总 OAV 中位数","低值替代数（MOA/EOA/MNA）"]
    body=[]
    for breed in order:
        rr=by[breed]
        vals=lambda c:[float(x[c]) for x in rr]
        ranges=["{:.2f}–{:.2f}".format(min(vals(c)),max(vals(c))) for c in
                ("MOA μg/g","EOA μg/g","MNA μg/g")]
        flags=[sum(str(x[c]).startswith("<") for x in rr) for c in
               ("MOA状态","EOA状态","MNA状态")]
        table_row=[breed,str(len(rr)),*ranges,f"{statistics.median(vals('总指数')):,.1f}",
                   "/".join(map(str,flags))]
        body.append("<tr>"+"".join("<td>"+html.escape(x)+"</td>" for x in table_row)+"</tr>")
    return "<table><thead><tr>"+ "".join("<th>"+x+"</th>" for x in headings)+\
           "</tr></thead><tbody>"+"".join(body)+"</tbody></table>"

markup = """<!doctype html><html lang="zh-CN"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>羊脂肪酸与气味指数｜交互分布图</title>
<style>
:root {font-family:system-ui,-apple-system,"Noto Sans CJK SC","Microsoft Yahei",sans-serif;color:#21354A;background:#F5F7F9}
body{margin:0 auto;max-width:1550px;padding:20px 24px 56px}
h1{font-size:28px;letter-spacing:.01em;margin:0 0 8px}
h2{font-size:20px;margin:20px 0 10px}
.lead{color:#64748B;margin:0 0 20px;line-height:1.65}
.tile{background:white;border:1px solid #E4E8EE;border-radius:15px;padding:20px;margin-bottom:18px}
.plot{width:100%;overflow-x:auto}
.note{font-size:14px;color:#506178;line-height:1.8}
.note strong{color:#21354A}
.legend{display:flex;flex-wrap:wrap;gap:16px;font-size:13px;color:#506178}
.legend span{white-space:nowrap}
table{width:100%;border-collapse:collapse;font-size:13px;text-align:left}
td,th{padding:11px 12px;border-bottom:1px solid #E7EDF3}
th{background:#F4F7FA;font-weight:650}
td:not(:first-child){font-variant-numeric:tabular-nums}
.info{background:#F8F9FB;border-left:3px solid #8797A8;padding:12px 14px;line-height:1.75;margin:10px 0}
@media(max-width:850px){body{padding:12px}.tile{padding:12px}table{min-width:1000px}}
</style></head><body>
<h1>不同品种的三种膻味脂肪酸与加权气味指数</h1>
<p class="lead">原始数据：98头动物、6个品种｜交互版：可悬停查看个体、缩放、平移、切换对数/线性坐标。</p>
<section class="tile"><h2>① 浓度范围、四分位区间与个体分布</h2>
<div class="legend">
<span>灰线 = 本组最小值至最大值</span>
<span>彩色粗线 = Q1–Q3（中间50%）</span>
<span>◆ = 中位数</span><span>● = 个体数据</span><span>△ = 低值替代涉及的个体</span>
<span>竖向虚线 = 最低校准点1.6667 μg/g（三酸图）</span>
</div>
<div class="plot">__FIG1__</div>
<p class="note">注意：三酸图中△代表该酸被标注为“&lt;1.6667 μg/g”且以0.83335 μg/g替代；总OAV图中△表示该动物至少有一种酸采用低值替代。替代值不代表真实检测浓度；因此绘制的最小值和IQR是<strong>用于分析的数据</strong>的分布，并非全部为实测浓度范围。</p>
</section>
<section class="tile"><h2>② 控制性别和月龄后的品种比较</h2>
<div class="plot">__FIG2__</div>
<p class="note">原始均值直接来自98头动物统计，校正后预测均值来自文件中“校正品种预测”工作表的模型结果，主模型有效n=96。桑塔因斯绵羊只有1头，未纳入主模型。原文件没有提供校正预测的95%置信区间，因此本图不虚构误差线。总体品种效应P≈0.00221；两两比较需以Holm校正结果为准。</p>
</section>
<section class="tile"><h2>③ 各品种原始范围与低值替代概况</h2>
<div style="overflow-x:auto">__TABLE__</div>
<div class="info note"><strong>图表解读：</strong>浓度图与OAV图使用不同的含义和数值范围，不能把三种酸的浓度简单相加解释气味强度。总OAV=MOA/0.02+EOA/0.006+MNA/0.65（基于Excel使用的阈值），仅为理论阈值加权指标，需确认阈值与样品基质的可比性。</div>
</section>
<section class="tile"><h2>统计与数据局限</h2>
<p class="note">1. 三酸浓度及总气味指数默认使用对数坐标，以免极端值压缩大多数样本。2. 桑塔因斯绵羊 n=1，不应对其作稳定的品种分布推断。3. EOA低值替代比例较高，尤其陕西关中奶山羊为33/55，图中大量重叠△不意味着都实测等于0.83335。4. 文件中将1.6667 μg/g描述为最低校准点而非经过验证的检出限。5. 品种、批次和动物来源可能混杂；文件中的校正结果仅控制性别、月龄，敏感性分析才进一步考虑批次。6. 陕西关中奶山羊是山羊，与绵羊相比还涉及物种差异。7. 面积/峰高选取规则和部分峰质量仍需方法学复核。</p>
</section></body></html>"""
a=pio.to_html(fig,include_plotlyjs=True,full_html=False,config={"displaylogo":False,"responsive":True})
b=pio.to_html(fig2,include_plotlyjs=False,full_html=False,config={"displaylogo":False,"responsive":True})
markup=markup.replace("__FIG1__",a).replace("__FIG2__",b).replace("__TABLE__",render_stats())
output=BASE/"sheep_fatty_acid_dashboard.html"
output.write_text(markup,encoding="utf-8")
print("Wrote",output,output.stat().st_size,"bytes")

# --- Static figures and independent high-resolution exports ---
# Keep the original 2x2 interactive dashboard while allowing every panel
# to be exported separately without taking screenshots or cropping the overview.
import copy
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, fontManager
from matplotlib.lines import Line2D


def choose_chinese_font() -> FontProperties:
    """Use a locally installed Chinese-capable font on Windows/macOS/Linux."""
    installed = {entry.name for entry in fontManager.ttflist}
    candidates = (
        "Microsoft YaHei", "SimHei", "PingFang SC", "Noto Sans CJK SC",
        "Noto Sans CJK JP", "WenQuanYi Micro Hei", "Arial Unicode MS",
    )
    for name in candidates:
        if name in installed:
            return FontProperties(family=name)
    print("警告：未发现支持中文的系统字体，导出图的中文可能变成方框。请安装微软雅黑或 Noto Sans CJK。")
    return FontProperties()


font = choose_chinese_font()
plt.rcParams.update({
    "font.family": font.get_family(),
    "axes.unicode_minus": False,
    "pdf.fonttype": 42,  # retain embeddable TrueType text in PDF
    "svg.fonttype": "path",  # avoids missing glyphs when viewing SVG elsewhere
    "savefig.facecolor": "white",
})


def percentile(values, q):
    """Linear interpolation matching the original Q1/Q3 display logic."""
    if len(values) == 1:
        return values[0]
    t = (len(values) - 1) * q
    lo, hi = math.floor(t), math.ceil(t)
    return values[lo] + (values[hi] - values[lo]) * (t - lo)


def draw_distribution(ax, column, title, key, statuscol, panel_index):
    """Shared renderer for the overview and four standalone panels."""
    color = COLORS[key]
    for yi, breed in enumerate(order):
        sample = by[breed]
        vals = sorted(float(z[column]) for z in sample)
        if len(vals) > 1:
            ax.plot([vals[0], vals[-1]], [yi, yi],
                    color="#ABB7C3", lw=1.8, zorder=1)
            ax.plot([percentile(vals, .25), percentile(vals, .75)],
                    [yi, yi], color=color, lw=9, alpha=.76,
                    solid_capstyle="butt", zorder=2)
        ax.scatter([statistics.median(vals)], [yi], marker="D", s=68,
                   c=[color], edgecolors="white", linewidths=1, zorder=5)
        rng = random.Random(98231 + yi * 761 + panel_index * 37)
        for z in sample:
            substituted = (
                str(z[statuscol]).startswith("<") if statuscol is not None
                else any(str(z[x]).startswith("<") for x in ("MOA状态", "EOA状态", "MNA状态"))
            )
            ax.scatter([float(z[column])], [yi + rng.uniform(-.23, .23)],
                       s=25 if not substituted else 42,
                       marker="o" if not substituted else "^",
                       facecolor=color if not substituted else "white",
                       alpha=.77, edgecolors=color, linewidths=.95, zorder=4)
    ax.set_xscale("log")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"{breed}  (n={len(by[breed])})" for breed in order],
                       fontproperties=font, fontsize=10.2)
    ax.invert_yaxis()
    ax.set_ylim(len(order) - .55, -.55)
    ax.grid(axis="x", alpha=.22, color="#8FA1AF")
    ax.set_axisbelow(True)
    ax.set_title(title, fontproperties=font, fontsize=15, pad=14)
    ax.set_xlabel("浓度 (μg/g)" if key != "OAV" else "总加权气味指数 (OAV)",
                  fontproperties=font, fontsize=11, labelpad=8)
    ax.tick_params(axis="x", labelsize=9.5)
    if key != "OAV":
        ax.axvline(1.6667, color="#C9914B", linewidth=1.25, linestyle="--", zorder=0)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#B4C1CE")


def panel_legend(color):
    """Legend uses the same symbols as the original Plotly dashboard."""
    return [
        Line2D([0], [0], color="#ABB7C3", lw=2, label="最小值–最大值"),
        Line2D([0], [0], color=color, lw=8, alpha=.76, label="Q1–Q3"),
        Line2D([0], [0], color=color, marker="D", markersize=7,
               markeredgecolor="white", lw=0, label="中位数"),
        Line2D([0], [0], color=color, marker="o", markersize=6,
               lw=0, label="个体值"),
        Line2D([0], [0], color=color, marker="^", markersize=7,
               markerfacecolor="white", lw=0, label="低值替代相关"),
    ]


def save_individual(fig, stem):
    """Write one plot in requested PNG densities + vector + optional TIFF."""
    for dpi in DPI_VALUES:
        dest = OUT / f"{stem}_{dpi}dpi.png"
        fig.savefig(dest, dpi=dpi, bbox_inches="tight", facecolor="white")
        print("导出", dest.name)
    if not args.no_tiff:
        dpi = max(DPI_VALUES)
        dest = OUT / f"{stem}_{dpi}dpi.tif"
        fig.savefig(dest, dpi=dpi, bbox_inches="tight", facecolor="white",
                    pil_kwargs={"compression": "tiff_lzw"})
        print("导出", dest.name)
    for ext in ("pdf", "svg"):
        dest = OUT / f"{stem}.{ext}"
        fig.savefig(dest, bbox_inches="tight", facecolor="white")
        print("导出", dest.name)


# 1. A clean 2x2 overview (300 dpi PNG, PDF, SVG).
figs, axes = plt.subplots(2, 2, figsize=(17, 10.5), layout="constrained")
for j, (column, title, key, statuscol) in enumerate(facets):
    draw_distribution(axes.flat[j], column, title, key, statuscol, j)
figs.suptitle("不同品种三种脂肪酸的范围与分布，以及加权气味指数",
              fontproperties=font, fontsize=21, fontweight="bold")
for output_file in (
    BASE / "sheep_fatty_acid_overview.png",
    OUT / "sheep_fatty_acid_overview_300dpi.png",
):
    figs.savefig(output_file, dpi=300, bbox_inches="tight", facecolor="white")
for ext in ("pdf", "svg"):
    figs.savefig(OUT / f"sheep_fatty_acid_overview.{ext}",
                 bbox_inches="tight", facecolor="white")
plt.close(figs)

# 2. Four independently rendered distribution panels.
for j, (column, title, key, statuscol) in enumerate(facets):
    standalone, ax = plt.subplots(figsize=(10.2, 5.7))
    draw_distribution(ax, column, title, key, statuscol, j)
    standalone.subplots_adjust(left=.22, right=.97, top=.87, bottom=.23)
    standalone.legend(handles=panel_legend(COLORS[key]), ncol=5,
                      loc="lower center", bbox_to_anchor=(.5, .025),
                      frameon=False, prop=font, fontsize=9)
    save_individual(standalone, f"panel_{key}")
    plt.close(standalone)

# 3. Adjusted breed comparison also gets independent static exports.
standalone, ax = plt.subplots(figsize=(10.2, 5.5))
for i, r in enumerate(adj):
    breed = r["品种"]
    orig = float(raw[breed]["总指数均值"])
    adjusted = float(r["标准化预测均值"])
    ax.plot([orig, adjusted], [i, i], color="#B5BFCA", linewidth=2, zorder=1)
    ax.scatter([orig], [i], s=74, marker="o", color="#ADB6C0", zorder=3)
    ax.scatter([adjusted], [i], s=92, marker="D", color="#266C95", zorder=4)
ax.set_yticks(range(len(adj)))
ax.set_yticklabels([f"{r['品种']}  (n={r['模型有效n']})" for r in adj],
                   fontproperties=font, fontsize=10.2)
ax.set_ylim(len(adj) - .55, -.55)
ax.set_xscale("log")
ax.set_xlabel("平均加权气味指数 (OAV，log 坐标)", fontproperties=font, fontsize=11)
ax.set_title("控制性别与月龄后的品种比较（模型有效 n=96）",
             fontproperties=font, fontsize=15, pad=14)
ax.grid(axis="x", alpha=.22, color="#8FA1AF")
ax.set_axisbelow(True)
for spine in ("top", "right", "left"):
    ax.spines[spine].set_visible(False)
standalone.subplots_adjust(left=.27, right=.96, top=.86, bottom=.23)
legend = [
    Line2D([0], [0], marker="o", lw=0, color="#ADB6C0", markersize=9, label="原始均值"),
    Line2D([0], [0], marker="D", lw=0, color="#266C95", markersize=9, label="校正后预测均值"),
]
standalone.legend(handles=legend, ncol=2, loc="lower center",
                  bbox_to_anchor=(.5, .05), frameon=False, prop=font)
save_individual(standalone, "panel_Adjusted_OAV")
plt.close(standalone)

# 4. Keep interactive output as before; additionally export independent HTML
# with hover labels to allow each individual plot to be inspected in isolation.
if not args.no_single_html:
    for j, (column, title, key, statuscol) in enumerate(facets):
        rr, cc = divmod(j, 2)
        standalone = go.Figure()
        for trace in fig.select_traces(row=rr + 1, col=cc + 1):
            trace_copy = copy.deepcopy(trace)
            trace_copy.update(xaxis="x", yaxis="y")
            standalone.add_trace(trace_copy)
        standalone.update_layout(
            template="plotly_white", width=1120, height=550,
            title=dict(text=title, x=.04, font=dict(size=21)),
            font=dict(family="Arial,'Noto Sans CJK SC','Microsoft Yahei',sans-serif", size=12),
            margin=dict(l=200, r=50, t=90, b=70),
            paper_bgcolor="white", plot_bgcolor="white",
            hoverlabel=dict(bgcolor="white", font_size=12),
            showlegend=False,
        )
        standalone.update_yaxes(tickvals=list(range(len(order))),
            ticktext=[f"{b}  n={len(by[b])}" for b in order],
            range=[len(order)-.40, -.57], showgrid=False)
        standalone.update_xaxes(type="log",
            title_text="浓度 (μg/g)" if j < 3 else "总气味指数 (OAV)",
            showgrid=True, gridcolor="#E6EBF0")
        if j < 3:
            standalone.add_vline(x=1.6667, line_color="#CF8D48",
                                 line_dash="dot", line_width=1.3)
        pio.write_html(standalone, OUT / f"panel_{key}_interactive.html",
                       include_plotlyjs=True, full_html=True,
                       config={"displaylogo": False, "responsive": True})
        print("导出", f"panel_{key}_interactive.html")
    pio.write_html(fig2, OUT / "panel_Adjusted_OAV_interactive.html",
                   include_plotlyjs=True, full_html=True,
                   config={"displaylogo": False, "responsive": True})
    print("导出 panel_Adjusted_OAV_interactive.html")

print("全部导出完成。高清单图目录：", OUT)
