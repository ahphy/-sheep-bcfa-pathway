"""Compact, self-contained BCFA mechanism: dot layout, RDKit structures, SVG editing.
Run with --lang zh (default), --lang en, or --lang both. Prior versions are untouched.
"""
from pathlib import Path
from importlib.metadata import version
import argparse, json, math, re, hashlib
import xml.etree.ElementTree as ET
import pygraphviz as pgv
import cairosvg
from PIL import ImageFont, Image
from rdkit import Chem
from rdkit.Geometry import Point2D
from rdkit.Chem import rdMolDescriptors, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
OUT=Path(__file__).resolve().parent
NS='http://www.w3.org/2000/svg';ET.register_namespace('',NS)
INK='#253844';MUT='#596d77';GREEN='#23776c';PURPLE='#8051a0';ORANGE='#b76122';BLUE='#246fa5'
COLORS=[BLUE,PURPLE,ORANGE];FONT='Noto Sans CJK SC'
FONT_PATH='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
FONT_BOLD='/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc'

def build(lang):
 def L(zh,en):return zh if lang=='zh' else en
 root=ET.Element('{'+NS+'}svg',{'width':'2000'});labels=[];structures=[]
 def rect(parent,x,y,w,h,fill='white',stroke='none',r=12):
  return ET.SubElement(parent,'{'+NS+'}rect',{'x':str(x),'y':str(y),'width':str(w),'height':str(h),'rx':str(r),'fill':fill,'stroke':stroke})
 def txt(parent,x,y,s,size=24,color=INK,bold=False,anchor='middle',track=True):
  attrs={'x':str(x),'y':str(y),'font-family':FONT,'font-size':str(size),'text-anchor':anchor,'fill':color}
  if bold:attrs['font-weight']='bold'
  el=ET.SubElement(parent,'{'+NS+'}text',attrs);el.text=s
  if track:
   f=ImageFont.truetype(FONT_BOLD if bold else FONT_PATH,size*4,index=2);bw=f.getlength(s)/4
   b=f.getbbox(s,anchor='ls');x0=x-bw/2 if anchor=='middle' else x if anchor=='start' else x-bw
   labels.append({'text':s,'bounds':[x0,y+b[1]/4,x0+bw,y+b[3]/4]})
  return el
 def arr(x1,y1,x2,y2,color=MUT,kind='solid'):
  attrs={'d':f'M{x1},{y1} L{x2},{y2}','stroke':color,'stroke-width':'2.5','fill':'none'}
  if kind!='solid':attrs['stroke-dasharray']='9 6' if kind=='candidate' else '2 6'
  ET.SubElement(root,'{'+NS+'}path',attrs)
  ang=math.atan2(y2-y1,x2-x1);pts=[(x2,y2),(x2-11*math.cos(ang)+4*math.sin(ang),y2-11*math.sin(ang)-4*math.cos(ang)),(x2-11*math.cos(ang)-4*math.sin(ang),y2-11*math.sin(ang)+4*math.cos(ang))]
  ET.SubElement(root,'{'+NS+'}polygon',{'points':' '.join(f'{x},{y}' for x,y in pts),'fill':color})
 rect(root,0,0,2000,2000)
 txt(root,1000,52,L('羊膻味脂肪酸的候选宿主合成机制','Candidate host synthesis of mutton odor-associated fatty acids'),36,INK,True)
 txt(root,1000,91,L('MOA / MNA / EOA：同一套延长逻辑，由起始引物与支链供体决定产物','MOA / MNA / EOA: primer choice + branched donor choice determine the product'),23,MUT)
 rect(root,30,115,675,355,'#f1f6f9','#cfdae0')
 rect(root,720,115,1250,355,'#f2f8f5','#cfddd7')
 txt(root,54,155,L('1  瘤胃发酵 → 上皮吸收','1  Rumen fermentation → absorption'),26,INK,True,'start')
 txt(root,745,155,L('2  脂肪细胞胞质：ACC 生成延长供体','2  Adipocyte cytosol: ACC supplies extension donors'),26,INK,True,'start')
 txt(root,143,190,L('饲料糖 → 微生物','Diet → microbes'),20,MUT)
 txt(root,282,190,L('上皮','Epithelium'),19,MUT)
 txt(root,475,190,L('宿主去向','Host handling'),20,MUT)
 txt(root,881,190,L('酰基-CoA 前体','Acyl-CoA precursors'),20,MUT)
 txt(root,1400,190,L('CoA 延长供体','CoA extension donors'),20,MUT)
 txt(root,1815,190,L('脱羧后保留','Retained after CO2 loss'),19,MUT)
 # Membrane is a compartment boundary, not a metabolic conversion.
 ET.SubElement(root,'{'+NS+'}path',{'d':'M274,209 L274,437 M285,209 L285,437','stroke':'#a2b6c2','stroke-width':'3','fill':'none'})
 acids=[L('乙酸 C2','Acetate C2'),L('丙酸 C3','Propionate C3'),L('丁酸 C4','Butyrate C4')]
 hosts=[(L('外周乙酸','Peripheral acetate'),L('可供脂肪组织合成','Available for adipose lipogenesis')),(L('门静脉丙酸','Portal propionate'),L('主要入肝 → 糖异生','Mainly liver → gluconeogenesis')),(L('吸收的丁酸','Absorbed butyrate'),L('上皮代谢 → β-羟丁酸','Epithelial metabolism → β-hydroxybutyrate'))]
 coas=[L('乙酰-CoA · C2','Acetyl-CoA · C2'),L('丙酰-CoA · C3','Propionyl-CoA · C3'),L('丁酰-CoA · C4','Butyryl-CoA · C4')]
 donors=[L('丙二酰-CoA · C3','Malonyl-CoA · C3'),L('甲基丙二酰-CoA · C4','Methylmalonyl-CoA · C4'),L('乙基丙二酰-CoA · C5','Ethylmalonyl-CoA · C5')]
 for i,(y,c,dc) in enumerate(zip([245,330,415],COLORS,[GREEN,PURPLE,ORANGE])):
  txt(root,142,y,acids[i],26,c,True);arr(242,y-8,348,y-8,c)
  txt(root,478,y-3,hosts[i][0],23,c,True);txt(root,478,y+24,hosts[i][1],17,MUT)
  arr(600,y-8,754,y-8,c,'solid' if i==0 else 'supply')
  txt(root,675,y-19,L('活化','Activation') if i==0 else L('供应待证','Supply unverified'),17,MUT)
  txt(root,878,y,coas[i],25 if lang=='zh' else 23,c,True)
  if i==0:txt(root,878,y+25,L('乙酰-CoA合成酶','Acetyl-CoA synthetase'),17,MUT)
  arr(1006,y-8,1190,y-8,dc,'solid' if i==0 else 'candidate')
  txt(root,1098,y-19,'ACC' if i==0 else L('ACC 旁反应*','ACC side reaction*'),18,dc,True)
  rect(root,1210,y-35,391,55,'#e1f0eb' if i==0 else '#eee7f4' if i==1 else '#f8ecdf')
  txt(root,1405,y,donors[i],24,dc,True)
  txt(root,1815,y,'+2 C' if i==0 else '+3 C' if i==1 else '+4 C',29,dc,True)
 txt(root,1000,504,L('下方 FASN 调用这些引物 / 供体；供体均为 CoA 形式，+2 / +3 / +4 为脱羧后净加碳数。','FASN uses these primers / CoA donors; +2 / +3 / +4 are net carbons retained after decarboxylation.'),23,MUT)
 txt(root,53,551,L('3  FASN 候选组装：三轮延长 → 支链酸','3  Candidate FASN assembly: three extensions → branched acids'),27,INK,True,'start')
 txt(root,1000,588,L('每轮：缩合 → 还原 → 脱水 → 还原（需 NADPH）  |  中间体结合 ACP；Cn 只计酰基碳数','Each cycle: condensation → reduction → dehydration → reduction (NADPH) | ACP-bound intermediates; Cn = acyl carbons'),21,MUT)
 # Automatic reaction-matrix layout; SVG editorial layer supplies shared headers.
 g=pgv.AGraph(directed=True,strict=False)
 g.graph_attr.update(rankdir='LR',bgcolor='transparent',pad='.1',nodesep='.10',ranksep='.35',splines='polyline',outputorder='edgesfirst',newrank='true')
 g.node_attr.update(shape='box',style='rounded,filled',fillcolor='white',fontname=FONT,fontsize='17',fontcolor=INK,margin='.10,.06',height='.70',penwidth='1.2')
 g.edge_attr.update(fontname=FONT,fontsize='13',style='dashed',penwidth='1.4',arrowsize='.65')
 rows=[('MOA',2,4,7,9,'methyl'),('MNA',3,5,8,10,'methyl'),('EOA',2,4,8,10,'ethyl')]
 for i,(name,p,first,bran,final,branch) in enumerate(rows):
  c=COLORS[i];bname=L('甲基','methyl') if branch=='methyl' else L('乙基','ethyl')
  labels_row=[name+'\n'+(coas[0] if p==2 else coas[1]),L('直链酰基','Straight acyl')+'\nC'+str(first),f'C{bran}\nC2 '+bname,f'C{final}\nC4 '+bname,'']
  for j,lab in enumerate(labels_row):
   g.add_node(f'n{i}_{j}',label=lab,color=c,group=name,width=['1.9','1.4','1.4','1.4','3.6'][j],height='1.20',fixedsize='true' if j==4 else 'false')
  for j,label in enumerate([L('+丙二酰\n+2C','+ malonyl\n+2C'),L('+甲基丙二酰\n+3C','+ methylmalonyl\n+3C') if branch=='methyl' else L('+乙基丙二酰\n+4C','+ ethylmalonyl\n+4C'),L('+丙二酰\n+2C','+ malonyl\n+2C'),L('提前水解*','Early hydrolysis*')]):
   g.add_edge(f'n{i}_{j}',f'n{i}_{j+1}',label=label,color=c,fontcolor=GREEN if j in (0,2) else PURPLE if j==1 and branch=='methyl' else ORANGE if j==1 else c)
  assert p+2+(3 if branch=='methyl' else 4)+2==final
 for j in range(5):
  g.add_subgraph([f'n{i}_{j}' for i in range(3)],rank='same')
  for i in range(2):g.add_edge(f'n{i}_{j}',f'n{i+1}_{j}',style='invis',weight='100')
 g.write(str(OUT/f'assembly_{lang}.dot'));g.layout(prog='dot');g.write(str(OUT/f'assembly_{lang}_positioned.dot'))
 boxes=[]
 for n in g.nodes():
  x,y=map(float,n.attr['pos'].split(','));w=float(n.attr['width'])*72;h=float(n.attr['height'])*72
  boxes.append((str(n),(x-w/2,y-h/2,x+w/2,y+h/2)))
 node_overlaps=[]
 for i,(na,a) in enumerate(boxes):
  for nb,b in boxes[i+1:]:
   if min(a[2],b[2])>max(a[0],b[0])+.1 and min(a[3],b[3])>max(a[1],b[1])+.1:node_overlaps.append([na,nb])
 assert not node_overlaps,node_overlaps
 svg=ET.fromstring(g.draw(format='svg').decode())
 graphgroup=next(z for z in svg if z.tag=='{'+NS+'}g')
 tr=graphgroup.get('transform');trans=re.search(r'translate\(([-.\d]+)[ ,]+([-.\d]+)\)',tr);tx,ty=map(float,trans.groups())
 vb=list(map(float,svg.get('viewBox').split()));scale=1940/vb[2]
 top=654;core_height=vb[3]*scale
 svg.set('x','30');svg.set('y',str(top));svg.set('width','1940');svg.set('height',str(core_height))
 # Shared headers align with the automatically calculated column positions.
 for j,head in enumerate([L('引物','Primer'),L('① 常规延长','① Extend'),L('② 引入支链','② Introduce branch'),L('③ 再延长','③ Extend again'),L('产物 / C4 支链','Products / C4 branch')]):
  x=float(g.get_node(f'n0_{j}').attr['pos'].split(',')[0]);txt(root,30+(x+tx)*scale,636,head,24,INK,True)
 products=[('MOA',L('4-甲基辛酸','4-Methyloctanoic acid'),'O=C(O)CCC(C)CCCC','C9H18O2',8,1),('MNA',L('4-甲基壬酸','4-Methylnonanoic acid'),'O=C(O)CCC(C)CCCCC','C10H20O2',9,1),('EOA',L('4-乙基辛酸','4-Ethyloctanoic acid'),'O=C(O)CCC(CC)CCCC','C10H20O2',8,2)]
 for group in svg.iter('{'+NS+'}g'):
  title=group.find('{'+NS+'}title')
  if title is None:continue
  if group.get('class')=='node' and re.fullmatch(r'n[0-2]_[0-3]',title.text or ''):
   # dot uses equal-height row cells to align all five columns. Render the
   # short intermediate cards within those cells, preserving center positions.
   nd=g.get_node(title.text);cx,cy=map(float,nd.attr['pos'].split(','));cy=-cy;nw=float(nd.attr['width'])*72
   for old in list(group):
    if old.tag in ('{'+NS+'}path','{'+NS+'}polygon'):group.remove(old)
   card=ET.Element('{'+NS+'}rect',{'x':str(cx-nw/2),'y':str(cy-26),'width':str(nw),'height':'52','rx':'10','fill':'white','stroke':COLORS[int(title.text[1])],'stroke-width':'1.2'})
   group.insert(1,card)
  if group.get('class')=='edge':
   # Dot positions edge labels; backgrounds protect glyphs from routed curves.
   for t in list(group):
    if t.tag!='{'+NS+'}text':continue
    fs=float(t.get('font-size'));f=ImageFont.truetype(FONT_PATH,round(fs*4),index=2);tw=f.getlength(t.text or '')/4
    x=float(t.get('x'));y=float(t.get('y'))
    bg=ET.Element('{'+NS+'}rect',{'x':str(x-tw/2-2),'y':str(y-fs*.94),'width':str(tw+4),'height':str(fs*1.18),'fill':'white','stroke':'none'})
    group.insert(list(group).index(t),bg)
  if not re.fullmatch(r'n[0-2]_4',title.text or ''):continue
  i=int(title.text[1]);name,full,sm,formula,parent,bl=products[i];c=COLORS[i]
  m=Chem.MolFromSmiles(sm);assert rdMolDescriptors.CalcMolFormula(m)==formula
  acid=m.GetSubstructMatch(Chem.MolFromSmarts('[CX3](=[OX1])[OX2H1]'))[0]
  cs=[a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum()==6]
  ends=[x for x in cs if x!=acid and sum(n.GetAtomicNum()==6 for n in m.GetAtomWithIdx(x).GetNeighbors())==1]
  path=max((Chem.GetShortestPath(m,acid,x) for x in ends),key=len)
  bp=[x for x in path if sum(n.GetAtomicNum()==6 for n in m.GetAtomWithIdx(x).GetNeighbors())==3]
  off=[x for x in cs if x not in path]
  assert len(path)==parent and len(cs)==parent+bl and len(off)==bl and len(bp)==1 and path.index(bp[0])+1==4
  rdDepictor.Compute2DCoords(m,canonOrient=False,coordMap={idx:Point2D(k*1.3,(k%2)*.60) for k,idx in enumerate(path)})
  d=rdMolDraw2D.MolDraw2DSVG(252,57);o=d.drawOptions();o.clearBackground=False;o.useBWAtomPalette();o.bondLineWidth=1.3;o.minFontSize=10;o.maxFontSize=13;o.atomLabels[bp[0]]='C4'
  rgb=tuple(int(c[z:z+2],16)/255 for z in (1,3,5));ids=bp+off
  d.DrawMolecule(m,highlightAtoms=ids,highlightAtomColors={x:rgb for x in ids});d.FinishDrawing()
  ms=d.GetDrawingText();molsvg=ET.fromstring(ms[ms.index('<svg'):]);cx,cy=map(float,g.get_node(title.text).attr['pos'].split(','));cy=-cy
  txt(group,cx,cy-31,name+' · '+full,13,c,True,track=False)
  molsvg.set('x',str(cx-126));molsvg.set('y',str(cy-26));group.append(molsvg)
  txt(group,cx,cy+36,formula,12,MUT,track=False)
  structures.append({'name':name,'smiles':sm,'formula':formula,'total_carbons':len(cs),'parent_chain':parent,'branch_carbons':bl,'branch_position':4,'stereochemistry':'unspecified'})
 root.append(svg)
 y=top+core_height+14
 rect(root,30,y,1940,67,'#edf5ef')
 txt(root,1000,y+42,L('为何都在 C4？末次延长在羧基端加入 2C，使已有支链从 C2 → C4。','Why C4? The final extension adds 2C at the carboxyl end, shifting the existing branch from C2 → C4.'),27,GREEN,True)
 y+=103
 for x,kind,label in [(55,'solid',L('已知生化 / 生理步骤','Established biochemistry / physiology')),(720,'candidate',L('候选反应 / 组装','Candidate reaction / assembly')),(1330,'supply',L('胞质供应未定','Unresolved cytosolic supply'))]:
  arr(x,y-7,x+55,y-7,MUT,kind);txt(root,x+72,y,label,21,MUT,anchor='start')
 y+=37
 txt(root,1000,y,L('* ACC 旁反应、FASN 供体接受与提前释放待验证；未证明羊体内完整通量，产物立体构型未指定。','* ACC side reactions, FASN donor acceptance and early release remain unverified in sheep; stereochemistry is unspecified.'),21,INK)
 y+=34
 txt(root,1000,y,L('ACC：乙酰-CoA羧化酶  ·  FASN：脂肪酸合酶  ·  CoA：辅酶A  ·  ACP：酰基载体蛋白','ACC, acetyl-CoA carboxylase · FASN, fatty acid synthase · CoA, coenzyme A · ACP, acyl carrier protein'),20,MUT)
 H=math.ceil(y+28);root.set('height',str(H));root.set('viewBox',f'0 0 2000 {H}')
 # Check the manually authored editorial text. Graph-native text and molecular
 # labels have their own layout and require preview inspection in addition.
 overlaps=[]
 for i,a in enumerate(labels):
  for b in labels[i+1:]:
   aa,bb=a['bounds'],b['bounds']
   if min(aa[2],bb[2])-max(aa[0],bb[0])>1 and min(aa[3],bb[3])-max(aa[1],bb[1])>1:overlaps.append([a['text'],b['text']])
 out_of_bounds=[a for a in labels if a['bounds'][0]<0 or a['bounds'][2]>2000 or a['bounds'][1]<0 or a['bounds'][3]>H]
 assert not overlaps,overlaps
 assert not out_of_bounds,out_of_bounds
 content=ET.tostring(root,encoding='utf-8',xml_declaration=True);base=OUT/f'BCFA_compact_v4_{lang}'
 base.with_suffix('.svg').write_bytes(content)
 cairosvg.svg2pdf(bytestring=content,write_to=str(base.with_suffix('.pdf')))
 cairosvg.svg2png(bytestring=content,write_to=str(base.with_suffix('.png')),output_width=6000,output_height=H*3)
 # Small display copy, not an extra mechanism version.
 im=Image.open(base.with_suffix('.png'));im.thumbnail((1600,1200));im.save(OUT/f'preview_{lang}.png')
 report={'language':lang,'structures':structures,'route_carbon_arithmetic':'passed','graph_nodes':len(g.nodes()),'graph_edges':sum(e.attr.get('style')!='invis' for e in g.edges()),'node_overlap_pairs':node_overlaps,'editorial_text_labels':len(labels),'editorial_text_overlap_pairs':overlaps,'png_size':[6000,H*3],'versions':{p:version(p) for p in ['pygraphviz','rdkit','CairoSVG','Pillow']},'sha256':{ext:hashlib.sha256(base.with_suffix('.'+ext).read_bytes()).hexdigest() for ext in ['svg','pdf','png']}}
 (OUT/f'validation_{lang}.json').write_text(json.dumps(report,indent=2,ensure_ascii=False));print(json.dumps(report,indent=2,ensure_ascii=False))

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--lang',choices=['zh','en','both'],default='zh');args=parser.parse_args()
 for lang in ['zh','en'] if args.lang=='both' else [args.lang]:build(lang)
