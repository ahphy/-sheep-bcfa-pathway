"""PyGraphviz dot layout + native RDKit SVG structures + SVG composition.
Run this file to export v3 SVG/PDF/PNG and editable layout DOT files.
"""
from pathlib import Path
import json, re, hashlib
from importlib.metadata import version
import xml.etree.ElementTree as ET
from PIL import ImageFont
import pygraphviz as pgv
import cairosvg
from rdkit import Chem
from rdkit.Chem import rdDepictor, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
OUT=Path(__file__).resolve().parent
NS='http://www.w3.org/2000/svg'; ET.register_namespace('',NS)
INK='#243746'; MUT='#566875'; GREEN='#26766b'; COLORS=['#256da6','#8052a0','#b96321']
checks=[]
def graph():
 g=pgv.AGraph(directed=True,strict=False)
 g.graph_attr.update(rankdir='TB',bgcolor='transparent',pad='.25',nodesep='.45',ranksep='.58',splines='polyline',outputorder='edgesfirst',fontname='DejaVu Sans',fontsize='18',compound='true',newrank='true')
 g.node_attr.update(shape='box',style='rounded,filled',fillcolor='white',fontname='DejaVu Sans',fontsize='16',fontcolor=INK,margin='.18,.12',width='4.7',height='.65',penwidth='1.1')
 g.edge_attr.update(fontname='DejaVu Sans',fontsize='13',color=MUT,fontcolor=MUT,arrowsize='.75',penwidth='1.3')
 return g
def cluster(g,name,label,fill):
 return g.add_subgraph(name='cluster_'+name,label=label,style='rounded,filled',fillcolor=fill,color='#cbd7dd',fontcolor=INK,fontsize='18',labeljust='l',margin='22')
def node(g,key,label,c=INK,**kw):g.add_node(key,label=label,color=c,**kw)
def edge(g,a,b,label='',c=MUT,style='solid',**kw):g.add_edge(a,b,label=label,color=c,fontcolor=c,style=style,**kw)
def rank(g,ids):g.add_subgraph(ids,rank='same')
def structure(smiles,c,formula,parent,branch):
 m=Chem.MolFromSmiles(smiles);assert m
 assert rdMolDescriptors.CalcMolFormula(m)==formula
 acid=m.GetSubstructMatch(Chem.MolFromSmarts('[CX3](=[OX1])[OX2H1]'))[0]
 cs=[a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum()==6]
 ends=[i for i in cs if i!=acid and sum(n.GetAtomicNum()==6 for n in m.GetAtomWithIdx(i).GetNeighbors())==1]
 path=max((Chem.GetShortestPath(m,acid,i) for i in ends),key=len)
 assert len(path)==parent and len(cs)==parent+branch
 off=[i for i in cs if i not in path]
 bp=[i for i in path if sum(n.GetAtomicNum()==6 for n in m.GetAtomWithIdx(i).GetNeighbors())==3]
 assert len(bp)==1 and path.index(bp[0])+1==4 and len(off)==branch
 rdDepictor.Compute2DCoords(m)
 d=rdMolDraw2D.MolDraw2DSVG(340,142);opt=d.drawOptions();opt.clearBackground=False;opt.useBWAtomPalette();opt.atomLabels[bp[0]]='C4';opt.bondLineWidth=1.8
 rgb=tuple(int(c[i:i+2],16)/255 for i in (1,3,5)); ids=off+bp
 d.DrawMolecule(m,highlightAtoms=ids,highlightAtomColors={i:rgb for i in ids});d.FinishDrawing()
 s=d.GetDrawingText();s=s[s.index('<svg'):]
 checks.append(dict(smiles=smiles,formula=formula,total_carbons=len(cs),parent_chain=parent,branch_carbons=branch,branch_position=4,stereochemistry='unspecified'))
 return ET.fromstring(s)
def svgtext(g,x,y,label,size=16,c=INK,bold=False):
 e=ET.SubElement(g,'{'+NS+'}text',{'x':str(x),'y':str(y),'text-anchor':'middle','font-family':'DejaVu Sans','font-size':str(size),'fill':c})
 if bold:e.set('font-weight','bold')
 e.text=label

g=graph()
r=cluster(g,'rumen','01  RUMEN · MICROBIAL FERMENTATION','#f1f7fa')
node(r,'diet','Dietary carbohydrates',width='4.7')
for i,(name,c) in enumerate(zip(['Acetate · C2','Propionate · C3','Butyrate · C4'],COLORS)):
 node(r,'sc'+str(i),name,c);edge(g,'diet','sc'+str(i),c=c)
rank(g,['sc0','sc1','sc2'])
e=cluster(g,'epithelium','02  RUMEN EPITHELIUM · ABSORPTION & HOST HANDLING','#f6f8fa')
for i,(label,c) in enumerate(zip(['Absorption → portal blood\n→ peripheral acetate supply','Absorption → portal blood\nMajor hepatic uptake / gluconeogenesis','Absorption + epithelial metabolism\nSubstantial β-hydroxybutyrate formation'],COLORS)):
 node(e,'ep'+str(i),label,c);edge(g,'sc'+str(i),'ep'+str(i),c=c)
rank(g,['ep0','ep1','ep2'])
a=cluster(g,'adipocyte','03  ADIPOCYTE CYTOSOL · ACC PRECURSOR NETWORK','#f1f8f4')
for i,(label,c) in enumerate(zip(['Acetate → acetyl-CoA\nAcetyl-CoA synthetase','Propionyl-CoA\nCytosolic supply unverified','Butyryl-CoA\nCytosolic supply unverified'],COLORS)):
 node(a,'coa'+str(i),label,c);edge(g,'ep'+str(i),'coa'+str(i),'Indirect supply / activation',c,'dotted')
rank(g,['coa0','coa1','coa2'])
for i,(label,c) in enumerate(zip(['Malonyl-CoA\nCanonical 2-carbon extension donor','Methylmalonyl-CoA\nCandidate methyl-branch donor','Ethylmalonyl-CoA\nCandidate ethyl-branch donor'],[GREEN,COLORS[1],COLORS[2]])):
 node(a,'don'+str(i),label,c)
 edge(g,'coa'+str(i),'don'+str(i),'ACC' if i==0 else 'ACC side reaction · candidate',c,'solid' if i==0 else 'dashed')
rank(g,['don0','don1','don2'])

f=graph();fc=cluster(f,'fasn','04  CANDIDATE FASN ASSEMBLY · THREE CARBON-BALANCED ROUTES','#fafbfd')
rows=[('MOA','Acetyl-CoA · C2','Butyryl-ACP · C4','2-Methylhexanoyl-ACP · C7\nBranch at C2','4-Methyloctanoyl-ACP · C9\nBranch at C4','methylmalonyl-CoA','3'),('MNA','Propionyl-CoA · C3','Pentanoyl-ACP · C5','2-Methylheptanoyl-ACP · C8\nBranch at C2','4-Methylnonanoyl-ACP · C10\nBranch at C4','methylmalonyl-CoA','3'),('EOA','Acetyl-CoA · C2','Butyryl-ACP · C4','2-Ethylhexanoyl-ACP · C8\nBranch at C2','4-Ethyloctanoyl-ACP · C10\nBranch at C4','ethylmalonyl-CoA','4')]
for i,(row,c) in enumerate(zip(rows,COLORS)):
 name,p,b,br,end,don,n=row
 for j,label in enumerate([name+'\n'+p,b,br,end]):node(fc,f'n{i}_{j}',label,c)
 for j,label in enumerate(['+ malonyl-CoA · +2 C retained','+ '+don+'\n+'+n+' C retained after CO₂ loss','+ malonyl-CoA · +2 C\nBranch C2 → C4']):edge(f,f'n{i}_{j}',f'n{i}_{j+1}',label,c,'dashed')
 node(fc,f'product{i}','',c,width='4.9',height='2.8',fixedsize='true')
 edge(f,f'n{i}_3',f'product{i}','Candidate early hydrolysis / release',c,'dashed')
 # Verify stated route carbon arithmetic independently of the product SMILES.
 pcount=[2,3,2][i]; assert pcount+2+int(n)+2==[9,10,10][i]
for j in range(4):rank(f,[f'n{i}_{j}' for i in range(3)])
rank(f,['product0','product1','product2'])
products=[('MOA · 4-Methyloctanoic acid','O=C(O)CCC(C)CCCC','C9H18O2',8,1),('MNA · 4-Methylnonanoic acid','O=C(O)CCC(C)CCCCC','C10H20O2',9,1),('EOA · 4-Ethyloctanoic acid','O=C(O)CCC(CC)CCCC','C10H20O2',8,2)]
roots=[]
layout_checks=[]

for name,gr in [('precursor_layout',g),('fasn_layout',f)]:
 gr.write(str(OUT/(name+'.dot')));gr.layout(prog='dot');gr.write(str(OUT/(name+'_positioned.dot')))
 root=ET.fromstring(gr.draw(format='svg').decode())
 # Verify actual dot node boxes do not overlap; this is independent of SVG scale.
 boxes=[]
 for n in gr.nodes():
  px,py=map(float,n.attr['pos'].split(','));nw=float(n.attr['width'])*72;nh=float(n.attr['height'])*72
  boxes.append((str(n),(px-nw/2,py-nh/2,px+nw/2,py+nh/2)))
 overlaps=[]
 for ii,(na,aa) in enumerate(boxes):
  for nb,bb in boxes[ii+1:]:
   if min(aa[2],bb[2])>max(aa[0],bb[0])+0.1 and min(aa[3],bb[3])>max(aa[1],bb[1])+0.1:overlaps.append([na,nb])
 assert not overlaps,overlaps
 layout_checks.append({'panel':name,'node_overlap_pairs':overlaps})
 # Mask lines behind edge labels / cluster titles using SVG text geometry.
 for group in root.iter('{'+NS+'}g'):
  cl=group.get('class')
  if cl not in ('edge','cluster'):continue
  title=group.find('{'+NS+'}title'); title=title.text if title is not None else ''
  fill={'cluster_rumen':'#f1f7fa','cluster_epithelium':'#f6f8fa','cluster_adipocyte':'#f1f8f4','cluster_fasn':'#fafbfd'}.get(title,'#ffffff')
  for t in list(group):
   if t.tag!='{'+NS+'}text':continue
   size=float(t.get('font-size','14'));font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',round(size*10))
   width=font.getlength(t.text or '')/10;tx=float(t.get('x'));ty=float(t.get('y'));anchor=t.get('text-anchor','start')
   left=tx-width/2 if anchor=='middle' else tx-width if anchor=='end' else tx
   bg=ET.Element('{'+NS+'}rect',{'x':str(left-2),'y':str(ty-size*.84),'width':str(width+4),'height':str(size*1.12),'fill':fill,'stroke':'none'})
   group.insert(list(group).index(t),bg)

 if name=='fasn_layout':
  for group in root.iter('{'+NS+'}g'):
   title=group.find('{'+NS+'}title')
   if title is None or not re.fullmatch('product[0-2]',title.text or ''):continue
   i=int(title.text[-1]);pn=gr.get_node(title.text);cx,cy=map(float,pn.attr['pos'].split(','));cy=-cy
   # Node positions and RDKit placement use dot's final geometry, not manual coordinates.
   titletext,sm,formula,pa,bl=products[i]
   svgtext(group,cx,cy-76,titletext,16,COLORS[i],True)
   mol=structure(sm,COLORS[i],formula,pa,bl);mol.set('x',str(cx-170));mol.set('y',str(cy-66));group.append(mol)
   svgtext(group,cx,cy+88,formula+' · C4 branch · configuration unspecified',12,MUT)
 roots.append(root)
# SVG editorial layer: consistent panel width, page title, context notes and legend.
W=1800;page=ET.Element('{'+NS+'}svg',{'width':str(W)})
ET.SubElement(page,'{'+NS+'}rect',{'width':'100%','height':'100%','fill':'white'})
svgtext(page,900,48,'Candidate host biosynthesis of mutton odor-associated fatty acids',31,INK,True)
svgtext(page,900,86,'MOA · MNA · EOA  |  Rumen supply, epithelial handling and adipocyte ACC–FASN synthesis',21,MUT)
y=112
for i,root in enumerate(roots):
 vb=list(map(float,root.attrib['viewBox'].split())); h=vb[3]*1720/vb[2]
 root.set('x','40');root.set('y',str(y));root.set('width','1720');root.set('height',str(h));page.append(root);y+=h+32
 if i==0:
  svgtext(page,900,y,'ACC uses ATP + HCO₃⁻. CoA pools are compartmentalized; adipose branched donor supply is unverified.',20,MUT);y+=38
  svgtext(page,900,y,'FASN model below selects primers and donors from the precursor network; each cycle requires NADPH.',20,MUT);y+=35
svgtext(page,900,y,'MITOCHONDRION: propionyl-CoA → methylmalonyl-CoA → succinyl-CoA  |  PCC · MCEE · MMUT / B12',19,MUT);y+=33
svgtext(page,900,y,'CYTOSOL: ECHDC1 converts ethylmalonyl-CoA → butyryl-CoA; methylmalonyl-CoA → propionyl-CoA (lower activity).',18,MUT);y+=33
svgtext(page,900,y,'Solid: established biochemical / physiological step   |   Dashed: candidate reaction   |   Dotted: indirect / uncertain supply',19,MUT);y+=38
svgtext(page,900,y,'Candidate model, not demonstrated sheep flux. FASN donor acceptance and early release require validation.',21,INK,True);y+=30
page.set('height',str(y));page.set('viewBox',f'0 0 {W} {y}')
svg=ET.tostring(page,encoding='utf-8',xml_declaration=True)
base=OUT/'BCFA_candidate_v3';base.with_suffix('.svg').write_bytes(svg)
cairosvg.svg2pdf(bytestring=svg,write_to=str(base.with_suffix('.pdf')))
cairosvg.svg2png(bytestring=svg,write_to=str(base.with_suffix('.png')),output_width=5400,output_height=round(y*3))
report={'structures':checks,'route_carbon_arithmetic':'passed','layout_checks':layout_checks,'layout_engine':'PyGraphviz / Graphviz dot','versions':{p:version(p) for p in ['pygraphviz','rdkit','CairoSVG','Pillow']},'graph_nodes':sum(len(x.nodes()) for x in [g,f]),'graph_edges':sum(len(x.edges()) for x in [g,f]),'png_size':[5400,round(y*3)],'sha256':{ext:hashlib.sha256(base.with_suffix('.'+ext).read_bytes()).hexdigest() for ext in ['svg','pdf','png']}}
(OUT/'validation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
