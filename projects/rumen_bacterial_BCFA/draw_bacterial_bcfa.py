"""One common FASN pathway; only R/X selections and final acids differ.
Run --lang both. All explanatory notes live in the separate figure caption.
"""
from pathlib import Path
import argparse,json,math,re,hashlib
from importlib.metadata import version
import xml.etree.ElementTree as ET
from PIL import Image,ImageFont
import pygraphviz as pgv
import cairosvg
from rdkit import Chem
from rdkit.Geometry import Point2D
from rdkit.Chem import rdDepictor,rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
OUT=Path(__file__).resolve().parent/'outputs';OUT.mkdir(exist_ok=True)
NS='http://www.w3.org/2000/svg';ET.register_namespace('',NS)
W,H=2200,1390
INK='#263b47';MUT='#567079';GREEN='#24776b';BLUE='#286fa3';PURPLE='#8054a0';ORANGE='#b96324';COLORS=[BLUE,PURPLE,ORANGE]
FONT='Noto Sans CJK SC';FP='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc';FB='/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc'
# Only these three R/X combinations are depicted; this is not a factorial model.
PRODUCTS=[('MOA','4-甲基辛酸','4-Methyloctanoic acid','O=C(O)CCC(C)CCCC','C9H18O2',8,1,'CH3','CH3'),('MNA','4-甲基壬酸','4-Methylnonanoic acid','O=C(O)CCC(C)CCCCC','C10H20O2',9,1,'CH2CH3','CH3'),('EOA','4-乙基辛酸','4-Ethyloctanoic acid','O=C(O)CCC(CC)CCCC','C10H20O2',8,2,'CH3','CH2CH3')]
GENERIC=['[1*]C(=O)S[2*]','[1*]CCC(=O)S[2*]','[1*]CCCC([3*])C(=O)S[2*]','[1*]CCCC([3*])CCC(=O)S[2*]','[1*]CCCC([3*])CCC(=O)O']

def build(lang):
 L=lambda zh,en:zh if lang=='zh' else en
 root=ET.Element('{'+NS+'}svg',{'width':str(W),'height':str(H),'viewBox':f'0 0 {W} {H}'});editorial=[];check=[]
 def rect(p,x,y,w,h,fill='white',stroke='none',r=12):
  return ET.SubElement(p,'{'+NS+'}rect',{'x':str(x),'y':str(y),'width':str(w),'height':str(h),'rx':str(r),'fill':fill,'stroke':stroke})
 def text(p,x,y,s,size=28,c=INK,bold=False,anchor='middle',record=True):
  e=ET.SubElement(p,'{'+NS+'}text',{'x':str(x),'y':str(y),'font-family':FONT,'font-size':str(size),'font-weight':'bold' if bold else 'normal','fill':c,'text-anchor':anchor});e.text=s
  if record:
   f=ImageFont.truetype(FB if bold else FP,size*4,index=2);width=f.getlength(s)/4;b=f.getbbox(s,anchor='ls');left=x if anchor=='start' else x-width if anchor=='end' else x-width/2
   editorial.append({'text':s,'bounds':[left,y+b[1]/4,left+width,y+b[3]/4]})
  return e
 def arrow(x1,y,x2,c=MUT,style='solid'):
  a={'d':f'M{x1},{y} L{x2},{y}','stroke':c,'stroke-width':'2.6','fill':'none'}
  if style!='solid':a['stroke-dasharray']='10 7' if style=='dashed' else '2 7'
  ET.SubElement(root,'{'+NS+'}path',a);ET.SubElement(root,'{'+NS+'}polygon',{'points':f'{x2},{y} {x2-12},{y-4.5} {x2-12},{y+4.5}','fill':c})
 def molecule(smiles,width,height,carrier=None,product=None):
  m=Chem.MolFromSmiles(smiles);assert m
  highlights={};labs={}
  if carrier:
   r=next(a.GetIdx() for a in m.GetAtoms() if a.GetIsotope()==1)
   dummies={a.GetIsotope():a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum()==0}
   for isotope,idx in dummies.items():labs[idx]={1:'R',2:carrier,3:'X'}[isotope]
   carbonyl=m.GetSubstructMatch(Chem.MolFromSmarts('[C](=O)[S,O]'))[0]
   backbone=Chem.GetShortestPath(m,r,carbonyl)
   if 3 in dummies:
    highlights[dummies[3]]=ORANGE
    bp=m.GetAtomWithIdx(dummies[3]).GetNeighbors()[0].GetIdx()
    assert len(Chem.GetShortestPath(m,carbonyl,bp))==(2 if GENERIC.index(smiles)==2 else 4)
   highlights[r]=BLUE
   # In the shared scheme, R and X are variable groups, not carbon atoms.
   nc=sum(a.GetAtomicNum()==6 for a in m.GetAtoms())
   expected=[1,3,5,7,7][GENERIC.index(smiles)];assert nc==expected
  else:
   name,_,_,_,formula,parent,bl,rv,xv=product
   assert rdMolDescriptors.CalcMolFormula(m)==formula
   carbonyl=m.GetSubstructMatch(Chem.MolFromSmarts('[CX3](=[OX1])[OX2H1]'))[0]
   cs=[a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum()==6]
   ends=[idx for idx in cs if idx!=carbonyl and sum(n.GetAtomicNum()==6 for n in m.GetAtomWithIdx(idx).GetNeighbors())==1]
   acid_path=max((Chem.GetShortestPath(m,carbonyl,i) for i in ends),key=len)
   off=[i for i in cs if i not in acid_path]
   bp=[i for i in acid_path if sum(n.GetAtomicNum()==6 for n in m.GetAtomWithIdx(i).GetNeighbors())==3]
   assert len(acid_path)==parent and len(cs)==parent+bl and len(off)==bl and len(bp)==1 and acid_path.index(bp[0])+1==4
   backbone=list(reversed(acid_path));r=backbone[0]
   c=COLORS[[p[0] for p in PRODUCTS].index(name)]
   highlights.update({i:c for i in off+bp});labs[bp[0]]='C4'
   check.append({'name':name,'formula':formula,'total_carbons':len(cs),'parent_chain':parent,'branch_carbons':bl,'branch_position':4,'R':rv,'X':xv,'stereochemistry':'unspecified'})
  coordmap={idx:Point2D(k*1.3,(k%2)*.62) for k,idx in enumerate(backbone)}
  if not carrier:
   # Keep all three real branch groups on the same side of their parent chain.
   bpoint=coordmap[bp[0]]
   branch_path=sorted(off,key=lambda i:len(Chem.GetShortestPath(m,bp[0],i)))
   for k,idx in enumerate(branch_path):coordmap[idx]=Point2D(bpoint.x-(1.20 if k else 0),bpoint.y-1.40-.65*k)
  rdDepictor.Compute2DCoords(m,canonOrient=False,coordMap=coordmap)
  d=rdMolDraw2D.MolDraw2DSVG(width,height);o=d.drawOptions();o.clearBackground=False;o.useBWAtomPalette();o.bondLineWidth=2.0;o.minFontSize=17;o.maxFontSize=24
  for i,lab in labs.items():o.atomLabels[i]=lab
  rgb={i:tuple(int(c[j:j+2],16)/255 for j in (1,3,5)) for i,c in highlights.items()}
  d.DrawMolecule(m,highlightAtoms=list(highlights),highlightAtomColors=rgb);d.FinishDrawing();s=d.GetDrawingText();return ET.fromstring(s[s.index('<svg'):])
 rect(root,0,0,W,H)
 text(root,1100,60,L('瘤胃细菌合成羊膻味脂肪酸的机制假说','Hypothesized rumen bacterial synthesis of odor-associated fatty acids'),42,INK,True)
 rect(root,30,105,2140,315,'#f2f7f5','#cfddd7')
 text(root,60,150,L('1  瘤胃细菌：引物与延长供体','1  Rumen bacteria: primers and extender supply'),28,INK,True,'start')
 text(root,350,225,L('乙酰-CoA / 丙酰-CoA','Acetyl-CoA / propionyl-CoA'),30,BLUE,True)
 text(root,350,325,L('起始引物 R','Starter group R'),28,BLUE)
 text(root,1050,220,L('丙二酰-CoA','Malonyl-CoA'),30,GREEN,True)
 text(root,1050,310,L('甲基丙二酰-CoA','Methylmalonyl-CoA'),30,PURPLE,True)
 text(root,1660,215,L('巴豆酰-CoA + CO2 + NADPH','Crotonyl-CoA + CO2 + NADPH'),27,ORANGE)
 text(root,1660,265,'CCR',28,ORANGE,True)
 ET.SubElement(root,'{'+NS+'}polygon',{'points':'1660,316 1655,304 1665,304','fill':ORANGE})
 text(root,1660,345,L('乙基丙二酰-CoA','Ethylmalonyl-CoA'),30,ORANGE,True)
 ET.SubElement(root,'{'+NS+'}path',{'d':'M1660,280 L1660,310','stroke':ORANGE,'stroke-width':'3'})
 text(root,60,475,L('2  PKS / 未鉴定合成系统：共享碳骨架模型','2  PKS / unidentified synthase: shared carbon assembly model'),30,INK,True,'start')
 # One directed path, laid out with dot. Variants never become separate routes.
 g=pgv.AGraph(directed=True)
 g.graph_attr.update(rankdir='LR',bgcolor='transparent',nodesep='.2',ranksep='.30',pad='.05',splines='line',outputorder='edgesfirst')
 g.node_attr.update(shape='box',label='',fixedsize='true',style='invis',height='1.48')
 g.edge_attr.update(style='dashed',color=INK,penwidth='1.6',arrowsize='.7',fontname=FONT,fontsize='17',fontcolor=GREEN)
 widths=[2.1,2.9,3.3,3.8,3.8]
 for j,w in enumerate(widths):g.add_node(f's{j}',width=str(w))
 edge_labels=[L('+ 丙二酰-CoA','+ malonyl-CoA'),L('+ 支链供体','+ branched donor'),L('+ 丙二酰-CoA','+ malonyl-CoA'),L('水解','Hydrolysis')]
 for j,label in enumerate(edge_labels):g.add_edge(f's{j}',f's{j+1}',label=label,fontcolor=ORANGE if j==1 else INK if j==3 else GREEN)
 g.write(str(OUT/f'shared_{lang}.dot'));g.layout(prog='dot');g.write(str(OUT/f'shared_{lang}_positioned.dot'))
 sroot=ET.fromstring(g.draw(format='svg').decode());vb=list(map(float,sroot.get('viewBox').split()));sc=2140/vb[2]
 gg=next(x for x in sroot if x.tag=='{'+NS+'}g');tt=re.search(r'translate\(([-.\d]+)[ ,]+([-.\d]+)\)',gg.get('transform'));tx,ty=map(float,tt.groups())
 top=550;actualh=vb[3]*sc
 sroot.set('x','30');sroot.set('y',str(top));sroot.set('width','2140');sroot.set('height',str(actualh))
 heads=[L('起始引物 R','Primer R'),L('模块1：延长','Module 1: extension'),L('模块2：C2支链','Module 2: C2 branch'),L('模块3：C4支链','Module 3: C4 branch'),L('游离脂肪酸','Free fatty acid')]
 coords=[]
 for j in range(5):
  nd=g.get_node(f's{j}');cx,cy=map(float,nd.attr['pos'].split(','));cy=-cy;nw=float(nd.attr['width'])*72
  coords.append((cx,cy,nw));text(root,30+(cx+tx)*sc,522,heads[j],27,INK,True)
  mol=molecule(GENERIC[j],round(nw),104,'CoA' if j==0 else 'ACP' if j<4 else 'OH')
  mol.set('x',str(cx-nw/2));mol.set('y',str(cy-52));gg.append(mol)
 root.append(sroot)
 # Three R/X selections, not three repeated pathways.
 cardtop=math.ceil(top+actualh+38);cardw=690;cardh=1210-cardtop-30
 for i,product in enumerate(PRODUCTS):
  name,zh,en,sm,formula,parent,bl,rv,xv=product;x=30+i*725;c=COLORS[i]
  rect(root,x,cardtop,cardw,cardh,'white',c,16)
  text(root,x+cardw/2,cardtop+45,name+' · '+L(zh,en),29,c,True)
  primer=L('乙酰-CoA','Acetyl-CoA') if rv=='CH3' else L('丙酰-CoA','Propionyl-CoA')
  donor=L('甲基丙二酰-CoA','Methylmalonyl-CoA') if xv=='CH3' else L('乙基丙二酰-CoA','Ethylmalonyl-CoA')
  text(root,x+cardw/2,cardtop+92,primer+'  /  R = '+rv,27,BLUE)
  text(root,x+cardw/2,cardtop+136,donor+'  /  X = '+xv,27,ORANGE)
  mol=molecule(sm,600,round(cardh-159),product=product);mol.set('x',str(x+45));mol.set('y',str(cardtop+149));root.append(mol)
 # Expand generic R/X choices to actual molecules and verify atom connectivity,
 # rather than relying on text labels alone.
 expansions=[]
 for name,zh,en,sm,formula,parent,bl,rv,xv in PRODUCTS:
  r='C' if rv=='CH3' else 'CC';x='C' if xv=='CH3' else 'CC'
  expanded=Chem.MolFromSmiles(r+'CCCC('+x+')CCC(=O)O')
  assert Chem.MolToSmiles(expanded)==Chem.MolToSmiles(Chem.MolFromSmiles(sm))
  expansions.append({'name':name,'R':rv,'X':xv,'canonical_smiles':Chem.MolToSmiles(expanded)})
 # Generic acyl carbon counts exclude R and X: 1 → 3 → 5 → 7 → 7.
 # R adds 1/2 carbons; X adds 1/2. Validate each stage for all combinations.
 for name,zh,en,sm,formula,parent,bl,rv,xv in PRODUCTS:
  nr=1 if rv=='CH3' else 2;nx=1 if xv=='CH3' else 2
  seq=[nr+1,nr+3,nr+5+nx,nr+7+nx]
  assert seq=={'MOA':[2,4,7,9],'MNA':[3,5,8,10],'EOA':[2,4,8,10]}[name]
 overlap=[]
 for i,a in enumerate(editorial):
  for b in editorial[i+1:]:
   aa,bb=a['bounds'],b['bounds']
   if min(aa[2],bb[2])-max(aa[0],bb[0])>1 and min(aa[3],bb[3])-max(aa[1],bb[1])>1:overlap.append([a['text'],b['text']])
 assert not overlap,overlap
 assert all(a['bounds'][0]>=0 and a['bounds'][2]<=W and a['bounds'][1]>=0 and a['bounds'][3]<=H for a in editorial)
 text(root,1100,1250,L('3  目标酸进入宿主的两条可能路径','3  Two possible routes of bacterial products into the host'),29,INK,True)
 text(root,570,1310,L('游离酸 → 瘤胃上皮吸收','Free acids → rumen epithelial uptake'),28,BLUE)
 text(root,1610,1310,L('微生物脂质 → 下游消化 / 吸收','Microbial lipids → downstream digestion / uptake'),28,PURPLE)
 text(root,1100,1360,L('→ 宿主利用与脂肪组织积累','→ Host utilization and adipose accumulation'),28,INK,True)
 content=ET.tostring(root,encoding='utf-8',xml_declaration=True);base=OUT/f'BCFA_bacterial_{lang}'
 base.with_suffix('.svg').write_bytes(content)
 cairosvg.svg2pdf(bytestring=content,write_to=str(base.with_suffix('.pdf')))
 cairosvg.svg2png(bytestring=content,write_to=str(base.with_suffix('.png')),output_width=6600,output_height=H*3)
 im=Image.open(base.with_suffix('.png'));im.thumbnail((1800,1100));im.save(OUT/f'preview_{lang}.png')
 report={'language':lang,'shared_reaction_path':{'nodes':5,'edges':4,'generic_smiles':GENERIC,'acyl_carbons_excluding_R_X':[1,3,5,7,7]},'products':check,'R_X_expansion_checks':expansions,'editorial_label_count':len(editorial),'editorial_overlap_pairs':overlap,'png_size':[6600,H*3],'versions':{p:version(p) for p in ['pygraphviz','rdkit','CairoSVG','Pillow']},'sha256':{e:hashlib.sha256(base.with_suffix('.'+e).read_bytes()).hexdigest() for e in ['svg','pdf','png']}}
 (OUT/f'validation_{lang}.json').write_text(json.dumps(report,ensure_ascii=False,indent=2));print(lang,cardtop,cardh,'validated and exported')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--lang',choices=['zh','en','both'],default='both');args=p.parse_args()
 for lang in ['zh','en'] if args.lang=='both' else [args.lang]:build(lang)
