"""RDKit-native, vector-first candidate pathway. Run: python draw_bcfa_v2.py.
No supplied original is executed or overwritten. Product stereochemistry unspecified.
"""
from pathlib import Path
import json, html, re, hashlib, platform
from importlib.metadata import version
import xml.etree.ElementTree as ET
from PIL import ImageFont
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
import cairosvg
OUT=Path(__file__).resolve().parent
W,H=1800,2340
INK='#243746'; MUT='#566875'; GREEN='#26766b'; BLUE='#256da6'; PURPLE='#8052a0'; ORANGE='#b96321'
parts=[]; labels=[]; results=[]
def add(s): parts.append(s)
def rect(x,y,w,h,fill='#ffffff',stroke='#ccd7dd',r=14):
 add(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" fill="{fill}" stroke="{stroke}"/>')
def text(x,y,s,size=22,color=INK,bold=False,anchor='middle',bg=None):
 s=str(s)
 if bg:
  f=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',size); tw=f.getlength(s)
  rect(x-tw/2-5,y-size,tw+10,size+5,bg,bg,0)
 add(f'<text x="{x}" y="{y}" fill="{color}" font-family="DejaVu Sans" font-size="{size}" font-weight="{"bold" if bold else "normal"}" text-anchor="{anchor}">{html.escape(s)}</text>')
 font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans'+('-Bold' if bold else '')+'.ttf',size)
 width=font.getlength(s); b=font.getbbox(s,anchor='ls'); left=x if anchor=='start' else x-width/2 if anchor=='middle' else x-width
 labels.append(dict(text=s,bounds=[left,y+b[1],left+width,y+b[3]]))
def arrow(x1,y1,x2,y2,color=MUT,kind='solid'):
 dash={'solid':'','candidate':'9 7','supply':'2 7'}[kind]
 add(f'<path d="M{x1},{y1} L{x2},{y2}" stroke="{color}" stroke-width="2.3" stroke-dasharray="{dash}" fill="none"/>')
 angle=__import__('math').atan2(y2-y1,x2-x1)
 import math
 pts=[(x2,y2),(x2-13*math.cos(angle)+5*math.sin(angle),y2-13*math.sin(angle)-5*math.cos(angle)),(x2-13*math.cos(angle)-5*math.sin(angle),y2-13*math.sin(angle)+5*math.cos(angle))]
 add('<polygon points="'+' '.join(f'{x},{y}' for x,y in pts)+f'" fill="{color}"/>')
def heading(y,n,title):
 text(65,y,n,24,GREEN,True,'start');text(110,y,title,25,INK,True,'start')
def node(cx,y,s,sub=None,color=INK):
 rect(cx-245,y-32,490,76 if sub else 60,stroke=color)
 text(cx,y+5,s,23,color,True)
 if sub:text(cx,y+32,sub,18,MUT)
def mol(smiles,x,y,color,expected,main,branchlen):
 m=Chem.MolFromSmiles(smiles);assert m is not None
 formula=rdMolDescriptors.CalcMolFormula(m); assert formula==expected
 # Identify carbonyl acid carbon by chemistry, never a hard-coded atom index.
 acid=m.GetSubstructMatch(Chem.MolFromSmarts('[CX3](=[OX1])[OX2H1]'))[0]
 carbons=[a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum()==6]
 ends=[i for i in carbons if sum(n.GetAtomicNum()==6 for n in m.GetAtomWithIdx(i).GetNeighbors())==1 and i!=acid]
 paths=[Chem.GetShortestPath(m,acid,i) for i in ends]
 chain=max(paths,key=len);assert len(chain)==main
 branches=[i for i in chain if sum(n.GetAtomicNum()==6 for n in m.GetAtomWithIdx(i).GetNeighbors())==3]
 assert len(branches)==1 and chain.index(branches[0])+1==4
 branch=[i for i in carbons if i not in chain];assert len(branch)==branchlen
 assert len(carbons)==main+branchlen
 rdDepictor.Compute2DCoords(m)
 d=rdMolDraw2D.MolDraw2DSVG(490,200);o=d.drawOptions();o.clearBackground=False;o.bondLineWidth=2.3;o.minFontSize=18;o.maxFontSize=24
 o.atomLabels[branches[0]]='C4';o.useBWAtomPalette()
 rgb=tuple(int(color[i:i+2],16)/255 for i in (1,3,5))
 d.DrawMolecule(m,highlightAtoms=[branches[0]]+branch,highlightAtomColors={i:rgb for i in [branches[0]]+branch});d.FinishDrawing()
 svg=d.GetDrawingText();svg=svg[svg.index('<svg'):];svg=re.sub(r'<svg\b',f'<svg x="{x}" y="{y}"',svg,count=1);add(svg)
 results.append(dict(smiles=smiles,canonical_smiles=Chem.MolToSmiles(m),formula=formula,total_carbons=len(carbons),parent_chain_carbons=main,branch_position=4,branch_carbons=branchlen,stereochemistry='unspecified'))

add(f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:rdkit="http://www.rdkit.org/xml" width="{W}" height="{H}" viewBox="0 0 {W} {H}">')
add('<defs><marker id="head" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 z" fill="context-stroke"/></marker></defs>')
rect(0,0,W,H,'#ffffff','#ffffff',0)
text(900,55,'Candidate host biosynthesis of mutton odor-associated fatty acids',34,INK,True)
text(900,94,'MOA · MNA · EOA  |  Rumen supply, epithelial handling and an adipocyte ACC–FASN model',23,MUT)
rect(40,125,1720,400,'#f3f8fa')
heading(166,'01','RUMEN FERMENTATION → EPITHELIAL ABSORPTION')
text(900,206,'Dietary carbohydrates → microbial fermentation → short-chain fatty acids',25,INK,True)
xs=[320,900,1480];colors=[BLUE,PURPLE,ORANGE]
for x,c,name in zip(xs,colors,['Acetate · C2','Propionate · C3','Butyrate · C4']):
 text(x,256,name,27,c,True);arrow(x,269,x,303,c)
text(320,332,'Epithelial absorption',23,BLUE,True)
text(320,369,'Portal → peripheral circulation',21)
text(900,332,'Epithelial absorption → portal blood',22,PURPLE,True)
text(900,369,'Major hepatic uptake / gluconeogenesis',21)
text(1480,332,'Absorption + epithelial metabolism',22,ORANGE,True)
text(1480,369,'Substantial conversion to β-hydroxybutyrate',20)
for x,c in zip(xs,colors):arrow(x,390,x,428,c,'supply')
text(320,458,'Acetate supply to adipose tissue',22,BLUE,True)
text(900,458,'Adipose propionyl-CoA supply uncertain',21,PURPLE,True)
text(1480,458,'Adipose butyryl-CoA supply uncertain',21,ORANGE,True)
text(900,498,'CoA esters are compartmentalized; no direct rumen-to-cytosol CoA transfer is implied.',20,MUT)
rect(40,550,1720,345,'#f4f9f6')
heading(593,'02','ADIPOCYTE CYTOSOL · ACC PRECURSOR NETWORK')
for x,c,s in zip(xs,colors,['Acetate → acetyl-CoA','Propionyl-CoA','Butyryl-CoA']): text(x,643,s,25,c,True)
text(320,674,'Acetyl-CoA synthetase',19,MUT)
for x,c,k in zip(xs,colors,['solid','candidate','candidate']):arrow(x,691,x,740,c,k)
for x,c,s in zip(xs,[GREEN,PURPLE,ORANGE],['Malonyl-CoA','Methylmalonyl-CoA','Ethylmalonyl-CoA']):text(x,775,s,25,c,True)
for x,s in zip(xs,['ACC: canonical carboxylation','ACC side reaction: candidate','ACC side reaction: candidate']):text(x,813,s,20,MUT)
text(900,861,'ATP + HCO₃⁻ → ADP + Pi   |   Branched donor availability in sheep adipocytes remains unverified.',21,MUT)
rect(40,920,1720,780,'#fcfcfd')
heading(963,'03','CANDIDATE FASN ASSEMBLY · THREE CARBON-BALANCED ROUTES')
text(900,1003,'Each extension includes condensation, reduction, dehydration and reduction; NADPH is required.',21,MUT)
rows=[('MOA',BLUE,'Acetyl-CoA · C2','Butyryl-ACP · C4','2-Methylhexanoyl-ACP · C7','4-Methyloctanoyl-ACP · C9','methylmalonyl-CoA','3'),('MNA',PURPLE,'Propionyl-CoA · C3','Pentanoyl-ACP · C5','2-Methylheptanoyl-ACP · C8','4-Methylnonanoyl-ACP · C10','methylmalonyl-CoA','3'),('EOA',ORANGE,'Acetyl-CoA · C2','Butyryl-ACP · C4','2-Ethylhexanoyl-ACP · C8','4-Ethyloctanoyl-ACP · C10','ethylmalonyl-CoA','4')]
for x,(name,c,primer,first,second,last,donor,retained) in zip(xs,rows):
 text(x,1052,name,30,c,True)
 node(x,1102,primer,'Primer selection / loading',c)
 arrow(x,1153,x,1213,c,'candidate');text(x,1188,'+ malonyl-CoA: +2 C retained',20,GREEN,bg='#fcfcfd')
 node(x,1255,first,None,c)
 arrow(x,1290,x,1361,c,'candidate');text(x,1320,'+ '+donor,20,c,bg='#fcfcfd');text(x,1346,'+'+retained+' C retained after CO₂ loss',18,MUT,bg='#fcfcfd')
 node(x,1405,second,'Branch at C2',c)
 arrow(x,1456,x,1516,c,'candidate');text(x,1491,'+ malonyl-CoA: branch C2 → C4',20,GREEN,bg='#fcfcfd')
 node(x,1559,last,'Branch at C4',c)
 arrow(x,1610,x,1658,c,'candidate')
 text(x,1684,'Candidate early hydrolysis / release',19,c)
# Product cards: native RDKit depictions with an explicitly labelled branch point.
for x,c,name,sm,formula,main,bl in [(320,BLUE,'MOA · 4-Methyloctanoic acid','O=C(O)CCC(C)CCCC','C9H18O2',8,1),(900,PURPLE,'MNA · 4-Methylnonanoic acid','O=C(O)CCC(C)CCCCC','C10H20O2',9,1),(1480,ORANGE,'EOA · 4-Ethyloctanoic acid','O=C(O)CCC(CC)CCCC','C10H20O2',8,2)]:
 rect(x-265,1725,530,275,stroke=c)
 text(x,1763,name,23,c,True);mol(sm,x-245,1770,c,formula,main,bl)
 text(x,1980,formula+' · C4 branch · stereochemistry unspecified',17,MUT)
rect(40,2025,840,180,'#f3f6f8')
rect(900,2025,860,180,'#f7f8f5')
text(460,2064,'MITOCHONDRION · PROPIONATE CATABOLISM',21,INK,True)
text(460,2104,'Propionyl-CoA → methylmalonyl-CoA → succinyl-CoA',20)
text(460,2141,'PCC → MCEE → MMUT / vitamin B12',20,MUT)
text(460,2178,'Separate pool; no direct transfer to cytosolic FASN shown.',18,MUT)
text(1330,2064,'CYTOSOL · ECHDC1 METABOLITE REPAIR',21,INK,True)
text(1330,2104,'Ethylmalonyl-CoA → butyryl-CoA',21)
text(1330,2141,'Methylmalonyl-CoA → propionyl-CoA (lower activity)',20)
text(1330,2178,'Competes with availability of branched extension donors.',18,MUT)
for x,k,s in [(80,'solid','Established biochemical / physiological step'),(700,'candidate','Candidate synthesis / side reaction'),(1260,'supply','Indirect or uncertain supply')]:
 arrow(x,2248,x+65,2248,MUT,k);text(x+80,2255,s,18,MUT,anchor='start')
text(900,2298,'Candidate model, not demonstrated sheep flux. FASN donor acceptance and short-chain release require validation.',22,INK,True)
add('</svg>')
# Check all authored text extents at export dimensions; molecular labels are
# rendered separately by RDKit and inspected in the raster/PDF previews.
collisions=[]
for i,a in enumerate(labels):
 for b in labels[i+1:]:
  p,q=a['bounds'],b['bounds']
  if min(p[2],q[2])-max(p[0],q[0])>1 and min(p[3],q[3])-max(p[1],q[1])>1: collisions.append([a['text'],b['text']])
assert not collisions,collisions
assert all(b['bounds'][0]>=0 and b['bounds'][2]<=W and b['bounds'][1]>=0 and b['bounds'][3]<=H for b in labels)
svg='\n'.join(parts); ET.fromstring(svg)
base=OUT/'BCFA_candidate_v2'
base.with_suffix('.svg').write_text(svg,encoding='utf-8')
cairosvg.svg2pdf(bytestring=svg.encode(),write_to=str(base.with_suffix('.pdf')))
cairosvg.svg2png(bytestring=svg.encode(),write_to=str(base.with_suffix('.png')),output_width=5400,output_height=7020)
report={'structures':results,'text_overlap_pairs':collisions,'text_labels_checked':len(labels),'outputs':{e:hashlib.sha256(base.with_suffix('.'+e).read_bytes()).hexdigest() for e in ['svg','pdf','png']},'versions':{p:version(p) for p in ['rdkit','cairosvg','Pillow']},'python':platform.python_version()}
(OUT/'validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
