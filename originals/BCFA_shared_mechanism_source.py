from pathlib import Path
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from rdkit import Chem
from rdkit.Chem import rdDepictor, rdMolDescriptors

OUT=Path(__file__).resolve().parent
mpl.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','pdf.fonttype':42})
fig,ax=plt.subplots(figsize=(20,12.5)); ax.set(xlim=(0,100),ylim=(0,100)); ax.axis('off')
ink='#263941'; muted='#5b6d75'; blue='#216ca3'; purple='#79509c'; orange='#bd641c'; green='#287b70'
def txt(x,y,s,n=12,c=ink,w='normal',ha='center',**kw): ax.text(x,y,s,fontsize=n,color=c,weight=w,ha=ha,va='center',**kw)
def arrow(a,b,c=ink,ls='-',lw=1.4): ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=12,color=c,lw=lw,linestyle=ls))
def rule(y): ax.plot([2,98],[y,y],lw=.7,color='#cbd5d8')
def structure(smiles,x,y,scale,color,formula):
    m=Chem.MolFromSmiles(smiles); rdDepictor.Compute2DCoords(m); assert rdMolDescriptors.CalcMolFormula(m)==formula
    assert len(Chem.GetShortestPath(m,1,5))==4 and m.GetAtomWithIdx(5).GetDegree()==3
    conf=m.GetConformer(); p=[(conf.GetAtomPosition(i).x,conf.GetAtomPosition(i).y) for i in range(m.GetNumAtoms())]
    mx=(min(a for a,b in p)+max(a for a,b in p))/2; my=(min(b for a,b in p)+max(b for a,b in p))/2
    p=[(x+(a-mx)*scale,y+(b-my)*scale) for a,b in p]
    for b in m.GetBonds():
        i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx(); x1,y1=p[i]; x2,y2=p[j]; pair={i,j}
        col=color if pair=={5,6} or (pair=={6,7} and 'CCC(CC)' in smiles) else ink
        ax.plot([x1,x2],[y1,y2],lw=2.1,color=col)
        if b.GetBondType()==Chem.BondType.DOUBLE:
            dx,dy=x2-x1,y2-y1; d=(dx*dx+dy*dy)**.5
            ax.plot([x1-dy/d*.22,x2-dy/d*.22],[y1+dx/d*.22,y2+dx/d*.22],lw=1.5,color=ink)
    for a in m.GetAtoms():
        if a.GetSymbol()=='O':
            u,v=p[a.GetIdx()]; txt(u,v,'OH' if a.GetTotalNumHs() else 'O',12,bbox=dict(facecolor='white',edgecolor='none',pad=.3))
    u,v=p[5]; ax.scatter([u],[v],s=38,c=color,zorder=5); txt(u-2 if 'CCC(CC)' not in smiles else u,v+1.7,'C4',10,color,'bold')

txt(50,97,'Proposed host synthesis of odor-associated branched-chain fatty acids',20,ink,'bold')
txt(2,91,'PRECURSOR SUPPLY',12,muted,'bold',ha='left')
txt(3,86,'Rumen fermentation',13,ink,'bold',ha='left')
txt(32,89,'Propionate',12); txt(32,82,'Butyrate',12)
arrow((18,86.8),(27,89)); arrow((18,85.2),(27,82))
txt(52,89,'Propionyl-CoA',12); txt(52,82,'Butyryl-CoA',12)
arrow((38,89),(45,89),muted,':'); arrow((38,82),(45,82),muted,':')
txt(79,89,'Methylmalonyl-CoA',12,purple,'bold'); txt(79,82,'Ethylmalonyl-CoA',12,orange,'bold')
arrow((59,89),(69,89),muted,'--'); arrow((59,82),(69,82),muted,'--')
txt(64,91.2,'ACC*',10,muted); txt(64,84.2,'ACC*',10,muted)
txt(50,76.5,'Acetyl-CoA  ── ACC ──▶  Malonyl-CoA   ·   shared two-carbon extension donor',12,green,'bold')
rule(72.8)
txt(2,69.8,'ONE SHARED FASN ASSEMBLY MODEL',12,muted,'bold',ha='left')
txt(12,65.6,'1  Primer',12,ink,'bold'); txt(36.5,65.6,'2  First extension',12,ink,'bold'); txt(62,65.6,'3  Branch introduction',12,ink,'bold'); txt(88,65.6,'4  Final extension',12,ink,'bold')
txt(12,59.8,'R–C(=O)–SCoA',14,ink,'bold'); txt(36.5,59.8,'R–CH₂–CH₂–C(=O)–SACP',12,ink,'bold')
txt(62,59.8,'R–(CH₂)₃–CH(X)–C(=O)–SACP',12,ink,'bold')
txt(88,59.8,'R–(CH₂)₃–CH(X)–CH₂–CH₂',11,ink,'bold'); txt(88,56.7,'–C(=O)–SACP',12,ink,'bold')
arrow((20,59.8),(25,59.8),muted,'--'); arrow((48,59.8),(51,59.8),muted,'--'); arrow((74,59.8),(77,59.8),muted,'--')
txt(22.5,62.5,'+ malonyl',10,green); txt(49.5,63,'+ branched donor',10,purple); txt(75.5,63,'+ malonyl',10,green)
txt(50,52.5,'C2 branch  →  two-carbon extension at the carboxyl end  →  C4 branch',14,ink,'bold')
txt(50,48.8,'Proposed thioester hydrolysis → R–(CH₂)₃–CH(X)–CH₂–CH₂–COOH',13,ink,'bold')
txt(50,45.8,'R selects the primer; X selects the branch. Only the three combinations below are proposed.',11,muted)
rule(43)
txt(2,40.5,'PRODUCTS AND THE ONLY THREE COMBINATIONS SHOWN',12,muted,'bold',ha='left')
items=[(17,blue,'MOA  ·  C9','R = CH₃   |   X = CH₃','acetyl primer + methylmalonyl','O=C(O)CCC(C)CCCC','C9H18O2'),(50,purple,'MNA  ·  C10','R = CH₂CH₃   |   X = CH₃','propionyl primer + methylmalonyl','O=C(O)CCC(C)CCCCC','C10H20O2'),(83,orange,'EOA  ·  C10','R = CH₃   |   X = CH₂CH₃','acetyl primer + ethylmalonyl','O=C(O)CCC(CC)CCCC','C10H20O2')]
for x,c,name,choice,desc,sm,formula in items:
    txt(x,37,name,15,c,'bold'); txt(x,33.8,choice,12,c,'bold'); txt(x,31.2,desc,10,muted)
    structure(sm,x,24,1.35,c,formula); txt(x,17.8,formula+'  ·  configuration unspecified',10,muted)
rule(15.7)
txt(2,13,'MITOCHONDRION',11,muted,'bold',ha='left')
txt(2,10.6,'Propionyl-CoA  —PCC→  (2S)-methylmalonyl-CoA  —MCEE→  (2R)-methylmalonyl-CoA  —MMUT / B12→  succinyl-CoA',10,ink,ha='left')
txt(2,7.9,'CYTOSOLIC REPAIR',11,muted,'bold',ha='left')
txt(2,5.6,'ECHDC1: ethylmalonyl-CoA → butyryl-CoA; methylmalonyl-CoA → propionyl-CoA (lower activity); limits branched donors.',10,ink,ha='left')
txt(2,2.4,'Solid: enzymatic basis   ·   Dashed: candidate assembly / ACC side reaction   ·   Dotted: indirect supply',10,muted,ha='left')
txt(98,2.4,'* Sheep cytosolic supply and early release remain unverified',9,muted,ha='right')
fig.subplots_adjust(left=.025,right=.975,top=.975,bottom=.025)
for ext in ['svg','pdf','png']: fig.savefig(OUT/f'BCFA_shared_mechanism.{ext}',dpi=300,facecolor='white')
print('RDKit formula and C4 checks passed')
