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
dark='#273943'; muted='#62747b'; blue='#216ca3'; purple='#79509c'; orange='#c66a20'
def t(x,y,s,fs=10,c=dark,weight='normal',ha='center',**kw): ax.text(x,y,s,fontsize=fs,color=c,weight=weight,ha=ha,va='center',**kw)
def line(x1,y1,x2,y2,c=dark,ls='-',lw=1.2,arrow=True):
    if arrow: ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=9,color=c,lw=lw,linestyle=ls))
    else: ax.plot([x1,x2],[y1,y2],color=c,lw=lw,ls=ls)
def structure(smiles,cx,cy,scale,c):
    m=Chem.MolFromSmiles(smiles); rdDepictor.Compute2DCoords(m)
    conf=m.GetConformer(); p=[(conf.GetAtomPosition(i).x,conf.GetAtomPosition(i).y) for i in range(m.GetNumAtoms())]
    mx=(min(x for x,y in p)+max(x for x,y in p))/2; my=(min(y for x,y in p)+max(y for x,y in p))/2
    p=[(cx+(x-mx)*scale,cy+(y-my)*scale) for x,y in p]
    # C1 is atom 1, C4 is atom 5 in these validated SMILES.
    branch=5
    assert len(Chem.GetShortestPath(m,1,branch))==4 and m.GetAtomWithIdx(branch).GetDegree()==3
    for b in m.GetBonds():
        i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx(); x1,y1=p[i]; x2,y2=p[j]; col=c if set((i,j)) in ({5,6},{6,7}) and ('CCC(CC)' in smiles or set((i,j))=={5,6}) else dark
        ax.plot([x1,x2],[y1,y2],lw=2,color=col)
        if b.GetBondType()==Chem.BondType.DOUBLE:
            dx,dy=x2-x1,y2-y1; d=(dx*dx+dy*dy)**.5
            ax.plot([x1-dy/d*.23,x2-dy/d*.23],[y1+dx/d*.23,y2+dx/d*.23],lw=1.4,color=dark)
    for atom in m.GetAtoms():
        if atom.GetSymbol()=='O':
            x,y=p[atom.GetIdx()]; t(x,y,'OH' if atom.GetTotalNumHs() else 'O',11,bbox={'facecolor':'white','edgecolor':'none','pad':.3})
    x,y=p[branch]; ax.scatter([x],[y],s=35,c=c,zorder=5); t(x+(0 if 'CCC(CC)' in smiles else -2),y+(2.3 if 'CCC(CC)' in smiles else 1.3),'C4',9,c,'bold')
    return rdMolDescriptors.CalcMolFormula(m)

t(50,97,'Proposed host biosynthetic pathways of mutton odor-associated branched-chain fatty acids',17,dark,'bold')
t(2,91,'PRECURSOR NETWORK',10,muted,'bold',ha='left')
t(12,85,'Rumen fermentation',10,dark,'bold'); t(33,89,'Propionate',10); t(33,81,'Butyrate',10)
line(20,86.7,27,89,muted); line(20,83.3,27,81,muted)
t(53,89,'Propionyl-CoA',10); t(53,81,'Butyryl-CoA',10)
line(38,89,46,89,muted,':'); line(38,81,46,81,muted,':')
t(78,89,'Methylmalonyl-CoA',10,purple,'bold'); t(78,81,'Ethylmalonyl-CoA',10,orange,'bold')
line(59,89,68,89,muted,'--'); line(59,81,68,81,muted,'--')
t(63.5,91.5,'ACC*',8,muted); t(63.5,83.5,'ACC*',8,muted)
t(91,89,'Acetyl-CoA',9); t(91,81,'Malonyl-CoA',9,'#31796e','bold'); line(91,87.5,91,83,muted); t(95,85,'ACC',8,muted)
t(50,76,'* Proposed cytosolic side reactions; dotted arrows denote uncertain supply / activation',8,muted)
line(2,73.5,98,73.5,'#d4dcdf',arrow=False)
t(2,71,'CANDIDATE FASN ASSEMBLY',10,muted,'bold',ha='left')
for x,s in [(16,'PRIMER'),(35,'FIRST EXTENSION'),(57,'BRANCH INTRODUCTION'),(82,'FINAL EXTENSION / RELEASE')]: t(x,68,s,9,muted,'bold')
rows=[(60,blue,'MOA','Acetyl-CoA','C2','Butyryl-ACP','C4','2-Methylhexanoyl-ACP','C7','4-Methyloctanoyl-ACP','C9','methylmalonyl-CoA'),(49,purple,'MNA','Propionyl-CoA','C3','Pentanoyl-ACP','C5','2-Methylheptanoyl-ACP','C8','4-Methylnonanoyl-ACP','C10','methylmalonyl-CoA'),(38,orange,'EOA','Acetyl-CoA','C2','Butyryl-ACP','C4','2-Ethylhexanoyl-ACP','C8','4-Ethyloctanoyl-ACP','C10','ethylmalonyl-CoA')]
for y,c,name,p,pc,f,fc,b,bc,end,ec,donor in rows:
    t(3,y,name,13,c,'bold',ha='left'); t(15,y,p,10,dark,'bold'); t(15,y-2.4,pc,9,c,'bold')
    t(35,y,f,10); t(35,y-2.4,fc,9,c,'bold'); t(59,y,b,10,dark,'bold'); t(59,y-2.4,bc+'  ·  C2 branch',9,c,'bold')
    t(84,y,end,10); t(84,y-2.4,ec+'  ·  C4 branch',9,c,'bold')
    for a,z in [(21,28),(41,50),(68,75)]: line(a,y,z,y,c,'--')
    t(25,y+2.5,'+ malonyl',8,muted); t(45.5,y+2.5,'+ '+donor,8,c,'bold'); t(71.5,y+2.5,'+ malonyl',8,muted)
    line(91,y,94.5,y,c,'--'); t(96.3,y,name,9,c,'bold'); t(93,y-3.2,'TE*',8,c)
    if y!=38: line(2,y-5.5,98,y-5.5,'#e5eaec',lw=.6,arrow=False)
t(50,31,'C2 branch  →  +2 carbons at the carboxyl end  →  C4 branch',11,dark,'bold')
line(2,28.5,98,28.5,'#d4dcdf',arrow=False)
structures=[('O=C(O)CCC(C)CCCC',18,blue,'MOA  ·  C9','C9H18O2'),('O=C(O)CCC(C)CCCCC',50,purple,'MNA  ·  C10','C10H20O2'),('O=C(O)CCC(CC)CCCC',82,orange,'EOA  ·  C10','C10H20O2')]
for sm,x,c,label,formula in structures:
    t(x,27,label,12,c,'bold'); assert structure(sm,x,20,1.55,c)==formula
    t(x,13.7,formula+'  ·  stereochemistry unspecified',8,muted)
t(2,10.5,'MITOCHONDRIAL CATABOLISM',9,muted,'bold',ha='left')
t(2,8.3,'Propionyl-CoA  —PCC→  (2S)-methylmalonyl-CoA  —MCEE→  (2R)-methylmalonyl-CoA  —MMUT / B12→  succinyl-CoA',8,dark,ha='left')
t(2,5.9,'CYTOSOLIC REPAIR',9,muted,'bold',ha='left')
t(2,3.8,'ECHDC1: ethylmalonyl-CoA → butyryl-CoA; methylmalonyl-CoA → propionyl-CoA (lower activity)  ·  limits branched donors',8,dark,ha='left')
t(2,1.5,'Solid: enzymatic basis   ·   Dashed: proposed assembly / side reaction   ·   Dotted: indirect supply   |   * TE identity and early release unverified',8,muted,ha='left')
fig.subplots_adjust(left=.025,right=.975,top=.975,bottom=.025)
for ext in ['svg','pdf','png']: fig.savefig(OUT/f'BCFA_mechanism_final.{ext}',dpi=300,facecolor='white')
print('formula and C4 connectivity validated')
