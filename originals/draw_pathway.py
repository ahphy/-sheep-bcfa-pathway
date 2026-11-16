from pathlib import Path
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from rdkit import Chem
from rdkit.Chem import rdDepictor, rdMolDescriptors

OUT = Path(__file__).resolve().parent
mpl.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','pdf.fonttype':42})
fig, ax = plt.subplots(figsize=(20,14)); ax.set(xlim=(0,100),ylim=(0,100)); ax.axis('off')
colors=['#2878b5','#8054a5','#dc8430']; grey='#40505a'
def text(x,y,s,size=10,color=grey,weight='normal',**kw):
    ax.text(x,y,s,ha='center',va='center',fontsize=size,color=color,weight=weight,**kw)
def box(x,y,w,h,s,edge='#bac5cb',fill='white',size=10,color=grey,weight='normal'):
    ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle='round,pad=0.2,rounding_size=0.9',ec=edge,fc=fill,lw=1.1))
    text(x,y,s,size,color,weight)
def arrow(x1,y1,x2,y2,color=grey,style='solid',lw=1.3):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=10,lw=lw,color=color,linestyle=style))
def molecule(smiles,cx,cy,scale=1.0):
    m=Chem.MolFromSmiles(smiles); rdDepictor.Compute2DCoords(m)
    p=m.GetConformer(); coords=[(p.GetAtomPosition(i).x,p.GetAtomPosition(i).y) for i in range(m.GetNumAtoms())]
    minx,maxx=min(a for a,b in coords),max(a for a,b in coords); miny,maxy=min(b for a,b in coords),max(b for a,b in coords)
    coords=[(cx+(a-(minx+maxx)/2)*scale,cy+(b-(miny+maxy)/2)*scale) for a,b in coords]
    for bond in m.GetBonds():
        i,j=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx(); x1,y1=coords[i]; x2,y2=coords[j]
        ax.plot([x1,x2],[y1,y2],color='#29363c',lw=1.5)
        if bond.GetBondType()==Chem.BondType.DOUBLE:
            dx,dy=x2-x1,y2-y1; d=(dx*dx+dy*dy)**.5
            ax.plot([x1-dy/d*.25,x2-dy/d*.25],[y1+dx/d*.25,y2+dx/d*.25],color='#29363c',lw=1)
    for a in m.GetAtoms():
        if a.GetSymbol()=='O':
            x,y=coords[a.GetIdx()]; text(x,y,'OH' if a.GetTotalNumHs() else 'O',10,'#29363c',bbox=dict(facecolor='white',edgecolor='none',pad=.4))
    return rdMolDescriptors.CalcMolFormula(m),m.GetNumAtoms()

text(50,98,'Candidate host pathways for 4-branched fatty acids',19,'#22323b','bold')
box(50,92.5,88,6,'RUMEN FERMENTATION  ·  microbial production',edge='#8fbab5',fill='#f4faf9',size=12,weight='bold')
box(30,86,19,3.8,'Propionate',edge='#8fbab5',fill='#f5fbfa',size=11)
box(70,86,19,3.8,'Butyrate',edge='#8fbab5',fill='#f5fbfa',size=11)
arrow(30,89.5,30,88); arrow(70,89.5,70,88)
text(50,82.5,'Absorption and host activation  ·  enzyme / compartment details require validation',9,color='#687980')
arrow(30,84,30,79.5,'#87999c',':'); arrow(70,84,70,79.5,'#87999c',':')
box(19,77,25,4,'Propionyl-CoA (C3)',size=10)
box(50,77,25,4,'Acetyl-CoA (C2)',size=10)
box(81,77,25,4,'Butyryl-CoA (C4)',size=10)
text(50,73,'CYTOSOL  ·  candidate FASN assembly',12,'#233b45','bold')
text(50,69.5,'Shared extension donor: Malonyl-CoA  ·  acetyl-CoA + CO₂ — ACC → malonyl-CoA',9,'#35786f')
xs=[17,50,83]
titles=['MOA  ·  methyl branch','MNA  ·  methyl branch','EOA  ·  ethyl branch']
primers=['Acetyl-CoA  |  primer C2','Propionyl-CoA  |  primer C3','Acetyl-CoA  |  primer C2']
first=['Butyryl-ACP  |  chain C4','Pentanoyl-ACP  |  chain C5','Butyryl-ACP  |  chain C4']
second=['2-Methylhexanoyl-ACP  |  C7','2-Methylheptanoyl-ACP  |  C8','2-Ethylhexanoyl-ACP  |  C8']
third=['4-Methyloctanoyl-ACP  |  C9','4-Methylnonanoyl-ACP  |  C10','4-Ethyloctanoyl-ACP  |  C10']
prod=['4-Methyloctanoic acid  ·  C9','4-Methylnonanoic acid  ·  C10','4-Ethyloctanoic acid  ·  C10']
smiles=['O=C(O)CCC(C)CCCC','O=C(O)CCC(C)CCCCC','O=C(O)CCC(CC)CCCC']
for i,x in enumerate(xs):
    c=colors[i]; text(x,66.5,titles[i],12,c,'bold')
    for y,s in [(63,primers[i]),(56,first[i]),(48,second[i]),(40,third[i])]: box(x,y,29,3.6,s,c,'white',9.2)
    for y1,y2 in [(61,58),(54,50),(46,42),(38,35)]: arrow(x,y1,x,y2,c,'--')
    labelbg=dict(facecolor='white',edgecolor='none',pad=1)
    text(x,59.5,'+ malonyl-CoA  ·  FASN',8,c,bbox=labelbg)
    text(x,52,('+ methylmalonyl-CoA  ·  +3 C retained' if i<2 else '+ ethylmalonyl-CoA  ·  +4 C retained'),8,c,bbox=labelbg)
    text(x,44,'+ malonyl-CoA  ·  branch C2 → C4',8,c,bbox=labelbg)
    text(x,36.5,'TE  ·  proposed early release',8,c,bbox=labelbg)
    ax.add_patch(FancyBboxPatch((x-14.5,22),29,12,boxstyle='round,pad=.25,rounding_size=1',ec=c,fc='#ffffff',lw=1.6))
    text(x,32,prod[i],10,c,'bold'); formula,n=molecule(smiles[i],x,27.5,.85)
    assert formula==('C9H18O2' if i==0 else 'C10H20O2')
    text(x,23.3,formula+'  ·  C4 branch  ·  stereochemistry unspecified',8,c)
text(50,20.7,'ACC side reactions (candidate cytosolic supply): propionyl-CoA → methylmalonyl-CoA; butyryl-CoA → ethylmalonyl-CoA',9,'#536d6c')
ax.add_patch(FancyBboxPatch((2,2),46,17,boxstyle='round,pad=.3,rounding_size=1.2',ec='#aab5ba',fc='#f4f6f7',lw=1))
text(25,17.4,'MITOCHONDRION  ·  propionate catabolism',10,'#34444c','bold')
text(25,14.7,'Propionyl-CoA  →  (2S)-Methylmalonyl-CoA',9)
text(25,12.5,'PCC  ·  biotin',8,'#667980')
text(25,10.3,'→  (2R)-Methylmalonyl-CoA  →  Succinyl-CoA  →  TCA',9)
text(25,8.1,'MCEE                               MMUT  ·  vitamin B12',8,'#667980')
text(25,4.4,'Separate compartment; no direct transfer to FASN implied',8,'#667980')
ax.add_patch(FancyBboxPatch((51,2),47,17,boxstyle='round,pad=.3,rounding_size=1.2',ec='#aab5ba',fc='white',lw=1))
text(74.5,17.4,'CYTOSOL  ·  metabolite repair',10,'#34444c','bold')
text(74.5,14.8,'Ethylmalonyl-CoA  →  Butyryl-CoA',9)
text(74.5,12.8,'Methylmalonyl-CoA  →  Propionyl-CoA  (lower activity)',9)
text(74.5,10.8,'ECHDC1  ·  limits branched donor availability  ·  distinct from MLYCD',8)
ax.plot([55,59],[7.8,7.8],color=grey,lw=1.3); text(65,7.8,'enzyme reaction',7.5)
ax.plot([72,76],[7.8,7.8],color=grey,lw=1.3,ls='--'); text(84,7.8,'proposed assembly',7.5)
ax.plot([55,59],[5,5],color=grey,lw=1.3,ls=':'); text(70,5,'indirect / uncertain supply',7.5)
text(89,5,'Evidence ≠ sheep flux',7.5)
fig.subplots_adjust(left=.02,right=.98,top=.98,bottom=.02)
for ext in ['svg','pdf','png']: fig.savefig(OUT/('BCFA_biosynthesis_pathway.'+ext),dpi=300,facecolor='white')
print('Exported; formula checks passed')
