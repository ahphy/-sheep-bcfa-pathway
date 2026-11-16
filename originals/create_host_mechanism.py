"""Build a host-transport version of the editable shared-pathway figure."""
from pathlib import Path

here = Path(__file__).resolve().parent
source = (here / 'BCFA_shared_mechanism_source.py').read_text(encoding='utf-8')
source = source.replace("fig,ax=plt.subplots(figsize=(20,12.5)); ax.set(xlim=(0,100),ylim=(0,100)); ax.axis('off')", "fig,ax=plt.subplots(figsize=(20,15)); ax.set(xlim=(0,100),ylim=(0,120)); ax.axis('off')")
start = source.index("txt(50,97,")
end = source.index("rule(72.8)", start)
top = '''txt(50,117,'Rumen–host candidate mechanism for odor-associated branched fatty acids',19,ink,'bold')
txt(2,112,'RUMEN FERMENTATION',11,muted,'bold',ha='left')
txt(50,111,'Dietary carbohydrates  →  microbial fermentation',13,ink,'bold')
for x,s,c in [(19,'Acetate  ·  C2',blue),(50,'Propionate  ·  C3',purple),(81,'Butyrate  ·  C4',orange)]:
    txt(x,106,s,13,c,'bold'); arrow((x,109.1),(x,107.3),muted)
txt(2,102,'RUMEN EPITHELIUM',11,muted,'bold',ha='left')
txt(19,99.5,'absorption',11,blue); txt(50,99.5,'absorption → portal blood',11,purple)
txt(81,99.5,'absorption + epithelial metabolism',10,orange)
for x,c in [(19,blue),(50,purple),(81,orange)]: arrow((x,104.6),(x,101),c)
txt(19,95.8,'peripheral circulation',11,blue)
txt(50,95.8,'liver uptake / gluconeogenesis',11,purple)
txt(81,95.8,'β-hydroxybutyrate',11,orange)
for x,c in [(19,blue),(50,purple),(81,orange)]: arrow((x,98),(x,97),c)
txt(50,91.8,'ADIPOSE CYTOSOL  ·  candidate precursor supply',12,ink,'bold')
txt(16,87.8,'Acetate → acetyl-CoA',11,blue)
txt(50,87.8,'Propionate ⋯→ propionyl-CoA',11,purple)
txt(83,87.8,'Butyrate ⋯→ butyryl-CoA',11,orange)
arrow((19,94.2),(16,89.5),blue,':')
txt(16,82.5,'acetyl-CoA → malonyl-CoA',11,green,'bold')
txt(50,82.5,'⇢ methylmalonyl-CoA',11,purple,'bold')
txt(83,82.5,'⇢ ethylmalonyl-CoA',11,orange,'bold')
txt(16,79.8,'ACC',10,muted); txt(50,79.8,'ACC side reaction*',10,muted); txt(83,79.8,'ACC side reaction*',10,muted)
txt(50,75.8,'Transport to adipose and branched precursor supply remain to be quantified',10,muted)
'''
source = source[:start] + top + source[end:]
source = source.replace("BCFA_shared_mechanism.{ext}", "BCFA_rumen_host_mechanism.{ext}")
(here / 'BCFA_rumen_host_mechanism_source.py').write_text(source, encoding='utf-8')
exec(compile(source, str(here / 'BCFA_rumen_host_mechanism_source.py'), 'exec'), {'__file__':str(here / 'BCFA_rumen_host_mechanism_source.py'), '__name__':'__main__'})
