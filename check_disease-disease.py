import pandas as pd
df = pd.read_table("data/disease_to_disease_ALL.tsv.gz", sep = "\t")

df[(df.GD_commonGeneIdSet > 100) & (df.GD_jaccard < 0.9)].sort_values(by = "GD_jaccard", ascending=False).head()

disease_a, disease_b = "C1306460", "C0242379"

from main import Solver

def get_something(disease_id):
    snps = set(Solver.snp2disease[Solver.snp2disease.diseaseId==disease_id].snpId)
    genes = set(Solver.snp2gene[Solver.snp2gene['snpId'].isin(snps)].geneSymbol)
    return genes

A = get_something(disease_a) 
B = get_something(disease_b) 
both = A & B

NA = A - both
NB = B - both

yes = 0
total = 0
for v in both:
    if not v in Solver.G.nodes: continue
    total += 1
    found = False
    neighbors = Solver.G.neighbors(v)
    for u in NA | NB:
        if u in neighbors:
            yes += 1
            break
    