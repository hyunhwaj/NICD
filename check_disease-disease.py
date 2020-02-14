import pandas as pd
from main import Solver
import random
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_table("data/disease_to_disease_ALL.tsv.gz", sep = "\t")

disease_ids = list(set(df[df.VD_variant1GeneIdSet > 100].diseaseId1.unique()) | set(df[df.VD_variant2GeneIdSet > 100].diseaseId2.unique()))

len(disease_ids)

def get_associated_genes(disease_id):
    snps = set(Solver.snp2disease[Solver.snp2disease.diseaseId==disease_id].snpId)
    genes = set(Solver.snp2gene[Solver.snp2gene['snpId'].isin(snps)].geneSymbol)
    
    return set([g for g in genes if g in Solver.G.nodes])


def calc_jaccard_similarity(A, B):
    return len(A&B) / len(A|B) 


def calc_intermediate_coverage(A, B):
    both = A & B

    NA = A - both
    NB = B - both

    yes = 0
    total = 0
    for v in both:
        if not v in Solver.G.nodes: continue
        total += 1
        neighbors = Solver.G.neighbors(v)
        for u in NA | NB:
            if u in neighbors:
                yes += 1
                break
    return (yes, total)


results = []

gene_sets = [ get_associated_genes(did) for did in disease_ids ]

for i in range(len(disease_ids)):
    disease_a = disease_ids[i]
    A = gene_sets[i]
    for j in range(i+1, len(disease_ids)):
        disease_b = disease_ids[j]
        B = gene_sets[j]
        if len(A & B) == 0: continue
        coverage = calc_intermediate_coverage(A, B)
        results.append(
            { 'a' : disease_a, 
            'b' : disease_b,
            'JS' : calc_jaccard_similarity(A, B), 
            'cover' : coverage[0] / coverage[1],
            'linked' : coverage[0],
            'total' : coverage[1]
            }
        )
    
df = pd.DataFrame(results)

plt.scatter(x=df.JS, y=df.cover, alpha=0.1)

plt.hist(df.cover)

plt.hist2d(x=df.JS, y=df.cover, normed=False, cmap='plasma')
cb = plt.colorbar()
plt.show()

plt.scatter(x=df.total, y=df.linked, alpha=0.1)

df2 = df[df.total >=10]

plt.hist(df2.cover, bins=50)
