deg = dict(Solver.G.degree)

df = Solver.snp2disease.merge(Solver.snp2gene, on = "snpId")

gene_snps = df.groupby(['geneSymbol'])['snpId'].count().to_dict()
gene_degs = dict(Solver.G.degree)


import pandas as pd

df = pd.read_table("./data/all_gene_disease_associations.tsv.gz")
gene_dis = df.groupby(['geneSymbol'])['diseaseId']
A = Solver.get_neighbors("C0520679")
B = Solver.get_neighbors("C0497327")

remained_A, remained_B = A['init'] - B['init'], B['init'] - A['init']

edges = Solver.format_neighbors_object(A, B)
E = len(edges)

used = [False for i in range(E)]

from math import sqrt, log
   
    
for _ in range(E):
    i, best_score = None, 0
    
    for j, v in enumerate(edges):
        if used[j]: continue
        penalty = 0.5 if not v['node'] in gene_snps else gene_snps[v['node']]
        con_score = sqrt(len(v['a']) * len(v['b']))
        if con_score == 0:
            if len(v['a']) > 0 or len(v['b']) > 0:
                con_score = 0.5
        score =  con_score / penalty
        if best_score < score:
            best_score = score
            i = j
   
    if i is None:
        print("Not found")
        break
    
    penalty = None if not edges[i]['node'] in gene_snps else gene_snps[edges[i]['node']]
    print(edges[i]['node'], best_score, len(edges[i]['a']), len(edges[i]['b']), penalty)

    used[i] = True
    remained_A -= edges[i]['a']
    remained_B -= edges[i]['b']    
    
    for j in range(E):
        if used[j]: continue
        edges[j]['a'] -= edges[i]['a']
        edges[j]['b'] -= edges[i]['b']            
        
    if len(remained_A) == 0 and len(remained_B) == 0:
        print("Found")
        break