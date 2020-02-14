from main import Solver
deg = dict(Solver.G.degree)

df = Solver.snp2disease.merge(Solver.snp2gene, on = "snpId")

gene_snps = df.groupby(['geneSymbol'])['snpId'].count().to_dict()
gene_degs = dict(Solver.G.degree)

import pandas as pd

df = pd.read_table("./data/all_gene_disease_associations.tsv.gz")
df = df[(df.diseaseType == "disease") | (df.diseaseType == "phenotype")]
gene_diss = df.groupby(['geneSymbol'])['diseaseId'].count().to_dict()

A = Solver.get_neighbors("C0520679")
B = Solver.get_neighbors("C0497327")

remained_A, remained_B = A['init'] - B['init'], B['init'] - A['init']

edges = Solver.format_neighbors_object(A, B)
E = len(edges)

used = [False for i in range(E)]

from math import sqrt, log   
    
def get_stats(v):
    n_snp = 0 if not v['node'] in gene_snps else gene_snps[v['node']]
    n_dis = 0 if not v['node'] in gene_diss else gene_diss[v['node']]
    con_a = len(v['a'])
    con_b = len(v['b'])
    g_deg = gene_degs[v['node']]

    return {
        'n_snp' : n_snp,
        'n_dis' : n_dis,
        'con_a' : con_a,
        'con_b' : con_b,
        'g_deg' : g_deg
    }

for _ in range(E):
    
    cur_status = []
    for j, v in enumerate(edges):
        stat_v = get_stats(v)
        stat_v.update(v)
        cur_status.append(stat_v)
        
    df = pd.DataFrame(cur_status)
    df['rank_con_a'] = df['con_a'].rank(ascending=False) 
    df['rank_con_b'] = df['con_b'].rank(ascending=False) 
    df['rank_n_snp'] = df['n_snp'].rank(ascending=False)  
    df['rank_n_dis'] = df['n_dis'].rank(ascending=False)  
    df['rank_g_deg'] = df['g_deg'].rank(ascending=False)  
    df['rank_agg'] = ((df['rank_con_a'] * df['rank_con_b']) ** (1/2) / (df['rank_n_snp'] * df['rank_n_dis'] * df['rank_g_deg']) ** (1/3))
    
    cur_status = df.to_dict('records')

    if _ == 0:
        df.sort_values('rank_agg').to_csv("results/init_score.csv")
        
    i, best_rank = None, 1
    for j, v in enumerate(cur_status):
        if len(v['a']) == 0 and len(v['b']) == 0: continue
        if not used[j] and best_rank > v['rank_agg']:
            best_rank = v['rank_agg']
            i = j

    if i is None:
        print("Not found")
        break
    
    used[i] = True
    remained_A -= edges[i]['a']
    remained_B -= edges[i]['b']    
   
    print(df.iloc[i].node, df.iloc[i].rank_agg, edges[i]['a'], edges[i]['b']) 
    for j in range(E):
        if used[j]: continue
        edges[j]['a'] -= edges[i]['a']
        edges[j]['b'] -= edges[i]['b']            
        
    edges[i]['a'].clear()
    edges[i]['b'].clear()
    if len(remained_A) == 0 and len(remained_B) == 0:
        print("Found")
        break
    