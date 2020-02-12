import networkx as nx
import pandas as pd


df = pd.read_table("data/PCNet_v1.3.tsv", sep = "\t")
G = nx.from_pandas_edgelist(df, 'from', 'to')

df_snp2disease = pd.read_table("data/all_variant_disease_associations.tsv.gz", sep = "\t")
df_snp2gene = pd.read_table("data/variant_to_gene_mappings.tsv.gz", sep = "\t")

df_snp2disease
df_snp2gene

len(df_snp2disease.diseaseId.value_counts())

def get_neighbors(disease_id):
    snps = set(df_snp2disease[df_snp2disease.diseaseId==disease_id].snpId)
    genes = set(df_snp2gene[df_snp2gene['snpId'].isin(snps)].geneSymbol)
    print(len(genes))
    edges = {}
    black_list = set()
    for n in genes:
        if not n in G.nodes: 
            black_list.add(n)
            continue
        for v in G.neighbors(n):
            if v in genes: continue
            if not v in edges: edges[v] = set()
            edges[v].add(n)

    return({
        'init' : genes - black_list,
        'neighbors': edges
    })



def format_neighbors_object(obj_A, obj_B):

    all_intermed_genes = set(obj_A['neighbors'].keys()) | set(obj_B['neighbors'].keys())
    edges = [ \
        {'node' : v, 
        'a' : obj_A['neighbors'][v] - obj_B['init'] if v in obj_A['neighbors'] else set() , \
        'b' : obj_B['neighbors'][v] - obj_A['init'] if v in obj_B['neighbors'] else set() } for v in all_intermed_genes ]
    return sorted(edges, key = lambda v : -len(v['a']) * len(v['b']))


A = get_neighbors("C0520679")
B = get_neighbors("C0497327")

ordered_intermed_genes = format_neighbors_object(A, B)

remained_A, remained_B = A['init'] - B['init'], B['init'] - A['init']
found = False
solution = []
solution_edges = { 'from' : [], 'to' : [] }

for v in ordered_intermed_genes:
    if len(remained_A) == 0 and len(remained_B) == 0:
        found = True
        break

    if len(remained_A & v['a']) == 0 and len(remained_B & v['b']) == 0:
        continue

    solution.append(v['node'])

    remained_A -= v['a']
    remained_B -= v['b']

    for n in v['a'] | v['b']:
        solution_edges['from'].append(v['node'])
        solution_edges['to'].append(n)

solution, found

remained_A, remained_B

pd.DataFrame(solution_edges).to_csv("results/test.csv", index=False)