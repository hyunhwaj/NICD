import networkx as nx
import pandas as pd
import logging

class Solver:
    df = pd.read_table("data/PCNet_v1.3.tsv", sep = "\t")
    G = nx.from_pandas_edgelist(df, 'from', 'to')
    G.remove_node("UBC")
    snp2disease = pd.read_table("data/all_variant_disease_associations.tsv.gz", sep = "\t")
    snp2gene = pd.read_table("data/variant_to_gene_mappings.tsv.gz", sep = "\t")


    @staticmethod
    def map_snp_to_genes(disease_id):
        snps = set(Solver.snp2disease[Solver.snp2disease.diseaseId==disease_id].snpId)
        genes = set(Solver.snp2gene[Solver.snp2gene['snpId'].isin(snps)].geneSymbol)
        return genes 
    
    
    @staticmethod
    def get_neighbors(disease_id):
        genes = Solver.map_snp_to_genes(disease_id)
        edges = {}
        black_list = set()
        for n in genes:
            if not n in Solver.G.nodes: 
                black_list.add(n)
                continue
            for v in Solver.G.neighbors(n):
                if v in genes: continue
                # if Solver.G.degree(v) >= 100: continue
                if not v in edges: edges[v] = set()
                edges[v].add(n)

        return({
            'init' : genes - black_list,
            'neighbors': edges
        })

    @staticmethod
    def format_neighbors_object(obj_A, obj_B):

        all_intermed_genes = set(obj_A['neighbors'].keys()) | set(obj_B['neighbors'].keys())
        edges = [ \
            {'node' : v, 
            'a' : obj_A['neighbors'][v] - obj_B['init'] if v in obj_A['neighbors'] else set() , \
            'b' : obj_B['neighbors'][v] - obj_A['init'] if v in obj_B['neighbors'] else set() } for v in all_intermed_genes ]
        return sorted(edges, key = lambda v : -len(v['a']) * len(v['b']))

    @staticmethod
    def solve_sequencial(obj_A, obj_B):
        ordered_intermed_genes = Solver.format_neighbors_object(obj_A, obj_B)
        remained_A, remained_B = A['init'] - B['init'], B['init'] - A['init']
        found = False
        solution = []
        solution_edges = { 'from' : [], 'to' : [], 'group': [] }

        for v in ordered_intermed_genes:
            if len(remained_A) == 0 and len(remained_B) == 0:
                found = True
                break

            if len(remained_A & v['a']) == 0 and len(remained_B & v['b']) == 0:
                continue

            solution.append(v['node'])

            remained_A -= v['a']
            remained_B -= v['b']

            for g in ['a', 'b']:
                for n in v[g]:
                    solution_edges['from'].append(v['node'])
                    solution_edges['to'].append(n)
                    solution_edges['group'].append(g)

        return {
            'node' : solution,
            'edges' : pd.DataFrame(solution_edges),
            'result' : found
        }


A = Solver.get_neighbors("C0520679")
B = Solver.get_neighbors("C0497327")
print(Solver.solve_sequencial(A, B))



