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
        if type(disease_id) is list:
            snps = set(Solver.snp2disease[Solver.snp2disease.diseaseId.isin(disease_id)].snpId)
        else:
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
        remained_A, remained_B = obj_A['init'] - obj_B['init'], obj_B['init'] - obj_A['init']
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
            'edges' : solution_edges,
            'result' : found
        }

    
    @staticmethod
    def solve_iterative_greedy(obj_A, obj_B):
        deg = dict(Solver.G.degree)
        df = Solver.snp2disease.merge(Solver.snp2gene, on = "snpId")
        gene_snps = df.groupby(['geneSymbol'])['snpId'].count().to_dict()
        gene_degs = dict(Solver.G.degree)

        df = pd.read_table("./data/all_gene_disease_associations.tsv.gz")
        df = df[(df.diseaseType == "disease") | (df.diseaseType == "phenotype")]
        gene_diss = df.groupby(['geneSymbol'])['diseaseId'].count().to_dict()

        remained_A, remained_B = obj_A['init'] - obj_B['init'], obj_B['init'] - obj_A['init']

        edges = Solver.format_neighbors_object(obj_A, obj_B)
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

        results = []
        valid = False
        solution_edges = { 'from' : [], 'to' : [], 'group': [] }
        
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

            i, best_rank = None, float('inf')
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
            results.append({
                'node' : df.iloc[i].node,
                'rank_agg' : df.iloc[i].rank_agg,
                'a' : list(edges[i]['a']),
                'b' : list(edges[i]['b']),
                'n_snp' : df.iloc[i].n_snp,
                'n_dis' : df.iloc[i].n_dis,
                'n_deg' : df.iloc[i].g_deg
            })

            for j in range(E):
                if used[j]: continue
                edges[j]['a'] -= edges[i]['a']
                edges[j]['b'] -= edges[i]['b']

            for g in ['a', 'b']:
                for n in edges[i][g]:
                    solution_edges['from'].append(edges[i]['node'])
                    solution_edges['to'].append(n)
                    solution_edges['group'].append(g)

            edges[i]['a'].clear()
            edges[i]['b'].clear()
            if len(remained_A) == 0 and len(remained_B) == 0:
                valid = True
                break
            
        print(remained_A) 
        print(remained_B)
            
        return {
            'result' : results,
            'edges' : solution_edges,
            'valid' : valid
        }