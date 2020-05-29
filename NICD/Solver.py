import networkx as nx
import pandas as pd
import logging
import os

class Solver:
    data_path = os.path.join(os.path.dirname(__file__), "data")
    df = pd.read_table(f"{data_path}/PCNet_v1.3.tsv", sep = "\t")
    
    G = nx.from_pandas_edgelist(df, 'from', 'to')
    G.remove_node("UBC")
    snp2disease = pd.read_table(f"{data_path}/all_variant_disease_associations.tsv.gz", sep = "\t")
    snp2gene = pd.read_table(f"{data_path}/variant_to_gene_mappings.tsv.gz", sep = "\t")
    
    df = pd.read_table(f"{data_path}/all_gene_disease_associations.tsv.gz")
    df = df[(df.diseaseType == "disease") | (df.diseaseType == "phenotype")]
    disease_ids = list(set(df.diseaseId))

    @staticmethod
    def map_snp_to_genes(disease_id):
        if type(disease_id) is list:
            snps = set(Solver.snp2disease[Solver.snp2disease.diseaseId.isin(disease_id)].snpId)
        else:
            snps = set(Solver.snp2disease[Solver.snp2disease.diseaseId==disease_id].snpId)
        genes = set(Solver.snp2gene[Solver.snp2gene['snpId'].isin(snps)].geneSymbol)
        return genes 
   
    @staticmethod
    def get_neighbor_genes(disease_id):
        init_genes = Solver.map_snp_to_genes(disease_id)
        return set().union(*[ Solver.G.neighbors(n) for n in init_genes if n in Solver.G.nodes ])
            
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
    def list_candidates(obj_A, obj_B):
        deg = dict(Solver.G.degree)
        df = Solver.snp2disease.merge(Solver.snp2gene, on = "snpId")
        gene_snps = df.groupby(['geneSymbol'])['snpId'].count().to_dict()
        gene_degs = dict(Solver.G.degree)

        data_path = os.path.join(os.path.dirname(__file__), "data")
        df = pd.read_table(f"{data_path}/all_gene_disease_associations.tsv.gz")
        df = df[(df.diseaseType == "disease") | (df.diseaseType == "phenotype")]
        gene_diss = df.groupby(['geneSymbol'])['diseaseId'].count().to_dict()

        remained_A, remained_B = obj_A['init'] - obj_B['init'], obj_B['init'] - obj_A['init']

        edges = Solver.format_neighbors_object(obj_A, obj_B)
        E = len(edges)

        used = [False for i in range(E)]

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

        cur_status = []
        for j, v in enumerate(edges):
            stat_v = get_stats(v)
            stat_v.update(v)
            cur_status.append(stat_v)

        df = pd.DataFrame(cur_status)
        
        return df
