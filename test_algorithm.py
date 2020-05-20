from Solver import Solver

dementia_ids = ["C0524851", "C0002736", "C0030567",
                "C0002395", "C0020179", "C0497327", 
                "C0338656", "C0011265", "C0038454",
                "C0242422", "C0233794", "C0026769",
                "C0338451", "C0752347"]

osa_ids = [ "C0520679", "C0037315", "C0520680"]


A = Solver.get_neighbors(dementia_ids)
B = Solver.get_neighbors(osa_ids)

ret_full = Solver.solve_iterative_greedy(A, B)

import random
random.seed(123)

def sample_neighbors(obj, rate = 0.75):
    obj['neighbors'] = { k : v for k, v in obj['neighbors'].items() if random.random() < rate }
    return obj

sampled_A = []
sampled_B = []
ret_sampled = [] 
for i in range(100):    
    sampled_A.append(sample_neighbors(A))
    sampled_B.append(sample_neighbors(B))
    ret_sampled.append(Solver.solve_iterative_greedy(sampled_A[-1], sampled_B[-1]))

genes_with_full = [ x['node'] for x in ret_full['result'] ]

from collections import defaultdict
genes_with_samples = defaultdict(int)

for sample in ret_sampled:
    for x in [ i['node'] for i in sample['result'] ]:
        genes_with_samples[x] += 1

for i, v in enumerate(genes_with_full):
    if genes_with_samples[v] > 0:
        print(v, ret_full['result'][i]['rank_agg'], genes_with_samples[v])
        
        
