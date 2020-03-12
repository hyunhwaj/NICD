import pandas as pd

from Solver import Solver

dementia_ids = ["C0524851", "C0002736", "C0030567",
                "C0002395", "C0020179", "C0497327", 
                "C0338656", "C0011265", "C0038454",
                "C0242422", "C0233794", "C0026769",
                "C0338451", "C0752347"]

osa_ids = [ "C0520679", "C0037315", "C0520680"]


A = Solver.get_neighbors(dementia_ids)
B = Solver.get_neighbors(osa_ids)

Solver.get_neighbor_genes(osa_ids)

ret = Solver.list_candidates(A, B)
ret = ret[(ret.con_a>0)&(ret.con_b>0)]

def time_now():
    from datetime import datetime
    return datetime.now().strftime("%H:%M:%S")
    
import random 
random.seed(123)
from collections import defaultdict
null_cnt = defaultdict(int)

NUM_ITER = len(ret.node) * 100

print(time_now(), "start")

null_common = []
for i in range(1, NUM_ITER+1):
    if i % 1000 == 0:
        print(time_now(), f"{i}th iteration...")
    null_ids = random.choices(Solver.disease_ids, k=len(dementia_ids) + len(osa_ids))
    null_a = Solver.get_neighbor_genes(null_ids[:len(dementia_ids)])
    null_b = Solver.get_neighbor_genes(null_ids[len(dementia_ids):])
    null_common.append(null_a & null_b)

for i in null_common:
    for v in i:
        null_cnt[v] += 1
                
print(time_now(), "end")

for v in ret.node:
    if null_cnt[v] / NUM_ITER  < 0.01:
        print(v, null_cnt[v], null_cnt[v] / NUM_ITER)
        
pval = [ null_cnt[v] / NUM_ITER for v in ret.node ]

import matplotlib.pyplot as plt
import numpy as np


import statsmodels

fdr = statsmodels.stats.multitest.multipletests(pval, method = "fdr_bh")
plt.hist(fdr[1])

sum(fdr[1]<0.05)

ret[fdr[1]<0.05]

fdr[1][fdr[1]<0.05]

plt.hist(np.array(pval) * NUM_ITER)

ret["pvalue"] = pval
ret["fdr"] = fdr[1]

ret.to_csv("results/March11-permutation-test.csv", index=False)

ret[ret.fdr<0.05].shape

len(set.union(*ret[ret.fdr<0.05].a))
len(set.union(*ret[ret.fdr<0.05].b))

valid_ret = ret[ret.fdr<0.05]
with open("results/edges.tsv", "w") as oup:
    print("from", "to", sep="\t", file=oup)
    for i in range(valid_ret.shape[0]):
        fr = valid_ret.iloc[i]['node']
        for to in valid_ret.iloc[i]['a']:
            print(fr, to, sep="\t", file=oup)
        for to in valid_ret.iloc[i]['b']:
            print(fr, to, sep="\t", file=oup)        
        
        
with open("results/node.tsv", "w") as oup:
    print("gene", "category", sep="\t", file=oup)
    for x in valid_ret.node:
        print(x, "common", sep="\t", file=oup)
    for a in set.union(*ret[ret.fdr<0.05].a):
        print(a, "dementia", sep="\t", file=oup)
    for b in set.union(*ret[ret.fdr<0.05].b):
        print(b, "OSA", sep="\t", file=oup)        