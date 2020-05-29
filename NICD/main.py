import sys, logging, random, os
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
from statsmodels.stats.multitest import multipletests        
import argparse

def main(disease_a, disease_b, 
         name_a="disease A", 
         name_b="disease_B",
         niter=100, cutoff=0.05, 
         random_seed=None, 
         outpath='results',
         verbose=False):
       
    outpath = os.path.abspath(outpath)
    FORMAT = '%(asctime)-15s %(name)s %(levelname)s %(message)s'
    
    if verbose:
        logging.basicConfig(level=logging.DEBUG, format=FORMAT)
    else:
        logging.basicConfig(level=logging.INFO, format=FORMAT)
            
    logging.info("Loading Solver module.")
    import pandas as pd
    from NICD.Solver import Solver

    logging.info("Starting NICD.")

    disease_a = disease_a.split(",")
    disease_b = disease_b.split(",")
    
    A = Solver.get_neighbors(disease_a)
    B = Solver.get_neighbors(disease_b)
    logging.debug(f"disease_a: {disease_a}")
    logging.debug(f"disease_b: {disease_b}")
    
    ret = Solver.list_candidates(A, B)
    logging.debug(
        f"{sum(ret.con_a>0)} and {sum(ret.con_a>0)} genes were found with both diseases."
    )
    
    ret = ret[(ret.con_a>0)&(ret.con_b>0)]

    def time_now():
        from datetime import datetime
        return datetime.now().strftime("%H:%M:%S")

    if random_seed is None:
        random_seed = random.randrange(sys.maxsize)
    random.seed(random_seed)
    
    logging.info(f"{random_seed} is the random seed.")
        
    from collections import defaultdict
    null_cnt = defaultdict(int)

    NUM_ITER = max(int(len(ret.node) * niter),1)

    logging.info(f"Starting a permutation Test. {NUM_ITER} iterations will be performed.")

    null_common = []
    for i in tqdm(range(1, NUM_ITER+1)):
        null_ids = random.choices(Solver.disease_ids, k=len(disease_a) + len(disease_b))
        null_a = Solver.get_neighbor_genes(null_ids[:len(disease_a)])
        null_b = Solver.get_neighbor_genes(null_ids[len(disease_a):])
        null_common.append(null_a & null_b)

    for i in null_common:
        for v in i:
            null_cnt[v] += 1

    logging.info("The end of the permutation Test.")

    pval = [null_cnt[v] / NUM_ITER for v in ret.node]


    logging.info("Calculating adjusted p-value of each candidate gene.")


    fdr = multipletests(pval, method = "fdr_bh")
    plt.hist(fdr[1])

    ret["pvalue"] = pval
    ret["padjust"] = fdr[1]

    if not os.path.exists(outpath):
        logging.info(f"NICD creates {outpath} to store outputs.")
        os.makedirs(outpath)

    logging.info(f"Generating outputs to {outpath}.")
    ret.to_csv(f"{outpath}/permutation-test.csv", index=False)

    logging.info("NICD found {} putative driver genes.".format(ret[ret.padjust<cutoff].shape[0]))


    valid_ret = ret[ret.padjust<cutoff]
    
    with open(f"{outpath}/edges.csv", "w") as oup:
        print("from", "to", sep=",", file=oup)
        for i in range(valid_ret.shape[0]):
            fr = valid_ret.iloc[i]['node']
            for to in valid_ret.iloc[i]['a']:
                print(fr, to, sep=",", file=oup)
            for to in valid_ret.iloc[i]['b']:
                print(fr, to, sep=",", file=oup)


    with open(f"{outpath}/nodes.csv", "w") as oup:
        print("gene", "category", sep=",", file=oup)
        for x in valid_ret.node:
            print(x, "common", sep=",", file=oup)
        for a in set.union(*ret[ret.padjust<cutoff].a):
            print(a, name_a, sep=",", file=oup)
        for b in set.union(*ret[ret.padjust<cutoff].b):
            print(b, name_b, sep=",", file=oup)

    logging.info("NICD has successfuly completed all tasks.")
