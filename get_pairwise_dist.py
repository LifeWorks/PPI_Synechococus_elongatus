import os, sys, re, csv
from pathlib import Path
home = str(Path.home())
import matplotlib.pyplot as plt
plt.style.use('default')
plt.rcParams['figure.facecolor'] = 'white'
import numpy as np
import scipy as sp
import pandas as pd
from scipy import optimize as opt
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from scipy.stats import pearsonr
import scipy.odr as odr
import networkx as nx
import seaborn as sns
import utils as ut


result_dir = home +"/workspace/data/PredPheno/cellbox/results/20230405_cyano_rna_5c08105256c016c79573735ab13502a2/"
analysis_dir = home + "/Dropbox/PNNL/PredPheno/SystemModeling/Modeling/S_elongatus/cellbox/RNA/20230408_analysis/"

result_folders = [f for f in os.listdir(result_dir) if os.path.isdir(os.path.join(result_dir, f))]
result_scores = []
for folder in result_folders:
    # print(folder)
    seed = int(re.findall(r'\d+\.\d+|\d+', folder)[0])
    w_files = [listed for listed in os.listdir(result_dir + folder) if "best.W" in listed]
    # print(w_files)
    result_scores += [[seed] + [float(num) if '.' in num else int(num) for num in re.findall(r'\d+\.\d+|\d+', file)] for file in w_files]
result_scores = pd.DataFrame(result_scores, columns=["seed", "rep", "score"])

scores_1000 = result_scores.sort_values(by="score", ascending=True).head(1000).copy()

scores = scores_1000.copy()
pairwise_dist = []
for i in range(scores.shape[0]):
    for j in range(i+1, scores.shape[0]):
        a_seed = '{:0>{}}'.format(str(scores.iloc[i, 0]), 3)
        b_seed = '{:0>{}}'.format(str(scores.iloc[j, 0]), 3)
        a_rep = str(scores.iloc[i, 1])
        b_rep = str(scores.iloc[j, 1])
        a_file = [listed for listed in os.listdir(result_dir + "seed_" + a_seed) if a_rep + "_best.W" in listed][0]
        b_file = [listed for listed in os.listdir(result_dir + "seed_" + b_seed) if b_rep + "_best.W" in listed][0]
        a_df = pd.read_csv(result_dir + "seed_" + a_seed + '/' + a_file, index_col=0)
        b_df = pd.read_csv(result_dir + "seed_" + b_seed + '/' + b_file, index_col=0)
        dist = [np.linalg.norm(a_df.to_numpy() - b_df.to_numpy())]
        if a_seed == b_seed:
            dist += ["Same run"]
        else:
            dist += ["Different run"]
        pairwise_dist += [dist]
pairwise_dist = pd.DataFrame(pairwise_dist, columns=["dist", "run"])
pairwise_dist.to_pickle(analysis_dir + "pairwise_distances_top1k.pkl")
