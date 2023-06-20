import pandas as pd
import numpy as np 
from scipy.cluster import hierarchy as hi
import matplotlib.pyplot as plt 
import seaborn as sns 
import glob 
import os 

cwd = os.getcwd()
if not os.path.exists("FIGS"):
    os.makedirs("FIGS")

os.chdir(cwd + "/InfMat")
infls = sorted(glob.glob("*_InfMat.csv"))

for infl in infls:
    INF = pd.read_csv(infl,index_col=0)
    df = pd.read_table(cwd+"/TeamsCalc/"+infl.replace("_InfMat.csv",".teams"),delimiter=":",names=["Team/Level","Nodes"])
    df.dropna(inplace=True,ignore_index=True)
    p = ",".join(df.loc[2:,"Nodes"].to_list())
    peri = p.split(",")
    INF.drop(index=peri,columns=peri,inplace=True)
    INF.sort_index(axis = 1, inplace = True)
    INF.sort_index(axis = 0, inplace = True)
    d = hi.distance.pdist(INF)
    L  = hi.linkage(d, method = "complete", optimal_ordering = True)
    clust = hi.cut_tree(L, n_clusters = 2)
    cluster = np.transpose(clust)
    t1 = INF.index[cluster[0] == 0]
    t2 = INF.index[cluster[0] == 1]
    tot = t1.append(t2)
    tot = tot.to_list()
    tot1 = tot[::-1]
    cl_inf = INF.loc[tot,:].T.loc[tot1,:]

    sns.set(rc={"figure.figsize":(10.5,8)}, font_scale = 1.5)
    sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":22})
    sns.set_context({"font.weight":"normal","font.size":22,"font.style":"normal","axes.labelsize":22,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":22,"ytick.labelsize":22,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":22,"legend.fontsize":22})

    ax = sns.heatmap(data = cl_inf, annot = False, cmap = "RdBu",vmin = -1, vmax = 1, cbar_kws = {"shrink": 0.8}, xticklabels = 1, yticklabels = 1)
    ax.axhline(0, color = 'k', linewidth = 1.5)
    ax.axhline(cl_inf.shape[0], color = 'k', linewidth = 1.5)
    ax.axvline(0, color = 'k', linewidth = 1.5)
    ax.axvline(cl_inf.shape[1], color = 'k', linewidth = 1.5)

    sns.set(rc={"figure.figsize":(10.5,8)}, font_scale = 1.5)
    sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":22})
    sns.set_context({"font.weight":"normal","font.size":22,"font.style":"normal","axes.labelsize":22,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":22,"ytick.labelsize":22,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":22,"legend.fontsize":22})

    plt.tight_layout()
    plt.savefig(cwd + "/FIGS/" + infl[:-4] + ".png", dpi = 400)
    plt.clf()


#SCLC HEATMAP

#ax = sns.heatmap(data = cl_inf, annot = False, cmap = sns.diverging_palette(175, 44,s = 100, l = 40, n = 8, center = "light", as_cmap = True))
