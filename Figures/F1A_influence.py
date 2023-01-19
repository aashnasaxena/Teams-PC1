import pandas as pd 
import numpy as np 
from scipy.cluster import hierarchy as hi 
import matplotlib.pyplot as plt 
import seaborn as sns 

inf = pd.read_csv()

INF = inf.iloc[:,:23]
genes = INF['Unnamed: 0'].to_list()
INF.columns = [" "] + genes
INF.set_index(" ", inplace = True)
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

sns.set(rc={"figure.figsize":(10.5,8)})
sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":18})
sns.set_context({"font.weight":"normal","font.size":18,"font.style":"normal","axes.labelsize":18,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":18,"ytick.labelsize":18,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":18,"legend.fontsize":18})

ax = sns.heatmap(data = cl_inf, annot = False, cmap = "RdBu", vmin = -1, vmax = 1, cbar_kws = {"ticks":[-1,-0.5,0,0.5,1], "shrink": 0.8}, xticklabels = 1, yticklabels = 1)
ax.axhline(0, color = 'k', linewidth = 1.5)
ax.axhline(cl_inf.shape[0], color = 'k', linewidth = 1.5)
ax.axvline(0, color = 'k', linewidth = 1.5)
ax.axvline(cl_inf.shape[1], color = 'k', linewidth = 1.5)

#ax = sns.clustermap(data = INF, method = 'complete', annot = False, row_linkage = L, col_linkage = L, tree_kws={"linewidths":0}, cmap = 'RdBu',vmin = -1, vmax = 1,cbar_kws ={"ticks":[-1,-0.5,0,0.5,1],"shrink":0.6})

plt.tight_layout()
plt.savefig("F1A_infl.png", dpi = 400, pad_inches = 0)
plt.savefig("F1A_infl.svg", dpi = 400)
