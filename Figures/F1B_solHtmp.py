import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
from sklearn.cluster import KMeans
from scipy.stats import zscore

# Get the current working directory
cwd = os.getcwd()

# List of g/k normalized files
normLi = sorted(glob.glob("*norm.csv"))

if not os.path.exists("FIGS"):
    os.makedirs("FIGS")
    
for fl in normLi:
    sol = pd.read_csv(fl)
    sol = sol.iloc[:,3:]
    df = pd.read_table(cwd+"/TeamsCalc/"+fl.replace("_norm.csv",".teams"), delimiter = ":", names = ["Team/Level","Nodes"])
    df.dropna(inplace=True,ignore_index=True)
    p = ",".join(df.loc[2:,"Nodes"].to_list())
    peri = p.split(",")
    t1 = df.loc[0,"Nodes"].split(",")
    t2 = df.loc[1,"Nodes"].split(",")
    Nodes = t1 + peri + t2
    print(Nodes)
    zn = sol.apply(zscore)
    #print(zn)
    zn.dropna(inplace=True,ignore_index=True)
    zn = zn.reindex(columns = Nodes)
    print(zn)

    kmeans = KMeans(n_clusters = 3).fit(zn)
    zn['Team'] = kmeans.predict(zn)

    zn.sort_values(by = "Team", inplace = True)
    zn.reset_index(inplace = True,drop=True)
    #zn.drop(columns = "index", inplace = True)

    for row in zn.index:
        if zn.loc[row,'Team'] == 1:
            t = row
            break

    for row in zn.index[t:]:
        if zn.loc[row,'Team'] == 2:
            tt = row
            break

    T1 = zn.iloc[:t,:]
    T2 = zn.iloc[t:tt,:]
    T3 = zn.iloc[tt:,:]

    T1.drop(columns = "Team", inplace = True)
    T2.drop(columns = "Team", inplace = True)
    T3.drop(columns = "Team", inplace = True)
    dfli = [T1,T2,T3]
    for a in range(len(dfli)):
        sns.set(rc={"figure.figsize":(10.5,8)})
        sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":14})
        sns.set_context({"font.weight":"normal","font.size":14,"font.style":"normal","axes.labelsize":14,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":14,"ytick.labelsize":14,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":14,"legend.fontsize":14})

        ax = sns.heatmap(data = dfli[a], annot = False, cmap = sns.diverging_palette(220, 20, as_cmap=True), vmin = -3, vmax = 3, xticklabels = 1, yticklabels = False)
        ax.axhline(0, color = 'k', linewidth = 1.5)
        ax.axhline(dfli[a].shape[0], color = 'k', linewidth = 1.5)
        ax.axvline(0, color = 'k', linewidth = 1.5)
        ax.axvline(dfli[a].shape[1], color = 'k', linewidth = 1.5)
        
        plt.tight_layout()
        plt.savefig(cwd+"/FIGS/"+fl.replace("_norm.csv","_")+str(a)+".png",dpi=400)
        plt.clf()

