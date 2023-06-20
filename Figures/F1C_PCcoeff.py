import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import os
import glob 

# Get the current working directory
cwd = os.getcwd()

# Get the list of .teams files
tmLi = sorted(glob.glob("TeamsCalc/*.teams"))
print(tmLi)

# List of g/k normalized files
normLi = sorted(glob.glob("*norm.csv"))

for fl in normLi:
    df = pd.read_table("TeamsCalc/"+fl.replace("_norm.csv",".teams"),delimiter=":",names=["Team/Level","Nodes"])
    df.dropna(inplace=True,ignore_index=True)
    t1 = df.loc[0,"Nodes"].split(",")
    t2 = df.loc[1,"Nodes"].split(",")
    p = ",".join(df.loc[2:,"Nodes"].to_list())
    peri = p.split(",")
    #coeff = pd.read_csv("PC_coefficients.csv")
    sol = pd.read_csv(fl)
    sol = sol.iloc[:,3:]
    net = fl.replace("_norm.csv","")
    pca = PCA(n_components=1)
    df_pca = pca.fit_transform(sol)
    coefficients = pca.components_
    print(coefficients)
    nodes = sol.columns.to_list()
    dt = pd.DataFrame(columns = ["Node",net])
    dt["Node"] = nodes
    dt[net] = coefficients[0]
    #numnodes = int(len(coeff.index)/3)
    #dt = coeff.loc[:numnodes-1,[net,"Node"]]
    dt.set_index("Node",inplace=True)
    dt.drop(index=peri,inplace=True)
    dt.sort_values(by=net,inplace=True)
    
    clrs = ['#023858' if (x in t2) else '#3690c0' for x in dt.index.to_list()]
    
    sns.set(rc={"figure.figsize":(6.5,9.5)})
    sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":18})
    sns.set_context({"font.weight":"normal","font.size":18,"font.style":"normal","axes.labelsize":18,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":18,"ytick.labelsize":18,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":18,"legend.fontsize":18})

    ax = sns.barplot(x = dt[net], y = dt.index, palette = clrs)
    ax.set_yticklabels(dt.index.to_list())
    plt.ylabel("",{"fontsize":24,"fontweight":"bold","fontfamily":"sans-serif"})
    plt.xlabel("PC1 Coefficient",{"fontsize":24,"fontweight":"bold","fontfamily":"sans-serif"})
    plt.tight_layout()
    plt.savefig(cwd+"/FIGS/"+net+"_barplot.png",dpi=400)
    plt.clf()
