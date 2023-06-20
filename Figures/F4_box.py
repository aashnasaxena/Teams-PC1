import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns 
from scipy.stats import f_oneway
from itertools import chain 

#To plot PC1 or ITS box plots
PC1 = True
ITS = True
num = True
if ITS == True:
    df = pd.read_csv("TeamsCalc/CompiledTS.csv")
    df1 = pd.DataFrame(columns = ["Team Strength","No. Nodes Removed"])
    df1['Team Strength'] = df['TS']
    df1['No. Nodes Removed'] = [0]+list(chain.from_iterable([list(i*np.ones(100)) for i in range(1,19)]))
    df1['No. Nodes Removed'] = df1['No. Nodes Removed'].astype(int)
    df1['Team Strength'] = df1['Team Strength'].astype(str)
    df2 = df1.iloc[1:,:]
    ind = []
    for i in range(int(len(df2.index)/100)):
        ind += np.arange(100).tolist()
    df2['index'] = ind
    piv = pd.pivot(df2,index=['index'],columns=['No. Nodes Removed'])
    piv.columns = [j for i,j in piv.columns]
    piv.reset_index(inplace=True,drop=True)
    piv[0] = df1.iloc[0,0]
    piv.sort_index(axis=1, inplace = True)
    print(piv)
    piv.to_csv("ITS_box.csv",index=False)
    df = pd.read_csv("ITS_box.csv")
    sns.set(rc={"figure.figsize":(10.5,8)})
    sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":18})
    sns.set_context({"font.weight":"normal","font.size":18,"font.style":"normal","axes.labelsize":18,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":18,"ytick.labelsize":18,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":18,"legend.fontsize":18})  
    ax = sns.boxplot(data=df,palette=sns.blend_palette(["#9f0245","#f6b4b0"],n_colors=19))
    #ax = sns.boxplot(data=df,palette = "Greens_r")
    #ax = sns.boxplot(data=df,palette = "Purples_r")
    #ax = sns.boxplot(data = df, palette = "Blues_r")
    plt.setp(ax.collections, edgecolor='k')
    plt.xlabel("No. of Nodes Removed",{"fontsize":24,"fontweight":"bold","fontfamily":"sans-serif"})
    plt.ylabel("Team Strength",{"fontsize":24,"fontweight":"bold","fontfamily":"sans-serif"})
    plt.tight_layout()
    plt.savefig("ITS_box.png",dpi=400)
    plt.clf()

    res = f_oneway(*[piv[col] for col in piv.columns])
    print(res)
    with open("ITS_anova.txt","w") as f:
        f.write("F-statistic "+str(res[0])+"\n"+"p-value "+str(res[1]))
 
if PC1 == True:
    df = pd.read_csv("PC1variance.csv")
    df1 = pd.DataFrame(columns = ["PC1 Variance","No. Nodes Removed"])
    df1['PC1 Variance'] = df.iloc[:,1:].mean(axis=1).to_numpy()
    df1['No. Nodes Removed'] = [0]+list(chain.from_iterable([list(i*np.ones(100)) for i in range(1,19)]))
    df1['No. Nodes Removed'] = df1['No. Nodes Removed'].astype(int)
    df1['PC1 Variance'] = df1['PC1 Variance'].astype(str)
    df2 = df1.iloc[1:,:]
    ind = []
    for i in range(int(len(df2.index)/100)):
        ind += np.arange(100).tolist()
    df2['index'] = ind
    piv = pd.pivot(df2,index=['index'],columns=['No. Nodes Removed'])
    piv.columns = [j for i,j in piv.columns]
    piv.reset_index(inplace=True,drop=True)
    piv[0] = df1.iloc[0,0]
    piv.sort_index(axis=1, inplace = True)
    print(piv)
    piv.to_csv("PC1_box.csv",index=False)
    df = pd.read_csv("PC1_box.csv")
    sns.set(rc={"figure.figsize":(10.5,8)})
    sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":18})
    sns.set_context({"font.weight":"normal","font.size":18,"font.style":"normal","axes.labelsize":18,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":18,"ytick.labelsize":18,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":18,"legend.fontsize":18})  
    ax = sns.boxplot(data=df,palette=sns.blend_palette(["#9f0245","#f6b4b0"],n_colors=19))
    #ax = sns.boxplot(data=df,palette = "Greens_r")
    #ax = sns.boxplot(data=df,palette = "Purples_r")
    #ax = sns.boxplot(data = df, palette = "Blues_r")
    plt.setp(ax.collections, edgecolor='k')
    plt.xlabel("No. of Nodes Removed",{"fontsize":24,"fontweight":"bold","fontfamily":"sans-serif"})
    plt.ylabel("PC1 Variance",{"fontsize":24,"fontweight":"bold","fontfamily":"sans-serif"})
    plt.tight_layout()
    plt.savefig("PC1_box.png",dpi=400)
    plt.clf()

    res = f_oneway(*[piv[col] for col in piv.columns])
    print(res)
    with open("PC1_anova.txt","w") as f:
        f.write("F-statistic "+str(res[0])+"\n"+"p-value "+str(res[1]))

if num == True:
    df = pd.read_csv("Num_PCs.csv")
    df1 = pd.DataFrame(columns = ["Number of PCs","No. Nodes Removed"])
    df1['Number of PCs'] = df.iloc[:,1:].mean(axis=1).to_numpy()
    df1['No. Nodes Removed'] = [0]+list(chain.from_iterable([list(i*np.ones(100)) for i in range(1,19)]))
    df1['No. Nodes Removed'] = df1['No. Nodes Removed'].astype(int)
    df1['Number of PCs'] = df1['Number of PCs'].astype(str)
    df2 = df1.iloc[1:,:]
    ind = []
    for i in range(int(len(df2.index)/100)):
        ind += np.arange(100).tolist()
    df2['index'] = ind
    piv = pd.pivot(df2,index=['index'],columns=['No. Nodes Removed'])
    piv.columns = [j for i,j in piv.columns]
    piv.reset_index(inplace=True,drop=True)
    piv[0] = df1.iloc[0,0]
    piv.sort_index(axis=1, inplace = True)
    print(piv)
    piv.to_csv("Num_box.csv",index=False)
    df = pd.read_csv("Num_box.csv")
    sns.set(rc={"figure.figsize":(10.5,8)})
    sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":18})
    sns.set_context({"font.weight":"normal","font.size":18,"font.style":"normal","axes.labelsize":18,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":18,"ytick.labelsize":18,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":18,"legend.fontsize":18})  
    ax = sns.boxplot(data=df,palette=sns.blend_palette(["#9f0245","#f6b4b0"],n_colors=19))
    #ax = sns.boxplot(data=df,palette = "Greens_r")
    #ax = sns.boxplot(data=df,palette = "Purples_r")
    #ax = sns.boxplot(data = df, palette = "Blues_r")
    plt.setp(ax.collections, edgecolor='k')
    plt.xlabel("No. of Nodes Removed",{"fontsize":24,"fontweight":"bold","fontfamily":"sans-serif"})
    plt.ylabel("Number of PCs",{"fontsize":24,"fontweight":"bold","fontfamily":"sans-serif"})
    plt.tight_layout()
    plt.savefig("Num_box.png",dpi=400)
    plt.clf()

    res = f_oneway(*[piv[col] for col in piv.columns])
    print(res)
    with open("Num_anova.txt","w") as f:
        f.write("F-statistic "+str(res[0])+"\n"+"p-value "+str(res[1]))

