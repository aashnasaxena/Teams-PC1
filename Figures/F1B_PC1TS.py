import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns 

inf = pd.read_csv(inf)
inf1 = inf.mean(axis = 1)
var = pd.read_csv(var)
var1 = var.mean(axis = 1)

df = pd.DataFrame(data = None, index = np.arange(101), columns = ['Network','Influence Team Strength','PC1 Variance'])
df['Network'] = "Random"
df.loc[0,'Network'] = "Biological"
df['Influence Team Strength'] = inf1.iloc[:].to_numpy()
df['PC1 Variance'] = var1.iloc[:].to_numpy()

sns.set(rc={"figure.figsize":(10.5,8)})
sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":18})
sns.set_context({"font.weight":"normal","font.size":18,"font.style":"normal","axes.labelsize":18,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":18,"ytick.labelsize":18,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":18,"legend.fontsize":18})

ax = sns.scatterplot(data = df1, x = "Influence Team Strength", y = "PC1 Variance", hue = "Network", palette = ['#0570b0','#74a9cf'],s = 60,edgecolor = "black")
plt.legend(bbox_to_anchor=(1.02,1),loc=2,borderaxespad=0)
plt.tight_layout()

sns.set(rc={"figure.figsize":(10.5,8)})
sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":18})
sns.set_context({"font.weight":"normal","font.size":18,"font.style":"normal","axes.labelsize":18,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":18,"ytick.labelsize":18,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":18,"legend.fontsize":18})

plt.savefig("F1B_PC1TS.png", dpi = 400)
plt.savefig("F1B_PC1TS.svg", dpi = 400)
