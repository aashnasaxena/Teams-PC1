import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns 

var = pd.read_csv()
var1 = var.mean(axis = 1)

bio = var1.iloc[0]

df = pd.DataFrame(data = None, index = np.arange(100), columns = ['Network', 'PC1 Variance'])
df['Network'] = var.iloc[1:,0].to_list()
df['PC1 Variance'] = var1.iloc[1:].to_numpy()

sns.set(rc={"figure.figsize":(10.5,8)})
sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":18})
sns.set_context({"font.weight":"normal","font.size":18,"font.style":"normal","axes.labelsize":18,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":18,"ytick.labelsize":18,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":18,"legend.fontsize":18})

sns.histplot(data = df, x = "PC1 Variance", stat = 'probability', bins = 'auto', element = 'step')
plt.axvline(bio, color = '#bd0026', linewidth = 3)
plt.xlim(0,1)
plt.savefig("F1A_PC1dom.png", dpi = 400, pad_inches = 0)
plt.savefig("F1A_PC1dom.svg", dpi = 400)

sns.set(rc={"figure.figsize":(10.5,8)})
sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":18})
sns.set_context({"font.weight":"normal","font.size":18,"font.style":"normal","axes.labelsize":18,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":18,"ytick.labelsize":18,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":18,"legend.fontsize":18})
