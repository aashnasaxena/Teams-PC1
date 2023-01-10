import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns 

num = pd.read_csv(num)

num1 = num.mean(axis = 1)
bio = num1.iloc[0]

df = pd.DataFrame(data = None, index = np.arnage(100), columns = ['Network', 'Number of PCs to Cover 90% Variance'])
df['Network'] = num.iloc[1:,0].to_list()
df['Number of PCs to Cover 90% Variance'] = num1.iloc[1:].to_numpy()

sns.set(rc={"figure.figsize":(10.5,8)})
sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":18})
sns.set_context({"font.weight":"normal","font.size":18,"font.style":"normal","axes.labelsize":18,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":18,"ytick.labelsize":18,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":18,"legend.fontsize":18})

sns.histplot(data = df, x = "Number of PCs to Cover 90% Variance", stat = 'probability', bins = 'auto', element = "step")
plt.axvline(bio, color = '#bd0026', linewidth = 3)
plt.tight_layout()

sns.set(rc={"figure.figsize":(10.5,8)})
sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":18})
sns.set_context({"font.weight":"normal","font.size":18,"font.style":"normal","axes.labelsize":18,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":18,"ytick.labelsize":18,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":18,"legend.fontsize":18})

plt.savefig("F1B_NumPC.png", dpi = 400)
plt.savefig("F1B_NumPC.svg", dpi = 400)
