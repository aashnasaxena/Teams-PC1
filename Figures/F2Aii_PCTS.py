import pandas as pd 
import numpy as np  
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns 

df = pd.read_csv("Consistency.csv")
df.sort_values(by = "Influence Team Strength",ascending = False, inplace = True, ignore_index = True)
print(df)

sns.set(rc={"figure.figsize":(10.5,8)})
sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":18})
sns.set_context({"font.weight":"normal","font.size":18,"font.style":"normal","axes.labelsize":18,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":18,"ytick.labelsize":18,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":18,"legend.fontsize":18})

ax = sns.scatterplot(data = df, x = df.loc[1:,'Influence Team Strength'], y = df.loc[1:,'PC1 Variance'], s = 55, color ='#74a9cf', alpha = 0.6, edgecolor = None)
bx = plt.scatter(x = df.loc[0,"Influence Team Strength"], y = df.loc[0,"PC1 Variance"], color = '#bd0026', s = 50, alpha = 1, edgecolor = None)
plt.yticks([0.2,0.4,0.6,0.8,1.0])

# arrow base x and y coordinates
abx = df.loc[0,"Influence Team Strength"]-0.015
aby = df.loc[0,"PC1 Variance"]-0.14
plt.arrow(abx,aby,0.012,0.11, head_width = 0.010, color = 'k', alpha = 1, width = 0.002)
plt.text(abx-0.02, aby-0.025,"Biological", fontsize = 16, fontweight = "bold",color = 'k', ha = 'center', va = 'center')
plt.yticks([0.2,0.4,0.6,0.8,1.0])
plt.tight_layout()

plt.savefig("PC1_TS.png",dpi=400)
plt.clf()

