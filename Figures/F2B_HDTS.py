import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns 

df = pd.read_csv("EMT_F1-F2.csv")
df100 = df.iloc[201:302,:]
df100.reset_index(inplace = True)
df100.drop(columns = "index", inplace = True)
final = pd.DataFrame(data = None, index = np.arange(101),columns = ['Network','Influence Team Strength','Hamming Distance'])
final['Network'] = "Random"
final.iloc[0,0] = "Biological"
final['Influence Team Strength'] = df100['Influence Team Strength']
final['Hamming Distance'] = df100['Hamming Distance']

sns.set(rc={"figure.figsize":(10.5,8)})
sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":18})
sns.set_context({"font.weight":"normal","font.size":18,"font.style":"normal","axes.labelsize":18,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":18,"ytick.labelsize":18,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":18,"legend.fontsize":18})

ax = sns.scatterplot(data = final, x = "Influence Team Strength", y = "Hamming Distance", hue = "Network", palette = ['#0570b0','#74a9cf'], s = 60, edgecolor = 'black')
plt.legend(bbox_to_anchor=(1.02,1), loc = 2, borderaxespad = 0)
plt.tight_layout()

sns.set(rc={"figure.figsize":(10.5,8)})
sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":18})
sns.set_context({"font.weight":"normal","font.size":18,"font.style":"normal","axes.labelsize":18,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":18,"ytick.labelsize":18,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":18,"legend.fontsize":18})

plt.savefig("F2B_HDTS.png", dpi = 400)
plt.savefig("F2B_HDTS.svg", dpi = 400)