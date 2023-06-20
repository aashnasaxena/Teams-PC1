import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns 

df = pd.read_csv("EMT.csv")
final = df.loc[:,["Network","Number of PCs","Hamming Distance"]
final['Network'] = "Random"
final.loc[0,"Network"] = "Biological"

sns.set(rc={"figure.figsize":(10.5,8)})
sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":18})
sns.set_context({"font.weight":"normal","font.size":18,"font.style":"normal","axes.labelsize":18,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":18,"ytick.labelsize":18,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":18,"legend.fontsize":18})

ax = sns.scatterplot(data = final, x = "Number of PCs", y = "Hamming Distance", hue = "Network", palette = ['#0570b0','#74a9cf'], s = 60, edgecolor = 'black')
plt.legend(bbox_to_anchor=(1.02,1), loc = 2, borderaxespad = 0)
plt.tight_layout()
plt.savefig("F2B_HDNum.png", dpi = 400)
plt.savefig("F2B_HDNum.svg", dpi = 400)
plt.clf()
