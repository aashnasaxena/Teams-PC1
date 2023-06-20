import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns 

df = pd.read_csv("EMT.csv")
bio = df.loc[0,'Hamming Distance']
final = df.loc[1:,"Hamming Distance"]

sns.set(rc={"figure.figsize":(10.5,8)})
sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":18})
sns.set_context({"font.weight":"normal","font.size":18,"font.style":"normal","axes.labelsize":18,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":18,"ytick.labelsize":18,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":18,"legend.fontsize":18})

sns.histplot(data = final, x = "Hamming Distance", stat = 'probability', bins = 'auto', element = 'step')
plt.axvline(bio, color = '#bd0026', linewidth = 3)
plt.tight_layout()
plt.savefig("F2B_HDhist.png", dpi = 400)
plt.savefig("F2B_HDhist.svg", dpi = 400)
plt.clf()
