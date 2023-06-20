import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os 
import glob 

# Get the current working directory
cwd = os.getcwd()

df = pd.read_csv("GON_box.csv")

sns.set(rc={"figure.figsize":(10.5,8)})
sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":18})
sns.set_context({"font.weight":"normal","font.size":18,"font.style":"normal","axes.labelsize":18,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":18,"ytick.labelsize":18,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":18,"legend.fontsize":18})

props = {"boxprops":{"facecolor":"#74a9cf","edgecolor":"black"},"medianprops":{"color":"black"},"whiskerprops":{"color":"black"}}
ax = sns.boxplot(data = df,color="#e8979e")
plt.setp(ax.collections, facecolor = "#e8979e",edgecolor = "black")
for i, artist in enumerate(ax.artists):
    artist.set_edgecolor("black")
    for j in range(i*6,i*6+6):
        line = ax.lines[j]
        line.set_color("black")
        line.set_mfc("black")
        line.set_mec("black")

plt.xlabel("Influence Team Strength",{"fontsize":24,"fontweight":"bold","fontfamily":"sans-serif"})
plt.ylabel("Number of PCs",{"fontsize":24,"fontweight":"bold","fontfamily":"sans-serif"})
ax.set_xticklabels(df.columns.to_list(),rotation=45)
plt.tight_layout()
plt.savefig(cwd+"/FIGS/"+"numPC_box.png",dpi=400)
