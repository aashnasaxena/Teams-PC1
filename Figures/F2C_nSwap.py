import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns 
import os
import glob 
from itertools import chain 

df = pd.read_csv("Consistency.csv")

#df["Number of Swaps"] = [0]+list(chain.from_iterable([list(i*np.ones(50)) for i in range(10,101,10)]))
#df["Number of Swaps"] = df["Number of Swaps"].astype(int)
df1 = df.iloc[1:,:]
ind = []
for i in range(int(len(df1.index)/50)):
    ind += np.arange(50).tolist()

df1["index"] = ind
pc1 = pd.pivot(df1[["index","Number of Swaps","PC1 Variance"]], index = ["index"], columns = ["Number of Swaps"])
pc1.columns = [j for i,j in pc1.columns]
pc1.reset_index(inplace=True,drop=True)
pc1[0] = df.loc[0,"PC1 Variance"]
pc1.sort_index(axis=1,inplace=True)
pcm = pc1.mean(axis=0).to_list()
pcml = [pcm[0]]
for i in pcm[1:]:
    pcml += [i]*50
df["Mean PC1"] = pcml

ts = pd.pivot(df1[["index","Number of Swaps","Influence Team Strength"]], index = ["index"], columns = ["Number of Swaps"])
ts.columns = [j for i,j in ts.columns]
ts.reset_index(inplace=True,drop=True)
ts[0] = df.loc[0,"Influence Team Strength"]
ts.sort_index(axis=1,inplace=True)
tsm = ts.mean(axis=0).to_list()
tsml = [tsm[0]]
for i in tsm[1:]:
    tsml += [i]*50
df["Mean TS"] = tsml

final = pd.DataFrame(columns = ["Number of Swaps","Influence Team Strength","PC1 Variance"])
final["Number of Swaps"] = df["Number of Swaps"].drop_duplicates().to_list()
final["Influence Team Strength"] = ts.mean(axis=0).to_list()
final["PC1 Variance"] = pc1.mean(axis=0).to_list()
dpc1 = tuple(pc1.std())
dts = tuple(ts.std())

sns.set(rc={"figure.figsize":(10.5,8)})
sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":18})
sns.set_context({"font.weight":"normal","font.size":18,"font.style":"normal","axes.labelsize":18,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":18,"ytick.labelsize":18,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":18,"legend.fontsize":18})

#clr = sns.cubehelix_palette(start=.5, rot=-.5,dark=0.1,light=0.7, as_cmap=True) #EMT
#clr = "Blues" #EMT box
#clr = sns.cubehelix_palette(start=2, rot=0, dark=0.1, light=.70, as_cmap=True) #SCLC
#clr = "Greens" #SCLC box
#clr = sns.cubehelix_palette(light = 0.7,as_cmap=True) #pluripotency
#clr = "Purples" #pluripotency box
clr = sns.blend_palette(["#f6b4b0","#9f0245"],as_cmap=True) #gonadal
ax = sns.boxplot(data = df, x = "Number of Swaps",y = "Influence Team Strength", hue = "Mean PC1", palette = sns.blend_palette(["#f6b4b0","#9f0245"],n_colors=11),dodge = False)
#ax = sns.scatterplot(data = final, x = "Number of Swaps",y = "Influence Team Strength",hue = "PC1 Variance", palette = clr, s = 80)
#plt.errorbar(x = final["Number of Swaps"], y = final["Influence Team Strength"], xerr = None, yerr = dts,fmt ="none",ecolor = "black", capsize =3.5,elinewidth=0.8,markeredgewidth=0.8)
#ax.set_xticks(final["Number of Swaps"].to_list())
#ax.set_ylim(0,1)
ax.get_legend().remove()
cb = plt.get_cmap(clr)
norm = plt.Normalize(min(df["Mean PC1"]), max(df["Mean PC1"])) #min-max of hue column
sm = plt.cm.ScalarMappable(norm = norm, cmap = cb)
sm.set_array([])
cbar = plt.colorbar(sm, ax = ax,aspect=20,shrink=1.0,pad=0.07)
#cbticks = np.linspace(round(np.ceil(min(final["Influence Team Strength"])*10)/10,1),round(np.floor(max(final["Influence Team Strength"])*10)/10,1),num=5,endpoint=True)
#cbar.set_ticks(cbticks)
cbar.set_label("PC1 Variance",labelpad = -40,y=1.05,rotation=0)
#cbar.outline.set_visible(False)
plt.xlabel("Number of Swaps",{"fontsize":24,"fontweight":"bold","fontfamily":"sans-serif"})
plt.ylabel("Team Strength",{"fontsize":24,"fontweight":"bold","fontfamily":"sans-serif"})
plt.tight_layout()
plt.savefig("ns_TSbox.png",dpi=400)
plt.clf()
