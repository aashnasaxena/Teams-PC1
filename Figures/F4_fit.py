import pandas as pd
import numpy as np 
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import os
import glob 

# Get the current working directory
cwd = os.getcwd()

# Scaled negative hill function
def nhill_sc(x,n,thresh,sc):
    return sc*(thresh**n/(x**n + thresh**n))

dfsa = pd.read_csv("N_repSA/PC1_box.csv")
dfsi = pd.read_csv("N_repSI/PC1_box.csv")

sa = pd.DataFrame(columns = ["No. Nodes Removed","PC1 Variance"])
si = pd.DataFrame(columns = ["No. Nodes Removed","PC1 Variance"])

sa["No. Nodes Removed"] = dfsa.columns.to_list()
sa["No. Nodes Removed"] = sa["No. Nodes Removed"].astype(int)
sa["PC1 Variance"] = dfsa.mean(axis=0).to_list()
dsa = tuple(dfsa.std())

si["No. Nodes Removed"] = dfsi.columns.to_list()
si["No. Nodes Removed"] = si["No. Nodes Removed"].astype(int)
si["PC1 Variance"] = dfsi.mean(axis=0).to_list()
dsi = tuple(dfsi.std())

paraSA,pcovSA = curve_fit(nhill_sc,sa["No. Nodes Removed"], sa["PC1 Variance"], bounds=([1.01,1,0],[np.inf,np.inf,1]))
paraSI,pcovSI = curve_fit(nhill_sc,si["No. Nodes Removed"], si["PC1 Variance"], bounds=([1.01,1,0],[np.inf,np.inf,1]))
dataSA = nhill_sc(sa["No. Nodes Removed"],*paraSA)
dataSI = nhill_sc(si["No. Nodes Removed"],*paraSI)
abserrSA = dataSA - sa["PC1 Variance"]
abserrSI = dataSI - si["PC1 Variance"]
seSA = np.square(abserrSA)
seSI = np.square(abserrSI)
mseSA = np.mean(seSA)
mseSI = np.mean(seSI)
rmseSA = np.sqrt(mseSA)
rmseSI = np.sqrt(mseSI)
RsqSA = 1 - (mseSA/np.var(sa["PC1 Variance"]))
RsqSI = 1 - (mseSI/np.var(si["PC1 Variance"]))

# Plot the data and fitted curve
sns.set(rc={"figure.figsize":(9.5,8)})
sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":18})
sns.set_context({"font.weight":"normal","font.size":18,"font.style":"normal","axes.labelsize":18,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":18,"ytick.labelsize":18,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":18,"legend.fontsize":18})

ax = sns.scatterplot(x = sa["No. Nodes Removed"], y = sa["PC1 Variance"], data = sa, color = "#5e81ac",s=50)
plt.errorbar(x = sa["No. Nodes Removed"], y = sa["PC1 Variance"], xerr = None, yerr = dsa,fmt ="none",ecolor = "#5e81ac", capsize =3.5,elinewidth=0.9,markeredgewidth=0.9)
ax = sns.lineplot(x = sa["No. Nodes Removed"], y = dataSA, color = "#5e81ac",linewidth=2.5)
ax = sns.scatterplot(x = si["No. Nodes Removed"], y = si["PC1 Variance"], data = si, color = "#bf616a",marker="v",s=75)
plt.errorbar(x = si["No. Nodes Removed"], y = si["PC1 Variance"], xerr = None, yerr = dsi, fmt = "none", ecolor = "#bf616a", capsize = 3.5, elinewidth=0.9,markeredgewidth=0.9)
ax = sns.lineplot(x = si["No. Nodes Removed"], y = dataSI, color = "#bf616a", linewidth=2.5, linestyle = "dashed")
le = np.arange(0,len(sa.index),3)
ax.set_xticks(le)
labels = [int(i) for i in le]
print(labels)
ax.set_xticklabels(labels)
plt.xlabel("No. Nodes Removed")
plt.ylabel("PC1 Variance")
#plt.legend(["SA Data","SA Hill fit","SI Data","SI Hill fit"],bbox_to_anchor=(1.02,1),loc=2,borderaxespad=0,fontsize=20)
leg_el = [Line2D([],[],marker="o",color="w",markerfacecolor="#5e81ac",label="SA Data",markersize=10),Line2D([],[],linestyle = "-",color = '#5e81ac', label = "SA Hill Fit", linewidth = 3),Line2D([],[],marker="v",color="w",markerfacecolor="#bf616a",label="SI Data",markersize=10), Line2D([],[],linestyle="--",color = '#bf616a', label = "SI Hill Fit", linewidth = 3)]
#leg = ax.legend(loc=2,bbox_to_anchor=(1.02,1),borderaxespad=0,handles=leg_el,frameon=True)
leg = ax.legend(loc = "upper center", bbox_to_anchor=(0.497,-0.15),ncol = 4, columnspacing = 1.2, handles = leg_el, frameon = True)
leg.get_frame().set_edgecolor('grey')
plt.tight_layout()
plt.savefig("PC1_SA_SIfit.png",dpi=400)
plt.clf()

dfw = pd.read_csv("N_repSA/PC1_box.csv")
dfr = pd.read_csv("N_repR/PC1_box.csv")

w = pd.DataFrame(columns = ["No. Nodes Removed","PC1 Variance"])
r = pd.DataFrame(columns = ["No. Nodes Removed","PC1 Variance"])

w["No. Nodes Removed"] = dfw.columns.to_list()
w["No. Nodes Removed"] = w["No. Nodes Removed"].astype(int)
w["PC1 Variance"] = dfw.mean(axis=0).to_list()
dw = tuple(dfw.std())

r["No. Nodes Removed"] = dfr.columns.to_list()
r["No. Nodes Removed"] = r["No. Nodes Removed"].astype(int)
r["PC1 Variance"] = dfr.mean(axis=0).to_list()
dr = tuple(dfr.std())

paraW,pcovW = curve_fit(nhill_sc,w["No. Nodes Removed"], w["PC1 Variance"], bounds=([1.01,1,0],[np.inf,np.inf,1]))
paraR,pcovR = curve_fit(nhill_sc,r["No. Nodes Removed"], r["PC1 Variance"], bounds=([1.01,1,0],[np.inf,np.inf,1]))
dataW = nhill_sc(w["No. Nodes Removed"],*paraW)
dataR = nhill_sc(r["No. Nodes Removed"],*paraR)
abserrW = dataW - w["PC1 Variance"]
abserrR = dataR - r["PC1 Variance"]
seW = np.square(abserrW)
seR = np.square(abserrR)
mseW = np.mean(seW)
mseR = np.mean(seR)
rmseW = np.sqrt(mseW)
rmseR = np.sqrt(mseR)
RsqW = 1 - (mseW/np.var(w["PC1 Variance"]))
RsqR = 1 - (mseR/np.var(r["PC1 Variance"]))

sns.set(rc={"figure.figsize":(9.5,8)})
sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":18})
sns.set_context({"font.weight":"normal","font.size":18,"font.style":"normal","axes.labelsize":18,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":18,"ytick.labelsize":18,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":18,"legend.fontsize":18})

ax = sns.scatterplot(x = w["No. Nodes Removed"], y = w["PC1 Variance"], data = w, color = "#5e81ac")
plt.errorbar(x = w["No. Nodes Removed"], y = w["PC1 Variance"], xerr = None, yerr = dw,fmt ="none",ecolor = "#5e81ac", capsize =3.5,elinewidth=0.9,markeredgewidth=0.9)
ax = sns.lineplot(x = w["No. Nodes Removed"], y = dataW, color = "#5e81ac",linewidth=2.5)
ax = sns.scatterplot(x = r["No. Nodes Removed"], y = r["PC1 Variance"], data = r, color = "#bf616a",marker="v",s=75)
plt.errorbar(x = r["No. Nodes Removed"], y = r["PC1 Variance"], xerr = None, yerr = dr, fmt = "none", ecolor = "#bf616a", capsize = 3.5, elinewidth=0.9,markeredgewidth=0.9)
ax = sns.lineplot(x = r["No. Nodes Removed"], y = dataR, color = "#bf616a", linewidth=2.5, linestyle = "dashed")
plt.xlabel("No. Nodes Removed")
le = np.arange(0,len(w.index),3)
ax.set_xticks(le)
labels = [int(i) for i in le]
print(labels)
ax.set_xticklabels(labels)
plt.ylabel("PC1 Variance")
leg_el = [Line2D([],[],marker="o",color="w",markerfacecolor="#5e81ac",label="W Data",markersize=10),Line2D([],[],linestyle = "-",color = '#5e81ac', label = "W Hill Fit", linewidth = 3),Line2D([],[],marker="v",color="w",markerfacecolor="#bf616a",label="R Data",markersize=10), Line2D([],[],linestyle="--",color = '#bf616a', label = "R Hill Fit", linewidth = 3)]
#leg = ax.legend(loc=2,bbox_to_anchor=(1.02,1),borderaxespad=0,handles=leg_el,frameon=True)
leg = ax.legend(loc = "upper center", bbox_to_anchor=(0.497,-0.15),ncol = 4, columnspacing = 1.2, handles = leg_el, frameon = True)
leg.get_frame().set_edgecolor('grey')
plt.tight_layout()
plt.savefig("PC1_W_Rfit.png",dpi=400)
plt.clf()

with open(cwd+"/SA_SI_fit.txt","w") as f:
    f.write("RMSE SA "+str(rmseSA)+"\n"+"Rsq. SA "+str(RsqSA)+"\n"+"Threshold SA "+str(paraSA[1])+"\n"+"n SA "+str(paraSA[0])+"\n"+"sc SA "+str(paraSA[2])+"\n"+"RMSE SI "+str(rmseSI)+"\n"+"Rsq. SI "+str(RsqSI)+"\n"+"Threshold SI "+str(paraSI[1])+"\n"+"n SI "+str(paraSI[0])+"\n"+"sc SI "+str(paraSI[2]))

with open(cwd+"/W_R_fit.txt","w") as f:
    f.write("RMSE W "+str(rmseW)+"\n"+"Rsq. W "+str(RsqW)+"\n"+"Threshold W "+str(paraW[1])+"\n"+"n W "+str(paraW[0])+"\n"+"sc W "+str(paraW[2])+"\n"+"RMSE R "+str(rmseR)+"\n"+"Rsq. R "+str(RsqR)+"\n"+"Threshold R "+str(paraR[1])+"\n"+"n R "+str(paraR[0])+"\n"+"sc R "+str(paraR[2]))



