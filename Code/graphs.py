import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns 
from scipy.stats import spearmanr

#BOXPLOT (for PC1 variance, Hamming Distance, Team Strength):

df = pd.read_csv(path)
df.sort_values(by = "Network", inplace = True)
df1 = df.mean(axis = 1)
df2 = pd.DataFrame(data = None, index = np.arange(10), columns = np.arange(1,101,10)) #depends on randomizations and no. of networks per randomization

for k in range(1,101,10):
	num = df1.iloc[k:k+10].to_numpy()
	df2[k] = num 

df2.columns = np.arange(10,101,10)

sns.boxplot(data = df2)
plt.xlabel('Number of Swaps', {"fontweight":"light","fontsize":10})
plt.ylabel('PC1 Variance',{"fontweight":"light","fontsize":10})
plt.tight_layout()
plt.savefig(path, dpi = 400, pad_inches = 0)
plt.clf()


#SCATTERPLOT (CTS vs. ITS, PC1 variance vs. Team Strength, PCA coefficients vs. Team Strength):

df_pc1 = pd.read_csv(path)
df_ts = = pd.read_csv(path)
df_pc1.sort_values(by = "Network", inplace = True)
df_ts.sort_values(by = "Network", inplace = True)
pc1 = df_pc1.mean(axis = 1)
ts = df_ts.mean(axis = 1)

df1 = pd.DataFrame(data = None, index = np.arange(10), columns = np.arange(1,101,10)) #index = no. of networks per randomization, cols = randomizations
df2 = pd.DataFrame(data = None, index = np.arange(10), columns = np.arange(1,101,10))

for k in range(1,101,10):
	num = pc1.iloc[k:k+10].to_numpy()
	df1[k] = num 

for k in range(1,101,10):
	num = ts.iloc[k:k+10].to_numpy()
	df2[k] = num 

d1 = [0] + list(df1.std())
c1 = tuple(d1)
d2 = [0] + list(df2.std())
c2 = tuple(d2)

Final = pd.DataFrame(data = None, index = np.arange(11), columns = ['Number of Swaps', 'PC1 Variance', 'Team Strength'])

Final['Number of Swaps'] = [0] + np.arange(10,101,10).tolist()
Final['PC1 Variance'] = [pc1.iloc[0]] + df1.mean(axis = 0).to_list()
Final['Team Strength'] = [ts.iloc[0]] + df2.mean(axis = 0).to_list()

sns.scatterplot(x = "PC1 Variance", y = "Team Strength", hue = "Number of Swaps", data = Final, palette = sns.color_palette("Spectral",11))
plt.errorbar(Final['PC1 Variance'], Final['Team Strength'], xerr = c1, yerr = c2, fmt = "none", marker = "none", elinewidth = 0.2)
plt.legend(bbox_to_anchor=(1.02,1), loc = 2, borderaxespad = 0)
plt.xlabel("PC1 Variance", {"fontweight":"light", "fontsize":10})
plt.ylabel("Team Strength", {"fontweight":"light", "fontsize":10})
plt.tight_layout()
plt.savefig(path, dpi = 400, pad_inches = 0)
plt.clf()

spearmanr(Final['PC1 Variance'], Final['Team Strength'])


#HISTOGRAM (PC1 variance, Number of PCs):

df = pd.read_csv(path)
df1 = df.mean(axis = 1)
ind = len(df1.index)
bio = df1.iloc[0]
df2. = pd.DataFrame(data = None, index = np.arange(ind), columns = ['Network', 'PC1 Variance'])
df2['Network'] = df.iloc[1:,0].to_list()
df2['PC1 Variance'] = df1.iloc[1:].to_list()

sns.histplot(data = df2, x = "PC1 Variance", y = None, stat = 'probability', bins = 'auto')
plt.axvline(bio, color = 'k', linewidth = 2)
plt.xlabel("PC1 Variance", {"fontweight":"light", "fontsize": 10})
plt.ylabel("Percentage", {"fontweight":"light", "fontsize": 10})
plt.tight_layout()
plt.savefig(path, dpi = 400, pad_inches = 0)
plt.clf()
