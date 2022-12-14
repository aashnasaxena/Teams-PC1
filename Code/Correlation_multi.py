import pandas as pd 
import numpy as np 
from scipy.cluster import hierarchy as hi 
import glob 
import os
import subprocess
from multiprocessing import Pool 
import matplotlib.pyplot as plt
import seaborn as sns 

cwd = os.getcwd()
tpdir = "/home/csb/Aashna/Range10/TOPO"
tlog = open("correlation_teams.txt","w")


os.chdir(tpdir)
results = (glob.glob("*.topo"))
tpfl = sorted(results)

os.chdir(cwd + "/Results")

#Get Correlation Matrix
def Correlation_Matrix(sol_path,cwd,t,n):
	df = pd.read_csv(sol_path)
	df1 = df.iloc[:,3:]
	df2 = df1.corr(method = 'spearman')
	#nodes = df1.columns
	#d = hi.distance.pdist(df2)
	#L = hi.linkage(d, method = 'complete')
	#clust = hi.cut_tree(L, n_clusters = 2)
	#cluster = np.transpose(clust)
	#t1 = nodes[cluster[0] == 0]
	#t2 = nodes[cluster[0] == 1]
	#tot = t1.append(t2)
	#df_clust = df2.loc[tot,:].T.loc[tot,:]
	df2.sort_index(axis = 0, inplace = True)
	df2.sort_index(axis = 1, inplace = True)
	df2.to_csv(cwd + "/Results/" + t[:-5] + "/" + str(n) + "/" + t[:-5] + "_correlation.csv")

#Get Correlation Figure
def Correlation_fig(correlation_path, cwd, t,tpfl,n):
	df_clust = pd.read_csv(correlation_path)
	df_clust.set_index("Unnamed: 0", inplace = True)
	d = hi.distance.pdist(df_clust)
	L = hi.linkage(d, method = "complete")
	sns.set(rc={'figure.figsize':(16,16)})
	sns.set_context("paper", rc={"font.weight":'bold',"legend.fontsize":8,"legend.title_fontsize":10,"font.size":8,"axes.titlesize":8,"axes.labelsize":8,"xtick.labelsize":10,"ytick.labelsize":10})
	ax = sns.clustermap(data = df_clust, method = "complete", annot = True, cmap = "coolwarm", row_linkage = L, col_linkage = L)
	plt.savefig(cwd + "/Results/" + tpfl[t][:-5] + "/" + str(n) + "/" + tpfl[t][:-5] + "_correlation.png", dpi = 400, pad_inches = 0)
	plt.clf()

#Get Team Strength
def Team_Strength(correlation_path,t,tpfl):
	df_clust = pd.read_csv(correlation_path)
	df_clust.set_index("Unnamed: 0", inplace = True)
	df_clust.drop(index = ['miR205','miR30c','miR9','VIM','CDH1','KLF8','TCF3'], columns = ['miR205','miR30c','miR9','VIM','CDH1','KLF8','TCF3'], inplace = True)
	nodes = df_clust.columns
	d = hi.distance.pdist(df_clust)
	L = hi.linkage(d, method = "complete")
	clust = hi.cut_tree(L, n_clusters = 2)
	cluster = np.transpose(clust)
	t1 = nodes[cluster[0] == 0]
	t2 = nodes[cluster[0] == 1]
	tot = t1.append(t2)
	df_cluster = df_clust.loc[tot,:].T.loc[tot,:]
	
	team1 = ""
	team2 = ""
	for g in t1:
		team1 = team1 + "," + g
	team1 = team1.replace(",","",1)
	for g in t2:
		team2 = team2 + "," + g
	team2 = team2.replace(",","",1)
	tlog.write(tpfl[t][:-5] + str(n) + " t1 " + str(len(t1)) + " " + team1 + "\n" + tpfl[t][:-5] + str(n) + " t2 " + str(len(t2)) + " " + team2 + "\n\n")

	df_t11 = df_cluster.loc[t1,t1]
	for g in t1:
		df_t11.loc[g,g] = 0
	df_t22 = df_cluster.loc[t2,t2]
	for g in t2:
		df_t22.loc[g,g] = 0
	df_t12 = df_cluster.loc[t1,t2]
	df_t21 = df_cluster.loc[t2,t1]

	num_t11 = df_t11.to_numpy()
	num_t22 = df_t22.to_numpy()
	num_t12 = df_t12.to_numpy()
	num_t21 = df_t21.to_numpy()

	t11 = abs(np.sum(num_t11, axis = None))/(len(num_t11)*(len(num_t11)-1))
	t22 = abs(np.sum(num_t22, axis = None))/(len(num_t22)*(len(num_t22)-1))
	t12 = abs(np.sum(num_t12, axis = None))/(len(num_t12)*len(num_t12[0]))
	t21 = abs(np.sum(num_t21, axis = None))/(len(num_t21)*len(num_t21[0]))

	ts = (t11 + t22 + t12 + t21)/4

	return ts 

pool1 = Pool(80)
pool1.starmap(Correlation_Matrix, [(glob.glob(cwd + "/Results/"+ t[:-5] + "/" + str(n) + "/*norm.csv")[0],cwd,t,n) for t in tpfl for n in range(1,4)])

pool1.close()
pool1.join()

df_cor = pd.DataFrame(data = None, index = np.arange(len(tpfl)), columns = ['Network',1,2,3])
names = []
for t in range(len(tpfl)):
	names.append(tpfl[t][:-5])
	for n in range(1,4):
		correlation_path = glob.glob(cwd + "/Results/" + tpfl[t][:-5] + "/" + str(n) + "/*correlation.csv")[0]
		ts = Team_Strength(correlation_path,t,tpfl)
		df_cor.iloc[t,n] = ts 

df_cor['Network'] = names
df_cor.sort_values(by = "Network", inplace = True)
df_cor.to_csv(cwd + "/correlation_ts.csv", index = False)







