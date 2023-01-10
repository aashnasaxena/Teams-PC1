import pandas as pd 
import numpy as np 
import os
import glob 
from multiprocessing import Pool
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans 
import math 

cwd = os.getcwd()

tlog = open("PCA12_teams.txt","w")

tpdir = "/home/csb/Aashna/Range10/TOPO"

os.chdir(tpdir)
results = glob.glob("*.topo")
tpfl = sorted(results)

def pc1var(sol_path):
	df = pd.read_csv(sol_path)
	df_sol = df.iloc[:,3:]
	l = len(df_sol.columns)
	pca = PCA(n_components = l)
	df_pca = pca.fit_transform(df_sol)
	variance = pca.explained_variance_ratio_
	return variance[0]

def pc2var(sol_path):
	df = pd.read_csv(sol_path)
	df_sol = df.iloc[:,3:]
	l = len(df_sol.columns)
	pca = PCA(n_components = l)
	df_pca = pca.fit_transform(df_sol)
	variance = pca.explained_variance_ratio_
	return variance[1]

def pcncomp(sol_path):
	df = pd.read_csv(sol_path)
	df_sol = df.iloc[:,3:]
	df_sol.drop(columns = ['miR205','miR9','miR30c','VIM','CDH1','KLF8','TCF3'], inplace = True)
	pca = PCA(0.9).fit(df_sol)
	df_pca = pca.transform(df_sol)
	num = pca.n_components_

	return num 

def cumulative(sol_path,t):
	df = pd.read_csv(sol_path)
	df_sol = df.iloc[:,3:]
	df_sol.drop(columns = ['miR205','miR9','miR30c','VIM','CDH1','KLF8','TCF3'], inplace = True)
	genes = df_sol.columns.to_list()
	pca = PCA(0.9).fit(df_sol)
	df_pca = pca.transform(df_sol)
	num = pca.n_components_
	feat = pca.components_
	e_val = pd.Series(np.sqrt(pca.explained_variance_))
	df2 = pd.DataFrame(data = feat, index = np.arange(num), columns = genes)

	for row in df2.index:
		df2.iloc[row,:] = df2.iloc[row,:]*e_val[row]

	df3 = df2.sum()
	values = df3.to_list()
	#df3.sort_index(axis = 0,inplace = True)
	#df3.columns = [t[:-5]]
	df4 = pd.DataFrame(data = None, index = np.arange(len(genes)), columns = ['Node',t[:-5]])
	df4['Node'] = genes
	df4[t[:-5]] = values
	df4.sort_values(by = "Node", inplace = True)
	return df4

def theta(sol_path,t):
	df = pd.read_csv(sol_path)
	df_sol = df.iloc[:,3:]
	df_sol.drop(columns = ['miR205','miR9','miR30c','VIM','CDH1','KLF8','TCF3'], inplace = True)
	genes = df_sol.columns.to_list()
	pca = PCA(0.96).fit(df_sol)
	df_pca = pca.transform(df_sol)
	n = pca.n_components_
	feat = pca.components_
	e_val = pd.Series(np.sqrt(pca.explained_variance_))
	df2 = pd.DataFrame(data = feat, index = np.arange(n), columns = genes)

	for row in df2.index:
		df2.iloc[row,:] = df2.iloc[row,:]*e_val[row]

	df2['Basis'] = 0
	df2.iloc[0,len(genes)] = 1 
	basis = df2.iloc[:,len(genes)]

	d1 = basis 
	vector = []
	angles = []

	for j in range(len(genes)):
		v = df2.iloc[:,j]
		vector.append(v)

		for k in range(len(vector)):
			angle = (np.arccos((np.dot(vector[k],d1))/((np.linalg.norm(vector[k]))*(np.linalg.norm(d1)))))*(180.0 / math.pi)

		angles.append(angle)

	df3 = pd.DataFrame(data = None, index = np.arange(len(genes)), columns = ['Node',t[:-5]])
	df3['Node'] = genes	
	df3[t[:-5]] = angles
	df3.sort_values(by = "Node", inplace = True)

	return df3

def theta12(sol_path,t):
	df = pd.read_csv(sol_path)
	df_sol = df.iloc[:,3:]
	df_sol.drop(columns = ['CDH1','KLF8','TCF3','VIM','miR205','miR30c','miR9'], inplace = True)
	genes = df_sol.columns.to_list()
	pca = PCA(n_components = 2)
	df_pca = pca.fit_transform(df_sol)
	n = pca.n_components_
	feat = pca.components_ 
	e_val = pd.Series(np.sqrt(pca.explained_variance_))
	df2 = pd.DataFrame(data = feat, index = np.arange(n), columns = genes)

	for row in df2.index:
		df2.iloc[row,:] = df2.iloc[row,:]*e_val[row]

	df2['Basis'] = 0
	df2.iloc[0,len(genes)] = 1
	basis = df2.iloc[:,len(genes)]

	d1 = basis 
	vector = []
	angles = []

	for j in range(len(genes)):
		v = df2.iloc[:,j]
		vector.append(v)

		for k in range(len(vector)):
			angle = (np.arccos((np.dot(vector[k],d1))/((np.linalg.norm(vector[k]))*(np.linalg.norm(d1)))))*(180.0/math.pi)

		angles.append(angle)

	df3 = pd.DataFrame(data = None, index = np.arange(len(genes)), columns = ['Node',t[:-5]])
	df3['Node'] = genes
	df3[t[:-5]] = angles
	df3.sort_values(by = "Node", inplace = True)

	return df3

pool = Pool(70)

PC1v = pool.map(pc1var,[(glob.glob(cwd + "/Results/" + t[:-5] + "/" + str(n) + "/*norm.csv")[0]) for n in range(1,4) for t in tpfl])
PC2v = pool.map(pc2var,[(glob.glob(cwd + "/Results/" + t[:-5] + "/" + str(n) + "/*norm.csv")[0]) for n in range(1,4) for t in tpfl])
PCnum = pool.map(pcncomp,[(glob.glob(cwd + "/Results/" + t[:-5] + "/" + str(n) + "/*norm.csv")[0]) for n in range(1,4) for t in tpfl])
PCcoeff = pool.starmap(cumulative,[(glob.glob(cwd + "/Results/" + t[:-5] + "/" + str(n) + "/*norm.csv")[0],t) for t in tpfl for n in range(1,4)])
PCtheta = pool.starmap(theta,[(glob.glob(cwd + "/Results/" + t[:-5] + "/" + str(n) + "/*norm.csv")[0],t) for t in tpfl for n in range(1,4)])
PCtheta12 = pool.starmap(theta,[(glob.glob(cwd + "/Results/" + t[:-5] + "/" + str(n) + "/*norm.csv")[0],t) for t in tpfl for n in range(1,4)])

pool.close()
pool.join()

var1 = pd.DataFrame(data = None, index = np.arange(len(tpfl)), columns = ['Network','1','2','3'])
var2 = pd.DataFrame(data = None, index = np.arange(len(tpfl)), columns = ['Network','1','2','3'])
nums = pd.DataFrame(data = None, index = np.arange(len(tpfl)), columns = ['Network','1','2','3'])

names = []
for t in tpfl:
	names.append(t[:-5])

var1['Network'] = names
var1['1'] = PC1v[:len(tpfl)]
var1['2'] = PC1v[len(tpfl):(2*len(tpfl))]
var1['3'] = PC1v[2*len(tpfl):]
var1.sort_values(by = "Network", inplace = True)

var2['Network'] = names
var2['1'] = PC2v[:len(tpfl)]
var2['2'] = PC2v[len(tpfl):(2*len(tpfl))]
var2['3'] = PC2v[2*len(tpfl):]
var2.sort_values(by = "Network", inplace = True)

nums['Network'] = names
nums['1'] = PCnum[:len(tpfl)]
nums['2'] = PCnum[len(tpfl):(2*len(tpfl))]
nums['3'] = PCnum[2*len(tpfl):]
nums.sort_values(by = "Network", inplace = True)

Cfin = []
for c in np.arange(0,(3*len(tpfl)),3):
	co = pd.concat(PCcoeff[c:c+3], axis = 0)
	Cfin.append(co)
Cframes = []
i = 1 
for c in Cfin:
	if i == 1:
		Cframes.append(c)
	else:
		Cframes.append(c.iloc[:,1])
	i += 1

coeff = pd.concat(Cframes, axis = 1)
coeff.sort_index(axis = 1, inplace = True)

Afin = []
for a in np.arange(0,(3*len(tpfl)),3):
	ang = pd.concat(PCtheta[a:a+3], axis = 0)
	Afin.append(ang)
Aframes = []
j = 1
for a in Afin:
	if j == 1:
		Aframes.append(a)
	else:
		Aframes.append(a.iloc[:,1])
	j += 1

PCangle = pd.concat(Aframes, axis = 1)
PCangle.sort_index(axis = 1, inplace = True)

Afin12 = []
for a in np.arange(0,(3*len(tpfl)),3):
	ang = pd.concat(PCtheta12[a:a+3], axis = 0)
	Afin12.append(ang)
Aframes12 = []
j = 1
for a in Afin12:
	if j == 1:
		Aframes12.append(a)
	else:
		Aframes12.append(a.iloc[:,1])
	j += 1

PCangle12 = pd.concat(Aframes12, axis = 1)
PCangle12.sort_index(axis = 1, inplace = True)

var1.to_csv(cwd + "/PC1variance.csv", index = False)
var2.to_csv(cwd + "/PC2variance.csv", index = False)
nums.to_csv(cwd + "/Num_PCs.csv", index = False)
coeff.to_csv(cwd + "/PC_coefficients.csv", index = False)
PCangle.to_csv(cwd + "/Theta.csv", index = False)
PCangle12.to_csv(cwd + "/Theta12.csv", index = False)

def PCA_teams(theta):
	df = pd.read_csv(theta)
	#df.drop(columns = "Unnamed: 0", inplace = True)
	df.sort_index(axis = 1, inplace = True)
	cols = df.columns[:-1].to_list() #if network name begins with a letter > "N", cols = df.columns[1:].to_list()
	n = int(len(df.iloc[:,0])/3)
	THETA = pd.DataFrame(data = None, index = np.arange(n), columns = cols)
	for k in cols:
		THETA[k] = (df.loc[0:n-1,k].to_numpy() + df.loc[n:(2*n)-1,k].to_numpy() + df.loc[(2*n):(3*n)-1,k].to_numpy())/3.0
	THETA['Node'] = df.loc[:n-1,'Node'].to_list()
	TS = []
	for col in cols:
		kmeans = KMeans(n_clusters = 2).fit(THETA[col].to_numpy().reshape(-1,1))
		cluster = kmeans.predict(THETA[col].to_numpy().reshape(-1,1))
		t1 = THETA[col][cluster == 0]
		T1 = []
		for i in t1.index:
			T1.append(THETA.loc[i,'Node'])
		t2 = THETA[col][cluster == 1]
		T2 = []
		for i in t2.index:
			T2.append(THETA.loc[i,'Node'])
		team1 = ""
		team2 = ""
		for g in T1:
			team1 = team1 + "," + g
		team1 = team1.replace(",","",1)
		for g in T2:
			team2 = team2 + "," + g
		team2 = team2.replace(",","",1)
		tlog.write(col + " t1 " + str(len(T1)) + " " + team1 + "\n" + col + " t2 " + str(len(T2)) + " " + team2 + "\n\n")

PCA_teams(cwd + "/Theta12.csv")
		
	

