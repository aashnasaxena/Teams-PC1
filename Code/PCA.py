import pandas as pd
import numpy as np
import os
import glob 
from multiprocessing import Pool
from sklearn.decomposition import PCA
import math 

cwd = os.getcwd()

tpdir = (cwd + "/TOPO")

os.chdir(tpdir)
tpfl = sorted(glob.glob("*.topo"))

def pc1var(sol_path):
	df = pd.read_csv(sol_path)
	df_sol = df.iloc[:,3:]
	pca = PCA(n_components = len(df_sol.columns))
	df_pca = pca.fit_transform(df_sol)
	variance = pca.explained_variance_ratio_
	return variance[0]

def pc2var(sol_path):
	df = pd.read_csv(sol_path)
	df_sol = df.iloc[:,3:]
	pca = PCA(n_components = len(df_sol.columns))
	df_pca = pca.fit_transform(df_sol)
	variance = pca.explained_variance_ratio_
	return variance[1]

def pcncomp(sol_path):
	df = pd.read_csv(sol_path)
	df_sol = df.iloc[:,3:]
	#with open(cwd+"/TeamsCalc/"+sol_path.split("/")[-3]+".teams",'r') as f:
	#	a = f.readlines()[2:]
	#P = []
	#for p in a:
	#	P = P + p.strip().split(":")[1].split(",")
	#P = [*set(P)]
	#if '' in P:
	#	P.remove('')
	#df_sol.drop(columns = P, inplace = True)
	pca = PCA(0.9).fit(df_sol)
	df_pca = pca.transform(df_sol)
	num = pca.n_components_

	return num 

def cumulative(sol_path,t):
	df = pd.read_csv(sol_path)
	df_sol = df.iloc[:,3:]
	with open(cwd+"/TeamsCalc/"+t[:-5]+".teams",'r') as f:
		a = f.readlines()[2:]
	P = []
	for p in a:
		P = P + p.strip().split(":")[1].split(",")
	P = [*set(P)]
	if '' in P:
		P.remove('')
	df_sol.drop(columns = P, inplace = True)
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
	df4 = pd.DataFrame(data = None, index = np.arange(len(genes)), columns = ['Node',t[:-5]])
	df4['Node'] = genes
	df4[t[:-5]] = values
	df4.sort_values(by = "Node", ignore_index = True, inplace = True)
	return df4

def theta(sol_path,t):
	df = pd.read_csv(sol_path)
	df_sol = df.iloc[:,3:]
	with open(cwd+"/TeamsCalc/"+t[:-5]+".teams",'r') as f:
		a = f.readlines()[2:]
	P = []
	for p in a:
		P = P + p.strip().split(":")[1].split(",")
	P = [*set(P)]
	if '' in P:
		P.remove('')
	df_sol.drop(columns = P, inplace = True)

	genes = df_sol.columns.to_list()
	pca = PCA(0.9).fit(df_sol)
	df_pca = pca.transform(df_sol)	
	n = pca.n_components_
	feat = pca.components_
	e_val = pd.Series(np.sqrt(pca.explained_variance_))
	df2 = pd.DataFrame(data = feat, index = np.arange(n), columns = genes)

	for row in df2.index:
		df2.iloc[row,:] = df2.iloc[row,:]*e_val[row]

	df2['Basis'] = 0
	df2.loc[0,'Basis'] = 1
	basis = df2.loc[:,'Basis']

	d1 = basis 
	vector = []
	angles = []

	for j in range(len(df2.columns)-1):
		v = df2.iloc[:,j]
		vector.append(v)

		for k in range(len(vector)):
			angle = (np.arccos((np.dot(vector[k],d1))/((np.linalg.norm(vector[k]))*(np.linalg.norm(d1)))))*(180.0 / math.pi)

		angles.append(angle)

	df3 = pd.DataFrame(data = None, index = np.arange(len(df2.columns)-1), columns = ['Node',t[:-5]])
	df3['Node'] = genes	
	df3[t[:-5]] = angles
	df3.sort_values(by = "Node", ignore_index = True, inplace = True)

	return df3

pool = Pool(48)

PC1v = pool.map(pc1var,[(glob.glob(cwd + "/Results/" + t[:-5] + "/*norm.csv")[0]) for t in tpfl])
PC2v = pool.map(pc2var,[(glob.glob(cwd + "/Results/" + t[:-5] + "/*norm.csv")[0]) for t in tpfl])
PCnum = pool.map(pcncomp,[(glob.glob(cwd + "/Results/" + t[:-5] + "/"+"/*norm.csv")[0]) for t in tpfl])
#PCcoeff = pool.starmap(cumulative,[(glob.glob(cwd + "/Results/" + t[:-5] + "/" + str(n) + "/*norm.csv")[0],t) for t in tpfl for n in range(1,4)])
#PCtheta = pool.starmap(theta,[(glob.glob(cwd + "/Results/" + t[:-5] + "/" + str(n) + "/*norm.csv")[0],t) for t in tpfl for n in range(1,4)])

pool.close()
pool.join()

var1 = pd.DataFrame(data = None, index = np.arange(len(tpfl)), columns = ['Network','PC1 Variance'])
var2 = pd.DataFrame(data = None, index = np.arange(len(tpfl)), columns = ['Network','PC2 Variance'])
nums = pd.DataFrame(data = None, index = np.arange(len(tpfl)), columns = ['Network','Number of PCs'])

names = []
for t in tpfl:
	names.append(t[:-5])

var1['Network'] = names
#var1['1'] = PC1v[:len(tpfl)]
#var1['2'] = PC1v[len(tpfl):(2*len(tpfl))]
#var1['3'] = PC1v[2*len(tpfl):]
var1['PC1 Variance'] = PC1v
var1.sort_values(by = "Network", inplace = True)

var2['Network'] = names
#var2['1'] = PC2v[:len(tpfl)]
#var2['2'] = PC2v[len(tpfl):(2*len(tpfl))]
#var2['3'] = PC2v[2*len(tpfl):]
var2['PC2 Variance'] = PC2v
var2.sort_values(by = "Network", inplace = True)

nums['Network'] = names
#nums['1'] = PCnum[:len(tpfl)]
#nums['2'] = PCnum[len(tpfl):(2*len(tpfl))]
#nums['3'] = PCnum[2*len(tpfl):]
nums["Number of PCs"] = PCnum
#nums.sort_values(by = "Network", inplace = True)

#Cfin = []
#for c in np.arange(0,(3*len(tpfl)),3):
#	co = pd.concat(PCcoeff[c:c+3], ignore_index = True, axis = 0)
#	Cfin.append(co)

#Cframes = []
#i = 1
#for c in Cfin:
#	if i == 1:
#		Cframes.append(c)
#	else:
#		Cframes.append(c.iloc[:,1])
#	i += 1
#print(Cframes)
#coeff = pd.concat(Cframes, axis = 1)
#coeff.sort_index(axis = 1, inplace = True)

#Afin = []
#for a in np.arange(0,(3*len(tpfl)),3):
#	ang = pd.concat(PCtheta[a:a+3], ignore_index = True, axis = 0)
#	Afin.append(ang)

#Aframes = []
#j = 1
#for a in Afin:
#	if j == 1:
#		Aframes.append(a)
#	else:
#		Aframes.append(a.iloc[:,1])
#	j += 1

#PCangle = pd.concat(Aframes, axis = 1)
#PCangle.sort_index(axis = 1, inplace = True)

var1.to_csv(cwd + "/PC1variance.csv", index = False)
var2.to_csv(cwd + "/PC2variance.csv", index = False)
nums.to_csv(cwd + "/Num_PCs.csv", index = False)
#coeff.to_csv(cwd + "/PC_coefficients.csv", index = False)
#PCangle.to_csv(cwd + "/Theta.csv", index = False)

def PCA_teams(coeff_path):
	df = pd.read_csv(coeff_path)
	df.sort_index(axis = 1, inplace = True)
	numnds = int(len(df.index)/3)
	cols = df.columns.to_list()
	cols.remove("Node")
	COEFF = pd.DataFrame(data = None, index = np.arange(numnds), columns = np.arange(len(cols)))
	for k in range(len(cols)):
		COEFF[k] = (df.iloc[0:numnds,k] + df.iloc[numnds:2*numnds,k] + df.iloc[2*numnds:3*numnds,k])/3.0
	COEFF.columns = cols 
	COEFF['Node'] = df['Node']
	TS = []
	for col in cols:
		COEFF.sort_values(by = [col], inplace = True)
		for row in COEFF.index:
			if COEFF.loc[row,col] >= 0:
				t = row
				break
		team1 = COEFF.loc[:t,col]
		team1.drop(index = t, inplace = True)
		team2 = COEFF.loc[t:,col]

		t1 = abs(team1.sum())
		t2 = abs(team2.sum())
		t_av = (t1+t2)/2.0
		TS.append(t_av)

	PCT = pd.DataFrame(data = None, index = np.arange(len(cols)), columns = ['Network', 'PCA_ts'])
	PCT['Network'] = cols 
	PCT['PCA_ts'] = TS
	PCT.to_csv(cwd + "/PCA_ts.csv", index = False)

#coeff_path = glob.glob(cwd + "/*PC_coefficients.csv")[0]
#PCA_teams(coeff_path)

