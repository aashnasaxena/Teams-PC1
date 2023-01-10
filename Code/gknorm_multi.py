import pandas as pd 
import numpy as np 
import glob
import os 
import math 
from multiprocessing import Pool 

cwd = os.getcwd()

tpdir = "/home/csb/Aashna/Range10/TOPO"
os.chdir(tpdir)
results = glob.glob("*.topo")
tpfl = sorted(results)

def rmvsol(cwd,t,n):
	for file in glob.glob(cwd + "/Results/" + t[:-5] + "/" + str(n) + "/*solution_*.dat"):
		os.remove(file)

def gknorm(cwd, t, n, prs_path, cfg_path, parameter_path, sol_path):
	#Naming the solution and parameter files:

	df_prs = pd.read_table(prs_path)
	names_cfg = ['1','2','3','4']
	df_cfg = pd.read_table(cfg_path, index_col = '1', names = names_cfg)
	a = df_cfg.loc['NumberOfGenes','2']
	gen = int(a)

	names_sol = ['model_index','ss_no','percentage']

	for i in range(gen):
		x = df_prs.iloc[i,0].replace("Prod_of_","")
		names_sol.append(x)

	solution = pd.read_table(sol_path, names = names_sol)
	solution.to_csv(cwd + "/Results/" + t[:-5] + "/" + str(n) + "/solution.csv", index = False)

	names_para = ['model_index','ss_no']
	para_names = df_prs.iloc[:,0].to_list()
	names = names_para + para_names

	para = pd.read_table(parameter_path, header = None, names = names)
	
	#g/k normalization:
	
	antilog = solution.copy()
	for col in solution.columns[3:]:
		antilog[col] = 2**solution[col]

	merged = antilog.merge(para, on = "model_index")
	merged.drop(columns = "ss_no_y", inplace = True)
	lis = merged.columns.to_list()
	merged.columns = np.arange(len(lis))

	norm = pd.DataFrame(data = None, index = np.arange(len(solution.index)), columns = np.arange(25))
	norm[0] = solution["model_index"]
	norm[1] = solution["ss_no"]
	norm[2] = solution["percentage"]

	for i in range(3,gen+3):
		norm[i] = np.log2(merged[i] / (merged[i+gen]/merged[i+(2*gen)]))

	norm.columns = solution.columns

	norm.to_csv(cwd + "/Results/" + t[:-5] + "/" + str(n) + "/norm.csv", index = False)

#for t in tpfl:
#	for n in range(1,4):
#		prs_path = glob.glob(cwd + "/Results/" + t[:-5] + "/" + str(n) + "/*.prs")[0]
#		parameter_path = glob.glob(cwd + "/Results/" + t[:-5] + "/" + str(n) + "/*parameters.dat")[0]
#		sol_path = glob.glob(cwd + "/Results/" + t[:-5] + "/" + str(n) + "/*solution.csv")[0]
#		gknorm(cwd,t,prs_path,parameter_path,sol_path) 

with Pool(60) as pool:
	pool.starmap(gknorm, [(cwd,t,n,glob.glob(cwd + "/Results/" + t[:-5] + "/" + str(n) + "/*.prs")[0], glob.glob(cwd + "/Results/" + t[:-5] + "/" + str(n) + "/*.cfg")[0], glob.glob(cwd + "/Results/" + t[:-5] + "/" + str(n) + "/*parameters.dat")[0],glob.glob(cwd + "/Results/" + t[:-5] + "/" + str(n) + "/*solution.dat")[0]) for t in tpfl for n in range(1,4)])

