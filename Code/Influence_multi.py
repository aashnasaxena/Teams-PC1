import pandas as pd 
import numpy as np 
import subprocess
import glob
import os 
from multiprocessing import Pool 
from scipy.cluster import hierarchy as hi 
import matplotlib.pyplot as plt
import seaborn as sns

cwd = os.getcwd()
tpdir = "/home/csb/Aashna/Range10/TOPO"
tlog = open("influence_teams.txt","w")


os.chdir(tpdir)
results = glob.glob("*.topo")
tpfl = sorted(results)

#Get Influence Matrix:
def Influence_Matrix(cwd,t,n):
        topo = pd.read_table(t, delimiter = ' ')
        n1 = topo['Source'].drop_duplicates()
        n2 = topo['Target'].drop_duplicates()
        N1 = n1.to_list()
        N2 = n2.to_list()
        N = N1 + N2
        Nodes = [*set(N)]

        #NOTE: index = source, colnames = target

        adj = pd.DataFrame(data = 0, index = Nodes, columns = Nodes)
        topo['Type'].replace(2,-1,inplace = True)

        for row in adj.index:
                for i in topo.index:
                        if topo.iloc[i][0] == row:
                                j = topo.iloc[i][1]
                                adj.loc[row, j] = topo.iloc[i][2]

        #print(adj)

        N = float(np.count_nonzero(adj))
        inf = adj.copy()
        M_max = adj.copy()
        M_max[M_max != 0] = 1.0
        for i in range(2,10):
                a = np.linalg.matrix_power(adj, i).astype(float)
                b = np.linalg.matrix_power(M_max, i).astype(float)
                inf = inf + np.divide(a, b, out=np.zeros_like(a), where=b!=0)

        inf = inf/10

        influence = pd.DataFrame(data=inf,index=Nodes,columns=Nodes)
        #nodes = influence.columns
        #d = hi.distance.pdist(influence)
        #L = hi.linkage(d, method = "complete")
        #clust = hi.cut_tree(L, n_clusters =2)
        #cluster = np.transpose(clust)
        #t1 = nodes[cluster[0] == 0]
        #t2 = nodes[cluster[0] == 1]
        #tot = t1.append(t2)
        #cl_influence = influence.loc[tot,:].T.loc[tot,:]
        influence.to_csv(cwd + "/Results/" + t[:-5] + "/" + str(n) + "/" + t[:-5] + "_influence.csv")

#Get influence matrix figure:
def Influence_fig(influence_path,cwd,t):
        influence = pd.read_csv(influence_path)
        influence.set_index("Unnamed: 0", inplace = True)
        d = hi.distance.pdist(influence)
        L = hi.linkage(d, method = 'complete')
        sns.set(rc={'figure.figsize':(16,16)})
        sns.set_context("paper", rc={"font.weight":'bold',"legend.fontsize":8,"legend.title_fontsize":10,"font.size":8,"axes.titlesize":8,"axes.labelsize":8,"xtick.labelsize":10,"ytick.labelsize":10})
        ax = sns.clustermap(data = influence, method = 'complete', annot = True, row_linkage = L, col_linkage = L, cmap = "coolwarm")
        plt.savefig(cwd + "/Results/" + tpfl[t][:-5] + "/" + tpfl[t][:-5] + "_influence.png", dpi = 400, pad_inches = 0)
        plt.clf()

#Get Team Strength:
def Team_Strength(influence_path,tpfl,t):
        influence = pd.read_csv(influence_path)
        influence.set_index("Unnamed: 0", inplace = True)
        influence.drop(index = ['miR9','miR30c','miR205','VIM','CDH1','KLF8','TCF3'], columns = ['miR9','miR30c','miR205','VIM','CDH1','KLF8','TCF3'], inplace = True)
        nodes = influence.columns
        d = hi.distance.pdist(influence)
        L = hi.linkage(d, method = 'complete')
        clust = hi.cut_tree(L, n_clusters = 2)
        cluster = np.transpose(clust)
        t1 = nodes[cluster[0] == 0]
        t2 = nodes[cluster[0] == 1]
        tot = t1.append(t2)
        df_clust = influence.loc[tot,:].T.loc[tot,:]

        team1 = ""
        team2 = ""
        for g in t1:
                team1 = team1 + "," + g
        team1 = team1.replace(",","",1)
        for g in t2:
                team2 = team2 + "," + g
        team2 = team2.replace(",","",1)
        tlog.write(tpfl[t][:-5] + " t1 " + str(len(t1)) + " " + team1 + "\n" + tpfl[t][:-5] + " t2 " + str(len(t2)) + " " + team2 + "\n\n")

        df_t11 = df_clust.loc[t1,t1]
        df_t22 = df_clust.loc[t2,t2]
        df_t12 = df_clust.loc[t1,t2]
        df_t21 = df_clust.loc[t2,t1]

        num_t11 = df_t11.to_numpy()
        num_t22 = df_t22.to_numpy()
        num_t12 = df_t12.to_numpy()
        num_t21 = df_t21.to_numpy()

        t11 = abs(np.sum(num_t11, axis = None))/(len(num_t11)*len(num_t11[0]))
        t22 = abs(np.sum(num_t22, axis = None))/(len(num_t22)*len(num_t22[0]))
        t12 = abs(np.sum(num_t12, axis = None))/(len(num_t12)*len(num_t12[0]))
        t21 = abs(np.sum(num_t21, axis = None))/(len(num_t21)*len(num_t21[0]))

        ts = (t11 + t22 + t12 + t21)/4 

        return ts 

pool = Pool(70)
pool.starmap(Influence_Matrix, [(cwd,t,n) for t in tpfl for n in range(1,4)])

pool.close()
pool.join()


df_inf = pd.DataFrame(data = None, index = np.arange(len(tpfl)), columns = ['Network',1,2,3])
names = []

for t in range(len(tpfl)):
        names.append(tpfl[t][:-5])
        for n in range(1,4):
                influence_path = glob.glob(cwd + "/Results/" + tpfl[t][:-5] + "/" + str(n) + "/*influence.csv")[0]
                ts = Team_Strength(influence_path,tpfl,t)
                df_inf.iloc[t,n] = ts
                

df_inf['Network']  = names
df_inf.sort_values(by = "Network", inplace = True)

df_inf.to_csv(cwd + "/influence_ts.csv", index = False)
