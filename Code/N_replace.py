import pandas as pd 
import numpy as np 
import networkx as nx 
import itertools as itls 
import random  
import os 
import glob
from math import floor, comb
from multiprocessing import Pool 

cwd = os.getcwd()

if not os.path.exists("TOPO"):
    os.makedirs("TOPO")

tp = glob.glob("*.topo")[0]
#Read the topo file as a dataframe
df = pd.read_table(tp, delim_whitespace = True)
df['Type'].replace(2,-1,inplace = True)
#Convert the topo to an adjacency matrix
tpG = nx.from_pandas_edgelist(df, source = "Source", target = "Target", edge_attr = "Type",create_using = nx.DiGraph())
adj = nx.to_pandas_adjacency(tpG, weight = "Type")
N = list(tpG.nodes)
#Number of networks to be generated for each deletion category
n_samples = 100
dt = pd.DataFrame(columns = ['Network','No. Nodes Removed','Nodes Removed'])
names = []
numrmv = []
ndrmv = []
#For loop to remove nodes from the network and replace them with disconnected nodes
for i in range(1,19):
    #number of nodes to remove from the network
    ndrop = i
    #combinations of nodes to remove  
    #m = list(itls.combinations(N,ndrop))
    m = itls.combinations(N,ndrop)
    print("comb gen")
    if ndrop == 1 or ndrop == len(N)-1:
        #ndel = random.choices(m, k=100)
        ndel = random.choices(list(itls.combinations(N,ndrop)),k=100)
    else:
        #ndel = random.sample(m,n_samples)
        ndel = list(itls.islice(m,0,None,floor(comb(len(N),ndrop)/100)))
        #ndel = list(itls.islice(m,n_samples))
        nd = len(ndel)-n_samples
        if nd != 0:
            random.shuffle(ndel)
            del ndel[:nd]
    for j in range(len(ndel)):
        adjd = adj.drop(index = list(ndel[j]), columns = list(ndel[j]))
        ntmp = nx.from_pandas_adjacency(adjd, create_using = nx.DiGraph())
        #replace the removed nodes with an equal number of disconnected nodes
        for r in range(i):
            ntmp.add_edge("R"+format(r+1,'02d'),"R"+format(r+1,'02d'), weight = 1)
        tpf = nx.to_pandas_edgelist(ntmp)
        tpf.columns = ["Source","Target","Type"]
        tpf['Type'] = tpf['Type'].astype(int)
        tpf['Type'].replace(-1,2,inplace=True)
        numrmv.append(i)
        ndrmv.append(ndel[j])
        names.append(tp[:-5].replace("0","").replace("_","")+"_"+format(i,'02d')+"_"+format(j+1,'03d'))
        tpf.to_csv(cwd+"/TOPO/"+tp[:-5].replace("0","").replace("_","")+"_"+format(i,'02d')+"_"+format(j+1,'03d')+".topo", sep = " ", index = False)

dt['Network'] = names
dt['No. Nodes Removed'] = numrmv
dt['Nodes Removed'] = ndrmv 
dt.to_csv(cwd+"/NrmvDF.csv", index = False)

    
