import pandas as pd
import numpy as np
import networkx as nx
import glob
import itertools
import subprocess as subp
import os
from multiprocessing import Pool 

# Create a directory to save all the influence matrices
if not os.path.exists("InfMat"):
    os.makedirs("InfMat")

# Create a list of all the top files in the directory
tpfoldr = "TOPO"
topoli = sorted(glob.glob(tpfoldr+"/*.topo"))
#print(topoli)

# Max path length to consider to generate influence matrix
path_length = 50

# Clean topo file (remove motnotinic.. lines)
def cleanTopo(topoPath):
    with open(topoPath, "r") as tf:
        tfl = tf.readlines()
        if "\t" in tfl[0]:
            sep = "\t"
        else:
            sep = " "
        tpl = [i.strip("\n").split(sep) for i in tfl if "not monotonic" not in i]
        tpl = [i for i in tpl if "not essential" not in i]
    return tpl

# Function to generate influence matrix from topofile 
def infmat(tp):
    # Clean topo file (remove motnotinic.. lines)
    tpl = cleanTopo(tp)
    # Read the topofile as a dataframe
    tdf = pd.DataFrame(tpl[1:], columns=tpl[0])
    tdf["Type"] = tdf["Type"].replace({"2":"-1"})
    # Convert the dataframe into a networkx graph (use of DiGraph constructor is important)
    tpG = nx.from_pandas_edgelist(tdf, source="Source", target="Target", edge_attr="Type", create_using=nx.DiGraph)
    # Convert to pandas dataframe if needed
    adjMat = nx.to_pandas_adjacency(tpG, weight="Type")
    #print(adjMat)
    ## Find the input nodes which have all the rows as zero
    #outptLi = list(adjMat.loc[:, (adjMat==0).all(axis=1)].columns)
    #notoutptLi = list(adjMat.loc[:, ~(adjMat==0).all(axis=1)].columns)
    ## Subset the dataframe to remove the output nodes(??)
    #adjMat = adjMat.loc[notoutptLi, notoutptLi]
    ## Find the output nodes which have all the columns as zero
    #inptLi = list(adjMat.loc[(adjMat==0).all(axis=0)].index)
    #notinptLi = list(adjMat.loc[~(adjMat==0).all(axis=0)].index)
    ## Subset the dataframe to remove the input nodes(??)
    #adjMat = adjMat.loc[notinptLi, notinptLi]
    node_li = adjMat.columns
    # Convert the topo file dataframe to  numpy matrix
    adjMat = adjMat.to_numpy()
    ################
    # MOdified IFL
    # Shallow copy the adjacency matrix
    infMat1 = adjMat.copy()
    maxMat = adjMat.copy()
    # Convert the non-zero elements in the matrix to 1 to get Max Matrix
    maxMat[maxMat != 0] = 1.0
    #path_length = 100
    path_length = 10
    #path_length = nx.diameter(tpG.to_undirected())
    # Take powers of the metix to get the influence matrix and update the Max Martix
    for i in range(2,path_length+1):
        a = np.linalg.matrix_power(adjMat, i).astype(float)
        b = np.linalg.matrix_power(maxMat, i).astype(float)
        infMat1 = infMat1 + np.divide(a, b, out=np.zeros_like(a), where=b!=0)/(i)
    # Normalise by the path length
    infMat1 = infMat1/sum(np.reciprocal(np.arange(1, path_length+1), dtype=float))
    #print(pd.DataFrame(infMat1))
    # Convert the matrix into a dataframe with row and column values
    infMat1df = pd.DataFrame(infMat1, index=node_li, columns=node_li)
    # Write to a csv file
    infMat1df.to_csv("InfMat/"+tp.replace(tpfoldr+"/", "").replace(".topo","")+"_InfMatMew.csv")
    #################
    #################
    # Shallow copy the adjacency matrix
    infMat = adjMat.copy()
    maxMat = adjMat.copy()
    # Convert the non-zero elements in the matrix to 1 to get Max Matrix
    maxMat[maxMat != 0] = 1.0
    path_length = 10
    # Take powers of the metix to get the influence matrix and update the Max Martix
    for i in range(2,path_length+1):
        a = np.linalg.matrix_power(adjMat, i).astype(float)
        b = np.linalg.matrix_power(maxMat, i).astype(float)
        #infMat = infMat + np.divide(a, b, out=np.zeros_like(a), where=b!=0)/(i)
        infMat = infMat + np.divide(a, b, out=np.zeros_like(a), where=b!=0)
    # Normalise by the path length
    infMat = infMat/(path_length)
    #################
    # Convert the matrix into a dataframe with row and column values
    infMatdf = pd.DataFrame(infMat, index=node_li, columns=node_li)
    # Write to a csv file
    infMatdf.to_csv("InfMat/"+tp.replace(tpfoldr+"/", "").replace(".topo","")+"_InfMat.csv")
    
pool = Pool(30)
pool.map(infmat,[tp for tp in topoli])

pool.close()
pool.join()
