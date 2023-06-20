import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import glob
import os, pathlib
from multiprocessing import Pool 
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster import hierarchy as hi
from itertools import compress
import networkx as nx 

# List out all the influence matrix csv files
#oinfLi = sorted(glob.glob("InfMat/*InfMat.csv"))
#minfLi = sorted(glob.glob("InfMat/*InfMatMew.csv"))
topoli = sorted(glob.glob("TOPO/*.topo"))
# Create a plots directory to save all the plots
if not os.path.exists("TeamsCalc"):
    os.makedirs("TeamsCalc")

# Function to calculate team strengths from influence matrices 
def team_strength(im,tp):
    # Read the influence matrix as a dataframe
    idf = pd.read_csv(im, index_col = 0)
    # Read the topo file as a dataframe and convert it to an adjacency matrix 
    topo = pd.read_table(tp, delim_whitespace = True)
    topo['Type'].replace(2,-1, inplace = True)
    tpG = nx.from_pandas_edgelist(topo,source = "Source", target = "Target", edge_attr = "Type", create_using = nx.DiGraph())
    adj = nx.to_pandas_adjacency(tpG, weight = "Type")
    # A white loop to iteratively find out all the peripheral nodes 
    # Initialise no in and out 
    NoInOut = False
    # Intialise the counter to track which level of peripheral node
    perLvl = 1
    # Create a list to store the peripheral nodes by level
    tfli = []
    # Create a list of all peripheral nodes to drop from the adjacency matrix
    peri = []
    # Run the loop until all the peripheral nodes are found
    while NoInOut == False:
        perLvl = perLvl + 1 
        # Find the input nodes which have all the rows as zero
        outptLi = adj.loc[(adj==0).all(axis=1)].index
        notoutptLi = adj.loc[~(adj==0).all(axis=1)].index
        # Subset the dataframe to remove the input node (??)
        adj = adj.loc[notoutptLi, notoutptLi]
        # Find the output nodes which have all the columns as zero
        inptLi = adj.loc[:, (adj==0).all(axis=0)].columns
        notinptLi = adj.loc[:, ~(adj==0).all(axis=0)].columns
        # Subset the dataframe to remove the output node (??)
        adj = adj.loc[notinptLi, notinptLi]
        # Check if there are not more peripheral nodes
        if not (list(inptLi) + list(outptLi)):
            NoInOut = True
            break
        # Append to the input and output nodes line list
        tfli.append("Input_"+str(perLvl)+":"+ ",".join(sorted(inptLi))+"\n")
        tfli.append("Output_"+str(perLvl)+":"+ ",".join(sorted(outptLi))+"\n")
        peri += sorted(inptLi) + sorted(outptLi)
    # Remove peripheral nodes from the influence matrix 
    idf.drop(index = peri, columns = peri, inplace = True)
    # Make a deep copy of idf to get the teams
    CpIdf = idf.copy()
    # Convert all -ve to -1 and +ve to 1
    CpIdf = CpIdf.where(CpIdf >= 0, -1)
    CpIdf = CpIdf.where(CpIdf <= 0, 1)
    # If the dataframe is filled with zeroes give all column names to one team (to avoid error in hclust) 
    if (CpIdf==0).sum().sum() == len(CpIdf.index)*len(CpIdf.columns):
        tmli = [sorted(list(idf.columns)), []]
    else:
        # Set the clustering parameters
        hclust =  AgglomerativeClustering(n_clusters=2, affinity='euclidean', linkage='ward')
        # Find the clusters/teams and get their labels as an array of [1,0...]
        labels = hclust.fit_predict(idf)
        ## New Clustering Method
        #labels = hi.cut_tree(hi.linkage(hi.distance.pdist(idf), method="complete"), n_clusters=2)
        # Find the team nodes names based on (0, 1) labels
        tmli = [sorted([nd for nd, sl in zip(idf.columns, labels) if sl == 1]), sorted([nd for nd, sl in zip(idf.columns, labels) if not sl == 1])]
        # Sort tmli to get consistent order of t1 an t2 everytime
        tmli.sort()
        ## Set the clustering parameters
        #hclust =  AgglomerativeClustering(n_clusters=2, affinity='euclidean', linkage='ward')
        ## Find the clusters/teams and get their labels as an array of [1,0...]
        #labels = hclust.fit_predict(CpIdf)
        ### New Clustering Method
        ##labels = hi.cut_tree(hi.linkage(hi.distance.pdist(idf), method="complete"), n_clusters=2)
        ## Find the team nodes names based on (0, 1) labels
        #t1 = sorted(list(compress(CpIdf.columns, labels)))
        #t2 = sorted([i for i in CpIdf.columns if i not in t1])
        
    # Write the teams identified into a .teams file
    os.chdir("TeamsCalc")
    with open(im.replace("InfMat/","").replace("_InfMat","").replace(".csv","")+".teams", "w") as tmf:
        tmf.write("T1:"+ ",".join(tmli[0])+"\n")
        tmf.write("T2:"+ ",".join(tmli[1])+"\n")
        for l in tfli:
            tmf.write(l)
    os.chdir("..")
    # Calcualtion of team strengths
    tsli = []
    for i in tmli:
        for j in tmli:
            tdf = idf.loc[i, j]
            gs = np.nanmean(tdf)
            # Convert NaN to 0 for cases with only one team
            if np.isnan(gs):
                gs = 0.0
            tsli.append(gs)
    ts = np.mean(np.absolute(tsli))
    T = [im.replace("InfMat/","").replace("_InfMat","").replace(".csv","")]+tsli+[ts]
    
    return T

pool = Pool(30)
T = pool.starmap(team_strength,[(glob.glob("InfMat/"+tp.replace("TOPO/","").replace(".topo","")+"_InfMat.csv")[0],tp) for tp in topoli])
TM = pool.starmap(team_strength,[(glob.glob("InfMat/"+tp.replace("TOPO/","").replace(".topo","")+"_InfMatMew.csv")[0],tp) for tp in topoli])

pool.close()
pool.join()

tsdf = pd.DataFrame(data = T, index = np.arange(len(topoli)), columns=["Network", "TS11", "TS12", "TS21", "TS22", "TS"])
mtsdf = pd.DataFrame(data = TM, index = np.arange(len(topoli)), columns = ["Network","TS11_Mew","TS12_Mew","TS21_Mew","TS22_Mew","TS_Mew"])
# Print the dataframe
os.chdir("TeamsCalc")
tsdf.to_csv("CompiledOldTS.csv", index = False)
print(tsdf)
# Save the new TS values in a df
mtsdf.to_csv("CompiledMewTS.csv", index = False)
print(mtsdf)
# Merge the two df's along rows 
mergedf = pd.concat([tsdf, mtsdf.iloc[:,1:]], axis=1)
# Save the compiled TS dataframe
mergedf.to_csv("CompiledTS.csv",index=False)
print(mergedf)
os.chdir("..")
