
import pandas as pd 
import numpy as np 
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt 
import seaborn as sns 
import os
import glob

# Get the current working directory
cwd = os.getcwd()

# list of g/k normalized files:
normLi = sorted(glob.glob("*norm.csv"))

for fl in normLi:
    df = pd.read_csv(fl)
    df_sol = df.iloc[:,3:]
    pca = PCA(n_components = 2)
    df_pca = pca.fit_transform(df_sol)
    variance = pca.explained_variance_ratio_
    v1 = variance[0]*100
    v2 = variance[1]*100
    final = pd.DataFrame(data = None, index = np.arange(len(df_pca)), columns = ['PC1','PC2','Team'])
    x = []
    y = []

    for i in range(len(df_pca)):
        x.append(df_pca[i][0])
        y.append(df_pca[i][1])

    final['PC1'] = x
    final['PC2'] = y

    kmeans = KMeans(n_clusters = 2).fit(df_pca)
    final['Team'] = kmeans.predict(df_pca)

    sns.set(rc={"figure.figsize":(10.5,8)})
    sns.set_style("ticks",rc={"axes.facecolor":"white","axes.edgecolor":"black","font.family":'sans-serif',"font.sans-serif":'Arial',"font.size":18})
    sns.set_context({"font.weight":"normal","font.size":18,"font.style":"normal","axes.labelsize":18,"axes.labelpad":8,"axes.labelweight":"bold","xtick.labelsize":18,"ytick.labelsize":18,"legend.edgecolor":"white","legend.facecolor":"white","legend.title_fontsize":18,"legend.fontsize":18})
    #sns.scatterplot(data = final, x = "PC1", y = "PC2", hue = "Team", palette = ['#509e90', '#1d6c8a'], s = 3)
    #['#0b559f','#aa1016']
    #plt.legend(bbox_to_anchor=(1.02,1), loc = 2, borderaxespad = 0)
    sns.scatterplot(data = final, x = "PC1", y = "PC2", s=3)
    plt.xlabel("PC1 "+"("+str(round(v1,1))+"% var.)",{"fontsize":24,"fontweight":"bold","fontfamily":"sans-serif"})
    plt.ylabel("PC2 "+"("+str(round(v2,1))+"% var.)",{"fontsize":24,"fontweight":"bold","fontfamily":"sans-serif"})
    plt.tight_layout()
    plt.savefig(cwd+"/FIGS/"+fl.replace("_norm.csv","_pca.png"), dpi = 400, pad_inches = 0)
    plt.clf()
