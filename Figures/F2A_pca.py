import pandas as pd 
import numpy as np 
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt 
import seaborn as sns 

df = pd.read_csv(sol)
df_sol = df.iloc[:,3:]
pca = PCA(n_components = 2)
df_pca = pca.fit_transform(df_sol)

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

sns.scatterplot(data = final, x = "PC1", y = "PC2", hue = "Team", palette = ["#a50f15","#3690c0"], s = 3)
#palette = "crest"
plt.legend(bbox_to_anchor=(1.02,1), loc = 2, borderaxespad = 0)
plt.tight_layout()
plt.savefig("F2A_pca.png", dpi = 400, pad_inches = 0)
plt.savefig("F2A_pca.svg", dpi = 400)
