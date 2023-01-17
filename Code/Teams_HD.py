import pandas as pd 
import numpy as np 
import os
import glob

cwd = os.getcwd()

def Team_Hamming(inft, pcat):
	inft.drop_duplicates(inplace = True)
	inft.reset_index(inplace = True)
	inft.drop(columns = "index", inplace = True)

	pcat.drop_duplicates(inplace = True)
	pcat.reset_index(inplace = True)
	pcat.drop(columns = "index", inplace = True)

	genes = inft.iloc[0,3] + "," + inft.iloc[1,3]
	Nodes = genes.split(sep = ",")

	final = pd.DataFrame(data = None, index = np.arange(int((len(inft.index))/2)), columns = ['Network','Inf T1','Inf T2','PC12 T1','PC12 T2','T11', 'T12','T21','T22','Hamming Distance','Inv. Hamming Distance'])
	nets = inft["Network"].drop_duplicates().to_list()
	final['Network'] = nets

	T11 = []
	T12 = []
	T21 = []
	T22 = []
	HD = []
	IHD = []
	infT1 = []
	infT2 = []
	PCAT1 = []
	PCAT2 = []

	for row in np.arange(0,(2*len(nets)),2):
		t11 = []
		t22 = []
		t12 = []
		t21 = []
		inft1 = ""
		inft2 = ""
		pcat1 = ""
		pcat2 = ""
		for node in Nodes:
			if node in inft.iloc[row,3]:
				inft1 = inft1 + ";" + node 
			else:
				inft2 = inft2 + ";" + node
		inft1 = inft1.replace(";","",1)
		inft2 = inft2.replace(";","",1)

		for node in Nodes:
			if node in pcat.iloc[0,3]:
				pcat1 = pcat1 + ";" + node 
			else:
				pcat2 = pcat2 + ";" + node 
		pcat1 = pcat1.replace(";","",1)
		pcat2 = pcat2.replace(";","",1)

		for node in Nodes:
			if node in inft.iloc[row,3] and node in pcat.iloc[row,3]:
				t11.append(node)
			elif node in inft.iloc[row+1,3] and node in pcat.iloc[row+1,3]:
				t22.append(node)
			elif node in inft.iloc[row,3] and node in pcat.iloc[row+1,3]:
				t12.append(node)
			else:
				t21.append(node)
		hd = (len(t12) + len(t21))/len(Nodes)
		ihd = (len(t11) + len(t22))/len(Nodes)
		HD.append(hd)
		IHD.append(ihd)
		T11.append(t11)
		T22.append(t22)
		T12.append(t12)
		T21.append(t21)
		infT1.append(inft1)
		infT2.append(inft2)
		PCAT1.append(pcat1)
		PCAT2.append(pcat2)

	final['T11'] = T11
	final['T22'] = T22
	final['T12'] = T12
	final['T21'] = T21
	final['Hamming Distance'] = HD
	final['Inv. Hamming Distance'] = IHD
	final['Inf T1'] = infT1
	final['Inf T2'] = infT2
	final['PC12 T1'] = PCAT1
	final['PC12 T2'] = PCAT2

	final.to_csv(cwd + "/Teams_HD.csv", index = False)


inft = pd.read_table(cwd + "/influence_teams.txt", delimiter = " ", names = ['Network', 'Type', 'No.', 'Team'])
pcat = pd.read_table(cwd + "/PCA_teams.txt", delimiter = " ", names = ['Network', 'Type','No.','Team'])
Team_Hamming(inft, pcat)











