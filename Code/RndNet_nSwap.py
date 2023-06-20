import random
import numpy as np 
import pandas as pd 
import subprocess 
import os 
import glob 

# get the current working directory
cwd = os.getcwd()

# change to the topo directory
tpdir = (cwd + "/TOPO")
os.chdir(tpdir)

results = glob.glob("*.topo")
tpfl = sorted(results) # path of the biological/parent topo file

def Rnd_Net_Range(bio, output,k):
	R = []
	Edge = []
	with open(bio,'r') as f:
		for line in f:
			i = line.strip().split(' ')
			R.append([i[0], i[1], i[2]])

	u = len(R) -1
	count = 0
	for j in range(1,len(R)):
		Edge.append(R[j][0] + "->" + R[j][1] + " " + R[j][2]) 
	for m in range(k):
		n1 = random.randint(1,u)
		n2 = random.randint(1,u)
		if n1 == n2:
			n2 = random.randint(1,u)
		type1 = R[n1][2]
		type2 = R[n2][2]
		target1 = R[n1][1]
		target2 = R[n2][1]
		source1 = R[n1][0]
		source2 = R[n2][0]
		edge1 = source1 + "->" + target1 + " " + type2
		edge2 = source2 + "->" + target2 + " " + type1
		#while edge1 in Edge or edge2 in Edge: #In type shuffling, this loop makes all swaps useful. In target shuffling, it only checks for duplicate edges. 
		#	n1 = random.randint(1,u)
		#	n2 = random.randint(1,u)
		#	if n1 == n2:
		#		n2 = random.randint(1,u)
		#	type1 = R[n1][2]
		#	type2 = R[n2][2]
		#	target1 = R[n1][1]
		#	target2 = R[n2][1]
		#	source1 = R[n1][0]
		#	source2 = R[n2][0]
		#	edge1 = source1 + "->" + target1 + " " + type2
		#	edge2 = source2 + "->" + target2 + " " + type1
		R[n1][2] = type2	
		R[n2][2] = type1
		Edge[n1-1] = source1 + "->" + target1 + " " + type2
		Edge[n2-1] = source2 + "->" + target2 + " " + type1
		if R[n1][2] != R[n2][2]:
			count +=1 

	for r in R:
		output.write(r[0] + ' ' + r[1] + ' ' + r[2] + "\n")

	return count


useful = []
names = []
for t in tpfl:
    bio = t
    for k in range(10,101,10): # range of swaps to cover
	    for n in range(1,51): # number of networks generated for each swap number
		    output = open(t[:3] + format(k, '03d') + "_" + format(n, '02d') + ".topo", 'w')
		    count = Rnd_Net_Range(bio, output, k)
		    useful.append(count)
		    names.append(t[:3] + format(k,'03d') + "_" + format(n,'02d'))

df = pd.DataFrame(data = None, index = np.arange(len(names)), columns = ['Network', 'Useful Swaps'])
df['Network'] = names
df['Useful Swaps'] = useful
df.to_csv(cwd + "/nSwaps.csv", index = False)




