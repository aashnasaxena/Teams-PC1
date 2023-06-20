import pandas as pd 
import numpy as np 
import os
import subprocess
import glob

cwd = os.getcwd()

rlog = open("racipe_log.txt", "w")

tpdir = cwd+"/TOPO"

os.chdir(tpdir)
results = (glob.glob("*.topo"))
tpfl = sorted(results)


rcpli = []
numC = 70

os.chdir(cwd)
subprocess.run("mkdir Results", shell=True)
os.chdir(cwd + "/Results")

for t in tpfl:
	subprocess.run("mkdir " + t[:-5], shell=True)
	repeats = "mkdir " + t[:-5] + "/1 " + t[:-5] + "/2 " + t[:-5] + "/3"
	subprocess.run(repeats,shell = True)
	for rep in range(1,4):
		tfc = "cp " + tpdir + "/" + t + " " + cwd + "/Results/" + t[:-5] + "/" + str(rep) + "/"
		subprocess.call(tfc, shell=True)
		rcpli.append("RACIPE " + cwd + "/Results/" + t[:-5] + "/" + str(rep) + "/*.topo -num_paras 10000 -num_ode 100 -threads " + str(numC))
	
for rr in rcpli:
    rnRCP = rr + " & wait"
    print(rnRCP)
    subprocess.run(rnRCP, shell=True)
    rlog.write(rnRCP + "\n")

   
