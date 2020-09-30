import os
from joblib import Parallel, delayed
 

files = os.listdir("./")
programs = []

for i in files:
	if ".lp" in i:
		programs.append(i)

Parallel(n_jobs=32)(delayed(os.system)("gurobi_cl ResultFile=" + i[:-3] + ".sol " + i) for i in programs)
        
#for i in programs:
#	os.system("gurobi_cl ResultFile=" + i[:-3] + ".sol " + i)

