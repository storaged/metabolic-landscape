import os
from joblib import Parallel, delayed

SOL_DIR = "./example-datasets/solutions/"
N_JOBS = 2

files = os.listdir(SOL_DIR)
programs = []
print files
for i in files:
	if ".lp" in i:
		programs.append(i)

#Parallel(n_jobs=N_JOBS)(delayed(os.system)("gurobi_cl ResultFile=" + i[:-3] + ".sol " + i) for i in programs)
        
print programs 
for i in programs:
        out_f = SOL_DIR + i[:-3] + ".sol"
        in_f = SOL_DIR + i
	os.system("gurobi_cl ResultFile=" + out_f + " " + in_f)

