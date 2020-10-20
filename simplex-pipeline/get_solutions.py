#!/usr/bin/env python

import os
import matplotlib.pyplot as plt
from define import get_reaction_list_irr, get_patients_dict, istrue
import numpy as np
import sys

def invert(n):
	if n == 0:
		return 1
	elif n == 1:
		return 0
	else:
		return None

def get_solution_files(directory):
	# get directory to .sol files
	files = os.listdir(directory)
	solutions = []

	for i in files:
	    if ".sol" in i:
		solutions.append(i)
	return solutions
	
def get_flow_dict(model_file, patient, genes_matrix, sol_directory):
	r = get_reaction_list_irr(model_file)#"../data/model_transform.sfba")
	patients_dict = get_patients_dict(genes_matrix)
        
        activity_dict = istrue(model_file, patient, patients_dict[patient])		

        reversal_dict = {}	
	for i in r[2]:
		activity = activity_dict[i]
		if activity == 1:
			val = True
		elif activity == 0:
			val = False
		elif activity == "NA":
			val = None
		for j in r[2][i]:
			reversal_dict[j] = val
		
	solution_dict = {}
	f = open(sol_directory + "/program" + patient + ".sol")
	f.readline()
	for i in f:
		if i[0] == "i":
			l = i.strip().split(" ")
			if reversal_dict[l[0]] == True:
				solution_dict[l[0]] = int(round(float(l[1])))
			elif reversal_dict[l[0]] == False:
				solution_dict[l[0]] = int(invert(round(float(l[1]))))
			elif reversal_dict[l[0]] == None:
				solution_dict[l[0]] = None
		elif "Objective" in i:
                        continue
                else:
                        print "Done."
			break
	return solution_dict

def calculate_patient(s, model_file, table_file, sol_directory, r):
        import csv
        patient = s[7:-4]
	flux_row = [patient]
	flow_dict = get_flow_dict(model_file, patient, table_file, sol_directory)
        #print flow_dict.keys()
        for i in r[1]:
	        if i in flow_dict:
		        flux_row.append(flow_dict[i])
		else:
			flux_row.append("NA")
	#matrix.append(flux_row)
        #print flux_row
        #print sol_directory+"tmp"+patient+".tmp"
        with open(sol_directory+"/tmp"+patient+".tmp", "w") as f:
                #f.write("\t".join([str(j) for j in flux_row]) + "\n")
                writer = csv.writer(f)
                writer.writerow(flux_row)
                #f.write(str(flux_row))


def create_matrix(model_file, table_file, sol_directory):	
	from joblib import Parallel, delayed
        import csv

        N_JOBS = 2
        
        r = get_reaction_list_irr(model_file)
	patients = get_patients_dict(table_file)
	matrix = []
	first_row = ["patient"]
	first_row.extend(r[1])
	matrix.append(first_row)

	solutions = get_solution_files(sol_directory)

        #Parallel(n_jobs=N_JOBS)(delayed(calculate_patient)(s, model_file, table_file, sol_directory, r) for s in solutions)

        for s in solutions:
                calculate_patient(s, model_file, table_file, sol_directory, r)

        
        for s in solutions:
                patient = s[7:-4]
                with open(sol_directory+"/tmp"+patient+".tmp", "r") as f:
                        #content = f.readlines()
                        reader = csv.reader(f)
                        i = False
                        for r in reader:
                                if i:
                                        print("ERRRRORRR!")
                                matrix.append(r)
                                i = True
                                
                        #print content[0].strip()
                        #matrix.append(content[0].strip())
        
	#for s in solutions:
	#	patient = s[7:-4]
	#	flux_row = [patient]
	#	flow_dict = get_flow_dict(model_file, patient, table_file, sol_directory)
        #        for i in r[1]:
	#		if i in flow_dict:
	#			flux_row.append(flow_dict[i])
	#		else:
	#			flux_row.append("NA")
	#	matrix.append(flux_row)

	return matrix
	
	
def reduce_matrix(m):
	first_row = m[0]
			
			
def count_NA(v):
	c = 0
	for i in v:
		if i == "NA":
			c += 1
	return c
	
def save_matrix(m, f):
	t_m = np.transpose(np.array(m))
	f = open(f, "w+")
	for i in t_m:
		f.write("\t".join([str(j) for j in i]) + "\n")
	f.close()
	
	
def integrate(v1, v2):
	if v1 == '0' and v2 == '0':
		return '0'
	if v1 == '0' and v2 == '1':
		return '1'
	if v1 == '0' and v2 == "NA":
		return '0'
	if v1 == "NA" and v2 == "0":
		return '0'
	if v1 == "NA" and v2 == "1":
		return '1'
	if v1 == '1' and v2 == '0':
		return '1'
	if v1 == '1' and v2 == "NA":
		return '1'
	if v1 == "NA" and v2 == "NA":
		return '0'
	if v1 == '1' and v2 == '1':
		return '1'

def filter1(infile, outfile):
	""" Wykonuje pierwsze filtrowanie.
	IN: matrix.csv
	OUT: filtered1.csv"""
	f = open(infile)
	g = open(outfile, "w+")
	
	head = f.readline().strip().split("\t")[1:]
	g.write("reaction" + "\t" + "\t".join(["X" + i for i in head]) + "\n")
	
	for i in f:
		l = i.strip().split("\t")
		if count_NA(l) < 0.8*(len(l)-1):
			g.write(i)
	g.close()
	f.close()

def filter2(infile, outfile):
	""" Wykonuje drugie filtrowanie.
	IN: filtered1.csv
	OUT: filtered2.csv"""
	f = open(infile)
	m = []
	n = []
	names = []
	for i in f:
		m.append(i.strip().split("\t"))
	for i in range(len(m)-1):
		name = m[i][0][:-1]
		row = []
		if m[i][0][:-1] == m[i+1][0][:-1]:
			row.append(name)
			for j in range(1, len(m[i])):
				res = integrate(m[i][j], m[i+1][j])
				if res == None:
					print m[i][j], m[i+1][j]
				row.append(res)
			if name not in names:
				n.append(row)
				names.append(name)
		else:
			if name not in names:
				row = [name]
				row.extend(m[i][1:])
				n.append(row)
				names.append(name)
	f.close()
		
	g = open(outfile, "w+")
	for i in n:
		g.write("\t".join([str(j) for j in i]) + "\n")
	g.close()

def hamming(x, y):
	result = 0
	for i in range(len(x)):
		if x[i] != y[i]:
			result += 1
	return result

def filter3(infile, outfile):
	""" Hamming filter """
	f = open(infile)
	g = open(outfile, "w+") 

	s = []

	g.write(f.readline())
	for i in f:
		l = i.strip().split("\t")
		#reactions = l[0:]
		reactions = l[1:]
	
		if len(s) == 0:
			s.append(reactions)
		else:
			c = 0
			for j in s:
				h = hamming(j, reactions)
				if h == 0:
					c += 1
			if c == 0:
				s.append(reactions)
				g.write(i)	
	g.close()
	f.close()

if __name__ == "__main__":
	metabolic_model = sys.argv[1]
	gene_matrix = sys.argv[2]
	sols = sys.argv[3]
	m = create_matrix(metabolic_model, gene_matrix, sols)
        save_matrix(m, "matrix_0.csv")
	filter1("matrix_0.csv", "matrix_1.csv")
	filter2("matrix_1.csv", "matrix_2.csv")
	#filter3("matrix_2.csv", "matrix_3.csv")
