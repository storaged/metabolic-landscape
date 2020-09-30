import numpy as np
import re
import random
import os
import sys

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
        
def get_patients_dict(table):
	""" Tworzy slownik pacjent-aktywnosc genu."""
	f = open(table)
	patients = f.readline().strip().split("\t")[1:]
		 
	patients_dict = {}
	for i in patients:
		patients_dict[i] = {}
		 
	for i in f:
		l = i.strip().split("\t")
		gene = l[0]
		for j in range(len(l[1:])):
			patients_dict[patients[j]][gene] = int(l[1:][j])
	return patients_dict


def istrue(m, name, gene_activity):
	""" Sprawdza wartosc logiczna formul genetycznych dla danego pacjenta"""
	
	activity_dict = {}
	
	f = open(m)
	for i in f:
		l = i.strip().split("\t")
		rule = l[4]
		if len(rule) > 0: # niepusta zasada
			dig = ""
			rule_string = ""
			started = False
			for i in rule:
				if i.isdigit():
					started = True
					dig += i
				else:
					if started:
						if dig not in gene_activity:
							rule_string += "0"
							rule_string += i
							dig = ""
							started = False
						else:
							rule_string += str(gene_activity[dig])
							rule_string += i
							dig = ""
							started = False
					else:
						rule_string += i
			if started:
				if dig not in gene_activity:
					rule_string += "0"
				else:
					rule_string += str(gene_activity[dig])
				
			try:
				phrase = rule_string.replace("OR", "|").replace("AND", "&")
				phrase_eval = eval(phrase)
				activity_dict[l[0]] = phrase_eval
			except:
				c += 1				
		else:
			activity_dict[l[0]] ="NA"
	return activity_dict
        
def get_metabolite_list(m):
	all_metabolites = []
	f = open(m)
	for i in f:
		l = i.strip().split("\t")
		equation = l[1] + " "
		pattern = "[x_M|M].*?[ ]"
		result = re.findall(pattern, equation)
		result_strip = [i.strip() for i in result]
		for j in result_strip:
			if j not in all_metabolites:
				all_metabolites.append(j)
	all_metabolites.sort()
	return all_metabolites
	
def get_reaction_list_irr(m):
	all_reactions = []
	states = []
	d = {}
	c = 1
	f = open(m)
	for i in f:
		l = i.strip().split("\t")
		reaction = l[0]
		if float(l[2]) < 0:
			states.append("i"+str(c)+"p")
			states.append("i"+str(c)+"m")
			d[l[0]] = ["i"+str(c)+"p", "i"+str(c)+"m"]
			all_reactions.append(l[0])
			c+=1
		else:
			states.append("i"+str(c)+"p")
			d[l[0]] = ["i"+str(c)+"p"]
			all_reactions.append(l[0])
			c+=1
	return all_reactions, states, d
	
def get_matrix_irr(m):
	matrix = []
	metabolite_list = get_metabolite_list(m)

	f = open(m)
	for i in f:
		l = i.strip().split("\t")
		reaction = l[0]
		equation = l[1]
		min_flow = float(l[2])
		max_flow = float(l[3])

		metabolites = equation.split("=")
		left = [i.strip() for i in metabolites[0].strip().split(" ")]
		right = [i.strip() for i in metabolites[1].strip().split(" ")]
			
		left_coefficients = []
		left_metabolites = []
		right_coefficients = []
		right_metabolites = []
		
		gotCoefficient = False
		for i in left:
			if i != "+":
				if is_number(i):
					left_coefficients.append(-float(i))
					gotCoefficient = True
				else:
					left_metabolites.append(i)
					if not gotCoefficient:
						left_coefficients.append(-1)
					else:
						gotCoefficient = False
	
		gotCoefficient = False
		for i in right:
			if i != "+":
				if is_number(i):
					right_coefficients.append(float(i))
					gotCoefficient = True
				else:
					right_metabolites.append(i)
					if not gotCoefficient:
						right_coefficients.append(1)
					else:
						gotCoefficient = False
	
		forward_row = [0]*len(metabolite_list)
		for i in zip(left_coefficients, left_metabolites):
			m_index = metabolite_list.index(i[1])
			forward_row[m_index] = i[0]
		for j in zip(right_coefficients, right_metabolites):
			m_index = metabolite_list.index(j[1])
			forward_row[m_index] = j[0]
		matrix.append(forward_row)
		if sum(forward_row) != 0:
			## print l[0]
			## print zip(left_coefficients, left_metabolites)
			## print zip(right_coefficients, right_metabolites)
			pass
			
	
	return np.transpose(np.array(matrix))


def get_fluxes_irr(m):
	reactions = []
	f = open(m)
	command = []
	num = 1
	for i in f:
		l = i.strip().split("\t")
		reaction = l[0]
		f_min = float(l[2])
		f_max = float(l[3])

		command.append(str(f_min) + " <= " + "v" + str(num) + " <= " + str(f_max))
		reactions.append(reaction + " (+)")		
		num += 1
	return reactions, command
	

def sv(matrix):
	""" Dla zadanej macierzy stechiometrycznej zwraca zasady Sv = 0"""
	rules = []
	c = 0
	for row in matrix:
		c +=1
		rule = ""
		subrules = []

		for j in range(len(row)):
			if row[j] == 0:
				pass
			else:
				if row[j] == 1:
					subrules.append("v" + str(j+1))

				elif row[j] == -1:
					subrules.append("- v" + str(j+1))

				else:
					if row[j].is_integer():
						subrules.append(str(int(row[j])) + " v" + str(j+1))
					else:
						subrules.append(str(row[j]) + " v" + str(j+1))
		
		rule = subrules[0]
		for s in subrules[1:]:
			if s[0] != "-":
				rule = rule + " + " + s
			else:
				rule = rule + " " + s	
		rules.append(rule + " = 0")
	return rules

def write_matrix(m):
	f = open("matrix.csv", "w+")
	f.write("\t".join(all_mets)+"\n")
	for i in matrix:
		f.write("\t".join([str(j) for j in i])+"\n")
	f.close()
	
def get_v_dict_irr(model):
	d = {}
	r = get_reaction_list_irr(model)[0]
	for i in range(1, len(r)+1):
		d[r[i-1]] = i
	return d
	
def get_vs_irr(model):
	reactions = []
	f = open(model)
	d = {}
	num = 1
	for i in f:
		l = i.strip().split("\t")
		reaction = l[0]
		f_min = float(l[2])
		f_max = float(l[3])
		d[reaction] = [f_min, f_max]
		num += 1
	return d
	
def get_rules_irr(metabolic_model, reacs, activity_dict):
	v_dict_irr = get_v_dict_irr(metabolic_model)
	vs_irr = get_vs_irr(metabolic_model)
	r = get_reaction_list_irr(metabolic_model)
	
	obj = []
	rules = []
	for i in reacs:
		if activity_dict[i] == 0: # reakcja nieaktywna
			plus = r[2][i][0]
			obj.append(plus)
			rule1 = str(vs_irr[i][0]) + " " + plus + " + v" + str(v_dict_irr[i]) + " >= "  + str(vs_irr[i][0])
			rule2 = "v" + str(v_dict_irr[i]) + " + " + str(vs_irr[i][1]) + " " + plus + " <= " + str(vs_irr[i][1])  
			rules.append(rule1)
			rules.append(rule2)
			
		elif activity_dict[i] == 1: # reakcja aktywna
			if len(r[2][i]) == 2:
				plus = r[2][i][0]
				minus = r[2][i][1]
				obj.append(plus)
				obj.append(minus)
				rules.append("v" + str(v_dict_irr[i]) + " + " + str(vs_irr[i][0] - 1) +  " " + plus + " >= " + str(vs_irr[i][0]))
				rules.append("v" + str(v_dict_irr[i]) + " + " + str(vs_irr[i][1] + 1) + " " + minus + " <= " + str(vs_irr[i][1]))
			else:
				plus = r[2][i][0]
				obj.append(plus)
				rules.append("v" + str(v_dict_irr[i]) + " + " + str(vs_irr[i][0] - 1) +  " " + plus + " >= " + str(vs_irr[i][0]))			
		elif activity_dict[i] == 'NA':
			pass
			
	obj_str = " + ".join(obj)
	return obj_str, rules



def writelines(l, fname):
	f = open(fname, "w+")
	c = 1
	for i in l:
		f.write("c" + str(c) + ": " + i + "\n")
		c+=1
	f.close()

	
def writematrix(m):
	f = open("matrix.csv", "w+")
	for i in m:
		f.write("\t".join([str(j) for j in i]) + "\n")
	f.close()
	
def write_program(metabolic_model, name, gene_activity):
	""" Zapisuje program liniowy jednego pacjenta """
	print name
	r = get_reaction_list_irr(metabolic_model)
	v_dict_irr = get_v_dict_irr(metabolic_model)
	bounds = get_vs_irr(metabolic_model)
	activity_dict = istrue(metabolic_model, name, gene_activity)
	rules = get_rules_irr(metabolic_model, r[0], activity_dict)	
	
	matrix = get_matrix_irr(metabolic_model)
	s = sv(matrix)	

	g = open("program" + name + ".lp", "w+")
	g.write("Max obj: ")
	g.write(rules[0])
	g.write("\n\n")
	g.write("Subject To\n")
	c = 1
	for subject in s:
		g.write("c" + str(c) + ": " + subject+"\n")
		c += 1
	for rule in rules[1]:
		g.write("c" + str(c) + ": " + rule.replace("+ -", "-") + "\n")
		c += 1

	g.write("\n")
	
	g.write("Bounds\n")
	for i in r[0]:
		g.write(str(bounds[i][0]) + " <= " + "v" + str(v_dict_irr[i]) + " <= " + str(bounds[i][1])+"\n")
	

	g.write("\n")
	g.write("Binary\n")
	for var in r[1]:
		g.write(var + "\n")
	g.write("\n")
	g.write("End")
	g.close()
		
	
def write_all(metabolic_model, infile):
	""" Zapisuje program liniowy dla wszystkich pacjentow w biezacym folderze"""
	from joblib import Parallel, delayed
        patients = get_patients_dict(infile)
	
        Parallel(n_jobs=48)(delayed(write_program)(metabolic_model, i, patients[i]) for i in patients)
        #for i in patients:
	#	print i
	#	write_program(metabolic_model, i, patients[i])

if __name__ == "__main__":
	metabolic_model = sys.argv[1]
	gene_matrix = sys.argv[2]
	write_all(metabolic_model, gene_matrix)
