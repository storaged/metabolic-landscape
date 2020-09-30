import glob, os, sys
import csv

def get_flux_matrix(solutions_dir):
    os.chdir(solutions_dir)

    matrix = []
    sol_files = glob.glob("*.sol")
    sol_files = [sol_file for sol_file in sol_files if os.stat(sol_file).st_size]
    for idx, file_name in enumerate(glob.glob("*.sol")):
        if idx == 0 :
            print "x"
            matrix.append([os.path.splitext(file_name)[0][len('program'):]])
            with open(file_name, 'r') as solfile:
                for line in solfile:
                    if line[0] == 'v':
                        matrix.append([x.strip() for x in line.split(" ")])
            print matrix
        else:
            matrix[0].append(os.path.splitext(file_name)[0][len('program'):])
            with open(file_name, 'r') as solfile:
                no = 1
                for line in solfile:
                    if line[0] == 'v':
                        matrix[no].append(line.split(" ")[1].strip())
                        no += 1
    
    with open("matrix_flux.csv", 'w') as csvfile:       
        writer = csv.writer(csvfile, delimiter="\t")
        [writer.writerow(r) for r in matrix]


    
if __name__ == "__main__":
    solutions_dir = sys.argv[1]
    get_flux_matrix(solutions_dir)
