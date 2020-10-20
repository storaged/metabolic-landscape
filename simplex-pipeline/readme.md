# LANDSCAPE DETERMINATION

The dictionary contains scipts and pipeline file that solves the linear problem using SIMPLEX algorithm. 
To perform the computations it is necessary to match the requirements

## Requirements: 
1. Prepare the environment: 
* `gurobi_cl` from Gurobi Command-Line Tool, see [Gurobi Docs](https://www.gurobi.com/documentation/9.0/refman/grb_command_line_tool.html) and [Gurobi License](https://www.gurobi.com/academia/academic-program-and-licenses/)
* `python2.7` libraries: `joblib`, `numpy`, `matplotlib`, `glob`

2. Set *paths* to needed files in the main pipeline script: **perform_full_metabolic_analysis.sh**
* **model** a path to the .sfba model file
* **soldir** a path to the dictionary in which resulting files will be stored 
* **act_mat** a path to the file with gene activity values for each sample in the study 

## Example usage:
`sh perform_full_metabolic_analysis.sh` 
