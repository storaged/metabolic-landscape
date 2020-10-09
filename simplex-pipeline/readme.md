# LANDSCAPE DETERMINATION

The dictionary contains scipts and pipeline file that solves the linear problem using SIMPLEX algorithm. 
To perform the computations it is necessary to match the requirements

Requirements: 
- gurobi_cl from Gurobi Command-Line Tool, see: https://www.gurobi.com/documentation/9.0/refman/grb_command_line_tool.html
- python 2.7 libraries: joblib, numpy, matplotlib, glob

Example usage:
sh perform_full_metabolic_analysis.sh ../example-datasets/random_activity_matrix.tab 
