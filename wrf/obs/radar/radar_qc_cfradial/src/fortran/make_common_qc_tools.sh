#This script compiles the fortran code and generates a python module

f2py  -c -lgomp --f90flags="-fopenmp -lgomp -O3" -m common_qc_tools common_qc_tools.f90


