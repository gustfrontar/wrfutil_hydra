#This script compiles the fortran code and generates a python module

f2py3  -c -lgomp --f90flags="-fcheck=all -fopenmp -lgomp -O3" -m common_qc_tools common_qc_tools.f90


