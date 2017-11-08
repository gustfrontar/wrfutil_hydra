#This script compiles the fortran code and generates a python module
export PATH=/share/anaconda3/bin/:$PATH
export LD_LIBRARY_PATH=/share/anaconda3/lib/:$LD_LIBRARY_PATH

export FC=gfortran
export F90=gfortran

f2py  -c -lgomp --f90flags="-fopenmp -lgomp " -m common_functions common_functions.f90

#f2py3  -c --f90flags="-g" -m common_functions common_functions.f90
