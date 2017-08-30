#This script compiles the fortran code and generates a python module
. $HOME/set_env_anaconda3.sh

export FC=ifort
export F90=ifort

f2py3  -c -lgomp --f90flags="-fopenmp -lgomp " -m common_functions common_functions.f90

#f2py3  -c --f90flags="-g" -m common_functions common_functions.f90
