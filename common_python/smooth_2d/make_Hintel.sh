#This script compiles the fortran code and generates a python module
 
#For paralel with openmp.
F90=ifort
FC=ifort 

f2py --f90flags="-fopenmp -lgomp " -m smooth_2d -c common_smooth2d.f90
