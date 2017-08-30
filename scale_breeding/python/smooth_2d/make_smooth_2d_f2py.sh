#This script compiles the fortran code and generates a python module
 
#For paralel with openmp.
F90=gfortran

f2py  -m smooth_2d -c common_smooth2d.f90
