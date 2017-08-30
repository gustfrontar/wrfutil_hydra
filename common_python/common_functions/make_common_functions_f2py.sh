#This script compiles the fortran code and generates a python module
. /data/opt/Env/Env_intel
source /etc/profile
#module load common/anaconda/4.2.0_python3.5  #Anaconda no funciona bien en hakshu.
module load common_intel/python/3.4.3

f2py3  -c -lgomp --f90flags="-fopenmp -lgomp " -m common_functions common_functions.f90

#f2py3  -c --f90flags="-g" -m common_functions common_functions.f90
