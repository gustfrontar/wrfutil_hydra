#This script compiles the fortran code and generates a python module
export PATH=/opt/anaconda3/bin/:$PATH
export LD_LIBRARY_PATH=/opt/anaconda3/lib/:$LD_LIBRARY_PATH

export FC=ifort
export F77=ifort
export F90=ifort
export F2PY=f2py
#export F2PY=f2py3

export FFLAGS='-fopenmp -lgomp -O3 -fPIC'
#export FFLAGS='-g -fbacktrace -fPIC'

rm *.o *.mod *.so

$F2PY -c -lgomp --f90flags="$FFLAGS" -m so_methods ./so_methods.f90
if ls so_methods.cpython*.so 1> /dev/null 2>&1
then
     echo "Compilacion EXITOSA !!!!!"
else
     echo "ERROR: Algo salio mal"
fi
