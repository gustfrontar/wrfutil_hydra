#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

###  Leyenda de USO

export WRFUTILDIR="/vol0301/data/hp150019/u10335/PREVENIR/wrfutil_hydra/POST/4km/r1"

export FILEIN="test_input/input1.nc,test_input/input2.nc,test_input/input3.nc"
export HORA=0
export MIEM=01
export PATHOUT="test_output"
export ANALISYS=3600
mkdir -p $PATHOUT 

### CONFIGURACION

##### FIN INICIALIZACION ######

ulimit -s unlimited

#conda activate myenv 
python $WRFUTILDIR/post_ens_wrfout_to_netcdf.py 

echo "Terminando"


