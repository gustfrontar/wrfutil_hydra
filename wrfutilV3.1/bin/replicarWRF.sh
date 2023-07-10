#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

###  Leyenda de USO

read -r -d USO << EOF
        Ud. deberia usar este escript de la siguiente manera:
                $0 <nuevo nombre entonro WRF > </path/al/directorion/instalacion/WRF/> <archivo configuracion librerias>
        Recuerde exportar la variable de entorno WRFUTILDIR=<directorio de instalacion del wrfutil>
EOF
: ${1:?"$USO"}  ${2:?"$USO"} ${3:?"$USO"} ${WRFUTILDIR:?"$USO"}
[ ! -f $3 ] && echo "archivo de configuracion de entorno inexistente" && exit 1

### PARAMETROS

WRFENV=$1
WRFPATH=$2
ENTORNO=$3

### CONFIGURACION

source $WRFUTILDIR/lib/errores.env
CONFIG=$WRFUTILDIR/config.env
#CONFIG=$(readlink -e ${CONFIG:-"$WRFUTILDIR"/config/config.env} 2> /dev/null)
[ -z "$CONFIG" ] && dispararError 2 "$CONFIG"
source $CONFIG


##### FIN INICIALIZACION ######

## Verifica que no exista una replica previa
[ -f "$WRFUTILDIR/WRFs/${WRFENV}.tar" ] && dispararError 3 "$WRFENV"

## Verifica que exista el path 
[ ! -d "$WRFPATH" ] && dispararError 4 "$WRFPATH"

### WRF 

mkdir -p $TMPDIR/WRF
[ $? -ne 0 ] && dispararError 2 "$TMPDIR/WRF"
cd $TMPDIR/WRF

ln -s $WRFPATH/main/wrf.exe
ln -s $WRFPATH/main/tc.exe
ln -s $WRFPATH/main/real.exe
ln -s $WRFPATH/main/ndown.exe
ln -s $WRFPATH/run/aerosol.formatted
ln -s $WRFPATH/run/aerosol_lat.formatted
ln -s $WRFPATH/run/aerosol_lon.formatted
ln -s $WRFPATH/run/aerosol_plev.formatted
ln -s $WRFPATH/run/bulkdens.asc_s_0_03_0_9
ln -s $WRFPATH/run/bulkradii.asc_s_0_03_0_9
ln -s $WRFPATH/run/CAM_ABS_DATA
ln -s $WRFPATH/run/CAM_AEROPT_DATA
ln -s $WRFPATH/run/CAMtr_volume_mixing_ratio.A1B
ln -s $WRFPATH/run/CAMtr_volume_mixing_ratio.A2
ln -s $WRFPATH/run/CAMtr_volume_mixing_ratio.RCP4.5
ln -s $WRFPATH/run/CAMtr_volume_mixing_ratio.RCP6
ln -s $WRFPATH/run/CAMtr_volume_mixing_ratio.RCP8.5
ln -s $WRFPATH/run/capacity.asc
ln -s $WRFPATH/run/CCN_ACTIVATE.BIN
ln -s $WRFPATH/run/CLM_ALB_ICE_DFS_DATA
ln -s $WRFPATH/run/CLM_ALB_ICE_DRC_DATA
ln -s $WRFPATH/run/CLM_ASM_ICE_DFS_DATA
ln -s $WRFPATH/run/CLM_ASM_ICE_DRC_DATA
ln -s $WRFPATH/run/CLM_DRDSDT0_DATA
ln -s $WRFPATH/run/CLM_EXT_ICE_DFS_DATA
ln -s $WRFPATH/run/CLM_EXT_ICE_DRC_DATA
ln -s $WRFPATH/run/CLM_KAPPA_DATA
ln -s $WRFPATH/run/CLM_TAU_DATA
ln -s $WRFPATH/run/co2_trans
ln -s $WRFPATH/run/coeff_p.asc
ln -s $WRFPATH/run/coeff_q.asc
ln -s $WRFPATH/run/constants.asc
ln -s $WRFPATH/run/ETAMPNEW_DATA
ln -s $WRFPATH/run/ETAMPNEW_DATA_DBL
ln -s $WRFPATH/run/ETAMPNEW_DATA.expanded_rain
ln -s $WRFPATH/run/ETAMPNEW_DATA.expanded_rain_DBL
ln -s $WRFPATH/run/GENPARM.TBL
ln -s $WRFPATH/run/grib2map.tbl
ln -s $WRFPATH/run/gribmap.txt
ln -s $WRFPATH/run/kernels.asc_s_0_03_0_9
ln -s $WRFPATH/run/kernels_z.asc
ln -s $WRFPATH/run/LANDUSE.TBL
ln -s $WRFPATH/run/masses.asc
ln -s $WRFPATH/run/MPTABLE.TBL
ln -s $WRFPATH/run/ozone.formatted
ln -s $WRFPATH/run/ozone_lat.formatted
ln -s $WRFPATH/run/ozone_plev.formatted
ln -s $WRFPATH/run/README.namelist
ln -s $WRFPATH/run/README.tslist
ln -s $WRFPATH/run/RRTM_DATA
ln -s $WRFPATH/run/RRTM_DATA_DBL
ln -s $WRFPATH/run/RRTMG_LW_DATA
ln -s $WRFPATH/run/RRTMG_LW_DATA_DBL
ln -s $WRFPATH/run/RRTMG_SW_DATA
ln -s $WRFPATH/run/RRTMG_SW_DATA_DBL
ln -s $WRFPATH/run/SOILPARM.TBL
ln -s $WRFPATH/run/termvels.asc
ln -s $WRFPATH/run/tr49t67
ln -s $WRFPATH/run/tr49t85
ln -s $WRFPATH/run/tr67t85
ln -s $WRFPATH/run/URBPARM.TBL
ln -s $WRFPATH/run/URBPARM_UZE.TBL
ln -s $WRFPATH/run/VEGPARM.TBL
ln -s $WRFPATH/run/wind-turbine-1.tbl
#cp  $WRFPATH/../envvars.sh .
cp  $ENTORNO .
tar cf ../${WRFENV}.tar  .
mv ../${WRFENV}.tar $WRFUTILDIR/WRFs/
cd ..
rm -r $TMPDIR/WRF

