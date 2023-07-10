#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

###  Leyenda de USO

read -r -d '' USO << EOF

        Ud. deberia usar este escript de la siguiente manera:
                $0 <nombre entorno >
        Recuerde exportar la variable de entorno WRFUTILDIR=<directorio de instalacion del wrfutil>
EOF
: ${1:?"$USO"}  ${WRFUTILDIR:?"$USO"}


### PARAMETROS

export NOMBRE=$1

### CONFIGURACION

[ ! -f $WRFUTILDIR/lib/errores.env ] && exit 1
source $WRFUTILDIR/lib/errores.env
CONFIG=$WRFUTILDIR/config.env
[ -f $WRFUTILDIR/RUNs/$NOMBRE/config.env ] && CONFIG=$WRFUTILDIR/RUNs/$NOMBRE/config.env
CONFIG=$(readlink -e ${CONFIG} 2> /dev/null)
[ -z "$CONFIG" ] && dispararError 4 "config.env"
source $CONFIG  $NOMBRE
[ ! -d "$BASEDIR" ] && dispararError 7 "$NOMBRE"
[ ! -f "$BASEDIR/$EXPCONF" ] && dispararError 4 "$BASEDIR/$EXPCONF"
source $BASEDIR/$EXPCONF
[ ! -f "$BASEDIR/$EXPDEP" ] && dispararError 4 "$BASEDIR/$EXPDEP"
source $BASEDIR/$EXPDEP
[ ! -f "$BASEDIR/$EXPASIM" ] && dispararError 4 "$BASEDIR/$EXPASIM"
source $BASEDIR/$EXPASIM


######### FIN VERIFICACION 

cd $BASEDIR/LETKF
#script de ejecucion
read -r -d '' QSCRIPTCMD << EOF
ulimit -s unlimited
ulimit -l unlimited
source envvars.sh
NSLOTS=$(($OBSWIN/$OBSFREC+1))

DIRANAL=$BASEDIR/HIST/ANA/$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours +$((10#$CICLO_MIN + $ANALISIS)) minutes" +"%Y%m%d_%H%M%S")/$NOMBRE/
DIRGUES=$BASEDIR/HIST/GUESS/$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours+$((10#$CICLO_MIN+$ANALISIS)) minutes" +"%Y%m%d_%H%M%S")/$NOMBRE/
###
# Post.procesamos los Analisis
###
source $CONDADIR/deactivate || $CONDADIR/conda deactivate
source $CONDADIR/activate wrfutil || $CONDADIR/conda activate wrfutil
export PATHOUT=$BASEDIR/HIST/POST/$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours +$((10#$CICLO_MIN + $ANALISIS)) minutes" +"%Y%m%d_%H%M%S")/$NOMBRE/
mkdir -p \$PATHOUT
export PATHOUT_PLOT=$BASEDIR/HIST/PLOT/$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours +$((10#$CICLO_MIN + $ANALISIS)) minutes" +"%Y%m%d_%H%M%S")/$NOMBRE/
export FILEIN_ANA=\$DIRANAL/anal\$(printf %05d $((10#$MIEMBRO_FIN+1)))
export FILEIN_GUES=\$DIRANAL/gues\$(printf %05d $((10#$MIEMBRO_FIN+1)))


### Loops sobre los archivos que se postprocesan
files=( \$FILEIN_ANA \$FILEIN_GUES )
for FILEIN in  \${files[@]}
do
	export FILEIN=\$FILEIN
	python $WRFUTILDIR/bin/python/write_post_smn_asim_oper_lev.py  
	python $WRFUTILDIR/bin/python/write_post_smn_asim_oper_sfc.py
done

python $WRFUTILDIR/bin/python/write_and_plot_smn_asim_rmsd_update.py
##
# Graficamos las observaciones asimiladas
##
PATH=$PATH:$CONDADIR

echo "Graficamos las observaciones asimiladas"
cd ${LETKFDIR}
${LETKFDIR}/obs_dat2txt.exe obs.dat > obs_assim.txt
FECHA_ASSIM="$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours +$((10#$CICLO_MIN + $ANALISIS)) minutes" +"%Y %m %d %H %M")"
source $CONDADIR/deactivate || $CONDADIR/conda deactivate
source $CONDADIR/activate wrfutil || $CONDADIR/conda activate wrfutil
python $WRFUTILDIR/bin/python/graficar_obs_desdetxt.py \$FECHA_ASSIM obs_assim.txt
# Guardamos las figuras: 
PATHOUT=$BASEDIR/HIST/PLOT/$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours +$((10#$CICLO_MIN + $ANALISIS)) minutes" +"%Y%m%d_%H%M%S")/$NOMBRE/
mkdir -p \$PATHOUT 
mv model*.jpg \$PATHOUT/

# Lo nuevo de Pau
export PATH_OBS=$DIROBSBIN
export PATH_PLOT=\$PATHOUT
python $WRFUTILDIR/bin/python/plot_obsdat_obscount.py \$FECHA_ASSIM \$NSLOTS
python $WRFUTILDIR/bin/python/write_obsdat_obsstats.py \$FECHA_ASSIM



EOF



# Parametros de encolamiento
QDEPEND_NAME=$DEPLETKFPOST
QDEPENDTYPE=$TYPELETKFPOST
QPROC_NAME=LETKFPOST_$1_$PASO
QPROC=$LETKFPOSTPROC
QTHREADS=$LETKFPOSTTHREADS
QARRAY=""
QWALLTIME=$LETKFPOSTWALLTIME
QEXCLU=1

# Encolar
queue

read -r -d '' QSCRIPTCMD << "EOF"
ulimit -s unlimited
ulimit -l unlimited
source envvars.sh
NSLOTS=$(($OBSWIN/$OBSFREC+1))

DIRANAL=$BASEDIR/HIST/ANA/$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours +$((10#$CICLO_MIN + $ANALISIS)) minutes" +"%Y%m%d_%H%M%S")/$NOMBRE/
###
# Post.procesamos el ensamble de los Analisis
###
source $CONDADIR/deactivate || $CONDADIR/conda deactivate
source $CONDADIR/activate wrfutil || $CONDADIR/conda activate wrfutil
export PATHOUT=$BASEDIR/HIST/POST/$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours +$((10#$CICLO_MIN + $ANALISIS)) minutes" +"%Y%m%d_%H%M%S")/$NOMBRE/
mkdir -p $PATHOUT
export PATHOUT_PLOT=$BASEDIR/HIST/PLOT//$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours +$((10#$CICLO_MIN + $ANALISIS)) minutes" +"%Y%m%d_%H%M%S")/$NOMBRE/


export FILEIN_ANA=$(ls $DIRANAL/anal* | grep -v $(printf %05d $((10#$MIEMBRO_FIN+1))) |grep -v sp )
files=( $FILEIN_ANA )
printf "%s\n" "${files[@]}" | xargs -P 10 -i sh -c 'export FILEIN={}; echo $FILEIN; python  $WRFUTILDIR/bin/python/write_post_smn_asim_ens_oper_ANA.py'

export ANAFILES=$(ls $PATHOUT/*ANA*)
export MIEMBROS=$((10#$MIEMBRO_FIN-10#$MIEMBRO_INI+1))
python $WRFUTILDIR/bin/python/write_spread_smn_asim_oper_ANA.py
EOF

# Parametros de encolamiento
QDEPEND_NAME=$DEPLETKFPOST
QDEPENDTYPE=$TYPELETKFPOST
QPROC_NAME=LETKFPOST_A_$1_$PASO
QPROC=$LETKFPOSTPROC
QTHREADS=$LETKFPOSTTHREADS
QARRAY=""
QWALLTIME=$LETKFPOSTWALLTIME
QEXCLU=1

# Encolar
queue

read -r -d '' QSCRIPTCMD << "EOF"
ulimit -s unlimited
ulimit -l unlimited
source envvars.sh
NSLOTS=$(($OBSWIN/$OBSFREC+1))

DIRGUES=$BASEDIR/HIST/GUESS/$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours+$((10#$CICLO_MIN+$ANALISIS)) minutes" +"%Y%m%d_%H%M%S")/$NOMBRE/
###
# Post.procesamos el ensamble de los First Guess
###
source $CONDADIR/deactivate || $CONDADIR/conda deactivate
source $CONDADIR/activate wrfutil || $CONDADIR/conda activate wrfutil
export PATHOUT=$BASEDIR/HIST/POST/$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours +$((10#$CICLO_MIN + $ANALISIS)) minutes" +"%Y%m%d_%H%M%S")/$NOMBRE/
mkdir -p $PATHOUT
export PATHOUT_PLOT=$BASEDIR/HIST/PLOT//$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours +$((10#$CICLO_MIN + $ANALISIS)) minutes" +"%Y%m%d_%H%M%S")/$NOMBRE/
export FILEIN_GUES=$DIRANAL/gues$(printf %05d $((10#$MIEMBRO_FIN+1)))


export FILEIN_GUES=$(ls $DIRGUES/*/wrf*)
files=( $FILEIN_GUES )
printf "%s\n" "${files[@]}" | xargs -P 10 -i sh -c 'export FILEIN={}; echo $FILEIN; python  $WRFUTILDIR/bin/python/write_post_smn_asim_ens_oper_GUES.py'

export GUESFILES=$(ls $PATHOUT/*GUE*)
export MIEMBROS=$((10#$MIEMBRO_FIN-10#$MIEMBRO_INI+1))
python $WRFUTILDIR/bin/python/write_spread_smn_asim_oper_GUES.py

EOF

# Parametros de encolamiento
QDEPEND_NAME=$DEPLETKFPOST
QDEPENDTYPE=$TYPELETKFPOST
QPROC_NAME=LETKFPOST_G_$1_$PASO
QPROC=$LETKFPOSTPROC
QTHREADS=$LETKFPOSTTHREADS
QARRAY=""
QWALLTIME=$LETKFPOSTWALLTIME
QEXCLU=1

# Encolar
queue

