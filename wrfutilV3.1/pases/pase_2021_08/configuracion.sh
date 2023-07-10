read -r -d '' USO << EOF

        Ud. deberia usar este escript SOLO UNA VEZ de la siguiente manera:
                $0 < desa | testing | prod >
        Recuerde exportar la variable de entorno WRFUTILDIR=<directorio de instalacion del wrfutil>
EOF
: ${1:?"$USO"}  ${WRFUTILDIR:?"$USO"}
AMBIENTE=$1
TIMESTAMP=$(date -d now +"%Y%m%d%H%M")

## Configuracion del config.env para establecer el ambiente
cp $WRFUTILDIR/config.env $WRFUTILDIR/config.env_$TIMESTAMP
echo "########  Ambiente de Ejecucion ###########" >> $WRFUTILDIR/config.env
echo "export AMBIENTE=$AMBIENTE" >> $WRFUTILDIR/config.env


## Configuracion del experimento.conf
cp $WRFUTILDIR/RUNs/deterministico/experimento.conf $WRFUTILDIR/RUNs/deterministico/experimento.conf_$TIMESTAMP
sed -i -e "s|export NIOT=6|export NIOT=10|g" $WRFUTILDIR/RUNs/deterministico/experimento.conf
sed -i -e "s|export NIOG=2|export NIOG=3|g" $WRFUTILDIR/RUNs/deterministico/experimento.conf
sed -i -e "s|export WRFPROC=1920|export WRFPROC=1938|g" $WRFUTILDIR/RUNs/deterministico/experimento.conf

## Agregados para correr el UPP
mkdir $WRFUTILDIR/RUNs/deterministico/POST/UPP
cp $WRFUTILDIR/templates/guardo_*.cmd $WRFUTILDIR/RUNs/deterministico/POST/UPP

## Agregar eel limpiado de los gribs en el limpiar 
echo 'MATCHFILE="WRFPRS_d01*"' >> $WRFUTILDIR/RUNs/deterministico/experimento.limp
echo 'TIMEPOGUARDO=10' >> $WRFUTILDIR/RUNs/deterministico/experimento.limp 
echo 'LIMPBASEDIR="POST"'>> $WRFUTILDIR/RUNs/deterministico/experimento.limp
echo "read -r -d '' PROCESAR <<< 'rm \$file'" >> $WRFUTILDIR/RUNs/deterministico/experimento.limp
echo 'limpiar '>> $WRFUTILDIR/RUNs/deterministico/experimento.limp
echo 'MATCHFILE="wrfprs_d01*"' >> $WRFUTILDIR/RUNs/deterministico/experimento.limp
echo 'TIMEPOGUARDO=10' >> $WRFUTILDIR/RUNs/deterministico/experimento.limp 
echo 'LIMPBASEDIR="POST"'>> $WRFUTILDIR/RUNs/deterministico/experimento.limp
echo "read -r -d '' PROCESAR <<< 'rm \$file'" >> $WRFUTILDIR/RUNs/deterministico/experimento.limp
echo 'limpiar '>> $WRFUTILDIR/RUNs/deterministico/experimento.limp

## Ejecutando la modificacion de los dataframes historicos
source activate wrfutil
python agrega_columnas.py


