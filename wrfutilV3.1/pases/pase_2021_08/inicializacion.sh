#!/bin/bash
read -r -d '' USO << EOF
        Recuerde exportar la variable de entorno WRFUTILDIR=<directorio de instalacion del wrfutil>
EOF
: ${WRFUTILDIR:?"$USO"}

BASEDIR=/data/calibracion/wrfutilV3.1/RUNs/
LOCALDIR=$WRFUTILDIR/pases/pase_2021_08/

cd $BASEDIR

echo "deterministico/CALIB/viento" > ${LOCALDIR}/filestar.txt
echo "ensamble/CALIB/viento">>${LOCALDIR}/filestar.txt
echo "CalibGFS/CALIB/viento">>${LOCALDIR}/filestar.txt
echo "CalibGEFS/CALIB/viento">>${LOCALDIR}/filestar.txt
echo "deterministico/PostDF/superficie">>${LOCALDIR}/filestar.txt
echo "ensamble/PostDF/superficie">>${LOCALDIR}/filestar.txt
echo "CalibGFS/PostDF/superficie">>${LOCALDIR}/filestar.txt
echo "CalibGEFS/PostDF/superficie">>${LOCALDIR}/filestar.txt

tar cf ${LOCALDIR}/operativo.tar --exclude='*.SLURM' -T  ${LOCALDIR}/filestar.txt 
cd $WRFUTILDIR/RUNs/
tar xvf ${LOCALDIR}/operativo.tar
