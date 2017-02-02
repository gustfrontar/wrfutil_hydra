#!/bin/bash
# Runing on LINUX

#
# DATA SOURCE
#
DATASRC=2 # 1:NOAA/NOMADS , 2:UCAR/DSS
export http_proxy=""

#export http_proxy="proxy.fcen.uba.ar:8080"

# Initial time
IDATE=20070101000000
# Final time
EDATE=20070102000000
# Frequency
DF=86400   #file frequency in seconds

EMAIL=gustfrontar@gmail.com
PASSWD=pr0n0st1c0
#if [ $# -ne 2 ] ; then
#  echo "USAGE: $0 EMAIL PASSWD"
#  exit
#fi

OBSNAME=PREPBUFREUROPE

OBSDIR=$HOME/share/DATA/OBS/prepbufr/$OBSNAME

source $HOME/share/LETKF_WRF/wrf/run/util.sh #Load time functions.


#
# WKDIR
#

mkdir -p $OBSDIR
cd $OBSDIR

#
# Settings for DATASRC
#


v=`wget -V |grep 'GNU Wget ' | cut -d ' ' -f 3`
a=`echo $v | cut -d '.' -f 1`
b=`echo $v | cut -d '.' -f 2`
tmp=`expr  100 \* $a  + $b `

if [ $tmp > 109 ] ; then
 opt='wget --no-check-certificate'
else
 opt='wget'
fi

DATAURL="http://rda.ucar.edu/data/ds337.0/tarfiles/"
WGET="wget --no-check-certificate"
OPT1='-O /dev/null --save-cookies auth.rda_ucar_edu --post-data'
OPT2="email=${EMAIL}&passwd=${PASSWD}&action=login"
$opt $OPT1="$OPT2" https://rda.ucar.edu/cgi-bin/login


OPT1="-N --load-cookies auth.rda_ucar_edu"
OPT2="$opt $OPT1 $DATAURL"

#
# Cycle run # MAIN LOOP #
#
CDATE=$IDATE

while [ $CDATE -le $EDATE ]
do
echo " >>>"
echo " >>> DOWNLOADING $CDATE "
echo " >>>"
#
# GET NCEP PREPBUFR
#

cy=`echo $CDATE | cut -c1-4`
cm=`echo $CDATE | cut -c5-6`
cd=`echo $CDATE | cut -c7-8`
ch=`echo $CDATE | cut -c9-10`
cn=`echo $CDATE | cut -c11-12`
cs=`echo $CDATE | cut -c13-14`

  echo ${OPT2}/${cy}/prepbufr.${cy}${cm}${cd}.wo40.tar.gz
  ${OPT2}/${cy}/prepbufr.${cy}${cm}${cd}.wo40.tar.gz
  #$WGET $OPT1 $DATAURL/${cy}/prepbufr.${cy}${cm}${cd}.wo40.tar.gz
  gunzip -f prepbufr.${cy}${cm}${cd}.wo40.tar.gz
  tar -xvf prepbufr.${cy}${cm}${cd}.wo40.tar
  rm -fr prepbufr.${cy}${cm}${cd}.wo40.tar

  mv ${cy}${cm}${cd}.wo40/prepbufr.gdas.${cy}${cm}${cd}00.wo40 ${cy}${cm}${cd}00.prepbufr.nr
  mv ${cy}${cm}${cd}.wo40/prepbufr.gdas.${cy}${cm}${cd}06.wo40 ${cy}${cm}${cd}06.prepbufr.nr
  mv ${cy}${cm}${cd}.wo40/prepbufr.gdas.${cy}${cm}${cd}12.wo40 ${cy}${cm}${cd}12.prepbufr.nr
  mv ${cy}${cm}${cd}.wo40/prepbufr.gdas.${cy}${cm}${cd}18.wo40 ${cy}${cm}${cd}18.prepbufr.nr

  rm -fr ${cy}${cm}${cd}.wo40

echo $CDATE $DF
CDATE=`date_edit2 $CDATE $DF `


done

rm -f auth.rda_ucar_edu
echo "NORMAL END"
