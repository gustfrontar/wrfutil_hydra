#!/bin/bash
# Runing on LINUX

#
# DATA SOURCE
#
DATASRC=2 # 1:NOAA/NOMADS , 2:UCAR/DSS
export http_proxy=""

#export http_proxy="proxy.fcen.uba.ar:8080"

# Initial time
IDATE=20160518000000
# Final time
EDATE=20160610000000
# Frequency
DF=86400   #file frequency in seconds

EMAIL=$1
PASSWD=$2
if [ $# -ne 2 ] ; then
  echo "USAGE: $0 EMAIL PASSWD"
  exit
fi

OBSDIR=$HOME/share/DATA/OBS/prepbufr/

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

DATAURL="http://rda.ucar.edu/data/ds337.0/prep48h/"
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

  ${OPT2}/${cy}/prepbufr.gdas.${cy}${cm}${cd}.t00z.nr.48h
  ${OPT2}/${cy}/prepbufr.gdas.${cy}${cm}${cd}.t06z.nr.48h
  ${OPT2}/${cy}/prepbufr.gdas.${cy}${cm}${cd}.t12z.nr.48h
  ${OPT2}/${cy}/prepbufr.gdas.${cy}${cm}${cd}.t18z.nr.48h

  mv prepbufr.gdas.${cy}${cm}${cd}.t00z.nr.48h ${cy}${cm}${cd}00.prepbufr.nr
  mv prepbufr.gdas.${cy}${cm}${cd}.t06z.nr.48h ${cy}${cm}${cd}06.prepbufr.nr
  mv prepbufr.gdas.${cy}${cm}${cd}.t12z.nr.48h ${cy}${cm}${cd}12.prepbufr.nr
  mv prepbufr.gdas.${cy}${cm}${cd}.t18z.nr.48h ${cy}${cm}${cd}18.prepbufr.nr

echo $CDATE $DF
CDATE=`date_edit2 $CDATE $DF `


done

rm -f auth.rda_ucar_edu
echo "NORMAL END"
