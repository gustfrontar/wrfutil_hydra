#!/bin/bash
# Runing on LINUX

#
# DATA SOURCE
#
DATASRC=2 # 1:NOAA/NOMADS , 2:UCAR/DSS
export http_proxy="proxy.fcen.uba.ar:8080"

# Initial time
IDATE=20100801000000
# Final time
EDATE=20101101000000
# Frequency
DF=21600   #file frequency in seconds

if [ $DATASRC -eq 2 ] ; then
 EMAIL=$1
 PASSWD=$2
 if [ $# -ne 2 ] ; then
  echo "USAGE: $0 EMAIL PASSWD"
  exit
 fi
fi

EMAIL=$1
PASSWD=$2

OBSNAME=PREPBUFRSA

OBSDIR=$HOME/share/OBS/$OBSNAME/prepbufr

source $HOME/share/LETKF_WRF/wrf/run/util.sh #Load time functions.


#
# WKDIR
#

mkdir -p $OBSDIR
cd $OBSDIR

#
# Settings for DATASRC
#
if test $DATASRC -eq 1
then
DATAURL="http://nomads.ncdc.noaa.gov/data/gdas"
WGET="wget"
else
DATAURL="http://rda.ucar.edu/data/ds337.0/tarfiles/"
WGET="wget --no-check-certificate"
rm -f auth.dss_ucar_edu
OPT1='-O /dev/null --save-cookies auth.rda_ucar_edu --post-data'
OPT2="email=${EMAIL}&passwd=${PASSWD}&action=login"
$WGET $OPT1="$OPT2" https://rda.ucar.edu/cgi-bin/login
OPT="-N --load-cookies auth.rda_ucar_edu"
fi

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

if [ ! -f $cy$cm$cd$ch.prepbufr.nr ] ; then
  if [ $DATASRC -eq 1 ] ; then
    $WGET $DATAURL/$IY$IM/$IY$IM$ID/gdas1.t${ch}z.prepbufr.nr
    mv gdas1.t${ch}z.prepbufr.nr $cy$cm$cd$ch.prepbufr.nr
  else
   if [ $ch -eq 00 ] ; then
     $WGET $OPT $DATAURL/${cy}/prepbufr.${cy}${cm}${cd}.nr.tar.gz
     gunzip -f prepbufr.${cy}${cm}${cd}.nr.tar.gz
     tar -xvf prepbufr.${cy}${cm}${cd}.nr.tar
     rm -fr prepbufr.${cy}${cm}${cd}.nr.tar
   fi
   mv ${cy}${cm}${cd}.nr/prepbufr.gdas.${cy}${cm}${cd}.t${ch}z.nr ${cy}${cm}${cd}${ch}.prepbufr.nr
 fi
fi

echo $CDATE $DF
CDATE=`date_edit2 $CDATE $DF `


done

rm -f auth.dss_ucar_edu
echo "NORMAL END"
