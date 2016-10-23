NFILES=`wc -l airslist_novdec2012.txt | awk '{print $1}'`
N=1
while test $N -le $NFILES
do
URL=`cat airslist_novdec2012.txt | head -n $N | tail -n 1`
FILE=`basename $URL`

echo $FILE | ./decoderconskip2
N=`expr $N + 1`
done


