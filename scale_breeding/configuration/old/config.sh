NPROCS=4   #Number of procs that will be used to run the model (should be consistente with scale namelist options)

DATABASEROOT="/home/jruiz/share/database_v170124/"

OUTPUTDIR="/home/jruiz/share/exp/${CONFIGURATION}/" #Output root directory where results will be stored.

TMPDIR="${SRCDIR}/../tmp/${CONFIGURATION}"   #TMP directory where computations will take place.

INIDATE=20130713000000   #Start date of the run.
ENDDATE=20130713010000   #End date of the run.


BDYDATAFREQ=30            #Boundary data frequency in seconds.
BDYDATAROOT="/home/jruiz/share/scale_input_data/wrf_osaka_1km_2013_07_13/"
BDYNFILES=5

FLENGTH=120     #Forecast length in seconds.
FOUTPUTFREQ=30  #Forecast output freq.





