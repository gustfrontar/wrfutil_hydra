#Generate a random number to save the current configuration of the working script.
script=$1

min=50
max=1000
RNAME=$((RANDOM % (max - min + 1) + min))

#Rename the main script to a random name.
cp $script ./work_$RNAME.sh 

#Submit the copy to the queue.
sbatch ./work_$RNAME.sh
