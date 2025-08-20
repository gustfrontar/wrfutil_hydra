NSUB=30

#job1_id=$(sbatch --parsable anal_post.sh)
count=1
while [ $count -le $NSUB ]; do
  echo "Count: $count"
#  job1_id=$(sbatch --dependency=afterany:$job1_id --parsable anal_post.sh)
  count=$((count + 1))
done

#job1_id=$(sbatch --parsable gues_post.sh)
count=1
while [ $count -le $NSUB ]; do
  echo "Count: $count"
#  job1_id=$(sbatch --dependency=afterany:$job1_id --parsable gues_post.sh)
  count=$((count + 1))
done

job1_id=$(sbatch --parsable dafcs_post.sh)
count=1
while [ $count -le $NSUB ]; do
  echo "Count: $count"
  job1_id=$(sbatch --dependency=afterany:$job1_id --parsable dafcs_post.sh)
  count=$((count + 1))
done




