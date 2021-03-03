ID=699492
for i in {1..3}
do
    ID=`bsub -w "ended($ID)" runall_batch.sh| awk '{print $2}'| sed  's/>//'| sed 's/<//' `
done
