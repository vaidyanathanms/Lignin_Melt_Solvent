ID=615462
for i in {1..4}
do
    ID=`bsub -w "ended($ID)" run_md.sh| awk '{print $2}'| sed  's/>//'| sed 's/<//' `
done
