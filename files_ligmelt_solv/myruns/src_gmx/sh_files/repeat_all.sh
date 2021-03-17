#!/bin/bash 

# submit the preprocessing job (run_preprocess.sh)
ID=$(bsub run_preprocess.sh | awk '{print $2}' | tr "<" " " | tr ">" " " | awk '{print $1}')

# submit the MD run job in a loop (run_md.sh)
for ((i=1;i<=3;i+=1))
do 
ID=$(bsub -w "ended("$ID")" run_md.sh | awk '{print $2}' | tr "<" " " | tr ">" " " | awk '{print $1}')
done
