#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "Usage:  $ grizli_extract_and_fit j023507-040202"
    echo ""
    echo "Available: "
    
    aws s3 ls s3://aws-grivam/Pipeline/ |grep -v footprint | awk '{print $2}'
    
    exit
fi

root=$1

echo "Running on root=${root}"

cd $HOME/GrizliExtract
rm *GrismFL* *npy ${root}*

aws s3 sync s3://aws-grivam/Pipeline/${root}/Extractions/ ./
aws s3 cp s3://aws-grivam/Pipeline/${root}_footprint.fits ./

grizli_extract_and_fit.py ${root} run

cpu_count=`python -c 'import multiprocessing; print(multiprocessing.cpu_count()//2)'`
mpiexec -n $cpu_count python -m mpi4py.futures $GRIZLICODE/grizli/pipeline/run_MPI.py

grizli_extract_and_fit.py ${root} summary

aws s3 sync --exclude "*" --include "${root}_*fits" ./ s3://aws-grivam/Pipeline/${root}/Extractions/
aws s3 sync --exclude "*" --include "${root}_*png" --include "*html" --acl public-read ./ s3://aws-grivam/Pipeline/${root}/Extractions/

