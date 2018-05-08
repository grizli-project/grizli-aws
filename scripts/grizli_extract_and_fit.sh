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

aws s3 sync s3://aws-grivam/Pipeline/${root}/Extractions/ ./
aws s3 cp s3://aws-grivam/Pipeline/${root}_footprint.fits ./

grizli_extract_and_fit.py ${root} run

mpiexec -n 16 python -m mpi4py.futures $GRIZLICODE/grizli/pipeline/run_MPI.py

grizli_extract_and_fit.py ${root} summary

aws s3 sync --exclude "*" --include "${root}_*fits" ./ s3://aws-grivam/Pipeline/${root}/Extractions/
aws s3 sync --exclude "*" --include "${root}_*png" --include "*html" --acl public-read ./ s3://aws-grivam/Pipeline/${root}/Extractions/

