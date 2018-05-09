#!/bin/bash

# Syncs from working directory on EC2 instance to S3 buket

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

# Clean up
cd $HOME/GrizliExtract
rm *GrismFL* *npy j[0-9]* *wcs*

# Copy from S3
aws s3 sync s3://aws-grivam/Pipeline/${root}/Extractions/ ./
aws s3 cp s3://aws-grivam/Pipeline/${root}_footprint.fits ./

## Extractions
grizli_extract_and_fit.py ${root} run

# Sync extractions
aws s3 sync --exclude "*" --include "${root}_*fits" ./ s3://aws-grivam/Pipeline/${root}/Extractions/
aws s3 sync --exclude "*" --include "${root}*png" --include "*html" --acl public-read ./ s3://aws-grivam/Pipeline/${root}/Extractions/

# Redshift fits
cpu_count=`python -c 'import multiprocessing; print(multiprocessing.cpu_count()//2)'`
mpiexec -n $cpu_count python -m mpi4py.futures $GRIZLICODE/grizli/pipeline/run_MPI.py

# Summary HTML & zhist figure
grizli_extract_and_fit.py ${root} summary

# Sync final
aws s3 sync --exclude "*" --include "${root}*fits" ./ s3://aws-grivam/Pipeline/${root}/Extractions/
aws s3 sync --exclude "*" --include "${root}*png" --include "${root}*reg" --include "*html" --acl public-read ./ s3://aws-grivam/Pipeline/${root}/Extractions/


