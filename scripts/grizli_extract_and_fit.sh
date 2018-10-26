#!/bin/bash

# Syncs from working directory on EC2 instance to S3 buket

BUCKET=aws-grivam
BUCKET=grizli-grism

if [ $# -eq 0 ]
  then
    echo "Usage:  $ grizli_extract_and_fit j023507-040202 [min_mag,max_max]"
    echo ""
    echo "Available: "
    
    aws s3 ls s3://${BUCKET}/Pipeline/ |grep -v footprint |grep PRE | awk '{print $2}'
    
    exit
fi

root=$1

if [[ -z $2 ]]; then
    maglim="16.5,26"
else
    maglim=$2
fi

echo "Running on root=${root} with maglim=${maglim}"

# Clean up
# sudo mkdir /GrizliImaging  (spectra)
# sudo chmod ugoa+rwx /GrizliImaging

rm -rf /GrizliSpectra/*
cd /GrizliSpectra

echo "Start:   `date`" > ${root}.log
aws s3 cp ${root}.log s3://${BUCKET}/Pipeline/Log/Start/

# Copy from S3
aws s3 sync s3://${BUCKET}/Pipeline/${root}/Extractions/ ./
aws s3 cp s3://${BUCKET}/Pipeline/${root}_footprint.fits ./

## Extractions
grizli_extract_and_fit.py ${root} run ${maglim}

# Sync extractions
aws s3 sync --exclude "*" --include "${root}_*fits" --acl public-read ./ s3://${BUCKET}/Pipeline/${root}/Extractions/
aws s3 sync --exclude "*" --include "${root}*png" --include "*html" --acl public-read ./ s3://${BUCKET}/Pipeline/${root}/Extractions/

# Log
echo "# Beams extracted:   `date`" >> ${root}.log
aws s3 cp ${root}.log s3://${BUCKET}/Pipeline/Log/Extract/ --acl public-read

dfits *beams.fits | fitsort RA DEC ID COUNT > ${root}.beams.log
aws s3 cp ${root}.beams.log s3://${BUCKET}/Pipeline/Log/Beams/ --acl public-read

# Check for corrupt full.fits files
grizli_check_bad_full.sh

num_beams=`ls *beams.fits | wc -l`
echo "${root} N=${num_beams}"

# First try with AWS Lambda
# fit_redshift_lambda.py ${root} False
# for iter in 1 2 3; do
#     fit_redshift_lambda.py ${root} False
# done
# 
# # Sync fits from lambda
# aws s3 sync s3://${BUCKET}/Pipeline/${root}/Extractions/ ./ --acl public-read

# Redshift fits, those that weren't fit by lambda
cpu_count=`python -c 'import multiprocessing; print(multiprocessing.cpu_count()//2)'`
mpiexec -n $cpu_count python -m mpi4py.futures $GRIZLICODE/grizli/pipeline/run_MPI.py


# Check for corrupt full.fits files
grizli_check_bad_full.sh

# Summary HTML & zhist figure
# root=`ls *phot.fits | sed "s/_/ /" | awk '{print $1}'`
grizli_extract_and_fit.py ${root} summary

#### Sync final
# PNG, HTML, REG - public-read
aws s3 sync --exclude "*" --include "${root}*png" --include "${root}*reg" --include "*html" --acl public-read ./ s3://${BUCKET}/Pipeline/${root}/Extractions/

# FITS files #- protected
aws s3 sync --exclude "*" --include "${root}*fits" --acl public-read ./ s3://${BUCKET}/Pipeline/${root}/Extractions/

# sudo shutdown now

## Log Finished
echo "Finished:   `date`" >> ${root}.log
aws s3 cp ${root}.log s3://${BUCKET}/Pipeline/Log/Finished/


