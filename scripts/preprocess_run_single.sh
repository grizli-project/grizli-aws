#!/bin/bash

# Syncs from working directory on EC2 instance to S3 buket

if [ $# -eq 0 ]
  then
    echo "Usage:  $ imaging_run_single j023507-040202"
    echo ""
    echo "Available: "
    
    aws s3 ls s3://aws-grivam/Pipeline/ |grep -v footprint |grep PRE | awk '{print $2}'
    
    exit
fi

root=$1

echo "Running on root=${root}"

# Initialize working directory
rm -rf /GrizliImaging/*
#sudo mkdir /GrizliImaging
#sudo chmod ugoa+rwx /GrizliImaging
cd /GrizliImaging

echo "Start:   `date`" > ${root}.log
aws s3 cp ${root}.log s3://grizli-preprocess/Pipeline/Log/Start/

# Copy from S3
aws s3 cp s3://grizli-preprocess/Pipeline/Fields/${root}_footprint.fits ./

## Extractions
preprocess_run_single.py ${root} 

cd /GrizliImaging

rm ${root}/Prep/*wcs-ref.fits
rm ${root}/Prep/*bkg.fits
rm ${root}/Prep/astrodrizzle.log

# Sync extractions
aws s3 sync --exclude "*" --include "${root}/Prep/*fl?.fits" \
                          --include "${root}/Prep/*.log" \
                          --include "${root}/Prep/*fail*" \
                          --include "${root}/Prep/*visits.npy" \
                          --include "${root}/Prep/*_fine*" \
                          --include "${root}/Prep/*png" \
                          --include "${root}/Prep/*_dr?_sci.fits" \
                          --include "${root}/Prep/*cat.fits" \
                          --include "${root}/Prep/*reg" \
                          --include "${root}/Prep/*radec" \
                          --acl public-read \
                          ./ s3://grizli-preprocess/Pipeline/

#if [ -e ${root}/Prep/${root}_fine.png ]; then 
failed=`ls ${root}/Prep/ |grep fail`
if [ -z "$failed" ] ; then 
    echo "${root}: Success" 
    echo "Finished:   `date`" > ${root}.log
    aws s3 cp ${root}.log s3://grizli-preprocess/Pipeline/Log/Finished/
    aws s3 rm s3://grizli-preprocess/Pipeline/Log/Failed/${root}.log
else
    echo "${root}: Fail..."
    echo "Failed:   `date`" > ${root}.log
    aws s3 cp ${root}.log s3://grizli-preprocess/Pipeline/Log/Failed/
    aws s3 rm s3://grizli-preprocess/Pipeline/Log/Finished/${root}.log
fi

# Done

cd /GrizliImaging
