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

# Clean up
cd $HOME/GrizliExtract
rm -rf $HOME/GrizliExtract/*

echo "Start:   `date`" > ${root}.log
aws s3 cp ${root}.log s3://grizli-imaging/Pipeline/Log/Start/

# Copy from S3
aws s3 cp s3://grizli-imaging/Pipeline/Fields/${root}_footprint.fits ./

## Extractions
imaging_run_single.py ${root} 

rm ${root}/Prep/*wcs-ref.fits
rm ${root}/Prep/astrodrizzle.log

# Sync extractions
aws s3 sync --exclude "*" --include "${root}/Prep/${root}*" --include "${root}/Prep/*.log" --include "${root}/Prep/*fine*" --acl public-read ./ s3://grizli-imaging/Pipeline/${root}/

# Done
echo "Finished:   `date`" > ${root}.log
aws s3 cp ${root}.log s3://grizli-imaging/Pipeline/Log/Finished/

cd $HOME/GrizliExtract
