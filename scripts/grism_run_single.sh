#!/bin/bash

# Syncs from working directory on EC2 instance to S3 buket

if [ $# -eq 0 ]
  then
    echo "Usage:  $ grism_run_single j023507-040202"
    echo ""
    echo "Available: "
    
    aws s3 ls s3://aws-grivam/Pipeline/ |grep -v footprint |grep PRE | awk '{print $2}'
    
    exit
fi

root=$1

echo "Running on root=${root}"

# Initialize working directory
rm -rf /GrizliImaging/*
#mkdir /GrizliImaging
#chmod ugoa+rwx /GrizliImaging
cd /GrizliImaging

echo "Start:   `date`" > ${root}.log
aws s3 cp ${root}.log s3://grizli-grism/Pipeline/Log/Start/
aws s3 rm s3://grizli-grism/Pipeline/Log/Failed/${root}.log

# Copy from S3
aws s3 cp s3://grizli-grism/Pipeline/Fields/${root}_footprint.fits ./
aws s3 cp s3://grizli-grism/Pipeline/Fields/${root}_master.radec ./

## Extractions
grism_run_single.py ${root} 

rm ${root}/Prep/*wcs-ref.fits
rm ${root}/Prep/*bkg.fits
rm ${root}/Prep/*ctx.fits
rm ${root}/Prep/astrodrizzle.log

# Sync extractions
#aws s3 sync --exclude "*" --include "${root}/Prep/${root}*" --include "${root}/Prep/*flt.fits" --include "${root}/Prep/*.log" --include "${root}/Prep/*fine*" --include "${root}/Extractions/*" --acl public-read ./ s3://grizli-grism/Pipeline/

aws s3 sync --exclude "*" --include "${root}/Prep/${root}*" --include "${root}/Prep/*flt.fits" --include "${root}/Prep/*.reg" --include "${root}/Prep/*.log" --include "${root}/Prep/*fine*" --include "${root}/Extractions/*" --exclude "${root}/Prep/FineBkup/*" --acl public-read ./ s3://grizli-grism/Pipeline/

if [ -e ${root}/Prep/${root}_phot.fits ]; then 
    echo "${root}: Success"
    echo "Finished:   `date`" > ${root}.log
    aws s3 cp ${root}.log s3://grizli-grism/Pipeline/Log/Finished/
else
    echo "${root}: Fail..."
    echo "Failed:   `date`" > ${root}.log
    aws s3 cp ${root}.log s3://grizli-grism/Pipeline/Log/Failed/
fi

# Done

cd /GrizliImaging
