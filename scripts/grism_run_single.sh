#!/bin/bash

# Syncs from working directory on EC2 instance to S3 buket

#BUCKET=aws-grivam
#BUCKET=grizli-grism
BUCKET=grizli

if [ $# -eq 0 ]
  then
    echo "Usage:  $ grism_run_single j023507-040202"
    echo ""
    echo "Available: "
    
    #aws s3 ls s3://aws-grivam/Pipeline/ |grep -v footprint |grep PRE | awk '{print $2}'
    
    exit
fi

root=$1

echo "Running on root=${root}"

# Initialize working directory
rm -rf /GrizliImaging/${root}*
#mkdir /GrizliImaging
#chmod ugoa+rwx /GrizliImaging
cd /GrizliImaging

echo "Start:   `date`" > ${root}.log
aws s3 cp ${root}.log s3://${BUCKET}/Pipeline/Log/Start/
aws s3 rm s3://${BUCKET}/Pipeline/Log/Failed/${root}.log

# Copy from S3
aws s3 cp s3://${BUCKET}/Pipeline/Fields/${root}_footprint.fits ./
aws s3 cp s3://${BUCKET}/Pipeline/Fields/${root}_master.radec ./
aws s3 cp s3://${BUCKET}/Pipeline/Fields/${root}_parent.radec ./

aws s3 sync s3://${BUCKET}/Pipeline/${root}/ ./${root}/
rm -rf ./${root}/Prep/*fail*
aws s3 rm --recursive --exclude "*" --include "*fail*" s3://${BUCKET}/Pipeline/${root}/Prep/

## Extractions
grism_run_single.py ${root} 

rm ${root}/Prep/*wcs-ref.fits
rm ${root}/Prep/*bkg.fits
rm ${root}/Prep/*ctx.fits
rm ${root}/Prep/astrodrizzle.log
rm -rf ${root}/Prep/FineBkup

cp ${root}.auto_script.log ./${root}/Prep/

# Sync extractions
#aws s3 sync --exclude "*" --include "${root}/Prep/${root}*" --include "${root}/Prep/*flt.fits" --include "${root}/Prep/*.log" --include "${root}/Prep/*fine*" --include "${root}/Extractions/*" --acl public-read ./ s3://${BUCKET}/Pipeline/

aws s3 sync --exclude "*" --include "Prep/${root}*"   \
                          --include "Prep/[ij]*_fl?.fits"  \
                          --include "Prep/u*_c??.fits" \
                          --include "Prep/*fail*"     \
                          --include "Prep/*.reg"      \
                          --include "Prep/*.log"      \
                          --include "Prep/*.png"      \
                          --include "Prep/*.radec"    \
                          --include "Prep/*fine*"     \
                          --include "RAW/*.png"       \
                          --include "Extractions/*"   \
                          --exclude "Prep/FineBkup/*" \
                          --acl public-read \
                          ./${root} s3://${BUCKET}/Pipeline/${root}

if [ -e ${root}/Prep/${root}_phot.fits ]; then 
    echo "${root}: Success"
    echo "Finished:   `date`" > ${root}.log
    aws s3 cp ${root}.log s3://${BUCKET}/Pipeline/Log/Finished/
    aws s3 rm s3://${BUCKET}/Pipeline/Log/Failed/${root}.log
    
    rm -rf /GrizliImaging/${root}*
else
    echo "${root}: Fail..."
    echo "Failed:   `date`" > ${root}.log
    aws s3 cp ${root}.log s3://${BUCKET}/Pipeline/Log/Failed/
    aws s3 rm s3://${BUCKET}/Pipeline/Log/Finished/${root}.log
fi

# Done

cd /GrizliImaging
