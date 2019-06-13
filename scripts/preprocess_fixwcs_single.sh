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

bucket=grizli-preprocess

echo "Running on root=${root}"

# Initialize working directory
rm -rf /GrizliImaging/${root}*
#sudo mkdir /GrizliImaging
#sudo chmod ugoa+rwx /GrizliImaging
cd /GrizliImaging

echo "Start:   `date`" > ${root}.log
aws s3 cp ${root}.log s3://${bucket}/FixWCS/Log/Start/

# Copy from S3
aws s3 cp s3://${bucket}/Pipeline/Fields/${root}_footprint.fits ./
aws s3 cp s3://${bucket}/Pipeline/Fields/${root}_master.radec ./
aws s3 cp s3://${bucket}/Pipeline/Fields/${root}_parent.radec ./

# For COSMOS
aws s3 cp s3://${bucket}/hsc-udeep-i25_corr_cosmos.radec ./
# For UDS
aws s3 cp s3://${bucket}/hsc-udeep-sxds_corr_uds.radec ./

# Sync entire directory
mkdir ${root}
aws s3 sync --exclude "Prep/FineBkup*" s3://${bucket}/Pipeline/${root} ./${root}
rm ${root}/Prep/*fail*
rm -rf ${root}/Prep/*fine.* ${root}/Prep/FineBk*

## Run the fix script
cd /GrizliImaging/${root}/Prep/
preprocess_fixwcs_single.py ${root} 

cd /GrizliImaging

rm ${root}/Prep/*wcs-ref.fits
rm ${root}/Prep/*bkg.fits
rm ${root}/Prep/astrodrizzle.log
rm ${root}/Prep/stwcs.log

# Remove old failed files
aws s3 rm --recursive --exclude "*" --include "*failed" s3://${bucket}/Pipeline/${root}/

#if [ -e ${root}/Prep/${root}_fine.png ]; then 
failed=`ls ${root}/Prep/ |grep fail`
if [ -z "$failed" ] ; then 
    # Sync extractions
    aws s3 sync --exclude "*" --include "Prep/[ij]*_fl?.fits" \
                              --include "Prep/u*_c??.fits" \
                              --include "Prep/*.log" \
                              --include "Prep/*wcs.fits" \
                              --include "Prep/*fail*" \
                              --include "Prep/*visits.npy" \
                              --include "Prep/*_fine*" \
                              --include "Prep/*png" \
                              --include "Prep/*_dr?_sci.fits" \
                              --include "Prep/*cat.fits" \
                              --include "Prep/*reg" \
                              --include "Prep/*radec" \
                              --acl public-read \
                              ./${root} s3://${bucket}/Pipeline/${root}

    echo "${root}: Success" 
    echo "Finished:   `date`" > ${root}.log
    aws s3 cp ${root}.log s3://${bucket}/FixWCS/Log/Finished/
    aws s3 rm s3://${bucket}/FixWCS/Log/Failed/${root}.log
    
    rm -rf /GrizliImaging/${root}*
    
else
    echo "${root}: Fail..."
    echo "Failed:   `date`" > ${root}.log
    aws s3 cp ${root}.log s3://${bucket}/FixWCS/Log/Failed/
    aws s3 rm s3://${bucket}/FixWCS/Log/Finished/${root}.log
    
    # Remove fine files
    aws s3 rm --recursive --exclude "*" --include "*fine.*" --include "FineBkup*" s3://${bucket}/Pipeline/${root}/Prep/
    
fi

# Done

cd /GrizliImaging
