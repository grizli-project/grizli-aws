#!/bin/bash

# grism_run_single.sh ${root} --run_extractions=True --extract_args.maglim=[17,21] --include_photometry_in_fit=True --noclean
# gunzip ${root}/*/*gz

# Syncs from working directory on EC2 instance to S3 buket

if [ $# -eq 0 ]; then
    echo "Usage:  $ grism_run_single j023507-040202"
    echo ""
    echo "Available: "
    
    #aws s3 ls s3://aws-grivam/Pipeline/ |grep -v footprint |grep PRE | awk '{print $2}'
    
    exit
fi

root=$1

echo "Running IRAC pipeline on root=${root}"

# Args
is_sync=0
clean=1

clean_full_fits=0

#BUCKET=aws-grivam
#BUCKET=grizli-grism
BUCKET=grizli-v1

for arg in "$@" ; do
    if [[ $arg == "--sync" ]] ; then
        is_sync=1
    elif [[ $arg == "--noclean" ]] ; then
        clean=0
    elif [[ ! -z `echo "${arg}" | grep "bucket="` ]]; then 
        BUCKET=`echo "${arg}" | sed "s/=/ /" | awk '{print $2}'`
    fi
done

echo "# ${root} is_sync=${is_sync}, clean=${clean}, BUCKET=${BUCKET}"

# Initialize working directory
if [ $clean -gt 0 ]; then
    rm -rf /GrizliImaging/${root}*
fi

#mkdir /GrizliImaging
#chmod ugoa+rwx /GrizliImaging
cd /GrizliImaging

fp_status=`aws s3 cp s3://${BUCKET}/Pipeline/IRAC/${root}_ipac.fits ./`
if [ -z "${fp_status}" ]; then 
    echo "Footprint s3://${BUCKET}/Pipeline/IRAC/${root}_ipac.fits not found"
    exit
fi

echo "Start:   `date`" > ${root}.log
aws s3 cp ${root}.log s3://${BUCKET}/IRAC/Log/Start/
aws s3 rm s3://${BUCKET}/IRAC/Log/Failed/${root}.log

# Copy from S3
aws s3 cp s3://${BUCKET}/IRAC/${root}_master.radec ./
aws s3 cp s3://${BUCKET}/IRAC/${root}_parent.radec ./
        
if [ $is_sync -gt 0 ]; then
    #aws s3 sync --exclude "*" --include "Prep/${root}*sci.fits.gz" --include 
    aws s3 sync s3://${BUCKET}/Pipeline/${root}/IRAC/ ./${root}/ --include "*"
    
    # Unzip zipped mosaics
    echo "gunzip files"
    gunzip --force ${root}/*fits.gz
        
fi

# aws s3 sync s3://${BUCKET}/Pipeline/${root} ./
# rm -rf ./${root}/Prep/*fail*
aws s3 rm --recursive --exclude "*" --include "*fail*" s3://${BUCKET}/Pipeline/${root}/IRAC/

## Run the pipeline
irac_run_single.py $@ #${root} 

if [ -e ${root}.success ]; then 
    echo "${root}: Success"
    echo "Finished:   `date`" > ${root}.log
    aws s3 cp ${root}.log s3://${BUCKET}/IRAC/Log/Finished/
    aws s3 rm s3://${BUCKET}/IRAC/Log/Failed/${root}.log
    
    if [ $clean -gt 0 ]; then
        rm -rf /GrizliImaging/${root}*
    fi
    
else
    echo "${root}: Fail..."
    echo "Failed:   `date`" > ${root}.log
    aws s3 cp ${root}.log s3://${BUCKET}/IRAC/Log/Failed/
    aws s3 rm s3://${BUCKET}/IRAC/Log/Finished/${root}.log
    
    # Sync what's there
    aws s3 sync /GrizliImaging/${root}/ s3://grizli-v1/Pipeline/${root}/IRAC/ --exclude "*" --include "${root}-ch*drz*fits*" --include "{root}.*png" --include "*-ch*psf*" --include "*log.fits" --include "*wcs.[lp]*" --include "*html" --include "*fail*" --acl public-read
    
fi

# Done

cd /GrizliImaging
