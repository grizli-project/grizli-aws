#!/bin/bash

# grism_run_single.sh ${root} --run_extractions=True --extract_args.maglim=[17,21] --include_photometry_in_fit=True --noclean
# gunzip ${root}/*/*gz

# Dummy loop for uncommented examples
if [ -n "${xxx}" ]; then

    # Initial processing
    grism_run_single.sh ${root} --run_extractions=True --extract_args.maglim=[17,21] --include_photometry_in_fit=True --noclean
    gunzip ${root}/*/*fits.gz

    # Rerun and extract more sources
    grism_run_single.sh ${root} --grism --run_extractions=True --extract_args.maglim=[16,26] --include_photometry_in_fit=True --noclean
    
    # Individual ID
    root=j022204m0412
    BUCKET=grizli
    fit_redshift_lambda.py ${root} --bucket_name=${BUCKET} --newfunc=False --skip_existing=True --sleep=True --ids=235 --zr=0.1,9
    
    root=j123624p6214
    grism_run_single.sh ${root} --grism --run_extractions=True --extract_args.maglim=[17,21] --extract_args.ids=[2728] --include_photometry_in_fit=True --noclean
    
fi

# Syncs from working directory on EC2 instance to S3 buket

#BUCKET=aws-grivam
#BUCKET=grizli-grism
BUCKET=grizli

if [ $# -eq 0 ]; then
    echo "Usage:  $ grism_run_single j023507-040202"
    echo ""
    echo "Available: "
    
    #aws s3 ls s3://aws-grivam/Pipeline/ |grep -v footprint |grep PRE | awk '{print $2}'
    
    exit
fi

root=$1

echo "Running on root=${root}"

clean=1
for i in "$@" ; do
    if [[ $i == "--grism" ]] ; then
        is_grism=1
    elif [[ $i == "--sync" ]] ; then
        is_sync=1
    elif [[ $i == "--noclean" ]] ; then
        clean=0
    fi
done

# Initialize working directory
if [ $clean -gt 0 ]; then
    rm -rf /GrizliImaging/${root}*
fi

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
        
if [ -n "${is_sync}" ] || [ -n "${is_grism}" ]; then
    #aws s3 sync --exclude "*" --include "Prep/${root}*sci.fits.gz" --include "Prep/${root}-ir*" --include "Prep/${root}*phot.fits" --include "Prep/${root}*psf.fits*" s3://${BUCKET}/Pipeline/${root}/ ./${root}/
    if [ -n "${is_sync}" ]; then 
        aws s3 sync s3://${BUCKET}/Pipeline/${root}/ ./${root}/ --include "*" --exclude "Extractions/*"
    else
        aws s3 sync s3://${BUCKET}/Pipeline/${root}/ ./${root}/
    fi
    
    # Make directories
    
    if [ ! -e "${root}/Persistence" ]; then
         mkdir ${root}/Persistence
    fi
        
    if [ ! -e "${root}/Extractions" ]; then
         mkdir ${root}/Extractions
    fi
    
    # If sync then clean previous files so that they'll be generated again
    if [ -n "${is_sync}" ]; then 
        rm ./${root}/Prep/${root}_phot.fits
        rm ./${root}/Prep/${root}-ir*
        rm ./${root}/Prep/${root}-*psf*
        rm ./${root}/Prep/*visits.npy
        rm ./${root}/Extractions/*
        
        aws s3 rm --recursive --exclude "*" --include "Prep/${root}_phot.fits" --include "Prep/${root}-ir*" --include "Prep/${root}-*psf*" --include "Prep/*visits.npy" --include "Extractions/*" s3://${BUCKET}/Pipeline/${root}
        
    fi 
    
    # Unzip zipped mosaics
    gunzip ${root}/Prep/*fits.gz
    
    # Symlinks to force skip already complete
    files=`ls ./${root}/Prep/*wcs.log | sed "s/_wcs.log/_drz_sci.fits/"`
    for file in $files; do 
        echo $file
        touch $file
    done
    
    files=`ls ./${root}/Prep/*_column.png | sed "s/_column.png/_drz_sci.fits/"`
    for file in $files; do 
        echo $file
        touch $file
    done
    
    # Make fake copies of flt-raw files
    cd ${root}/Prep/
    files=`ls *_flt.fits`
    cd ../RAW/
    for file in $files; do 
        out=`echo $file | sed "s/_flt/_raw/"`
        echo $file $out
        ln -s ../Prep/$file $out
        cp ../Prep/${file} .
    done
    cp ../Prep/*flc.fits .
    cd ../../
    
fi

# aws s3 sync s3://${BUCKET}/Pipeline/${root} ./
# rm -rf ./${root}/Prep/*fail*
aws s3 rm --recursive --exclude "*" --include "*fail*" s3://${BUCKET}/Pipeline/${root}/Prep/

## Extractions
grism_run_single.py $@ #${root} 

rm ${root}/Prep/*wcs-ref.fits
rm ${root}/Prep/*bkg.fits
rm ${root}/Prep/*ctx.fits
rm ${root}/Prep/astrodrizzle.log
rm -rf ${root}/Prep/FineBkup

cp ${root}.auto_script.log ./${root}/Prep/
cp ${root}*yml ./${root}/Prep/

echo "gzip mosaics"

gzip ${root}/Prep/${root}*_dr?_*fits
gzip ${root}/Prep/${root}*_seg.fits
gzip ${root}/Extractions/*grism*fits

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
                          --include "RAW/*info"       \
                          --include "Extractions/*"   \
                          --exclude "Prep/FineBkup/*" \
                          --acl public-read \
                          ./${root} s3://${BUCKET}/Pipeline/${root}

# Do we have beams files and/or redshift fit outputs?
nbeams=`ls ${root}/Extractions/ |grep -c beams.fits`
nfull=`ls ${root}/Extractions/ |grep -c full.fits`

if [ $nbeams -gt 0 ]; then
        
    # Run the lambda function
    if [ $nbeams -ne $nfull ]; then 
        echo "Run redshift fit (nbeams=${nbeams}, nfull=${nfull})"
        fit_redshift_lambda.py ${root} --bucket_name=${BUCKET} --newfunc=False --skip_existing=True --sleep=True
    else
        echo "Make redshift catalog (nbeams=${nbeams}, nfull=${nfull})"
    fi
    
    # Generate catalog
    cd ${root}/Extractions/
    
    aws s3 sync --exclude "*" --include "${root}*full.fits" --include "${root}*stack.png" --include "*cat.fits" --include "${root}*info.fits" --acl public-read s3://${BUCKET}/Pipeline/${root}/Extractions/ ./
    
    grizli_extract_and_fit.py ${root} summary
    
    aws s3 sync --exclude "*" --include "${root}*fits" ./ s3://${BUCKET}/Pipeline/${root}/Extractions/ --acl public-read
    aws s3 sync --exclude "*" --include "${root}*png" --include "${root}*reg" --include "*html" --acl public-read ./ s3://${BUCKET}/Pipeline/${root}/Extractions/
    
    # Back to initial directory
    cd ../../
    
fi

if [ -e /tmp/${root}.success ]; then 
    echo "${root}: Success"
    echo "Finished:   `date`" > ${root}.log
    aws s3 cp ${root}.log s3://${BUCKET}/Pipeline/Log/Finished/
    aws s3 rm s3://${BUCKET}/Pipeline/Log/Failed/${root}.log
    
    if [ $clean -gt 0 ]; then
        rm -rf /GrizliImaging/${root}*
    fi
    
else
    echo "${root}: Fail..."
    echo "Failed:   `date`" > ${root}.log
    aws s3 cp ${root}.log s3://${BUCKET}/Pipeline/Log/Failed/
    aws s3 rm s3://${BUCKET}/Pipeline/Log/Finished/${root}.log
fi

# Done

cd /GrizliImaging
