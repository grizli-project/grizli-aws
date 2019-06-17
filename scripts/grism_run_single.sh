#!/bin/bash

# grism_run_single.sh ${root} --run_extractions=True --extract_args.maglim=[17,21] --include_photometry_in_fit=True --noclean
# gunzip ${root}/*/*gz

# Dummy loop for uncommented examples
if [ -n "${xxx}" ]; then

    # Initial processing
    grism_run_single.sh ${root} --run_extractions=True --extract_args.maglim=[17,21] --include_photometry_in_fit=False --noclean
    gunzip ${root}/*/*fits.gz

    # Refine based on initial fits
    grism_run_single.sh ${root} --run_extractions=True --extract_args.maglim=[17,21] --include_photometry_in_fit=False --noclean --refine_with_fits=True --grism
    gunzip ${root}/*/*fits.gz
    
    # Rerun and extract more sources
    grism_run_single.sh ${root} --grism --run_extractions=True --extract_args.maglim=[16,25] --include_photometry_in_fit=False --noclean
    
    # Quasar templates
    ids=303
    BUCKET=grizli
    fit_redshift_lambda.py ${root} --bucket_name=${BUCKET} --newfunc=False --skip_existing=False --sleep=True --ids=${ids} --quasar_fit=True --check_wcs=True --output_path=self --use_psf=False --verbose=True
    
    # Individual ID
    root=j022204m0412
    BUCKET=grizli
    ids=235
    fit_redshift_lambda.py ${root} --bucket_name=${BUCKET} --newfunc=False --skip_existing=True --sleep=True --ids=${ids} --zr=0.1,9

    BUCKET=grizli
    zr=0.1,12
    zr=4,8
    fit_redshift_lambda.py ${root} --bucket_name=${BUCKET} --newfunc=False --skip_existing=True --sleep=True --ids=${ids} --zr=${zr} --verbose=True --check_wcs=True
    aws s3 ls s3://${BUCKET}/Pipeline/${root}/Extractions/ |grep full.fits | sort -k 1 -k 2
    
    root=j123624p6214
    grism_run_single.sh ${root} --grism --run_extractions=True --extract_args.maglim=[17,21] --include_photometry_in_fit=True --noclean --parse_visits_args.combine_minexp=1 --extract_args.ids=[2728]

    # Redo with F160W as reference
    grism_run_single.sh ${root} --sync --run_extractions=True --extract_args.maglim=[17,21] --include_photometry_in_fit=True --noclean --parse_visits_args.combine_minexp=1 --extract_args.ids=[2728] --grism_prep_args.gris_ref_filters.G141=[F160W]

    grism_run_single.sh ${root} --grism --run_extractions=True --extract_args.maglim=[17,23] --include_photometry_in_fit=True --noclean --parse_visits_args.combine_minexp=1  --grism_prep_args.gris_ref_filters.G141=[F160W]
    
fi

# Syncs from working directory on EC2 instance to S3 buket

if [ $# -eq 0 ]; then
    echo "Usage:  $ grism_run_single j023507-040202"
    echo ""
    echo "Available: "
    
    #aws s3 ls s3://aws-grivam/Pipeline/ |grep -v footprint |grep PRE | awk '{print $2}'
    
    exit
fi

root=$1

echo "Running on root=${root}"

# Args
is_grism=0
is_sync=0
clean=1
lambda_verbose="False"

clean_full_fits=0

#BUCKET=aws-grivam
#BUCKET=grizli-grism
BUCKET=grizli

for arg in "$@" ; do
    if [[ $arg == "--grism" ]] ; then
        is_grism=1
    elif [[ $arg == "--sync" ]] ; then
        is_sync=1
    elif [[ $arg == "--noclean" ]] ; then
        clean=0
    elif [[ $arg == "--lambda_verbose" ]] ; then
        lambda_verbose="True"
    elif [[ $arg == "--refine_with_fits=True" ]] ; then
        clean_full_fits=1
    elif [[ ! -z `echo "${arg}" | grep "bucket="` ]]; then 
        BUCKET=`echo "${arg}" | sed "s/=/ /" | awk '{print $2}'`
    fi
done

echo "# ${root} is_grism=${is_grism}, is_sync=${is_sync}, clean=${clean}, BUCKET=${BUCKET}"

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
        
if [ $is_sync -gt 0 ] || [ $is_grism -gt 0 ]; then
    #aws s3 sync --exclude "*" --include "Prep/${root}*sci.fits.gz" --include "Prep/${root}-ir*" --include "Prep/${root}*phot.fits" --include "Prep/${root}*psf.fits*" s3://${BUCKET}/Pipeline/${root}/ ./${root}/
    if [ $is_sync -gt 0 ]; then 
        aws s3 sync s3://${BUCKET}/Pipeline/${root}/ ./${root}/ --include "*" --exclude "Extractions/*" --exclude "Thumbnails/*fits"
    else
        aws s3 sync s3://${BUCKET}/Pipeline/${root}/ ./${root}/ --include "*" --exclude "Extractions/*1D*" --exclude "Extractions/*.R30*" --exclude "Extractions/*stack.fits" --exclude "Thumbnails/*fits"
    fi
    
    # Make directories
    
    if [ ! -e "${root}/Persistence" ]; then
         mkdir ${root}/Persistence
    fi
        
    if [ ! -e "${root}/Extractions" ]; then
         mkdir ${root}/Extractions
    fi
    
    # If sync then clean previous files so that they'll be generated again
    if [ $is_sync -gt 0 ]; then 
        rm ./${root}/Prep/${root}_phot.fits
        rm ./${root}/Prep/${root}-ir*
        rm ./${root}/Prep/${root}-*psf*
        rm ./${root}/Prep/*visits.npy
        rm ./${root}/Extractions/*
        
        aws s3 rm --recursive --exclude "*" --include "Prep/${root}_phot.fits" --include "Prep/${root}-ir*" --include "Prep/${root}-*psf*" --include "Prep/*visits.npy" --include "Extractions/*" s3://${BUCKET}/Pipeline/${root}        
    fi 
    
    # Unzip zipped mosaics
    echo "gunzip files"
    gunzip --force ${root}/*/*fits.gz
    
    # Symlinks to force skip already complete
    files=`ls ./${root}/Prep/*wcs.log | sed "s/_wcs.log/_drz_sci.fits/"`
    for file in $files; do 
        echo $file
        if [ ! -e ${file} ]; then
            touch $file
        fi
    done
    
    files=`ls ./${root}/Prep/*_column.png | sed "s/_column.png/_drz_sci.fits/"`
    for file in $files; do 
        echo $file
        if [ ! -e ${file} ]; then
            touch $file
        fi
    done
    
    # Make fake copies of flt-raw files
    cd ${root}/Prep/
    files=`ls *_flt.fits`
    cd ../RAW/
    for file in $files; do 
        out=`echo $file | sed "s/_flt/_raw/"`
        echo $file $out
        
        if [ ! -e ${out} ]; then 
            ln -s ../Prep/$file $out
        fi
        
        if [ ! -e "${file}" ]; then 
            cp ../Prep/${file} .
        fi
        
    done
    
    echo "Copy FLC"
    cp ../Prep/*flc.fits ../RAW/
    
    # Back to root
    cd ../../
    
fi

# aws s3 sync s3://${BUCKET}/Pipeline/${root} ./
# rm -rf ./${root}/Prep/*fail*
aws s3 rm --recursive --exclude "*" --include "*fail*" s3://${BUCKET}/Pipeline/${root}/Prep/

# If refining the model based on the fits, remove "beams" locally 
# and remove "full" on remote
if [ $clean_full_fits -gt 0 ]; then
    rm -rf /GrizliImaging/${root}/Extractions/*beams.fits
    aws s3 rm --recursive --exclude "*" --include "*full.fits" s3://${BUCKET}/Pipeline/${root}/Extractions/
fi

## Extractions
grism_run_single.py $@ #${root} 

# Now remove local full.fits if necessary so they won't be synced
if [ $clean_full_fits -gt 0 ]; then
    rm -rf /GrizliImaging/${root}/Extractions/*full.fits
fi

rm ${root}/Prep/*wcs-ref.fits
rm ${root}/Prep/*bkg.fits
rm ${root}/Prep/*ctx.fits
rm ${root}/Prep/astrodrizzle.log
rm -rf ${root}/Prep/FineBkup

cp ${root}.auto_script.log.txt ./${root}/Prep/
cp ${root}/RAW/*.txt ./${root}/Prep/
cp ${root}*yml ./${root}/Prep/

echo "gzip mosaics"

gzip --force ${root}/Prep/${root}*_dr?_*fits
gzip --force ${root}/Prep/${root}*_seg.fits
gzip --force ${root}/Extractions/*_grism*fits

# echo "gzip exposures"
# gzip --force ${root}/Prep/*_fl?.fits

# Copies of segmentation image
rm ${root}/Extractions/${root}*seg.fits
cp ${root}/Prep/${root}*seg.fits.gz ${root}/Extractions/

# Sync extractions
#aws s3 sync --exclude "*" --include "${root}/Prep/${root}*" --include "${root}/Prep/*flt.fits" --include "${root}/Prep/*.log" --include "${root}/Prep/*fine*" --include "${root}/Extractions/*" --acl public-read ./ s3://${BUCKET}/Pipeline/

aws s3 sync --exclude "*" --include "Prep/${root}*"   \
                          --include "Prep/[ij]*_fl?.*fits*"  \
                          --include "Prep/u*_c??.fits" \
                          --include "Prep/*fail*"     \
                          --include "Prep/*.reg"      \
                          --include "Prep/*.cat.fits"      \
                          --include "Prep/*.log"      \
                          --include "Prep/*.txt"      \
                          --include "Prep/*.png"      \
                          --include "Prep/*.radec"    \
                          --include "Prep/*fine*"     \
                          --include "RAW/*.png"       \
                          --include "RAW/*info"       \
                          --include "Thumbnails/*"   \
                          --include "Extractions/*"   \
                          --exclude "Prep/*_[0-9]*.thumb*"   \
                          --exclude "Prep/*_[0-9]*.rgb*"   \
                          --exclude "Extractions/templates*"   \
                          --exclude "Extractions/FILTER*"   \
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
        echo "fit_redshift_lambda.py ${root} --bucket_name=${BUCKET} --newfunc=False --skip_existing=True --sleep=True --verbose=${lambda_verbose}"
        fit_redshift_lambda.py ${root} --bucket_name=${BUCKET} --newfunc=False --skip_existing=True --sleep=True --verbose=${lambda_verbose} --check_wcs=True
    else
        echo "Make redshift catalog (nbeams=${nbeams}, nfull=${nfull})"
    fi
    
    # Generate catalog
    cd ${root}/Extractions/
    
    aws s3 sync --exclude "*" --include "${root}*full.fits" --include "${root}*stack.png" --include "*cat.fits" --include "${root}*info.fits" --acl public-read s3://${BUCKET}/Pipeline/${root}/Extractions/ ./
    
    #BUCKET=grizli
    grizli_extract_and_fit.py ${root} summary
    
    aws s3 sync --exclude "*" --include "${root}*fits" ./ s3://${BUCKET}/Pipeline/${root}/Extractions/ --acl public-read
    aws s3 sync --exclude "*" --include "${root}*png" --include "${root}*reg" --include "*html" --include "*json" --acl public-read ./ s3://${BUCKET}/Pipeline/${root}/Extractions/
    
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
