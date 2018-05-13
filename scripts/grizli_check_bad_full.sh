#!/bin/bash 

# Check for corrupt "full.fits" files that will break the summary script
#

files=`ls *full.fits`
for file in $files; do
    log=`dfits -x 0 $file | fitsort EXTNAME | grep DSCI`
    if [[ -z $log ]]; then
        rmfile=`echo $file |sed "s/.fits/*/"`
        echo "rm $rmfile    # Missing extensions"
        rm $rmfile
    fi
done

    