BUCKET=aws-grivam
BUCKET=grizli-grism

roots=`cat GrizliPrepPassed.log`
for root in $roots; do
    
    # Sync
    cat=`aws s3 ls s3://${BUCKET}/Pipeline/Log/Sync/${root}.log`
    
    # info.fits
    #cat=`aws s3 ls s3://${BUCKET}/Pipeline/${root}/Extractions/${root}.info.fits`
    
    if [[ -z $cat ]]; then
        echo "! Run $root" 
        
        #fi; done
        
        aws s3 cp ${root}_footprint.fits s3://${BUCKET}/Pipeline/
        aws s3 cp ${root}_footprint.pdf s3://${BUCKET}/Pipeline/
        convert -density 150  ${root}_footprint.pdf ${root}_footprint.png
        aws s3 cp ${root}_footprint.png s3://${BUCKET}/Pipeline/ --acl public-read
        
        cd ${root}/Extractions/
        sync_extractions_TO_s3 ${root}
        cd ../../
        
    else
        echo $cat
    fi
done

for root in $roots; do
    cd ${root}/Extractions/
    aws s3 sync --exclude "*" --include "*npy" \
                          --follow-symlinks ./ s3://${BUCKET}/Pipeline/${root}/Extractions/ --acl public-read
    
    cd ../../
done
    

    
###################
#
# Below are various checks meant to be run manually
#

# Make empty sync files for completed
roots=`aws s3 ls s3://${BUCKET}/Pipeline/Log/Sync/ | awk '{print $4}' | sed "s/.log//"`
for root in $roots; do
    cat=`aws s3 ls s3://${BUCKET}/Pipeline/${root}/Extractions/${root}.info.fits`
    if [[ ! -z $cat ]]; then
        sync=`aws s3 ls s3://${BUCKET}/Pipeline/Log/Start/${root}.log`
        if [[ -z $sync ]]; then
            echo $cat
            date > /tmp/${root}.log
            aws s3 cp /tmp/${root}.log s3://${BUCKET}/Pipeline/Log/Start/
            aws s3 cp /tmp/${root}.log s3://${BUCKET}/Pipeline/Log/Finished/
        fi
    fi
done

roots=`aws s3 ls s3://${BUCKET}/Pipeline/ | grep pdf |  awk '{print $4}' | sed "s/.pdf/.fits/"`
for root in $roots; do
    aws s3 cp ${root}  s3://${BUCKET}/Pipeline/
done

# Missing
roots=`ls *footp*pdf | sed "s/_footprint.pdf//"`
for root in $roots; do
    cat=`aws s3 ls s3://${BUCKET}/Pipeline/Log/Sync/${root}.log`
    if [[ -z $cat ]]; then
        echo "! Run $root"
    fi
done
