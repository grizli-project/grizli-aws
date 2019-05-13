# Generate catalogs for fields that don't have them

roots=`python -c "from grizli_aws import catalogs; roots, dates = catalogs.get_roots(verbose=False); print('\n'.join(roots[::-1]))"`

roots="j0226m0355 j0331m2843"

BUCKET=aws-grivam
BUCKET=grizli-grism

for root in $roots; do
    
    cat=`aws s3 ls s3://${BUCKET}/Pipeline/${root}/Extractions/${root}.info.fits`

    cat="" # Force
    
    cd /GrizliSpectra/
    rm -rf *
    
    if [[ -z $cat ]]; then
        echo $root  
    
        aws s3 sync --exclude "*" --include "${root}*full.fits" --include "*cat.fits" --acl public-read s3://${BUCKET}/Pipeline/${root}/Extractions/ ./
         
        grizli_extract_and_fit.py ${root} summary
        
        aws s3 sync --exclude "*" --include "${root}*fits" ./ s3://${BUCKET}/Pipeline/${root}/Extractions/ --acl public-read
        aws s3 sync --exclude "*" --include "${root}*png" --include "${root}*reg" --include "*html" --acl public-read ./ s3://${BUCKET}/Pipeline/${root}/Extractions/
        
        rm ${root}*
        
    else 
        echo "${root} file found"
    fi    
    
done

    