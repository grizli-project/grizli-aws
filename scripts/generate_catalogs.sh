
roots=`python -c "from grizli_aws import catalogs; roots, dates = catalogs.get_roots(verbose=False); print('\n'.join(roots[::-1]))"`

for root in $roots; do
    
    cat=`aws s3 ls s3://aws-grivam/Pipeline/${root}/Extractions/${root}.info.fits`

    cd ~/GrizliExtract/
    
    if [[ -z $cat ]]; then
        echo $root  
    
        aws s3 sync --exclude "*" --include "${root}*full.fits" --include "*cat.fits" --acl public-read s3://aws-grivam/Pipeline/${root}/Extractions/ ./
         
        grizli_extract_and_fit.py ${root} summary
        
        aws s3 sync --exclude "*" --include "${root}*fits" ./ s3://aws-grivam/Pipeline/${root}/Extractions/
        aws s3 sync --exclude "*" --include "${root}*png" --include "${root}*reg" --include "*html" --acl public-read ./ s3://aws-grivam/Pipeline/${root}/Extractions/
        
        rm ${root}*
        
    else 
        echo "${root} file found"
    fi    
    
done

    