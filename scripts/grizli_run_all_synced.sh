###
### Run all roots found in the synced directory
###

roots=`aws s3 ls s3://aws-grivam/Pipeline/Log/Sync/ | awk '{print $4}' | sed "s/.log//"`

roots=`aws s3 ls s3://aws-grivam/Pipeline/Log/Sync/ | awk '{print $4}' | sed "s/.log//" | grep -v j033` # Skip UDF

if [[ -z $1 ]]; then
    maglim="16.5,26"
else
    maglim=$2
fi

for root in $roots; do
    start=`aws s3 ls s3://aws-grivam/Pipeline/Log/Start/${root}.log | awk '{print $4}'`
    stop=`aws s3 ls s3://aws-grivam/Pipeline/Log/Finished/${root}.log| awk '{print $4}'`
    
    if [[ -z $start && -z $stop ]]; then
        echo "Run ${root}"
        grizli_extract_and_fit.sh ${root} ${maglim}
    else
        aws s3 ls s3://aws-grivam/Pipeline/${root}/Extractions/ > /tmp/ext
        beams=`grep beams.fits /tmp/ext | wc -l`    
        full=`grep full.fits /tmp/ext | wc -l`    
    
        echo "Skip ${root} (start=$start  stop=$stop) ${beams} / ${full}"
    fi
done


#source deactivate ; source activate grizli-dev; cd ~/python/grizli-aws/; git pull; python setup.py install; cd ~/python/grizli; git checkout grizli/version.py; git pull; python setup.py install; cd ~
# 