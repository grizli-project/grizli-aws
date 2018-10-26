###
### Run all roots found in the synced directory
###

BUCKET=aws-grivam
BUCKET=grizli-grism

roots=`aws s3 ls s3://${BUCKET}/Pipeline/Log/Sync/ | awk '{print $4}' | sed "s/.log//"`

roots=`aws s3 ls s3://${BUCKET}/Pipeline/Log/Sync/ | awk '{print $4}' | sed "s/.log//" | grep -v j033` # Skip UDF


aws s3 ls s3://${BUCKET}/Pipeline/Log/Start/ > /tmp/grizli_started.log
aws s3 ls s3://${BUCKET}/Pipeline/Log/Extract/ > /tmp/grizli_extracted.log
aws s3 ls s3://${BUCKET}/Pipeline/Log/Finished/ > /tmp/grizli_finished.log

if [[ -z $1 ]]; then
    maglim="16.5,26"
else
    maglim=$2
fi

date > /tmp/grizli_run_all_synced.log

for root in $roots; do
    
    # Check from existing
    ostart=`grep ${root}.log /tmp/grizli_started.log | awk '{print $4}'`
    ostop=`grep ${root}.log /tmp/grizli_finished.log | awk '{print $4}'`
    #if [[ -z $ostart && -z $ostop ]]; then echo $root; fi; done
    
    if [[ -n $ostart || -n $ostop ]]; then
        echo "Skip ${root} (start=$ostart --- stop=$ostop)"
        echo "Skip ${root} (start=$ostart --- stop=$ostop)" >> /tmp/grizli_run_all_synced.log
        continue
    fi
    
    # Check again in case done by another process
    start=`aws s3 ls s3://${BUCKET}/Pipeline/Log/Start/${root}.log | awk '{print $4}'`
    extract=`aws s3 ls s3://${BUCKET}/Pipeline/Log/Extract/${root}.log | awk '{print $4}'`
    stop=`aws s3 ls s3://${BUCKET}/Pipeline/Log/Finished/${root}.log| awk '{print $4}'`
    
    if [[ -z $start && -z $stop ]]; then
        echo "Run ${root}" 
        echo "Run ${root}" >> /tmp/grizli_run_all_synced.log
        grizli_extract_only.sh ${root} ${maglim}
    else
        aws s3 ls s3://${BUCKET}/Pipeline/${root}/Extractions/ > /tmp/ext
        beams=`grep beams.fits /tmp/ext | wc -l`    
        full=`grep full.fits /tmp/ext | wc -l`    
    
        echo "Skip ${root} (start=$start extract=$extract  stop=$stop) ${beams} / ${full}"
        echo "Skip ${root} (start=$start extract=$extract  stop=$stop) ${beams} / ${full}" >> /tmp/grizli_run_all_synced.log
    fi
done


#source deactivate ; source activate grizli-dev; cd ~/python/grizli-aws/; git pull; python setup.py install; cd ~/python/grizli; git checkout grizli/version.py; git pull; python setup.py install; cd ~
# 