###
### Run all roots found in the synced directory
###

BUCKET=aws-grivam
BUCKET=grizli-grism

roots=`aws s3 ls s3://${BUCKET}/Pipeline/Log/Extract/ | awk '{print $4}' | sed "s/.log//"`

roots=`aws s3 ls s3://${BUCKET}/Pipeline/Log/Extract/ | awk '{print $4}' | sed "s/.log//" | grep -v j033 | grep -v beams` # Skip UDF

aws s3 ls s3://${BUCKET}/Pipeline/Log/ExtractStart/ > /tmp/grizli_started.log
aws s3 ls s3://${BUCKET}/Pipeline/Log/Extract/ > /tmp/grizli_extracted.log
aws s3 ls s3://${BUCKET}/Pipeline/Log/Finished/ > /tmp/grizli_finished.log

if [[ -z $1 ]]; then
    maglim="16.5,26"
else
    maglim=$2
fi

date > /tmp/grizli_run_all_extracted.log

for root in $roots; do
    # Check from existing
    ostart=`grep ${root}.log /tmp/grizli_started.log | awk '{print $4}'`
    ostop=`grep ${root}.log /tmp/grizli_finished.log | awk '{print $4}'`
    #if [[ -z $ostart && -z $ostop ]]; then echo $root; fi; done #else echo "Skip ${root}"; fi; done
    
    if [[ -n $ostart || -n $ostop ]]; then
        echo "Skip ${root} (start=$ostart --- stop=$ostop)"
        echo "Skip ${root} (start=$ostart --- stop=$ostop)" >> /tmp/grizli_run_all_extracted.log
        continue
    fi
    
    # Check again in case done by another process
    start=`aws s3 ls s3://${BUCKET}/Pipeline/Log/ExtractStart/${root}.log | awk '{print $4}'`
    extract=`aws s3 ls s3://${BUCKET}/Pipeline/Log/Extract/${root}.log | awk '{print $4}'`
    stop=`aws s3 ls s3://${BUCKET}/Pipeline/Log/Finished/${root}.log| awk '{print $4}'`
    
    if [[ -z $start && -z $stop ]]; then
        echo "Run ${root}"
        date > ${root}.log
        aws s3 cp ${root}.log s3://${BUCKET}/Pipeline/Log/ExtractStart/
        grizli_extract_and_fit.sh ${root} ${maglim}
    else
        aws s3 ls s3://${BUCKET}/Pipeline/${root}/Extractions/ > /tmp/ext
        beams=`grep beams.fits /tmp/ext | wc -l`    
        full=`grep full.fits /tmp/ext | wc -l`    
    
        echo "Skip ${root} (start=$start extract=$extract  stop=$stop) ${beams} / ${full}"
    fi
done


#source deactivate ; source activate grizli-dev; cd ~/python/grizli-aws/; git pull; python setup.py install; cd ~/python/grizli; git checkout grizli/version.py; git pull; python setup.py install; cd ~
# 