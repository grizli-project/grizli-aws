###
### Run all roots found in the synced directory
###

#BUCKET=aws-grivam
#BUCKET=grizli-grism
BUCKET=grizli-v1

roots=`aws s3 ls s3://${BUCKET}/IRAC/ | grep ipac.fits | awk '{print $4}' | sed "s/_ipac.fits//" | shuf`

# Existing
aws s3 ls s3://${BUCKET}/IRAC/Log/Start/ > /tmp/grizli_started.log
aws s3 ls s3://${BUCKET}/IRAC/Log/Finished/ > /tmp/grizli_finished.log

date > /tmp/irac_run_all.log

for root in $roots; do
    
    ostart=`grep ${root}.log /tmp/grizli_started.log | awk '{print $4}'`
    ostop=`grep ${root}.log /tmp/grizli_finished.log | awk '{print $4}'`
    #if [[ -z $ostart && -z $ostop ]]; then echo $root; fi; done
    
    if [[ -n $ostart || -n $ostop ]]; then
        echo "Skip ${root} (start=$ostart --- stop=$ostop)"
        echo "Skip ${root} (start=$ostart --- stop=$ostop)" >> /tmp/irac_run_all.log
        continue
    fi
    
    start=`aws s3 ls s3://${BUCKET}/IRAC/Log/Start/${root}.log | awk '{print $4}'`
    stop=`aws s3 ls s3://${BUCKET}/IRAC/Log/Finished/${root}.log| awk '{print $4}'`
    
    if [[ -z $start && -z $stop ]]; then
        echo "Run ${root}"
        echo "Run ${root}" >> /tmp/irac_run_all.log
        irac_run_single.sh ${root} 
    else
        echo "Skip ${root} (start=$start stop=$stop)"
        echo "Skip ${root} (start=$start stop=$stop)" >> /tmp/irac_run_all.log
    fi
done


#source deactivate ; source activate grizli-dev; cd ~/python/grizli-aws/; git pull; python setup.py install; cd ~/python/grizli; git checkout grizli/version.py; git pull; python setup.py install; cd ~
#