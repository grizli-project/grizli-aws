###
### Run all roots found in the synced directory
###

BUCKET=grizli-imaging

for arg in "$@" ; do
    if [[ ! -z `echo "${arg}" | grep "bucket="` ]]; then 
        BUCKET=`echo "${arg}" | sed "s/=/ /" | awk '{print $2}'`
    fi
done

roots=`aws s3 ls s3://${BUCKET}/Pipeline/Fields/ | grep footprint.fits | awk '{print $4}' | sed "s/_footprint.fits//" | shuf`

# Existing
aws s3 ls s3://${BUCKET}/Pipeline/Log/Start/ > /tmp/grizli_started.log
aws s3 ls s3://${BUCKET}/Pipeline/Log/Finished/ > /tmp/grizli_finished.log

date > /tmp/imaging_run_all.log

for root in $roots; do
    
    ostart=`grep ${root}.log /tmp/grizli_started.log | awk '{print $4}'`
    ostop=`grep ${root}.log /tmp/grizli_finished.log | awk '{print $4}'`
    #if [[ -z $ostart && -z $ostop ]]; then echo $root; fi; done
    
    if [[ -n $ostart || -n $ostop ]]; then
        echo "Skip ${root} (start=$ostart --- stop=$ostop)"
        echo "Skip ${root} (start=$ostart --- stop=$ostop)" >> /tmp/imaging_run_all.log
        continue
    fi
    
    start=`aws s3 ls s3://${BUCKET}/Pipeline/Log/Start/${root}.log | awk '{print $4}'`
    stop=`aws s3 ls s3://${BUCKET}/Pipeline/Log/Finished/${root}.log| awk '{print $4}'`
    
    if [[ -z $start && -z $stop ]]; then
        echo "Run ${root}"
        echo "Run ${root}" >> /tmp/imaging_run_all.log
        grism_run_single.sh ${root} $@
        rm -rf /GrizliImaging/${root}*
    else
        echo "Skip ${root} (start=$start stop=$stop)"
        echo "Skip ${root} (start=$start stop=$stop)" >> /tmp/imaging_run_all.log
    fi
done


#source deactivate ; source activate grizli-dev; cd ~/python/grizli-aws/; git pull; python setup.py install; cd ~/python/grizli; git checkout grizli/version.py; git pull; python setup.py install; cd ~
# 