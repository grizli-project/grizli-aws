###
### Run all roots found in the synced directory
###

roots=`aws s3 ls s3://grizli-imaging/Pipeline/Fields/ | grep footprint.fits | awk '{print $4}' | sed "s/_footprint.fits//"`


for root in $roots; do
    start=`aws s3 ls s3://grizli-imaging/Pipeline/Log/Start/${root}.log | awk '{print $4}'`
    stop=`aws s3 ls s3://grizli-imaging/Pipeline/Log/Finished/${root}.log| awk '{print $4}'`
    
    if [[ -z $start && -z $stop ]]; then
        echo "Run ${root}"
        imaging_run_single.sh ${root} 
    else
        echo "Skip ${root} (start=$start stop=$stop)"
    fi
done


#source deactivate ; source activate grizli-dev; cd ~/python/grizli-aws/; git pull; python setup.py install; cd ~/python/grizli; git checkout grizli/version.py; git pull; python setup.py install; cd ~
# 