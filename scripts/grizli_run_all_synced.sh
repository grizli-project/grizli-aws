###
### Run all roots found in the synced directory
###

roots=`aws s3 ls s3://aws-grivam/Pipeline/Log/Sync/ | awk '{print $4}' | sed "s/.log//"`

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
        echo "Skip ${root} (start=$start  stop=$stop)"
    fi
done
