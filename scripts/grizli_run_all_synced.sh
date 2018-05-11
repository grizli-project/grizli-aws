"""
Run all roots found in the synched directory
"""

roots=`aws s3 ls s3://aws-grivam/Pipeline/Log/Sync/ | awk '{print $4}' | sed "s/.log//"`

for root in $roots; do
    start=`aws s3 ls s3://aws-grivam/Pipeline/Log/Start/${root}.log | awk '{print $4}'`
    stop=`aws s3 ls s3://aws-grivam/Pipeline/Log/Finished/${root}.log| awk '{print $4}'`
    
    if [[ -z $start && -z $stop ]]; then
        echo "Run ${root}"
        grizli_extract_and_fit.sh ${root}
    else
        echo "Skip ${root} (start=$start  stop=$stop)"
    fi
done
