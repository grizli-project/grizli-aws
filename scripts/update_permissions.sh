### Global update permissions, takes forever
roots=`aws s3 ls s3://aws-grivam/Pipeline/Log/Finished/ | awk '{print $4}' | sed "s/.log//"`

for root in $roots; do 

    echo ""; echo ""; echo "======="
    echo $root
    echo "======="; echo ""; echo ""; 

    aws s3 ls s3://aws-grivam/Pipeline/${root}/Extractions/ > ${root}.files
    files=`grep ${root} ${root}.files | grep fits | awk '{print $4}'`
    
    for file in $files; do 
        echo $file
        aws s3api put-object-acl --bucket aws-grivam --key Pipeline/${root}/Extractions/${file} --grant-read uri=http://acs.amazonaws.com/groups/global/AllUsers
    
    done
done