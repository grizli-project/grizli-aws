aws ec2  describe-instances --filters "Name=instance-state-name,Values=running" --query "Reservations[*].Instances[*].[InstanceId,InstanceType]" | sort

ids=`aws ec2  describe-instances --filters "Name=instance-state-name,Values=running" --query "Reservations[*].Instances[*].[InstanceId,InstanceType,PublicDnsName]" | sort -k 1 | awk '{print $1}' `

ids=`aws ec2  describe-instances --filters "Name=instance-state-name,Values=running" --query "Reservations[*].Instances[*].[InstanceId,InstanceType,PublicDnsName]" | sort -k 1 |  grep c5 | awk '{print $1}' `


### Drizzler
id=i-0a8bae4cfe82dbce5
aws ec2 start-instances --instance-id $id --profile oliveraws
aws ec2 describe-instances --instance-id $id --profile oliveraws --query "Reservations[*].Instances[*].[InstanceId,InstanceType,PublicDnsName]"


ids=`aws ec2  describe-instances --filters "Name=instance-state-name,Values=running" --query "Reservations[*].Instances[*].[InstanceType,InstanceId,PublicDnsName]" | sort -k 2 | awk '{print $2}'`

# aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters '{"commands":["auto_run_preprocess_single cos-j1001p0217-f140w-022"],"executionTimeout":["172000"]}' --timeout-seconds 600 --region us-east-1

###### Check status
for ext in Start Finished Failed; do aws s3 ls s3://grizli-v1/Pipeline/Log/${ext}/ > /tmp/${ext}.log; echo ""; echo $ext `wc -l /tmp/${ext}.log`; echo ""; cat /tmp/${ext}.log | sort -k 1 -k 2 | tail; done

## CLEAN LOGS
bucket=grizli-v1

aws s3 rm --recursive s3://${bucket}/Pipeline/Log/Start/
aws s3 sync s3://${bucket}/Pipeline/Log/Finished/ s3://${bucket}/Pipeline/Log/Start/ --acl public-read
aws s3 sync s3://${bucket}/Pipeline/Log/Failed/ s3://${bucket}/Pipeline/Log/Start/ --acl public-read

# Refresh all repos
for id in $ids; do     
    echo $id
    
    aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters '{"commands":["refresh_repositories > /tmp/refresh.log; aws s3 sync s3://grizli-preprocess/Scripts/ /usr/local/bin/ >> /tmp/refresh.log; chmod +x /usr/local/bin/auto_run*"],"executionTimeout":["1800"]}' --timeout-seconds 600 --region us-east-1

done

# INit, don't kill
#init="\"refresh_repositories > /tmp/refresh.log; aws s3 sync s3://grizli-preprocess/Scripts/ /usr/local/bin/ >> /tmp/refresh.log; chmod +x /usr/local/bin/auto_run*\","
init="\"refresh_repositories > /tmp/refresh.log; chmod +x /usr/local/bin/auto_run*\","
halt=""

# Don't init
init=""
halt=",\"halt\""

# Halt if all done
halt=",\"hasfiles=`ls /GrizliImaging/ | grep log`; if [ -z \$hasfiles ]; then halt; fi\""

for id in $ids; do
    # Halt if done
    aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters "{\"commands\":[\"refresh_repositories; halt_if_grizli_done\"],\"executionTimeout\":[\"172000\"]}" --timeout-seconds 600 --region us-east-1
done

for id in $ids; do     
    echo $id
    
    #aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters '{"commands":["auto_run_spectral_fits","halt"],"executionTimeout":["172000"]}' --timeout-seconds 600 --region us-east-1
    
    ########### Run pipeline on each sub-field
    #aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters '{"commands":["auto_run_preprocess","halt"],"executionTimeout":["172000"]}' --timeout-seconds 600 --region us-east-1
    
    ########### new preprocessing for CANDELS fields
    aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters "{\"commands\":[${init}\"auto_run_imaging --only_preprocess=True --make_mosaics=False --make_thumbnails=False --make_phot=False --bucket=grizli-v1\"${halt}],\"executionTimeout\":[\"172000\"]}" --timeout-seconds 600 --region us-east-1
        
    ############# Run grism preprep
    #aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters "{\"commands\":[${init}\"auto_run_grism\"${halt}],\"executionTimeout\":[\"172000\"]}" --timeout-seconds 600 --region us-east-1
    
    # ############# Run fixwcs
    # aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters "{\"commands\":[${init}\"auto_run_preprocess_fixwcs\",\"halt\"],\"executionTimeout\":[\"172000\"]}" --timeout-seconds 600 --region us-east-1
    
    #aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters '{"commands":["auto_run_grism","halt"],"executionTimeout":["172000"]}' --timeout-seconds 600 --region us-east-1
    #aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters '{"commands":["auto_run_preprocess_single j1000p0206-f606w-005"],"executionTimeout":["14400"]}' --timeout-seconds 600 --region us-east-1

    #aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters '{"commands":["auto_fit_redshift_local j0332m2743","halt"],"executionTimeout":["21600"]}' --timeout-seconds 600 --region us-east-1


    sleep 8 # Force asynchronous

done


# By DNS
aws ec2  describe-instances --filters "Name=instance-state-name,Values=running" --query "Reservations[*].Instances[*].[InstanceId,PublicDnsName]" | sort -k 1

dns_names=`aws ec2  describe-instances --filters "Name=instance-state-name,Values=running" --query "Reservations[*].Instances[*].[InstanceId,PublicDnsName]" | sort -k 1 | awk '{print $2}'`
for dns in $dns_names; do
    rsync -avz --progress --dry-run -e "ssh -o StrictHostKeyChecking=no -i $HOME/OwnMacbook.pem" ec2-user@${dns}:/GrizliSpectra/ ./ > /tmp/dns.files

    log=`grep beams.log /tmp/dns.files`
    nfull=`grep full.fits /tmp/dns.files | wc -l`
    nbeams=`grep beams.fits /tmp/dns.files | wc -l`
    echo "$dns  $log    $nbeams  $nfull"
done