aws ec2  describe-instances --filters "Name=instance-state-name,Values=running" --query "Reservations[*].Instances[*].[InstanceId,InstanceType]" | sort

ids=`aws ec2  describe-instances --filters "Name=instance-state-name,Values=running" --query "Reservations[*].Instances[*].[InstanceId,InstanceType,PublicDnsName]" | grep large| awk '{print $1}' | sort`

ids=`aws ec2  describe-instances --filters "Name=instance-state-name,Values=running" --query "Reservations[*].Instances[*].[InstanceId]" | awk '{print $1}' | sort `

# aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters '{"commands":["auto_run_preprocess_single cos-j1001p0217-f140w-022"],"executionTimeout":["172000"]}' --timeout-seconds 600 --region us-east-1

for id in $ids; do     
    echo $id
    
    aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters '{"commands":["refresh_repositories > /tmp/refresh.log; aws s3 sync s3://grizli-preprocess/Scripts/ /usr/local/bin/ >> /tmp/refresh.log; chmod +x /usr/local/bin/auto_run*"],"executionTimeout":["1800"]}' --timeout-seconds 600 --region us-east-1

done

for id in $ids; do     
    echo $id
    
    #aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters '{"commands":["auto_run_spectral_fits","halt"],"executionTimeout":["172000"]}' --timeout-seconds 600 --region us-east-1
    #aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters '{"commands":["auto_run_preprocess","halt"],"executionTimeout":["172000"]}' --timeout-seconds 600 --region us-east-1
    aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters '{"commands":["auto_run_preprocess_fixwcs","halt"],"executionTimeout":["172000"]}' --timeout-seconds 600 --region us-east-1
    #aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters '{"commands":["auto_run_grism","halt"],"executionTimeout":["172000"]}' --timeout-seconds 600 --region us-east-1
    #aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters '{"commands":["auto_run_preprocess_single j1000p0206-f606w-005"],"executionTimeout":["14400"]}' --timeout-seconds 600 --region us-east-1

    #aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters '{"commands":["auto_fit_redshift_local j0332m2743","halt"],"executionTimeout":["21600"]}' --timeout-seconds 600 --region us-east-1


    sleep 5 # Force asynchronous

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