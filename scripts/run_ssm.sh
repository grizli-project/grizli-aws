aws ec2  describe-instances --filters "Name=instance-state-name,Values=running" --query "Reservations[*].Instances[*].[InstanceId,InstanceType]" | sort

ids=`aws ec2  describe-instances --filters "Name=instance-state-name,Values=running" --query "Reservations[*].Instances[*].[InstanceId,InstanceType,PublicDnsName]" | sort -k 1 | awk '{print $1}' `

### c5 instances
ec2_type="r5d.large" 
ids=`aws ec2  describe-instances --filters "Name=instance-state-name,Values=running" --query "Reservations[*].Instances[*].[InstanceId,InstanceType,PublicDnsName]" | sort -k 1 |  grep ${ec2_type} | awk '{print $1}' `
echo `echo $ids | wc -w` "${ec2_type} instances"

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
init="\"refresh_repositories > /tmp/refresh.log\",\"chmod +x /usr/local/bin/auto_run*\",\"mount_nvme > /tmp/mount.log\","
init="\"refresh_repositories\","
#halt=""
halt=",\"halt_if_grizli_done\""

# Don't init
init=""
halt=",\"halt\""

# Halt if all done
#halt=",\"hasfiles=`ls /GrizliImaging/ | grep log`; if [ -z \$hasfiles ]; then halt; fi\""

# for id in $ids; do
#     # Halt if done
#     aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters "{\"commands\":[\"refresh_repositories; halt_if_grizli_done\"],\"executionTimeout\":[\"172000\"]}" --timeout-seconds 600 --region us-east-1
# done

for id in $ids; do     
    echo $id

    ####### Halt
    #aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters '{"commands":["halt"],"executionTimeout":["172000"]}' --timeout-seconds 600 --region us-east-1

    ####### Init only
    # aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters "{\"commands\":[${init}],\"executionTimeout\":[\"172000\"]}" --timeout-seconds 600 --region us-east-1


    #aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters '{"commands":["auto_run_spectral_fits","halt"],"executionTimeout":["172000"]}' --timeout-seconds 600 --region us-east-1
    
    ########### Run pipeline on each sub-field
    #aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters '{"commands":["auto_run_preprocess","halt"],"executionTimeout":["172000"]}' --timeout-seconds 600 --region us-east-1
    
    ################################################################
    
    ########### new preprocessing for CANDELS fields
    # aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters "{\"commands\":[${init}\"auto_run_imaging --only_preprocess=True --make_mosaics=False --make_thumbnails=False --make_phot=False --is_parallel_field=False --extra_filters=g800l  --bucket=grizli-v1\"${halt}],\"executionTimeout\":[\"172000\"]}" --timeout-seconds 600 --region us-east-1

    ########### new preprocessing for COSMOS fields
    aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters "{\"commands\":[${init}\"auto_run_imaging --only_preprocess=True --make_mosaics=False --make_thumbnails=False --make_phot=False --is_parallel_field=False --extra_filters=g800l  --bucket=grizli-cosmos-v2\"${halt}],\"executionTimeout\":[\"172000\"]}" --timeout-seconds 600 --region us-east-1
    
    ########### new preprocessing for CANDELS blue filters
    #aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters "{\"commands\":[${init}\"auto_run_imaging --only_preprocess=True --make_mosaics=False --make_thumbnails=False --make_phot=False --is_parallel_field=False --extra_filters=f435w,f275w,f336w --mosaic_args.optical_filters=f435w,f275w,f336w --visit_prep_args.max_err_percentile=95 --bucket=grizli-v1\"${halt}],\"executionTimeout\":[\"172000\"]}" --timeout-seconds 600 --region us-east-1
    
    ### Redo mosaics with fixed persistence and mask_spikes
    aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters "{\"commands\":[${init}\"auto_run_imaging --run_fine_alignment=2 --bucket=grizli-v1 --run_parse_visits=False --preprocess_args.skip_single_optical_visits=True --extra_filters=g800l --redo_persistence_mask=True --mask_spikes=True --persistence_args.err_threshold=0.1 --persistence_args.reset=False --sync\"${halt}],\"executionTimeout\":[\"172000\"]}" --timeout-seconds 600 --region us-east-1
    
    # Reset Persistence
    aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters "{\"commands\":[${init}\"auto_run_imaging --run_fine_alignment=2 --bucket=grizli-v1 --run_parse_visits=False --preprocess_args.skip_single_optical_visits=True --extra_filters=g800l --redo_persistence_mask=True --mask_spikes=True --persistence_args.err_threshold=0.6 --persistence_args.reset=True --sync\"${halt}],\"executionTimeout\":[\"172000\"]}" --timeout-seconds 600 --region us-east-1
    
    # Redo from start
    ##  --run_fine_alignment=True --extra_filters=g800l  --bucket=grizli-v1 --is_parallel_field=False --preprocess_args.skip_single_optical_visits=True --mask_spikes=True --persistence_args.err_threshold=0.2
    
    ## --visit_prep_args.reference_catalogs=PS1,DES,DSC,SDSS,GAIA,WISE
    ## --visit_prep_args.align_mag_limits=20,23
    ## --is_parallel_field=True --parse_visits_args.combine_same_pa=1 --parse_visits_args.max_dt=0.3 --preprocess_args.skip_single_optical_visits=True

    # Parallel fields
    ## aws s3 rm s3://grizli-v1/Pipeline/${root}/ --recursive ;
    ## grism_run_single.sh ${root} --is_parallel_field=True --parse_visits_args.combine_same_pa=1 --parse_visits_args.max_dt=0.4 --run_fine_alignment=True --extra_filters=g800l --bucket=grizli-v1 --preprocess_args.skip_single_optical_visits=True --mask_spikes=True --persistence_args.err_threshold=0.5
    
    ### Frontier fields
    # gunzip ${root}/Prep/${root}*gz ${root}/Extractions/${root}*gz
    # grism_run_single.sh ${root} --run_fine_alignment=True --bucket=grizli-v1 --filters=F814W --is_parallel_field=False --noclean
    
    # grism_run_single.sh ${root} --run_fine_alignment=True --bucket=grizli-v1 --extra_filters=f435w,f475w,g800l --is_parallel_field=False --noclean
    
    ########### new preprocessing for CHARGE fields
    #aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters "{\"commands\":[${init}\"auto_run_imaging --run_fine_alignment=True --extra_filters=g800l  --bucket=grizli-v1\"${halt}],\"executionTimeout\":[\"172000\"]}" --timeout-seconds 600 --region us-east-1
    aws ssm send-command --document-name "AWS-RunShellScript" --instance-ids "${id}" --parameters "{\"commands\":[${init}\"auto_run_imaging --run_fine_alignment=True --extra_filters=g800l --bucket=grizli-v1 --preprocess_args.skip_single_optical_visits=True --mask_spikes=True --persistence_args.err_threshold=0.5 --visit_prep_args.single_image_CRs=False \"${halt}],\"executionTimeout\":[\"172000\"]}" --timeout-seconds 600 --region us-east-1
    
    #################################################################
    
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