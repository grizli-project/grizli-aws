#!/usr/bin/env python
def fit_lambda(root='j100025+021706', newfunc=True, bucket_name='aws-grivam', skip_existing=True):
    import time
    import os
    import numpy as np
    import boto3
    import json

    beams, files = get_needed_paths(root, bucket_name=bucket_name, skip_existing=skip_existing)
    if len(beams) == 0:
        print('{0}: No beams to fit'.format(root))
        
        return False
    # Auth to create a Lambda function (credentials are picked up from above .aws/credentials)
    session = boto3.Session()

    # Make sure Lambda is running in the same region as the HST public dataset
    client = session.client('lambda', region_name='us-east-1')

    if newfunc:
        func = 'GrizliZ'+root.replace('+','p').replace('-','m')
        
        try:
            # Create the function
            response = client.create_function(
                FunctionName=func,
                Runtime='python3.6',
                Role='arn:aws:iam::521547342229:role/grivam_lambda',
                Handler='process.handler',
                Code={
                    'S3Bucket': 'aws-grivam', 
                    'S3Key': 'venv.zip'
                },
                Description='Lambda redshift fits: '+root,
                Timeout=300,
                MemorySize=1024,
                Publish=True
            )
        except:
            pass
    else:
        func = 'GrizliTestFunction'
     
    print ('Lambda function: {0} (newfunc={1})'.format(func, newfunc))
               
    # Auth to create a Lambda function 
    session = boto3.Session()
    client = session.client('lambda', region_name='us-east-1')

    for obj in beams:
        print(obj)
        event = {
              's3_object_path': obj,
              'verbose':      "False",
              'skip_started': "True",
              'check_wcs' :   "False",
              'bucket': bucket_name,
            }

        # Invoke Lambda function
        response = client.invoke(
            FunctionName=func,
            InvocationType='Event',
            LogType='Tail',
            Payload=json.dumps(event))
    
    sleep_time = 303*np.ceil(len(beams)/950)
    print('{0}: sleep {1}'.format(time.ctime(), sleep_time))
    
    time.sleep(sleep_time)
    
    # Status again to check products
    beams, files = get_needed_paths(root, bucket_name=bucket_name, skip_existing=False)
    
def get_needed_paths(root, get_string=False, bucket_name='aws-grivam', skip_existing=True):
    """
    Get the S3 paths of the "beams.fits" files that still need to be fit.
    """
    import boto3
    import time
    
    s3 = boto3.resource('s3')
    s3_client = boto3.client('s3')
    bkt = s3.Bucket(bucket_name)
    
    files = [obj.key for obj in bkt.objects.filter(Prefix='Pipeline/{0}/Extractions/'.format(root))]

    beams = []
    logs = []
    full = []
    start = []
    
    for file in files:
        if 'beams.fits' in file:
            beams.append(file)

        if 'log_par' in file:
            logs.append(file)
        
        if 'start.log' in file:
            start.append(file)
        
        if 'full.fits' in file:
            full.append(file)
    
    label = '{0} / {1} / Nbeams: {2}, Nfull: {3}, Nlog: {4}, Nstart: {5}'.format(root, time.ctime(), len(beams), len(full), len(logs), len(start))
    if get_string:
        return label
        
    print(label)
    
    for i in range(len(beams))[::-1]:
        test = (beams[i].replace('.beams.fits', '.full.fits') in full)
        test |= (beams[i].replace('.beams.fits', '.start.log') in start)
        if test & skip_existing:
            beams.pop(i)

        
    return beams, files
    
if __name__ == "__main__":
    import sys
    import numpy as np
    from grizli import utils
    from grizli.pipeline import auto_script
    utils.set_warnings()
    
    if len(sys.argv) < 2:
        print('Usage: fit_redshift_lambda.py {field}')
        exit 
    
    root = sys.argv[1]
    if len(sys.argv) == 3:
        newfunc = bool(sys.argv[2].lower() == 'true')
    else:
        newfunc = False
        
    fit_lambda(root=root, newfunc=newfunc, bucket_name='grizli-grism', skip_existing=False)
    