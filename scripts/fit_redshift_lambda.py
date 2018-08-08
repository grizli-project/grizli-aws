#!/usr/bin/env python
def fit_lambda(root='j100025+021706', newfunc=True):
    import time
    import os
    import numpy as np
    import boto3
    import json

    beams, files = get_needed_paths(root)
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
              'verbose': "False",
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
    beams, files = get_needed_paths(root)
    
def get_needed_paths(root):
    """
    Get the S3 paths of the "beams.fits" files that still need to be fit.
    """
    import boto3
    
    s3 = boto3.resource('s3')
    s3_client = boto3.client('s3')
    bkt = s3.Bucket('aws-grivam')
    
    files = [obj.key for obj in bkt.objects.filter(Prefix='Pipeline/{0}/Extractions/'.format(root))]

    beams = []
    logs = []
    full = []

    for file in files:
        if 'beams.fits' in file:
            beams.append(file)

        if 'log_par' in file:
            logs.append(file)

        if 'full.fits' in file:
            full.append(file)
    
    print('{3} / Nbeams: {0}, Nfull: {1}, Nlog: {2}'.format(len(beams), len(full), len(logs), root))
    
    for i in range(len(beams))[::-1]:
        if  (beams[i].replace('.beams.fits', '.full.fits') in full):
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
        newfunc = True
        
    fit_lambda(root=root, newfunc=newfunc)
    