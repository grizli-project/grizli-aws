#!/usr/bin/env python
def fit_lambda(root='j100025+021706', newfunc=False, bucket_name='aws-grivam', skip_existing=True, check_wcs=False, use_psf=False, verbose=False, skip_started=True, quasar_fit=False, **kwargs):
    """
    check_wcs=True for ACS
    """
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
        func = 'GrizliLambda-0-12-0-41'
        
    print ('Lambda function: {0} (newfunc={1})'.format(func, newfunc))
               
    # Auth to create a Lambda function 
    session = boto3.Session()
    client = session.client('lambda', region_name='us-east-1')
    
    # s3 object = s3://aws-grivam/{obj}
    # obj = 'Pipeline/sparcs0034/Extractions/sparcs0034_00441.beams.fits'
     
    for obj in beams:
        print(obj)
        if False:
            # Defaults
            check_wcs=False # Set to true for ACS
            use_psf=False   # Use point source profile
            verbose=False   # Logging
            skip_started=True # SKip objects already started
            quasar_fit=False  # Fit with quasar templates
            
        event = {
              "s3_object_path": obj,
              "bucket": bucket_name,
              "verbose":      str(verbose),
              "skip_started": str(skip_started),
              "check_wcs" :   str(check_wcs),
              "quasar_fit" : str(quasar_fit),
              "use_psf"   : str(use_psf)
            }
        
        # output_path
        
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
    beams, files = get_needed_paths(root, bucket_name=bucket_name, skip_existing=True)
    
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

    #bucket_name = 'grizli-grism'
    bucket_name = 'aws-grivam'
    skip_existing = True
    newfunc = False
    
    kwargs = {'root':root, 
              'bucket_name':bucket_name, 
              'skip_existing':True, 
              'newfunc':False}
    
    # Args passed to the lambda event
    kwargs['check_wcs'] = False # Set to true for ACS
    kwargs['use_psf'] = False   # Use point source profile
    kwargs['verbose'] = False   # Logging
    kwargs['skip_started'] = True # SKip objects already started
    kwargs['quasar_fit'] = False  # Fit with quasar templates
    
    if len(sys.argv) > 2:
        for args in sys.argv[2:]:
            keypair = args.strip('--').split('=')
            if keypair[0] in ['newfunc','skip_existing']:
                if len(keypair) == 1:
                    kwargs[keypair[0]] = True
                else:
                    kwargs[keypair[0]] = keypair[1].lower() in ['true']
            else:
                if keypair[0] in kwargs:
                    kwargs[keypair[0]] = keypair[1]
            
    fit_lambda(**kwargs) #newfunc=newfunc, bucket_name=bucket_name, skip_existing=skip_existing)
    