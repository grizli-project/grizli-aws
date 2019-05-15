#!/usr/bin/env python
def fit_lambda(root='j100025+021706', beams=[], newfunc=False, bucket_name='aws-grivam', skip_existing=True, sleep=True, check_wcs=False, use_psf=False, verbose=False, skip_started=True, quasar_fit=False, zr=None, output_path=None, show_event=False, **kwargs):
    """
    check_wcs=True for ACS
    """
    import time
    import os
    import yaml
    
    import numpy as np
    import boto3
    import json

    if len(beams) == 0:
        beams, files = get_needed_paths(root, bucket_name=bucket_name, skip_existing=skip_existing)
    else:
        sleep = False
        
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
    
    # s3 object = s3://{bucket_name}/{s3_object_path}
    # e.g., 'Pipeline/sparcs0034/Extractions/sparcs0034_00441.beams.fits'
     
    NB = len(beams)
    for i, s3_object_path in enumerate(beams):
        print('{0:>5} / {1:>5} : {2}'.format(i+1, NB, s3_object_path))
        
        if False:
            # Defaults
            check_wcs=False # Set to true for ACS
            use_psf=False   # Use point source profile
            verbose=False   # Logging
            skip_started=True # SKip objects already started
            quasar_fit=False  # Fit with quasar templates
            
        event = {
              "s3_object_path": s3_object_path,
              "bucket": bucket_name,
              "verbose":      str(verbose),
              "skip_started": str(skip_started),
              "check_wcs" :   str(check_wcs),
              "quasar_fit" : str(quasar_fit),
              "use_psf"   : str(use_psf)
            }
        
        if output_path is not None:
            if output_path == 'self':
                event['output_path'] = os.path.dirname(s3_object_path)
            else:
                event['output_path'] = output_path
        
        if zr is not None:
            event['zr'] = zr.split(',')
            
        if show_event:
            print('Event: \n\n', event.__str__().replace('\'', '"').replace(',', ',\n'))
            
        # Invoke Lambda function
        response = client.invoke(
            FunctionName=func,
            InvocationType='Event', 
            LogType='Tail',
            Payload=json.dumps(event))
    
    # Sleep for status check    
    if sleep:
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
    import yaml
    
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
              'newfunc':False,
              'sleep':True,
              'beams':[],
              'output_path':None,
              'zr':None,
              'show_event':False}
    
    # Args passed to the lambda event
    kwargs['check_wcs'] = False # Set to true for ACS
    kwargs['use_psf'] = False   # Use point source profile
    kwargs['verbose'] = False   # Logging
    kwargs['skip_started'] = True # SKip objects already started
    kwargs['quasar_fit'] = False  # Fit with quasar templates
    
    dryrun=False
    
    if len(sys.argv) > 2:
        for args in sys.argv[2:]:
            keypair = args.strip('--').split('=')
            
            # Booleans
            if keypair[0] in ['newfunc','skip_existing','sleep','checkwcs','use_psf','verbose','skip_started','quasar_fit', 'show_event']:
                if len(keypair) == 1:
                    kwargs[keypair[0]] = True
                else:
                    kwargs[keypair[0]] = keypair[1].lower() in ['true']
            
            # List of s3 beams.fits paths
            elif keypair[0] == 'beams':
                # Fit a single object
                kwargs['beams'] = keypair[1].split(',')
            
            # List of ids associated with {root}
            elif keypair[0] == 'ids':
                kwargs['beams'] = ['Pipeline/{0}/Extractions/{0}_{1:05d}.beams.fits'.format(root, int(id)) for id in keypair[1].split(',')]
            
            # don't run the script
            elif keypair[0] == 'dryrun':
                dryrun=True
            
            # Everything else
            else:
                if keypair[0] in kwargs:
                    kwargs[keypair[0]] = keypair[1]
    
    print('Arguments: \n\n', '  '+yaml.dump(kwargs).replace('\n', '\n   '))
    
    if not dryrun:
        fit_lambda(**kwargs) #newfunc=newfunc, bucket_name=bucket_name, skip_existing=skip_existing)
    