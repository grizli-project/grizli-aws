#!/usr/bin/env python
def fit_lambda(root='j100025+021706', newfunc=True):
    import time
    import os
    import numpy as np
    import boto3
    import json

    from grizli_aws.fit_redshift_single import run_grizli_fit
    
    
    beams, files = get_needed_paths(root)
    if len(beams) == 0:
        print('{0}: No beams to fit'.format(root))
        
    for obj in beams:
        print(obj)
        event = {
              's3_object_path': obj,
              'verbose':      "True",
              'skip_started': "True",
              'check_wcs' :   "False",
            }
        
        try:
            run_grizli_fit(event)
        except:
            pass
            
    # Status again to check products
    beams, files = get_needed_paths(root)
    
def get_needed_paths(root, get_string=False):
    """
    Get the S3 paths of the "beams.fits" files that still need to be fit.
    """
    import boto3
    import time
    
    s3 = boto3.resource('s3')
    s3_client = boto3.client('s3')
    bkt = s3.Bucket('aws-grivam')
    
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
        if test:
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
        
    fit_lambda(root=root, newfunc=newfunc)
    