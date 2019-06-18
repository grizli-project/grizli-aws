"""
Move exposures to a general directory like


"""

def move_all():
    import grizli_aws.utils
    roots=glob.glob('*wfc3ir*')
    for i, root in enumerate(roots):
        print('\n\n\n# {0} {1}\n\n'.format(i, root))
        
        grizli_aws.utils.move_exposures(root=root, bucket='grizli-v1', base_path='Exposures', verbose=True)
        os.system('mv {0}_visits.npy {0}/Prep/'.format(root))
        
def move_exposures(root='j021732m0512_1716_bfu_marshall_wfc3ir_f125w-f160w-g141', bucket='grizli-v1', base_path='Exposures', verbose=True):
    """
    Move FLT/FLC files from visit directories to a common exposure bucket
     
    Paths like s3://stpubdata/hst/public/icwb/icwb1iu5q/icwb1iu5q_raw.fits   
    
    But s3://[bucket]/[[base_path]]/icwb/icwb1iu5q
        
    """
    import boto3
    import os
    import numpy as np
    
    s3 = boto3.resource('s3')
    s3_client = boto3.client('s3')
    bkt = s3.Bucket(bucket)
    
    files = [obj.key for obj in bkt.objects.filter(Prefix='Pipeline/{0}/Prep/'.format(root))]
    files += [obj.key for obj in bkt.objects.filter(Prefix='Pipeline/{0}/RAW/'.format(root))]
    
    extensions = ['flc.fits', 'flc.wcs.log.fits', 'flt.fits', 'flt.wcslog.fits', 'ramp.png', 'trails.png']
    
    copied_files = []
    
    for file in files:
        ext = file.split('_')[-1]
        if (ext in extensions) | ('mask.reg' in file):
            file_root = os.path.basename(file)
            prog = file_root[:4]
            dataset = file_root[:9]
            
            newfile = os.path.join(base_path, prog, dataset, file_root)
            
            if verbose:
                print('Rename {0} -> {1} on s3://{2}'.format(file, newfile, bucket))
            s3_client.copy_object(Bucket=bucket,
                                  CopySource=os.path.join(bucket, file), 
                                  Key=newfile, ACL='public-read')
                                  
            s3_client.delete_object(Bucket=bucket, Key=file)
            copied_files.append(file_root)
            
    # Visits file
    files = [obj.key for obj in bkt.objects.filter(Prefix='Pipeline/{0}/Prep/{0}_visits.npy'.format(root))]
    for file in files:
        print('Update awspath in {0}'.format(file))
        bkt.download_file(file, os.path.basename(file), 
                          ExtraArgs={"RequestPayer": "requester"})
        
        visits, groups, info = np.load(os.path.basename(file))  
        for v in visits:
            v['awspath'] = []
            for i, fi in enumerate(v['files']):
                file_root = os.path.basename(fi)
                prog = file_root[:4]
                dataset = file_root[:9]

                aws = os.path.join(bucket, base_path, prog, dataset)
                v['awspath'].append(aws)
        
        np.save(os.path.basename(file), [visits, groups, info])
        bkt.upload_file(os.path.basename(file), file, 
                          ExtraArgs={'RequestPayer': 'requester',
                                     'ACL': 'public-read'})
          
