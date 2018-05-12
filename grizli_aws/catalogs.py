#!/usr/bin/env python
# encoding: utf-8
"""
make_summary.py

Created by Gabriel Brammer on 2018-05-09.
Copyright (c) 2018 __MyCompanyName__. All rights reserved.
"""

def get_roots(verbose=True):
    
    import os
    from grizli import utils
    import boto3
    
    s3 = boto3.resource('s3')
    
    bucket = s3.Bucket('aws-grivam')
    
    roots = []
    dates = []
    
    footprint_keys = bucket.objects.filter(Prefix='Pipeline/j', Delimiter='/')
    for obj in footprint_keys:
        if 'footprint.fits' not in obj.key:
            continue
            
        #print(obj.key)
        root=obj.key[9:23]
        zhist_list = bucket.objects.filter(Prefix='Pipeline/{0}/Extractions/{0}_zhist.png'.format(root))
        for zh_obj in zhist_list:
            if zh_obj.key:
                date = zh_obj.last_modified.isoformat().split('+00')[0]
                if verbose:
                    print(zh_obj.key, root, date)
                    
                roots.append(root)
                dates.append(date)
    
    return roots, dates
    
def make_summary_html(output_table='grizli_aws.html', verbose=True):
    import os
    from grizli import utils
    import boto3
    
    s3 = boto3.resource('s3')
    
    bucket = s3.Bucket('aws-grivam')

    roots, dates = get_roots(verbose=verbose)
    
    # Table
    tab = utils.GTable()
    tab['root'] = roots
    tab['modified'] = dates
    
    tab['Full'] = ['<a href=https://s3.amazonaws.com/aws-grivam/Pipeline/{0}/Extractions/{0}-fit.html>Full</a>'.format(root.replace('+','%2B')) for root in roots]
    
    tab['Best'] = ['<a href=https://s3.amazonaws.com/aws-grivam/Pipeline/{0}/Extractions/{0}-fit.zq.html>Best</a>'.format(root.replace('+','%2B')) for root in roots]
    
    #tab['footprint'] = ['<a href=https://s3.amazonaws.com/aws-grivam/Pipeline/{0}_footprint.pdf> <img width=400 src=https://s3.amazonaws.com/aws-grivam/Pipeline/{0}_footprint.pdf></a>'.format(root.replace('+','%2B')) for root in roots]
    
    tab['zhist'] = ['<a href=https://s3.amazonaws.com/aws-grivam/Pipeline/{0}/Extractions/{0}_zhist.png> <img width=400 src=https://s3.amazonaws.com/aws-grivam/Pipeline/{0}/Extractions/{0}_zhist.png></a>'.format(root.replace('+','%2B')) for root in roots]
    
    tab.write_sortable_html(output_table, localhost=False)
    print("""
Done: {0}

Sync to s3: aws s3 cp {0} s3://aws-grivam/Pipeline/ --acl public-read
""".format(output_table))
    
    return tab
    
def fix_html_tables(verbose=True):
    """
    Replace + with %2B in html files
    """
    
    roots, dates = get_roots(verbose=verbose)
    
    import os
    from grizli import utils
    import boto3
    
    s3 = boto3.resource('s3')
    
    bucket = s3.Bucket('aws-grivam')
    
    if False:
        for root in roots:
            if '+' not in root:
                continue

            print("aws s3 cp {0}-fit.html s3://aws-grivam/Pipeline/{0}/Extractions/ --acl public-read".format(root))
            print("aws s3 cp {0}-fit.zq.html s3://aws-grivam/Pipeline/{0}/Extractions/ --acl public-read".format(root))
    
    for root in roots:
        if '+' not in root:
            continue
        
        for ext in ['fit', 'fit.zq']:
            html_file = './{0}-{1}.html'.format(root, ext)
            Key = 'Pipeline/{0}/Extractions/{0}-{1}.html'.format(root, ext)
            #print(html_file, Key)
            
            s3.meta.client.download_file('aws-grivam', Key, html_file)
        
            lines = open(html_file).readlines()
            new_lines = [line.replace(root, root.replace('+','%2B')) for line in lines] 
            fp = open(html_file,'w')
            fp.writelines(new_lines)
            fp.close()
            
            print("aws s3 cp {0} s3://aws-grivam/Pipeline/{1}/Extractions/ --acl public-read".format(html_file, root))
            
            #s3.meta.client.upload_file(html_file, 'aws-grivam', Key, ExtraArgs={'ACL':'public-read'})
            
             
            
    