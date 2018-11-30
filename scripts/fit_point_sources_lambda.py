#!/usr/bin/env python
def fit_quasars_lambda(root='j100025+021706', newfunc=True):
    import time
    import os
    import numpy as np
    import boto3
    import json
    from astropy import table
    
    from grizli import utils
    
    s3 = boto3.resource('s3')
    bkt = s3.Bucket('aws-grivam')
    #bkt.download_file('/Pipeline/wisps-aug10.fits', 'wisps-aug10', ExtraArgs={"RequestPayer": "requester"})
    
    fit = utils.read_catalog('wisps-aug10.fits')
    
    fit = table.vstack([utils.read_catalog('grizli-3dhst-18.08.05.fits'), utils.read_catalog('grizli-18.05.17-full.fits')])

    # All 
    fit = table.vstack([utils.read_catalog('wisps-aug10.fits'), utils.read_catalog('grizli-3dhst-18.08.05.fits'), utils.read_catalog('grizli-18.05.17-full.fits')])
    
    fit = utils.read_catalog('grizli-3dhst-18.08.05.fits')
    fit = utils.read_catalog('grizli-18.05.17-full.fits')
    
    roots=np.unique(fit['root'])
        
    pa = np.polyfit([16, 21, 24, 25.5], [4.5, 2.6, 2.3, 2.1], 2)
    py = np.polyval(pa, fit['mag_auto'])
    py = np.interp(fit['mag_auto'], [16, 18, 19, 21, 24, 25.5], [7, 7, 5, 2.6, 2.3, 2.1])
    point_source = (fit['flux_radius'] < py) #& (fit['mag_auto'] < 23)
    fit['is_point'] = point_source*2-1
    
    bad = (fit['mag_auto'] > 22) & (fit['flux_radius'] < 1.2)
    fit['too_small'] = bad*2-1
    
    stars = (fit['is_point'] > 0) & (fit['too_small'] < 0) & (fit['mag_auto'] < 24.5)
    stars &= fit['ninput'] < 40
    
    beams = ['Pipeline/{0}/Extractions/{0}_{1:05d}.beams.fits'.format(root, id) for root, id in zip(fit['root'][stars], fit['id'][stars])]
    ninput = fit['ninput'][stars]
    
    # Auth to create a Lambda function (credentials are picked up from above .aws/credentials)
    session = boto3.Session()
        
    # Make sure Lambda is running in the same region as the HST public dataset
    client = session.client('lambda', region_name='us-east-1')

    func = 'GrizliFitQuasar'
    
    print ('Lambda function: {0}'.format(func))
               
    # Auth to create a Lambda function 
    session = boto3.Session()
    client = session.client('lambda', region_name='us-east-1')

    done_files = [os.path.basename(obj.key) for obj in bkt.objects.filter(Prefix='Pipeline/QuasarFit/')]
    
    count = 0
    
    for beam, ninp in zip(beams, ninput):        
        # if count > 400:
        #     break
        
        if os.path.basename(beam).replace('.beams','.full') in done_files:
            #print('{0} {1} {2}'.format(count, 'Skip', beam))
            continue
        else:
            count += 1
            print('{0} {1} {2} ({3})'.format(count, '    ', beam, ninp))
            
        event = {
              's3_object_path': beam,
              'verbose': "True",
              "use_psf": str(ninp < 15),
            }

        # Invoke Lambda function
        response = client.invoke(
            FunctionName=func,
            InvocationType='Event',
            LogType='Tail',
            Payload=json.dumps(event))
    
    files=glob.glob('*full.fits')
    roots=np.unique([file.split('_')[0] for file in files])
    for root in roots:
        if os.path.exists('{0}.info.fits'.format(root)):
            continue
        
        try:
            auto_script.summary_catalog(field_root=root, dzbin=0.01, use_localhost=False, filter_bandpasses=None)
        except:
            pass
    
    #
    #fit = utils.read_catalog('wisps-aug10.fits')
    
    from astropy import table
    files = glob.glob('*info.fits')
    roots=np.unique([file.split('.info')[0] for file in files])
    
    tabs = []
    for root in roots:
        print(root)
        file='{0}.info.fits'.format(root)
        if os.path.exists(file):
            tabs.append(utils.read_catalog(file))
        
    tab = table.vstack(tabs)
    
    tab['bic_diff_spl'] = tab['bic_spl'] - tab['bic_temp']
    tab['bic_diff_spl'].format = '.1f'

    tab['chinusp'] = tab['chi2spl']/tab['dof']
    tab['chinusp'].format = '.1f'

    tab['star14'] = tab['splf03'] / np.maximum(tab['splf04'], tab['sple04']) 
    tab['star14'].format = '.1f'

    tab['star14e'] = tab['splf03']/np.sqrt(tab['sple03']**2+tab['sple04']**2)
    tab['star14e'].format = '.2f'
    
    tab['spl14'] = tab['splf03']/tab['sple03']
    tab['spl14'].format = '.0f'
    
    cols = ['root', 'idx','ra', 'dec', 't_g102', 't_g141', 'mag_auto', 'is_point', 'z_map', 'chinu', 'bic_diff', 'chinusp', 'bic_diff_spl', 'star14', 'star14e', 'spl14', 'zwidth1', 'a_image', 'sn_SIII', 'sn_Ha', 'sn_OIII', 'sn_Hb', 'sn_OII', 'log_mass', 'aws_png_stack', 'png_full', 'png_line']
    
    newcol = {}
    for col in ['png_full','png_line']:
        newcol[col] = [item.replace('+','%2B') for item in tab[col]]
        tab.remove_column(col)
        tab[col] = newcol[col]
        
    idx, dr = fit.match_to_catalog_sky(tab)
    for c in cols:
        if c not in tab.colnames:
            tab[c] = fit[c][idx]
    
    clip = tab['bic_diff'] > 0
    
    outroot='point_sources_aug15'
    
    tab.write(outroot+'.fits', overwrite=True)
    
    tab[cols].filled(fill_value=-1).write_sortable_html(outroot+'.html', replace_braces=True, localhost=False, max_lines=50000, table_id=None, table_class='display compact', css=None, filter_columns=['mag_auto', 'z_map', 'z02', 'z97', 'bic_diff', 'chinusp', 'star14', 'star14e', 'spl14', 'bic_diff_spl', 'chinu', 'a_image', 'flux_radius', 'zwidth1', 'is_point', 'sn_SIII', 'sn_Ha', 'sn_OIII', 'sn_Hb', 'sn_OII'], use_json=True)
            
    sleep_time = 303*np.ceil(len(beams)/950)
    print('{0}: sleep {1}'.format(time.ctime(), sleep_time))
    
    time.sleep(sleep_time)
    
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

    for file in files:
        if 'beams.fits' in file:
            beams.append(file)

        if 'log_par' in file:
            logs.append(file)

        if 'full.fits' in file:
            full.append(file)
    
    label = '{4} / {3} / Nbeams: {0}, Nfull: {1}, Nlog: {2}'.format(len(beams), len(full), len(logs), root, time.ctime())
    if get_string:
        return label
        
    print(label)
    
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
    