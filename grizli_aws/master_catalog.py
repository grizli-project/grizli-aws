def summary():
    
    from hsaquery import overlaps
    
    overlaps.summary_table(output='pointing_summary')
    
    tab = utils.GTable.gread('pointing_summary.fits')
    
    roots = tab['NAME']
    
    tab['Full'] = ['<a href=https://s3.amazonaws.com/aws-grivam/Pipeline/{0}/Extractions/{0}-full.html?chinu_max=2&bic_diff_min=30&zwidth1_max=0.01>Full</a>'.format(root.replace('+','%2B')) for root in roots]
    
    tab['fp'] = ['<a href=https://s3.amazonaws.com/aws-grivam/Pipeline/{0}_footprint.png> <img height=200 src=https://s3.amazonaws.com/aws-grivam/Pipeline/{0}_footprint.png></a>'.format(root.replace('+','%2B')) for root in roots]
    
    tab['zhist'] = ['<a href=https://s3.amazonaws.com/aws-grivam/Pipeline/{0}/Extractions/{0}_zhist.png> <img height=200 src=https://s3.amazonaws.com/aws-grivam/Pipeline/{0}/Extractions/{0}_zhist.png></a>'.format(root.replace('+','%2B')) for root in roots]
    
    cols = ['NAME', 'RA', 'DEC', 'E(B-V)', 'GalLat', 'GalLon', 'NFILT', 'filter', 'target', 'target_description', 'proposal_id', 'pi_name', 'TexpG102', 'PAG102', 'TexpG141', 'PAG141', 'MAST', 'Full', 'fp', 'zhist']
    
    tab[cols].write_sortable_html('summary-18.05.17.html', replace_braces=True, localhost=False, max_lines=50000, table_id=None, table_class='display compact', css=None, filter_columns=['mag_auto', 'z_map', 'bic_diff', 'chinu', 'zwidth1', 'is_point', 'sn_SIII', 'sn_Ha', 'sn_OIII', 'sn_Hb', 'sn_OII'], use_json=False)
    
    print('aws s3 cp summary-18.05.17.html s3://aws-grivam/Pipeline/ --acl public-read')
    
def regenerate_webpages():
    """

    roots=`python -c "from grizli_aws import catalogs; roots, dates = catalogs.get_roots(verbose=False); print('\n'.join(roots))"`
    
    # ACS
    roots=`ls ../*footprint.fits | sed "s/\// /g" | sed "s/_foot/ /" | awk '{print $2}'`
    
    for root in $roots; do
        aws s3 cp s3://aws-grivam/Pipeline/${root}/Extractions/${root}.info.fits ./
    done

    """
    import glob
    import numpy as np
    
    from grizli import utils
    
    ## retrieve
    
    roots = [file.split('.info')[0] for file in glob.glob('*info.fits')]
    
    fpi = open('sync_to_s3.sh', 'w')
    fpi.write('# \n')
    fpi.close()
        
    for root in roots:
        print(root)
        fit = utils.GTable.gread(root+'.info.fits')
        
        pa = np.polyfit([16, 21, 24, 25.5], [4.5, 2.6, 2.3, 2.1], 2)
        py = np.polyval(pa, fit['mag_auto'])

        # Point source
        point_source = (fit['flux_radius'] < py) & (fit['mag_auto'] < 23)
        fit['is_point'] = point_source*1
        
        fit['log_mass'] = np.log10(fit['stellar_mass'])
        fit['log_mass'].format = '.2f' 
        
        cols = ['root', 'idx','ra', 'dec', 't_g102', 't_g141', 'mag_auto', 'is_point', 'z_map', 'chinu', 'bic_diff', 'zwidth1', 'a_image', 'sn_SIII', 'sn_Ha', 'sn_OIII', 'sn_Hb', 'sn_OII', 'log_mass', 'png_stack', 'png_full', 'png_line']

        # ACS
        if 't_g800l' in fit.colnames:
            cols = ['root', 'idx','ra', 'dec', 't_g800l', 'mag_auto', 'is_point', 'z_map', 'z02', 'z97', 'chinu', 'bic_diff', 'zwidth1', 'a_image', 'sn_SIII', 'sn_Ha', 'sn_OIII', 'sn_Hb', 'sn_OII', 'log_mass', 'png_stack', 'png_full', 'png_line']
        else:
            cols = ['root', 'idx','ra', 'dec', 't_g102', 't_g141', 'mag_auto', 'is_point', 'z_map', 'z02', '97', 'chinu', 'bic_diff', 'zwidth1', 'a_image', 'sn_SIII', 'sn_Ha', 'sn_OIII', 'sn_Hb', 'sn_OII', 'log_mass', 'png_stack', 'png_full', 'png_line']
            
        fit['ra'].format = '.4f'
        fit['dec'].format = '.4f'
        fit['z02'].format = '.2f'
        fit['z97'].format = '.2f'
        fit['mag_auto'].format = '.2f'
        fit['t_g800l'].format = '.0f'
        fit['t_g102'].format = '.0f'
        fit['t_g141'].format = '.0f'
        fit['zq'].format = '.1f'
        fit['zwidth1'].format = '.3f'
        fit['bic_diff'].format = '.0f'
        fit['a_image'].format = '.1f'
        for l in ['Ha','OIII','Hb','OII','SIII']:
            fit['sn_'+l].format = '.1f'
            
        fit[cols].write_sortable_html(root+'-full.html', replace_braces=True, localhost=False, max_lines=50000, table_id=None, table_class='display compact', css=None, filter_columns=['mag_auto', 'z_map', 'z02', 'z97', 'bic_diff', 'chinu', 'zwidth1', 'is_point', 'sn_SIII', 'sn_Ha', 'sn_OIII', 'sn_Hb', 'sn_OII'], use_json=True)
        
        if '+' in root:
            ext = 'html'
            lines = open(root+'-full.'+ext).readlines()
            for i, line in enumerate(lines):
                if root in line:
                    lines[i] = line.replace(root, root.replace('+', '%2B'))
        
            fp = open(root+'-full.'+ext,'w')
            fp.writelines(lines)
            fp.close()
            
            ext = 'json'
            lines = open(root+'-full.'+ext).readlines()
            for i, line in enumerate(lines):
                if (root in line) & ('href' in line):
                    lines[i] = line.replace(root, root.replace('+', '%2B'))
        
            fp = open(root+'-full.'+ext,'w')
            fp.writelines(lines)
            fp.close()
            
            
        fpi = open('sync_to_s3.sh', 'a')
        fpi.write('aws s3 cp {0}-full.html s3://aws-grivam/Pipeline/{0}/Extractions/ --acl public-read\n'.format(root))
        fpi.write('aws s3 cp {0}-full.json s3://aws-grivam/Pipeline/{0}/Extractions/ --acl public-read\n'.format(root))
        fpi.close()
        
    fpi.close()
    print ('bash sync_to_s3.sh')

def master_catalog(outroot='grizli-18.05.17-full'):
    
    import glob
    import numpy as np
    
    from grizli import utils
    
    ## retrieve
    
    roots = [file.split('.info')[0] for file in glob.glob('*info.fits')]
    
    tabs = [utils.GTable.gread(root+'.info.fits') for root in roots]
    fit = utils.GTable(table.vstack(tabs))
    
    N = len(fit)
    aws_col = {}
    for c in ['png_stack', 'png_full', 'png_line']:
        aws_col[c] = []
        
    for i in range(N):
        root = fit['root'][i]
        
        aws = 'https://s3.amazonaws.com/aws-grivam/Pipeline/{0}/Extractions/{0}'.format(root)
        
        for c in ['png_stack', 'png_full', 'png_line']:
            pre = fit[c][i]
            aws_col[c].append(pre.replace(root, aws).replace(root, root.replace('+','%2B')))
            
    for c in ['png_stack', 'png_full', 'png_line']:
        fit['aws_'+c] = aws_col[c]
    
    
    pa = np.polyfit([16, 21, 24, 25.5], [4.5, 2.6, 2.3, 2.1], 2)
    py = np.polyval(pa, fit['mag_auto'])

    # Point source
    point_source = (fit['flux_radius'] < py) & (fit['mag_auto'] < 23)
    point_source = (fit['flux_radius'] < py) & (fit['mag_auto'] < 24) # ACS
    fit['is_point'] = point_source*1
    
    fit['log_mass'] = np.log10(fit['stellar_mass'])
    fit['log_mass'].format = '.2f' 
    
    fit['ra'].format = '.4f'
    fit['dec'].format = '.4f'
    fit['mag_auto'].format = '.2f'
    fit['t_g102'].format = '.0f'
    fit['t_g141'].format = '.0f'
    fit['zq'].format = '.1f'
    fit['zwidth1'].format = '.3f'
    fit['bic_diff'].format = '.0f'
    fit['a_image'].format = '.1f'
    for l in ['Ha','OIII','Hb','OII','SIII']:
        fit['sn_'+l].format = '.1f'
        
    if 't_g800l' in fit.colnames:
        cols = ['root', 'idx','ra', 'dec', 't_g800l', 'mag_auto', 'is_point', 'z_map', 'chinu', 'bic_diff', 'zwidth1', 'a_image', 'sn_SIII', 'sn_Ha', 'sn_OIII', 'sn_Hb', 'sn_OII', 'log_mass', 'aws_png_stack', 'aws_png_full', 'aws_png_line']
    else:
        cols = ['root', 'idx','ra', 'dec', 't_g102', 't_g141', 'mag_auto', 'is_point', 'z_map', 'chinu', 'bic_diff', 'zwidth1', 'a_image', 'sn_SIII', 'sn_Ha', 'sn_OIII', 'sn_Hb', 'sn_OII', 'log_mass', 'aws_png_stack', 'aws_png_full', 'aws_png_line']
        
    
    if False:
        clip = (fit['bic_diff'] > 20) & (fit['chinu'] < 2) & (fit['zwidth1'] < 0.01)
    else:
        clip = fit['mag_auto'] > 0
        
    fit[cols][clip].filled(fill_value=-1).write_sortable_html(outroot+'.html', replace_braces=True, localhost=False, max_lines=50000, table_id=None, table_class='display compact', css=None, filter_columns=['mag_auto', 'z_map', 'bic_diff', 'chinu', 'zwidth1', 'is_point', 'sn_SIII', 'sn_Ha', 'sn_OIII', 'sn_Hb', 'sn_OII', 'log_mass', 'a_image'], use_json=True)
    
    new = utils.GTable()
    for col in fit.colnames:
        new[col] = fit[col][clip]
        
    new.write(outroot+'.fits', overwrite=True)
    
    print('aws s3 cp {0}.fits s3://aws-grivam/Pipeline/ --acl public-read\n'.format(outroot))
    print('aws s3 cp {0}.html s3://aws-grivam/Pipeline/ --acl public-read\n'.format(outroot))
    print('aws s3 cp {0}.json s3://aws-grivam/Pipeline/ --acl public-read\n'.format(outroot))

def summary_table():
    from hsaquery import overlaps
    
    output_table = 'summary_glass-acs-2018.05.21'
    overlaps.summary_table(output=output_table)
    
    tab = utils.GTable.gread(output_table+'.fits')
    
    tab['Browse'] = ['<a href=https://s3.amazonaws.com/aws-grivam/Pipeline/{0}/Extractions/{0}-full.html>Browse</a>'.format(root.replace('+', '%2B')) for root in tab['NAME']]
    
    pixscale = np.array([0.06]*len(tab))
    pixscale[9:] = 0.03
    tab['driz_scale'] = pixscale
    
    cols = ['NAME', 'RA', 'DEC', 'E(B-V)', 'filter', 'proposal_id', 'pi_name', 'driz_scale', 'MAST', 'Browse']
    tab[cols].write_sortable_html(output_table+'.html', use_json=False, localhost=False, max_lines=1000)
    
    print('aws s3 cp {0}.html s3://aws-grivam/Pipeline/ --acl public-read'.format(output_table))
    

    