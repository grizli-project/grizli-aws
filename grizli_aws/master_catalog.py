def summary(outroot='summary-18.05.17'):
    
    from hsaquery import overlaps
    from grizli import utils
    
    overlaps.summary_table(output='pointing_summary')
    
    tab = utils.GTable.gread('pointing_summary.fits')
    
    roots = tab['NAME']
    
    tab['Full'] = ['<a href=https://s3.amazonaws.com/aws-grivam/Pipeline/{0}/Extractions/{0}-full.html?chinu_max=2&bic_diff_min=30&zwidth1_max=0.01>Full</a>'.format(root.replace('+','%2B')) for root in roots]
    
    tab['fp'] = ['<a href=https://s3.amazonaws.com/aws-grivam/Pipeline/{0}_footprint.png> <img height=200 src=https://s3.amazonaws.com/aws-grivam/Pipeline/{0}_footprint.png></a>'.format(root.replace('+','%2B')) for root in roots]
    
    tab['zhist'] = ['<a href=https://s3.amazonaws.com/aws-grivam/Pipeline/{0}/Extractions/{0}_zhist.png> <img height=200 src=https://s3.amazonaws.com/aws-grivam/Pipeline/{0}/Extractions/{0}_zhist.png></a>'.format(root.replace('+','%2B')) for root in roots]
    
    #cols = ['NAME', 'RA', 'DEC', 'E(B-V)', 'GalLat', 'GalLon', 'NFILT', 'filter', 'target', 'target_description', 'proposal_id', 'pi_name', 'TexpG102', 'PAG102', 'TexpG141', 'PAG141', 'MAST', 'Full', 'fp', 'zhist']
    cols = ['NAME', 'RA', 'DEC', 'E(B-V)', 'GalLat', 'GalLon', 'NFILT', 'filter', 'target', 'proposal_id', 'pi_name', 'TexpG102', 'PAG102', 'TexpG141', 'PAG141', 'MAST', 'Full', 'fp', 'zhist']
    
    tab[cols].write_sortable_html(outroot+'.html', replace_braces=True, localhost=False, max_lines=50000, table_id=None, table_class='display compact', css=None, filter_columns=['mag_auto', 'z_map', 'bic_diff', 'chinu', 'zwidth1', 'is_point', 'sn_SIII', 'sn_Ha', 'sn_OIII', 'sn_Hb', 'sn_OII'], use_json=False)
    
    print('aws s3 cp {0}.html s3://aws-grivam/Pipeline/ --acl public-read'.format(outroot))
    
def regenerate_webpages(outroot='master'):
    """

    roots=`python -c "from grizli_aws import catalogs; roots, dates = catalogs.get_roots(verbose=False); print('\n'.join(roots))"`

    roots=`ls *phot.fits | sed "s/_phot/ /" | awk '{print $1}'`
    
    # ACS
    BUCKET=grizli-grism
    
    roots=`ls ../*footprint.fits | sed "s/\// /g" | sed "s/_foot/ /" | awk '{print $2}'`
    
    for root in $roots; do
        if [ ! -e ${root}.info.fits ]; then 
            aws s3 cp s3://${BUCKET}/Pipeline/${root}/Extractions/${root}.info.fits ./
            #aws s3 cp s3://${BUCKET}/Pipeline/${root}/Extractions/${root}_phot.fits ./
        else
            echo $root
        fi
    done
    
    for root in $roots; do
        aws s3 ls s3://aws-grivam/Pipeline/${root}/Extractions/ > /tmp/x
        echo $root
        grep beams.fits /tmp/x | wc
        grep full.fits /tmp/x | wc
    done
    

    """
    import glob
    import numpy as np
    
    from grizli import utils
    from astropy import table
    ## retrieve
    
    roots = [file.split('.info')[0] for file in glob.glob('*info.fits')]
    
    fpi = open('sync_to_s3.sh', 'w')
    fpi.write('# \n')
    fpi.close()
    
    tables = []
    
    for root in roots:
        fit = utils.GTable.gread(root+'.info.fits')
        phot = utils.GTable.gread(root+'_phot.fits')
        ix, dr = phot.match_to_catalog_sky(fit)
        fit['phot_dr'] = dr
        for c in phot.colnames:
            if c not in fit.colnames:
                fit[c] = phot[ix][c]
            else:
                pass
                #print(c)
                
        print(root, len(fit))
        
        # # DQ on mag differences
        # apcorr = fit['flux_auto']/fit['flux_aper_1']
        # m140 = 23.9-2.5*np.log10(fit['f140w_flux_aper_1']*apcorr)
        # dm140 = m140-fit['mag_wfc3,ir,f140w']
        # mask = fit['{0}_mask_aper_{1}'.format('f140w', 1)]
        # bad =  mask > 0.2*np.percentile(mask[np.isfinite(mask)], 99.5)
        # dm140[bad] = np.nan
        # 
        # m160 = 23.9-2.5*np.log10(fit['f160w_flux_aper_1']*apcorr)
        # dm160 = m160-fit['mag_wfc3,ir,f160w']
        # mask = fit['{0}_mask_aper_{1}'.format('f160w', 1)]
        # bad =  mask > 0.2*np.percentile(mask[np.isfinite(mask)], 99.5)
        # dm160[bad] = np.nan
        # 
        # dm = dm140
        # dm[(~np.isfinite(dm140)) & np.isfinite(dm160)] = dm160[(~np.isfinite(dm140)) & np.isfinite(dm160)]
        # dm[~np.isfinite(dm)] = -99
        # 
        # fit['fit_dmag'] = dm
        # fit['fit_dmag'].format = '4.1f'
        
        # Point source
        pa = np.polyfit([16, 21, 24, 25.5], [4.5, 2.6, 2.3, 2.1], 2)
        py = np.polyval(pa, fit['mag_auto'])
        point_source = (fit['flux_radius'] < py) #& (fit['mag_auto'] < 23)
        fit['is_point'] = point_source*2-1
        
        bad = (fit['mag_auto'] > 22) & (fit['flux_radius'] < 1.2)
        fit['too_small'] = bad*2-1
        
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
        
        
        fit['log_mass'] = np.log10(fit['stellar_mass'])
        fit['log_mass'].format = '.2f' 
        
        cols = ['root', 'idx','ra', 'dec', 't_g102', 't_g141', 'mag_auto', 'is_point', 'z_map', 'chinu', 'bic_diff', 'zwidth1', 'a_image', 'sn_SIII', 'sn_Ha', 'sn_OIII', 'sn_Hb', 'sn_OII', 'log_mass', 'png_stack', 'png_full', 'png_line']
        
        # Check grisms
        cols = ['root', 'idx','ra', 'dec', 'mag_auto', 'is_point', 'z_map', 'z02', 'z97', 'chinu', 'bic_diff', 'zwidth1', 'a_image', 'flux_radius', 'sn_SIII', 'sn_Ha', 'sn_OIII', 'sn_Hb', 'sn_OII', 'log_mass', 'aws_png_stack', 'aws_png_full', 'aws_png_line']
        for grism in ['g800l', 'g102', 'g141'][::-1]:
            if np.isfinite(fit['t_'+grism]).sum() > 0:
                cols.insert(4, 't_'+grism)
                            
        fit['ra'].format = '.4f'
        fit['dec'].format = '.4f'
        fit['z02'].format = '.2f'
        fit['z97'].format = '.2f'
        fit['flux_radius'].format = '.1f'
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
            
        fit[cols].write_sortable_html(root+'-full.html', replace_braces=True, localhost=False, max_lines=50000, table_id=None, table_class='display compact', css=None, filter_columns=['mag_auto', 'z_map', 'z02', 'z97', 'bic_diff', 'chinu', 'a_image', 'flux_radius', 'zwidth1', 'is_point', 'sn_SIII', 'sn_Ha', 'sn_OIII', 'sn_Hb', 'sn_OII'], use_json=True)
        
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
        
        tables.append(fit)

    master = table.vstack(tables)
    
    cols = ['root', 'idx','ra', 'dec', 'mag_auto', 'is_point', 'z_map', 'z02', 'z97', 'chinu', 'bic_diff', 'zwidth1', 'a_image', 'flux_radius', 'sn_SIII', 'sn_Ha', 'sn_OIII', 'sn_Hb', 'sn_OII', 'log_mass', 'aws_png_stack', 'aws_png_full', 'aws_png_line']
    for grism in ['g800l', 'g102', 'g141'][::-1]:
        if np.isfinite(master['t_'+grism]).sum() > 0:
            cols.insert(4, 't_'+grism)
    
            
    master[cols].write_sortable_html(outroot+'.html', replace_braces=True, localhost=False, max_lines=500000, table_id=None, table_class='display compact', css=None, filter_columns=['mag_auto', 't_g102', 't_g141', 'z_map', 'z02', 'z97', 'bic_diff', 'chinu', 'a_image', 'flux_radius', 'zwidth1', 'is_point', 'sn_SIII', 'sn_Ha', 'sn_OIII', 'sn_Hb', 'sn_OII', 'log_mass'], use_json=True)
    
    #sel = master['bic_diff'] > 30
    #master[sel][cols].write_sortable_html(outroot+'.html', replace_braces=True, localhost=False, max_lines=500000, table_id=None, table_class='display compact', css=None, filter_columns=['mag_auto', 't_g102', 't_g141', 'z_map', 'z02', 'z97', 'bic_diff', 'chinu', 'a_image', 'flux_radius', 'zwidth1', 'is_point', 'sn_SIII', 'sn_Ha', 'sn_OIII', 'sn_Hb', 'sn_OII', 'log_mass'], use_json=True)
    
    new = utils.GTable()
    for col in fit.colnames:
        new[col] = master[col]
        
    new.write(outroot+'.fits', overwrite=True)
    
    print('aws s3 cp {0}.html s3://aws-grivam/Pipeline/ --acl public-read\n'.format(outroot))
    print('aws s3 cp {0}.json s3://aws-grivam/Pipeline/ --acl public-read\n'.format(outroot))
    print('aws s3 cp {0}.fits s3://aws-grivam/Pipeline/ --acl public-read\n'.format(outroot))
    
    
    print ('bash sync_to_s3.sh')

def master_catalog(outroot='grizli-18.05.17-full'):
    
    import glob
    import numpy as np
    
    from grizli import utils
    from astropy import table
    ## retrieve
    
    roots = [file.split('.info')[0] for file in glob.glob('*info.fits')]
    
    tabs = [utils.GTable.gread(root+'.info.fits') for root in roots]
    fit = utils.GTable(table.vstack(tabs))
    
    fit['flux_radius'].format = '.1f'
    
    N = len(fit)
    aws_col = {}
    for c in ['png_stack', 'png_full', 'png_line']:
        aws_col[c] = []
    
    bucket='aws-grivam'
    bucket='grizli-grism'
        
    for i in range(N):
        root = fit['root'][i]
        
        aws = 'https://s3.amazonaws.com/{0}/Pipeline/{1}/Extractions/{1}'.format(bucket, root)
        
        for c in ['png_stack', 'png_full', 'png_line']:
            pre = fit[c][i]
            aws_col[c].append(pre.replace(root, aws).replace(root, root.replace('+','%2B')))
            
    for c in ['png_stack', 'png_full', 'png_line']:
        fit['aws_'+c] = aws_col[c]
    
    fit['aws_png_sed'] = [l.replace('full.png','sed.png') for l in fit['aws_png_full']]
    
    psx = [16, 18.5, 21, 24, 25.5]
    psy = [8, 4., 2.6, 2.3, 2.1]
    
    # New pixel scale
    psy = np.array(psy)*0.06/0.1
    
    # ACS
    if False:
        psx = [16, 21, 24, 25.5, 27]
        psy = [4.5, 2.9, 2.9, 2.6, 2.6]
        
    pa = np.polyfit(psx, psy, 2)
    py = np.interp(fit['mag_auto'], psx, psy)

    # Point source
    point_source = (fit['flux_radius'] < py) & (fit['mag_auto'] < 25.5)
    #point_source = (fit['flux_radius'] < py) & (fit['mag_auto'] < 24) # ACS
    fit['is_point'] = point_source*1
    
    fit['log_mass'] = np.log10(fit['stellar_mass'])
    fit['log_mass'].format = '.2f' 
    
    ############
    # Warnings for ambiguous line IDs
    dz_line = np.array([5007, 6563])*(1./5007.-1/6563.)
    col = 'ambiguous_HaOIII'

    dz_line = np.array([3727, 6563])*(1./3727.-1/6563.)
    col = 'ambiguous_HaOII'
    
    zw1 = fit['zwidth1']/(1+fit['z_map'])
    zw2 = fit['zwidth2']/(1+fit['z_map'])
    fit['zw1'] = zw1
    fit['zw1'].format = '.3f'
    fit['zw2'] = zw2
    fit['zw2'].format = '.3f'
    
    for dz_line, col in zip([np.array([5007, 6563])*(1./5007.-1/6563.), np.array([3727, 6563])*(1./3727.-1/6563.)], ['ambiguous_HaOIII', 'ambiguous_HaOII']):
        if col in fit.colnames: fit.remove_column(col)

        for dz_i in dz_line:
            if col not in fit.colnames:
                fit[col] = np.abs(zw1 - dz_i) < 0.008
            else:
                fit[col] |= np.abs(zw1 - dz_i) < 0.008

    ############
    # Reliable redshifts, based on 3D-HST COSMOS
    ambiguous = (fit['ambiguous_HaOIII'] | fit['ambiguous_HaOII']) & (zw2 / zw1 < 1.1)
    fit['ok_width'] = (fit['zw1'] < 0.02) | ambiguous
    fit['ok_width'].format = 'd'
    
    # Overall data quality
    fit['use_spec'] = (fit['ok_width']) & (fit['chinu'] < 10) & (fit['bic_diff'] > 0) & (fit['flux_radius'] > 1)

    # ACS GLASS
    if False:
        fit['use_spec'] &= (fit['flux_radius'] > 1.9)
    
    fit['use_spec'].format = 'd'
    
    fit['z_map'].format = '.4f'
    fit['chinu'].format = '.2f'
    
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
    
    cols = ['root', 'idx','ra', 'dec', 't_g800l', 't_g102', 't_g141', 'mag_auto', 'use_spec', 'is_point', 'z_map', 'chinu', 'bic_diff', 'zw1', 'ok_width', 'flux_radius', 'sn_SIII', 'sn_Ha', 'sn_OIII', 'sn_Hb', 'sn_OII', 'log_mass', 'aws_png_stack', 'aws_png_full', 'aws_png_line']
    
    for g in ['g102', 'g141', 'g800l']:
        if 't_'+g in fit.colnames:
            bad = ~np.isfinite(fit['t_'+g])
            fit['t_'+g][bad] = 0.
            
            if fit['t_'+g].max() == 0:
                pop = True
            else:
                pop = False
        else:
            pop = True
        
        if pop:
            cols.pop(cols.index('t_'+g))
    
    if 'sparcs' in outroot:
        cols += ['aws_png_sed']
        
    filter_cols = ['mag_auto', 'z_map', 'z02', 'z97', 'bic_diff', 'chinu', 'flux_radius', 'zw1', 'use_spec', 'is_point', 'sn_SIII', 'sn_Ha', 'sn_OIII', 'sn_Hb', 'sn_OII', 't_g800l', 't_g102', 't_g141']
    
    if False:
        clip = (fit['bic_diff'] > 20) & (fit['chinu'] < 2) & (fit['zwidth1'] < 0.01)
    else:
        clip = fit['mag_auto'] > 0
        
    fit[cols][clip].filled(fill_value=-1).write_sortable_html(outroot+'.html', replace_braces=True, localhost=False, max_lines=50000, table_id=None, table_class='display compact', css=None, filter_columns=filter_cols, use_json=True)

    print('aws s3 cp {0}.html s3://aws-grivam/Pipeline/ --acl public-read\n'.format(outroot))
    print('aws s3 cp {0}.json s3://aws-grivam/Pipeline/ --acl public-read\n'.format(outroot))
        
    new = utils.GTable()
    for col in fit.colnames:
        new[col] = fit[col][clip]
        
    new.write(outroot+'.fits', overwrite=True)
    
    print('aws s3 cp {0}.fits s3://aws-grivam/Pipeline/ --acl public-read\n'.format(outroot))
    
    ######################
    # use_spec versions
    outroot += '.use_spec'
    clip &= (fit['use_spec'])

def summary_table(output_table='summary_glass-acs-2018.05.21'
):
    from hsaquery import overlaps
    from grizli import utils
    
    overlaps.summary_table(output=output_table)
    
    tab = utils.GTable.gread(output_table+'.fits')
    
    tab['Browse'] = ['<a href=https://s3.amazonaws.com/aws-grivam/Pipeline/{0}/Extractions/{0}-full.html>Browse</a>'.format(root.replace('+', '%2B')) for root in tab['NAME']]
    
    pixscale = np.array([0.06]*len(tab))
    pixscale[9:] = 0.03
    tab['driz_scale'] = pixscale
    
    url = "http://archive.stsci.edu/hst/search.php?action=Search&RA={0}&DEC={1}&radius=5.&sci_aec=S"
    root_url = [('< a href='+url+'>{2}</a>').format(line['RA'], line['DEC'], line['NAME']) for line in tab]
    
    tab['root'] = root_url
    
    cols = ['NAME', 'RA', 'DEC', 'MW_EBV', 'filter', 'proposal_id', 'proposal_pi', 'driz_scale', 'MAST', 'Browse']
    tab[cols].write_sortable_html(output_table+'.html', use_json=False, localhost=False, max_lines=1000)
    
    print('aws s3 cp {0}.html s3://aws-grivam/Pipeline/ --acl public-read'.format(output_table))
    

    