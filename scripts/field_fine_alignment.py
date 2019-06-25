import os
import glob
import numpy as np

from grizli import prep, utils
from grizli.pipeline import auto_script
import astropy.units as u

shell = """

# Misc shell things

# clean out ACS
aws s3 rm --recursive --exclude "*" --include "*wfc3uvis*" --include "*acswfc*" --exclude "*wfc3ir*" s3://grizli-v1/Pipeline/

"""

def go():

    os.chdir('/home/ec2-user/Mosaics')
    
    root = 'j123656p6215'

    root = 'j021732m0512'
    
    root = 'j141956p5255'
    
    root = 'j033236m2748'
    
    # Sync needed files for GOODSN
    if root == 'j123656p6215':
        # os.system('aws s3 sync --exclude "*" --include "{0}*/Prep/*expflag*" --include "{0}*/Prep/*fail*" --include "{0}*/Prep/*cat.fits" --include "{0}*/Prep/{0}_*dr*fits.gz" --include "{0}*/Prep/*visits.npy" --include "{0}*/Prep/*[._]wcs*" --include "{0}*/Prep/*shifts.*" s3://grizli/Pipeline/ .'.format(root))
        # os.system('aws s3 sync --exclude "*" --include "{0}*/Prep/*wcs.log" --include "{0}*/Prep/*shifts.*" s3://grizli/Pipeline/ .'.format(root))

        os.system('aws s3 sync --exclude "*" --include "{0}*/Prep/*expflag*" --include "{0}*/Prep/*fail*" --include "{0}*/Prep/*cat.fits" --include "{0}*/Prep/*visits.npy" --include "{0}*/Prep/*[._]wcs*" --include "{0}*/Prep/*shifts.*" s3://grizli/Pipeline/ .'.format(root))

    else:
        #os.system('aws s3 sync --exclude "*" --include "{0}*/Prep/*expflag*" --include "{0}*/Prep/*fail*" --include "{0}*/Prep/*cat.fits" --include "{0}*/Prep/*visits.npy" --include "{0}*/Prep/*[._]wcs*" --include "{0}*/Prep/*_fl?.*" --include "{0}*/Prep/*shifts.*" s3://grizli-v1/Pipeline/ .'.format(root))
        os.system('aws s3 sync --exclude "*" --include "{0}*/Prep/*expflag*" --include "{0}*/Prep/*fail*" --include "{0}*/Prep/*cat.fits" --include "{0}*/Prep/*visits.npy" --include "{0}*/Prep/*[._]wcs*" --include "{0}*/Prep/*shifts.*" s3://grizli-v1/Pipeline/ .'.format(root))
    
    #os.system('files=`find . |grep "dr[cz]_" |grep fits.gz |grep -v "\-ir_dr"`; for file in $files; do echo $file; gunzip -f $file; done')
    
    # Remake catalogs
    mos_files = glob.glob('{0}*/Prep/*_dr*sci.fits'.format(root))
    mos_files.sort()
    
    thresh = 10
    force=False
    
    for i, file in enumerate(mos_files):
        root = file.split('_drc_sci')[0].split('_drz_sci')[0]
        print(i, root)
        
        cat_file = root+'.cat.fits'
        if (not os.path.exists(cat_file)) | (force):
            prep.make_SEP_catalog(root=root, threshold=thresh, get_background=True, phot_apertures=[1*u.arcsec])
    
    os.chdir('Fine/')
    os.system('rsync -avz  ../{0}*/Prep/*cat.fits .'.format(root))
    os.system('ln -s ../{0}*/Prep/*_dr?_sci.fits .'.format(root))
    os.system('rsync -avz  ../{0}*/Prep/*visits.npy .'.format(root))
    os.system('rsync -avz  ../{0}*/Prep/*fail* .'.format(root))
    
    visit_files = glob.glob('*visits.npy')

    #visit_files = glob.glob('*-tile-*visits.npy')

    visit_files.sort()
    
    filt = 'f850lp'
    
    all_visits = []
    for f in visit_files:
        root_i = f.split('_visits.npy')[0]
        cat_file = root_i+'-'+filt+'.cat.fits'
        if not os.path.exists(cat_file):
            continue
            
        visits, _, _ = np.load(f, allow_pickle=True)
        has_failed=False
        for v in visits:
            has_failed |= os.path.exists(v['product']+'.failed')
        
        if has_failed:
            continue
        
        all_visits.append({'product':root_i+'-'+filt})
        
    # Fine alignment
    radec, ref_catalog = prep.get_radec_catalog(ra=189.2316435, 
            dec=62.249739865958, 
            product='xxx',
            reference_catalogs=['PS1'], radius=30)
    
    auto_script.fine_alignment(field_root='xxx', all_visits=all_visits, gaia_by_date=False, catalogs=['PS1'], redrizzle=False, shift_only=False, radec=radec, tol=1.e-3)
    
    ##### Done
    import numpy as np
    import os
    import glob
    from grizli import prep, utils
    from drizzlepac import updatehdr
    import astropy.io.fits as pyfits
    import astropy.wcs as pywcs
    
    fine_visits, fine_fit = np.load('xxx_fine.npy', allow_pickle=True)
    
    N = len(fine_visits)

    trans = np.reshape(fine_fit.x, (N,-1))#/10.
    sh = trans.shape
    if sh[1] == 2:
        pscl = np.array([10.,10.])
        trans = np.hstack([trans/pscl, np.zeros((N,1)), np.ones((N,1))])
    elif sh[1] == 3:
        pscl = np.array([10.,10.,100])
        trans = np.hstack([trans/pscl, np.ones((N,1))])
    elif sh[1] == 4:
        pscl = np.array([10.,10.,100,100])
        trans = trans/pscl
    
    #  Update WCS
    # Update direct WCS
    for ix, direct in enumerate(fine_visits):
        root_i = direct['product'][:-7]
        print('\n\n\n\n##### ', ix, root_i, '#######\n\n\n')
        
        # Sync products
        os.system('aws s3 sync s3://grizli/Pipeline/{0}/Prep/ ./{0}/Prep/ --exclude "*" --include "{0}*drc_sci.fits.gz" --include "*_flc.*fits"'.format(root_i))

        #os.system('aws s3 ls s3://grizli/Pipeline/{0}/Prep/'.format(root_i))

        # Just visits and fail
        os.system('aws s3 sync s3://grizli/Pipeline/{0}/Prep/ ./{0}/Prep/ --exclude "*" --include "*visits.npy" --include "*fail*"'.format(root_i))
        
        #direct = visits[ix]
        out_shift, out_rot = trans[ix,:2], trans[ix,2]
        out_scale = trans[ix,3]
        
        xyscale = trans[ix,:4]
        wcs_ref_file = str('{0}/Prep/{0}-f850lp_drc_sci.fits.gz'.format(root_i))
        wcs_ref = pywcs.WCS(pyfits.open(wcs_ref_file)[0].header, 
                            relax=True)
        
        
        for file in glob.glob('{0}/Prep/*flc.fits'.format(root_i)):
            prep.update_wcs_fits_log(file, wcs_ref,
                                xyscale=xyscale,
                                initialize=False,
                                replace=('.fits', '.wcslog.fits'), 
                                wcsname='FINE')
        
            updatehdr.updatewcs_with_shift(file, 
                                  wcs_ref_file,
                                  xsh=out_shift[0], ysh=out_shift[1],
                                  rot=out_rot, scale=out_scale,
                                  wcsname='FINE', force=True,
                                  reusename=True, verbose=True,
                                  sciext='SCI')

            ### Bug in astrodrizzle? Dies if the FLT files don't have MJD-OBS
            ### keywords
            im = pyfits.open(file, mode='update')
            im[0].header['MJD-OBS'] = im[0].header['EXPSTART']
            im.flush()
    
    ## Fields with catalogs
    visit_files = glob.glob('j*/Prep/*f814w*visits.npy')
    visit_files.sort()
    
    ## Mosaic now that they're aligned
    
    # Shift parameters
    if False:
        os.system('aws s3 sync --exclude "*" --include "{0}*/Prep/*expflag*" --include "{0}*/Prep/*fail*" --include "{0}*/Prep/*cat.fits" --include "{0}*/Prep/*visits.npy" --include "{0}*/Prep/*[._]wcs*" --include "{0}*/Prep/*shifts.*" --include "{0}*/Prep/*expflag*" s3://grizli-v1/Pipeline/ .'.format(root))
        
        # No catalogs 
        os.system('aws s3 sync --exclude "*" --include "{0}*/Prep/*expflag*" --include "{0}*/Prep/*fail*" --include "{0}*/Prep/*visits.npy" --include "{0}*/Prep/*[._]wcs*" --include "{0}*/Prep/*shifts.*" s3://grizli-v1/Pipeline/ .'.format(root))
        
        # ACS
        os.system('aws s3 sync --exclude "*" --include "{0}*acswfc*/Prep/*expflag*" --include "{0}*acswfc*/Prep/*fail*" --include "{0}*acswfc*/Prep/*visits.npy" --include "{0}*acswfc*/Prep/*[._]wcs*" --include "{0}*acswfc*/Prep/*shifts.*" s3://grizli-v1/Pipeline/ .'.format(root))
    
    
    # Shifts
    os.system('echo "# root visit exp xsh ysh rot scl n xrms yrms" > shifts_results.log')
    os.system('grep "_fl"  j*/Prep/*shifts.log | grep -v "\#" | sed "s/:/ /" | sed "s/.Prep./ /g" | sed "s/_shifts.log//" >> shifts_results.log')
    sh = utils.read_catalog('shifts_results.log')
    badsh = (sh['xrms'] > 1) | (sh['yrms'] > 1) #| (sh['n'] < 10)
    shift_skip_visits = list(np.unique(sh['visit'][badsh]))
    
    # Flag
    os.system('echo "group,x,expm,flag" > expflag_results.csv')
    os.system('grep "fits" j*/Prep/*expflag*txt | sed "s/:/,/g" | sed "s/.txt/ /" | sed "s/.Prep//" | sed "s/\//,/g" | sed "s/_raw/_flt/" >> expflag_results.csv')
    expf = utils.read_catalog('expflag_results.csv')
    bad_expf = expf['flag'] != 'NORMAL'
    exp_visits = np.unique(expf['group'][bad_expf])
    exp_visits = []
    
    # wcs
    os.system('echo "# dir visit x xs ys rot scl rms n" > shift_wcs.log')
    os.system('grep " 0 " j*/Prep/*wcs.log | sed "s/.Prep./ /" | sed "s/_wcs.log://" >> shift_wcs.log')
    shifts = utils.read_catalog('shift_wcs.log')
    skip = (shifts['n'] < 7)
    skip |= (shifts['xs'] == 0.) | (np.abs(shifts['rot']) > 0.2)
    skip_visits = ['{0}_{1}'.format('_'.join(d.split('_')[:2]), v) for d, v in zip(shifts['dir'][skip], shifts['visit'][skip])]
        
    skip_visits += list(exp_visits)
    skip_visits += shift_skip_visits
    skip_visits = list(np.unique(skip_visits))
    
    skip_keys = ['uds-12-bhm']
    skip_keys += ['94s-245.0'] # bad shifts.txt
    skip_keys += ['sn2002zx-'] # bad scale
    
    visit_files = glob.glob('j*/Prep/*visits.npy')
    #visit_files = glob.glob('*shizels*/Prep/*visits.npy')
    
    all_visits = []
    all_groups = []
    all_info = []
    
    for file in visit_files:
        root_i = os.path.basename(file).split('_visits')[0]
        assoc = '_'.join(os.path.basename(file).split('_')[:2])
        
        visits, groups, info = np.load(file)
        info['keep'] = False
        
        ic = 0
        
        for v in visits:
            failed = glob.glob('{0}/Prep/{1}*fail*'.format(root_i, v['product']))
            if len(failed) > 0:
                print(failed)
                continue
            
            for k in skip_keys:
                if k in v['product']:
                    print('SKIP', v['product'])
                    continue
                    
            v['product'] = assoc+'_'+v['product']
            if v['product'] in skip_visits:
                print('SKIP ', v['product'])
                continue
            
            v['awspath'] = ['grizli-v1/Pipeline/{0}/Prep'.format(root_i)]*len(v['files'])
            all_visits.append(v)
            for f in v['files']:
                info['keep'][info['FILE'] == f] = True
        
        all_info.append(info[info['keep'] == True])
        
        for g in groups:
            failed = glob.glob('{0}/Prep/{1}*fail*'.format(root_i, g['direct']['product']))
            if len(failed) > 0:
                print(failed)
                continue
        
            # for ext in ['direct', 'grism']:
            #     g[ext]['product'] = assoc+'_'+g[ext]['product']
                
            all_groups.append(g)
                
    import astropy.table
    all_info = astropy.table.vstack(all_info)
    
    if root == 'j123656p6215':    
        out_root = 'gdn-'+root
        
    elif root == 'j021732m0512':
        out_root = 'uds-'+root
    elif root == 'j141956p5255':
        out_root = 'egs-'+root
    elif root == 'j033236m2748':
        out_root = 'gds-'+root
    else:
        out_root = 'xxx-'+root

    np.save('{0}_visits.npy'.format(out_root), [all_visits, all_groups, all_info])
    
    if True:
        os.system('aws s3 cp {0}_visits.npy s3://grizli-v1/Mosaics/ --acl public-read'.format(out_root))
        
    all_visits, all_groups, all_info = np.load('{0}_visits.npy'.format(out_root))
    
    ########
    # Remake catalogs and drizzled images
    from drizzlepac.astrodrizzle import AstroDrizzle
    import numpy as np
    import glob
    import os
    from grizli import prep
    
    all_visits, all_groups, all_info = np.load('{0}_visits.npy'.format(out_root))
    
    filters = ['f814w','f140w','f160w','f105w','f098m','f850lp']
    #filters = ['f850lp']
    count, products = 0, []
    for ii, direct in enumerate(all_visits[::-1]):
        filt = direct['product'].split('-')[-1]
        
        prod_files = glob.glob(direct['product']+'*sci.fits')
        
        if (filt not in filters):
            continue
        else:
            products.append(direct['product'])
            if (len(prod_files) > 0):
                print('Skip ', direct['product'])
                continue
            else:
                count+=1
                print('\n\n\n\n===========\n\n',ii, count, direct['product'], '\n\n\n========\n')
                #break
               
        isACS = '_flc' in direct['files'][0]
        if isACS:
            bits = 64+32+256
            driz_cr_snr = '3.5 3.0'
            driz_cr_scale = '1.2 0.7'
        else:
            bits = 576+256
            driz_cr_snr = '8.0 5.0'
            driz_cr_scale = '2.5 0.7'
        
        for aws, file in zip(direct['awspath'], direct['files']):
            if not os.path.exists(file):
                os.system('aws s3 cp s3://{0}/{1} ./'.format(aws, file))
                
        cr_corr = False
        AstroDrizzle(direct['files'], output=direct['product'], clean=True,
                     context=False, preserve=False, skysub=True,
                     driz_separate=cr_corr, driz_sep_wcs=cr_corr, 
                     median=cr_corr, 
                     blot=cr_corr, driz_cr=cr_corr, 
                     driz_cr_corr=cr_corr, 
                     driz_cr_snr=driz_cr_snr, driz_cr_scale=driz_cr_scale, 
                     driz_combine=True, final_bits=bits, coeffs=True, 
                     resetbits=4096*cr_corr, build=False, 
                     final_kernel='point', 
                     final_wht_type='IVM')        
                
        # Remake catalog
        cat = prep.make_SEP_catalog(root=direct['product'], threshold=5)
        os.system('rm {0}*_[wbs]?[gt].fits'.format(direct['product']))
        
        # Remove mosaic
        os.system('rm {0}*_sci.fits'.format(direct['product']))

        # Remove FLC
        if isACS:
            for f in direct['files']:
                os.remove(f)
        
        
    
    # Make mosaics
    kwargs = auto_script.get_yml_parameters()
    mos_args = {'mosaic_args':kwargs['mosaic_args'],
                'fix_stars':kwargs['visit_prep_args']['fix_stars'], 
                'mask_spikes':kwargs['mask_spikes'], 'skip_single_optical_visits':kwargs['preprocess_args']['skip_single_optical_visits']}    
    
    #mos_args['mosaic_args']['ir_filters'] = ['F140W']
    #mos_args['mosaic_args']['optical_filters'] = ['F814W']
    mos_args['mosaic_args']['fill_mosaics'] = False
    mos_args['mosaic_args']['half_optical_pixscale'] = True

    mos_args['mosaic_args']['kernel'] = 'square'
    mos_args['mosaic_args']['pixfrac'] = 0.33
    mos_args['mosaic_args']['wcs_params']['pixel_scale'] = 0.1
    
    if root == 'j123656p6215':
        mos_args['mosaic_args']['kernel'] = 'point'
        mos_args['mosaic_args']['pixfrac'] = 0.33
        mos_args['mosaic_args']['wcs_params']['pixel_scale'] = 0.1
        
    #mos_args['mosaic_args']['wcs_params']['filters'] = ['F140W']
    
    os.system('ln -s j*/Prep/*_fl?.fits .')
    
    auto_script.make_combined_mosaics(out_root, **mos_args) 
    
    # Grism models for subsets of visits....
    