#!/usr/bin/env python
def auto_run(root='j023507-040202'):
    
    import os
    import glob
    import numpy as np
    
    import astropy.io.fits as pyfits
    import astropy.wcs as pywcs
    
    from drizzlepac import updatehdr
    from stwcs import updatewcs
    
    from grizli import utils, prep
    from grizli.pipeline import auto_script
    
    visit_file = '{0}_visits.npy'.format(root)
    visits, all_groups, info = np.load(visit_file)
    
    # Something wrong with some files with bad shifts, reset wcs
    for visit in visits:
        for file in visit['files']:
            utils.fetch_hst_calibs(file, calib_types=['IDCTAB','NPOLFILE','IMPHTTAB'])
            updatewcs.updatewcs(file, verbose=True, use_db=False)
        
        # Apply shifts
        shift_log = '{0}_shifts.log'.format(visit['product'])
        if os.path.exists(shift_log):
            sh = utils.read_catalog(shift_log)
            flt0 = pyfits.open(sh['flt'][0])
            wcs_ref = pywcs.WCS(flt0['SCI',1].header, fobj=flt0, relax=True)
            shift_dict = {}
            for i in range(len(sh)):
                shift_dict[sh['flt'][i]] = [sh['xshift'][i], sh['yshift'][i]]
            
            prep.apply_tweak_shifts(wcs_ref, shift_dict, grism_matches={},
                                    verbose=False)
                                  
    # Redrizzle mosaics
    prep.drizzle_overlaps(visits, check_overlaps=False, skysub=False,
                          static=False, pixfrac=0.8, scale=None, 
                          final_wcs=False, fetch_flats=False, final_rot=None)
    
    ####### Alignment
    os.system('rm *wcs.*')
    
    # Radec
    master_radec = '{0}/../../{1}_master.radec'.format(os.getcwd(), root)
    hsc = '{0}/../../{1}'.format(os.getcwd(), 'hsc-udeep-i25_corr_cosmos.radec')
    
    if not os.path.exists(master_radec):
        if root.startswith('cos-') & os.path.exists(hsc):
            master_radec = hsc
        else:
            master_radec = None
    
    parent_radec = '{0}/../../{1}_parent.radec'.format(os.getcwd(), root)
    if not os.path.exists(parent_radec):
        parent_radec = None
    
    if master_radec is not None:
        radec = master_radec
    elif parent_radec is not None:
        radec = parent_radec
    else:
        radec = hsc
     
    print('master RADEC file: ', radec)
       
    thresh = 2.5
    for visit in visits:
        # Remake catalogs
        cat = prep.make_SEP_catalog(root=visit['product'], threshold=thresh)
        
        # Redo alignment
        try:
            result = prep.align_drizzled_image(root=visit['product'], radec=radec, mag_limits=[19,23], simple=False, max_err_percentile=80, clip=120, outlier_threshold=5)
        except:
            continue
            
        orig_wcs, drz_wcs, out_shift, out_rot, out_scale = result
        
        # Propagate shifts
        for file in visit['files']:
            updatehdr.updatewcs_with_shift(file, 
                                str('{0}_wcs.fits'.format(visit['product'])),
                                      xsh=out_shift[0], ysh=out_shift[1],
                                      rot=out_rot, scale=out_scale,
                                      wcsname='HSC', force=True,
                                      reusename=True, verbose=True,
                                      sciext='SCI')
        
            ### Bug in astrodrizzle? Dies if the FLT files don't have MJD-OBS
            ### keywords
            im = pyfits.open(file, mode='update')
            im[0].header['MJD-OBS'] = im[0].header['EXPSTART']
            im.flush()
    
    # Redrizzle mosaics again including new shifts
    prep.drizzle_overlaps(visits, check_overlaps=False, skysub=False,
                          static=False, pixfrac=0.8, scale=None, 
                          final_wcs=False, fetch_flats=False, final_rot=None)
    # Remake catalogs
    thresh = 2.5
    for visit in visits:
        # Remake catalogs
        cat = prep.make_SEP_catalog(root=visit['product'], threshold=thresh)
    
    # Update visits file
    v = auto_script.get_visit_exposure_footprints(visit_file=visit_file,   
                             check_paths=['./', '../RAW'], simplify=1.e-6)
    
    if False:
        # Mosaic
        auto_script.drizzle_overlaps(root, filters=['F160W'], 
                                     min_nexp=1, pixfrac=0.8, scale=0.1,
                                     make_combined=False,
                                     ref_image=None, 
                                     static=False) 
         
    # import os
    # import glob
    # 
    # import matplotlib.pyplot as plt
    # 
    # from grizli import utils, prep
    # from grizli.pipeline import auto_script, photoz
    # utils.set_warnings()
    # 
    # tab = utils.GTable.gread('{0}_footprint.fits'.format(root))
    # 
    # HOME_PATH = os.getcwd()
    # 
    # auto_script.VALID_FILTERS = ['F098M', 'F105W', 'F110W', 'F125W', 'F127M', 'F139M', 'F140W', 'F153M', 'F160W', 'F410M', 'F435W', 'F438W', 'F439W', 'F450W', 'F467M', 'F475W', 'F475X', 'F547M', 'F550M', 'F555W', 'F569W', 'F600LP', 'F606W', 'F621M', 'F622W', 'F625W', 'F675W', 'F689M', 'F702W', 'F763M', 'F775W', 'F791W', 'F814W', 'F845M', 'F850LP', 'F350LP']
    # 
    # IS_PARALLEL = utils.column_string_operation(tab['proposal_pi'], 'alkan', method='count', logical='or').sum() > 0
    #     
    # IS_DASH = list(np.unique(np.cast[int](tab['proposal_id']))) == [14114]
    # 
    # master_radec = '{0}/{1}_master.radec'.format(os.getcwd(), root)
    # if not os.path.exists(master_radec):
    #     if root.startswith('cos-') & os.path.exists('hsc-udeep-i25_corr_cosmos.radec'):
    #         master_radec = '{0}/{1}'.format(os.getcwd(), 'hsc-udeep-i25_corr_cosmos.radec')
    #     else:
    #         master_radec = None
    # 
    # parent_radec = '{0}/{1}_parent.radec'.format(os.getcwd(), root)
    # if not os.path.exists(parent_radec):
    #     parent_radec = None
    # 
    #     
    # BKG_PARAMS = {'bw': 1024, 'bh': 1024, 'fw': 3, 'fh': 3}
    # 
    # auto_script.go(root=root, maglim=[19, 23], HOME_PATH=HOME_PATH, inspect_ramps=False, manual_alignment=False, is_parallel_field=IS_PARALLEL, reprocess_parallel=True, only_preprocess=True, run_extractions=True, run_fit=False, s3_sync='cp', fine_radec=None, combine_all_filters=False, gaia_by_date=True, align_simple=False, align_clip=100, master_radec=master_radec, parent_radec=parent_radec, is_dash=IS_DASH, run_parse_visits=True, reference_wcs_filters=['F160W','F140W','F125W','F105W','F110W','F098M','F814W','F850LP', 'F606W','F435W'], make_phot=False, make_mosaics=False, align_rms_limit=4, align_min_overlap=2, imaging_bkg_params=BKG_PARAMS)
    #         
    
    # Make footprints
    
    # files = glob.glob('*dr?_wht.fits')
    # if len(files) > 0:
    #     print('Make footprint')
    #     fp = open('{0}.fp.reg'.format(root),'w')
    #     fp.write('fk5\n')
    # 
    #     for weight_image in files:
    #         root_i = '_dr'.join(weight_image.split('_dr')[:-1])
    #         reg = prep.drizzle_footprint(weight_image, shrink=10, ext=0,
    #                                      outfile=None, label=root_i) 
    #         fp.write(reg+'\n')
    #     
    #     fp.close()
    
        
if __name__ == "__main__":
    import sys
    import time
    
    import numpy as np
    from grizli import utils
    from grizli.pipeline import auto_script
    utils.set_warnings()
            
    root = sys.argv[1]
    #auto_run(root=root) 
    
    try:
        auto_run(root=root) 
    except:
        fp = open('{0}.failed'.format(root),'w')
        fp.write(time.ctime()+'\n')
        fp.close()
        