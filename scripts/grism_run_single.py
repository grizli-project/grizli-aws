#!/usr/bin/env python

def auto_run(root='j023507-040202'):
    import os
    from grizli import utils
    utils.set_warnings()
    from grizli.pipeline import auto_script
    
    params = {
    # Preprocessing flags         
    's3_sync': True,               # Fetch data from AWS
    'inspect_ramps': False,        # Visual inspection on ramps
    'reprocess_parallel': True,    # Reprocess ramps in parallel
    'is_dash': False,              # Flag if visits are in DASH mode
    'filters': auto_script.VALID_FILTERS,   #  Only use these filters
    'run_parse_visits': True,      # Parse visits to _visits.npy
    'is_parallel_field': False,    # If parallels, change visit parsing
    'skip_single_optical_visits': False, # Skip preprocess if single exp
    'fix_stars': True,             # Fill DQ flagged centers of stars
    'imaging_bkg_params': {'bh': 128, 'bw': 128, 'fh': 3, 'fw': 3,
                       'pixel_scale': 0.06}, # Imaging background 

    # Alignment
    'manual_alignment': False,    # Interactive alignment
    'align_clip': 120,            # Parameter for initial alignment pairs
    'align_mag_limits': [14, 24], # Magnitude range for alignment sources
    'align_min_overlap': 0.2,     # Overlap fraction to use HST catalogs
    'align_outlier_threshold': 4, # Parameter for triangle matches
    'align_rms_limit': 2,    # Force shifts = 0 if rms > rms limit
    'align_simple': False,   # ??
    'catalogs': ['PS1','DES','NSC','SDSS','GAIA','WISE'], # Ref cat order
    'gaia_by_date': True,   # Get GAIA catalogs by visit epoch
    'master_radec': None,    # Force reference
    'parent_radec': None,    # Reference to use if no HST overlaps

    'run_fine_alignment': False, # Run Fine alignment script
    'fine_radec': None,         # Fine alingment reference,  else GAIA

    # Image Mosaics                   
    'make_mosaics': True,             # Make full-field mosaics
    'fill_mosaics': 'grism',          # Fill empty parts of IR mosaics
    'combine_all_filters': False,      # Combine ACS & IR
    'mosaic_pixel_scale': None,       # Pixel scale of mosaics (0.06)
    'half_optical_pixscale': False,   # Make smaller ACS/UVIS pixel scale
    'mosaic_pixfrac': 0.75,           # Pixfrac of mosaic images
    'reference_wcs_filters': auto_script.VALID_FILTERS, # These define WCS if they exist
    'make_phot': True,                # Make photometric catalog

    # Spectral Extractions
    'only_preprocess': False, # Make grism models
    'run_extractions': False, # Don't extract grism spectra
    'maglim': [17, 26],       # Magnitude range of grism spectra 
    'run_fit': False          # Fit grism spectra
    }
      
    tab = utils.GTable.gread('{0}_footprint.fits'.format(root))
    
    IS_PARALLEL = utils.column_string_operation(tab['proposal_pi'], ['Malkan', 'Trenti'], method='count', logical='or').sum() > 0
    
    master_radec = '{0}/{1}_master.radec'.format(os.getcwd(), root)
    if not os.path.exists(master_radec):
        master_radec = None
    
    parent_radec = '{0}/{1}_parent.radec'.format(os.getcwd(), root)
    if not os.path.exists(parent_radec):
        parent_radec = None
    
    params['parent_radec'] = parent_radec
    params['master_radec'] = master_radec
    
    # Limited filters
    params['only_preprocess'] = False
    params['filters'] = auto_script.IR_W_FILTERS + auto_script.IR_GRISMS
    #params['filters'] += ['F814W', 'F850LP', 'F777W']
    
    params['is_parallel_field'] = IS_PARALLEL
    params['mosaic_pixel_scale'] = 0.06+0.02*params['is_parallel_field']
    params['mosaic_pixfrac'] = params['mosaic_pixel_scale']/0.12
    
    auto_script.go(root=root, HOME_PATH=os.getcwd(), **params)
    
    #auto_script.go(root=root, maglim=[19, 23], HOME_PATH=HOME_PATH, inspect_ramps=False, manual_alignment=False, is_parallel_field=IS_PARALLEL, reprocess_parallel=True, only_preprocess=False, run_extractions=False, run_fit=False, s3_sync='cp', fine_radec=None, combine_all_filters=False, gaia_by_date=True, align_simple=False, align_clip=100, master_radec=master_radec, parent_radec=parent_radec, is_dash=False, run_parse_visits=True, reference_wcs_filters=['F160W','F140W','F125W','F105W','F110W','F098M','F814W','F850LP', 'F606W','F435W'])
    
    os.chdir('../Prep/')
    auto_script.make_report(root)
    
    # plt.ioff()
    # #fig = auto_script.field_rgb(root=root, HOME_PATH=HOME_PATH, xsize=18)
    # plt.close(fig)
    
    # Photo-z
    # try:
    #     out = photoz.eazy_photoz(root, object_only=False, force=True,
    #                           aper_ix=1, sys_err=0.05, apply_prior=False,
    #                           beta_prior=True, 
    #                           external_limits=3, external_sys_err=0.3)
    # except:
    #     pass
    
if __name__ == "__main__":
    import sys
    import numpy as np
    from grizli import utils
    from grizli.pipeline import auto_script
    utils.set_warnings()
    
    if len(sys.argv) < 3:
        print('Usage: aws.py {field} {init/run/summary} [min_mag,max_mag]')
        exit 
        
    root = sys.argv[1]
    auto_run(root=root) 