#!/usr/bin/env python

import sys
import os
import time
import json

import golfir.model
import golfir.utils

def run(root, argv=[], ds9=None):
    
    run_dir = '/GrizliImaging/{0}'.format(root)
    if os.path.exists(run_dir):
        print('directory {0} exists'.format(run_dir))
        return True

    if os.path.exists('/tmp/{0}.finished.txt'.format(root)):
        print('/tmp/{0}.finished.txt'.format(root))
        return True
    
    #ds9 = None
    
    if ds9 is 'connect':
        import grizli.ds9
        ds9 = grizli.ds9.DS9()
        
    defaults = {'ds9': ds9, 
              'patch_arcmin': 1.0,      # Size of patch to fit
              'patch_overlap': 0.2,     # Overlap of automatic patches
              'mag_limit': [24, 27],    # Two-pass modeling.  Fit sources brighter than mag_limit in HST catalog
              'run_alignment': True,    # Run fine alignment between IRAC-HST, in between two steps of `mag_limit`
              'galfit_flux_limit': 70,  # Brightness limit (uJy) of objects to fit with GALFIT
              'refine_brightest': True, # Refine masked bright objects with galfit
              'any_limit': 16,          # Brightness limit below which to mask *any* sources
              'point_limit': 16,        # Brightness limit below which to mask point-like sources
              'bright_sn': 30,          # S/N threshold for masked pixels of bright object
              'bkg_kwargs': {'order_npix': 32},          # Arguments to the local background routine
              'channels': ['ch1', 'ch2'],  # Channels to try
              'psf_only': False,        
              'use_saved_components': False, # Use models from a "components" file if found
              'window': None,             # PSF-match windowing
              'fetch': True, 
              'PATH': '/GrizliImaging/', 
              'use_patches': True, 
              'sync_results': True}
          

    defaults['patch_arcmin'] = -1
    args, kwargs = golfir.utils.argv_to_dict(argv, defaults=defaults)
    
    golfir.model.run_all_patches(root, **kwargs)
        
    os.chdir('/GrizliImaging/')
    os.system(f'rm -rf /GrizliImaging/{root}')
    
    fp = open(f'/tmp/{root}.finished.txt','w')
    fp.write(time.ctime())
    fp.close()
    
    return True

if __name__ == '__main__':
    root = sys.argv[1]
             
    run(root, argv=sys.argv[1:])
