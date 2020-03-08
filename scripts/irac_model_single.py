#!/usr/bin/env python

import sys
import os
import time

import golfir.model

root = sys.argv[1]

start_file = '/GrizliImaging/{0}.start.txt'.format(root)
if os.path.exists(start_file):
    print('start file {0} exists'.format(start_file))
    return True
    
ds9 = None

kwargs = {'ds9': ds9, 
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
          'window': None            # PSF-match windowing
          }
          

kwargs['patch_arcmin'] = -1
golfir.model.run_all_patches(root, PATH='/GrizliImaging/', use_patches=True, sync_results=True, **kwargs)

fp = open(f'/tmp/{root}.finished.txt','w')
fp.write(time.ctime())
fp.close()


