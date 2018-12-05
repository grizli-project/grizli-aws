#!/usr/bin/env python
def auto_run(root='j023507-040202'):
    
    import os
    import glob
    
    import matplotlib.pyplot as plt
    
    from grizli import utils, prep
    from grizli.pipeline import auto_script, photoz
    utils.set_warnings()
    
    tab = utils.GTable.gread('{0}_footprint.fits'.format(root))
    
    HOME_PATH = os.getcwd()
    
    auto_script.VALID_FILTERS = ['F098M', 'F105W', 'F110W', 'F125W', 'F127M', 'F139M', 'F140W', 'F153M', 'F160W', 'F410M', 'F435W', 'F438W', 'F439W', 'F450W', 'F467M', 'F475W', 'F475X', 'F547M', 'F550M', 'F555W', 'F569W', 'F600LP', 'F606W', 'F621M', 'F622W', 'F625W', 'F675W', 'F689M', 'F702W', 'F763M', 'F775W', 'F791W', 'F814W', 'F845M', 'F850LP', 'F350LP']
    
    IS_PARALLEL = utils.column_string_operation(tab['proposal_pi'], 'alkan', method='count', logical='or').sum() > 0
        
    IS_DASH = list(np.unique(np.cast[int](tab['proposal_id']))) == [14114]
    
    master_radec = '{0}/{1}_master.radec'.format(os.getcwd(), root)
    if not os.path.exists(master_radec):
        if root.startswith('cos-') & os.path.exists('hsc-udeep-i25_corr_cosmos.radec'):
            master_radec = '{0}/{1}'.format(os.getcwd(), 'hsc-udeep-i25_corr_cosmos.radec')
        else:
            master_radec = None
    
    parent_radec = '{0}/{1}_parent.radec'.format(os.getcwd(), root)
    if not os.path.exists(parent_radec):
        parent_radec = None
       
    BKG_PARAMS = {'bw': 1024, 'bh': 1024, 'fw': 3, 'fh': 3}
    
    catalogs = ['PS1','DES','NSC', 'SDSS','GAIA','WISE']
    
    # MACS0717
    if 'macs0717' in root:
        catalogs = ['GAIA', 'PS1']
    
    run_fine_alignment=False
        
    auto_script.go(root=root, maglim=[19, 23], HOME_PATH=HOME_PATH, 
                   inspect_ramps=False, manual_alignment=False, 
                   is_parallel_field=IS_PARALLEL, reprocess_parallel=True, 
                   only_preprocess=True, run_extractions=True, run_fit=False, 
                   s3_sync='cp', catalogs=catalogs,
                   fine_radec=None, run_fine_alignment=run_fine_alignment, 
                   combine_all_filters=False, 
                   gaia_by_date=True, align_simple=False, align_clip=120, 
                   master_radec=master_radec, parent_radec=parent_radec, 
                   is_dash=IS_DASH, run_parse_visits=True, 
                   reference_wcs_filters=['F160W','F140W','F125W','F105W','F110W','F098M','F814W','F850LP', 'F606W','F435W'], 
                   make_phot=False, make_mosaics=False, align_rms_limit=4, 
                   align_min_overlap=2, imaging_bkg_params=BKG_PARAMS)
        
    
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
    import numpy as np
    from grizli import utils
    from grizli.pipeline import auto_script
    utils.set_warnings()
            
    root = sys.argv[1]
    auto_run(root=root) 