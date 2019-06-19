
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import glob
from grizli import prep, utils

def test_egs():
    
    size=(80,20)
    theta = 49.699
    dy = -1./60
    
    fields = {}
    fields['egs-rot'] = {'field':'egs-rot', 'ra':214.8288, 'dec':52.8234+dy, 'size':size, 'tile_size':6, 'overlap':0.3, 'pixscale':0.05, 'theta':theta}
    
    key = 'egs-rot'
    tiles = define_tiles(**fields[key])
    os.system('xpaset -p ds9 regions delete all; ds9_reg egs_50mas_tile_wcs_magenta.reg ; ds9_reg egs-rot_50mas_tile_wcs.reg')
    
IR_FILTERS = ['f098m', 'f110w', 'f140w', 'f105w', 'f160w','f125w' ]
OPTICAL_FILTERS = ['f606w', 'f775w', 'f814w', 'f850lp']
    
def make_candels_tiles(key='gdn', filts=OPTICAL_FILTERS, make_tile_catalogs=False):
    
    import numpy as np
    import matplotlib.pyplot as plt
    import astropy.io.fits as pyfits
    import glob
    from grizli import prep, utils
    
    try:
        test = key
    except:
        key = 'gdn'
        
    # Dimensions
    if 0:
        # Compute size from visit footprints
        xlim = [1e30, -1e30]
        ylim = [1e30, -1e30]
        for visit in all_visits:
            xp, yp = visit['footprint'].buffer(1./60).boundary.xy
            xlim[0] = np.minimum(xlim[0], np.min(xp))
            xlim[1] = np.maximum(xlim[1], np.max(xp))
        
            ylim[0] = np.minimum(ylim[0], np.min(yp))
            ylim[1] = np.maximum(ylim[1], np.max(yp))
    
        ra, dec = np.mean(xlim), np.mean(ylim)
        cosd = np.cos(dec/180*np.pi)
        dx = np.diff(xlim)[0]*cosd*60
        dy = np.diff(ylim)[0]*60
        print('\'ra\':{0:.4f}, \'dec\':{1:.4f}, \'size\':({2:.1f}, {3:.1f})'.format(ra, dec, dx, dy))
        
    fields = {}
    fields['gdn'] = {'field':'gdn', 'ra':189.236, 'dec':62.257, 'size':(24,24), 'tile_size':6, 'overlap':0.3, 'pixscale':0.05}
    fields['uds'] = {'field':'uds', 'ra':34.361, 'dec':-5.163, 'size':(36,24), 'tile_size':6, 'overlap':0.3, 'pixscale':0.05}
    fields['egs'] = {'field':'egs', 'ra':214.8288, 'dec':52.8234, 'size':(54, 60), 'tile_size':6, 'overlap':0.3, 'pixscale':0.05}
    fields['gds'] = {'field':'gds', 'ra':53.1180-0.5/60, 'dec':-27.8216-0.5/60, 'size':(32, 32), 'tile_size':6, 'overlap':0.3, 'pixscale':0.05}
    
    # Rotated EGS along strip
    size=(80,20)
    theta = 49.699
    dy = -1./60
    fields['egs'] = {'field':'egs', 'ra':214.8288, 'dec':52.8234+dy, 'size':size, 'tile_size':6, 'overlap':0.3, 'pixscale':0.05, 'theta':theta}
    
    from grizli_aws.field_tiles import define_tiles, drizzle_tiles
     
    tiles = define_tiles(**fields[key])
    
    if False:
        # test single tile
        t = '01.01'
        one_tile = {}
        one_tile[t] = tiles[t]
        drizzle_tiles(all_visits, one_tile, filts=filts, prefix=key, pixfrac=0.33, output_bucket='s3://grizli-v1/Mosaics/')
        
    #all_visits, all_groups, all_info = np.load('goodss-j033236m2748_visits.npy')
    
    #filts = ['f098m', 'f110w', 'f140w', 'f105w', 'f160w','f125w' ]
    
    #filts = ['f606w', 'f775w', 'f814w', 'f850lp']
    
    visit_file = glob.glob('{0}-j*visits.npy'.format(key))
    if len(visit_file) == 0:
        os.system('aws s3 sync s3://grizli-v1/Mosaics/ ./ --exclude "*" --include "{0}*npy"'.format(key))
        visit_file = glob.glob('{0}-j*visits.npy'.format(key))
    
    all_visits, all_groups, info = np.load(visit_file[0])
    
    drizzle_tiles(all_visits, tiles, filts=filts, prefix=key, pixfrac=0.33, output_bucket='s3://grizli-v1/Mosaics/')
    
    ##########
    # Combine final mosaics
    
    filt_i ='f814w'
    
    # Compute tiles needed for full mosaic extent
    init = True
    use_ref = 'all'
    
    if use_ref == 'wfc3ir':
        # wfc3/ir
        files = glob.glob(key+'*-[0-9][0-9].[0-9][0-9]*100mas*-*sci.fits.gz')
    elif use_ref == 'acswfc':
        # acs
        files = glob.glob(key+'*-[0-9][0-9].[0-9][0-9]*050mas*-*sci.fits.gz')
    else:
        files = glob.glob(key+'*-[0-9][0-9].[0-9][0-9]*-*sci.fits.gz')

    files.sort()

    tile_pos = np.array([int(file.split('-')[1].replace('.','')) for file in files])
    xp = tile_pos // 100
    yp = tile_pos % 100
    
    left, bot = xp.min(), yp.min()

    nx = (xp.max()-xp.min())+1
    ny = (yp.max()-yp.min())+1
    print('Define mosaic for {0}: {1} x {2}'.format(key, nx, ny))
    
    # WCS defined in lower left corner
    ll = '{0:02d}.{1:02d}'.format(left, bot)
    tiles = np.load('{0}_50mas_tile_wcs.npy'.format(key))[0]
    wcs = tiles[ll]

    ## Loop over filters
    #filts = ['f098m','f110w','f140w','f105w','f160w','f125w']
    for filt_i in filts:
    
        files = glob.glob('*-[0-9][0-9].[0-9][0-9]*-'+filt_i+'*sci.fits.gz')
        files.sort()
        if len(files) == 0:
            continue
            
        tile_pos = np.array([int(file.split('-')[1].replace('.','')) for file in files])
        xp = tile_pos // 100
        yp = tile_pos % 100
        
        ref = pyfits.open(files[0])[0]
        sh = ref.data.shape
        the_filter = utils.get_hst_filter(ref.header)
    
        is_fine = '50mas' in files[0]
        if is_fine:
            wh = wcs._header
        else:  
            wh = utils.to_header(tiles[ll])
            for k in ['NAXIS1','NAXIS2']:
                wh[k] //= 2
        
            wh['CRPIX1'] = (wh['CRPIX1']+0.5)/2.
            wh['CRPIX2'] = (wh['CRPIX2']+0.5)/2.
        
            for k in ['CD1_1','CD1_2', 'CD2_1','CD2_2']:
                if k in wh:
                    wh[k] *= 2
    
        for k in ['INSTRUME', 'DETECTOR', 'PHOTFNU', 'PHOTFLAM', 'PHOTBW','PHOTZPT', 'PHOTMODE', 'PHOTPLAM', 'FILTER', 'FILTER1', 'FILTER2']:
            if k in ref.header:
                wh[k] = ref.header[k]
    
        # Compute size of overlap in pixels (tiles generated with 0.3' overlap)
        coarse_pix = 0.1
        overlap_arcmin = 0.3    
        olap = int(overlap_arcmin*60/coarse_pix*(1+is_fine))
    
        for ext in ['sci','wht']:
            data = np.zeros((sh[0]*ny-olap*(ny-1), sh[1]*nx-olap*(nx-1)), dtype=np.float32)
            data_sh = data.shape
            for i, file in enumerate(files):
                x0 = (xp[i]-left)*(sh[1]-olap)
                y0 = (yp[i]-bot)*(sh[0]-olap)
                slx = slice(x0,x0+sh[1])
                sly = slice(y0,y0+sh[0])

                if (slx.stop > data_sh[1]) | (sly.stop > data_sh[0]):
                    print('skip', file.replace('_sci', '_'+ext))
                    continue
                else:
                    print(file.replace('_sci', '_'+ext))
                    
                sci_i = pyfits.open(file.replace('_sci', '_'+ext))
                data[sly, slx] = sci_i[0].data
        
            mos = pyfits.PrimaryHDU(data=data, header=wh)
            mos.writeto('{0}-{1:03d}mas-{2}_drz_{3}.fits'.format(key, 100//(1+is_fine), the_filter.lower(), ext), overwrite=True)
        
            del(data)
    
    # Gzip and sync
    os.system('gzip --force {0}-???mas*fits'.format(key))
    os.system('aws s3 sync --exclude "*" --include "{0}-???mas*gz" ./ s3://grizli-v1/Mosaics/ --acl public-read'.format(key))
    
    if make_tile_catalogs:
        # Make catalogs
        import os
        from grizli import prep
        import glob
        for filt_i in filts:
            files = glob.glob('*-[0-9][0-9].[0-9][0-9]*-'+filt_i+'*sci.fits.gz')
            files.sort()
            for i, file in enumerate(files):
                
                print('\n\n#######\n{0}\n#######\n'.format(froot))
                
                froot = file.split('_dr')[0]
                if os.path.exists('{0}.cat.fits'.format(froot)):
                    continue
                
                prep.make_SEP_catalog(froot, threshold=3,
                                      rescale_weight=False) 
                os.system('rm {0}_bkg.fits'.format(froot))
                #os.system('rm {0}_seg.fits'.format(froot))
                
def grism_model():
    
    import numpy as np
    from grizli import utils, multifit
    from grizli.pipeline import auto_script
    import astropy.io.fits as pyfits
    import os
    import glob
    
    from grizli.pipeline import auto_script
    kwargs = auto_script.get_yml_parameters()
    
    if key == 'uds':
        field_root = root = 'uds-grism-j021732m0512'
    elif key == 'egs':
        field_root = root = 'egs-grism-j141956p5255'
    elif key == 'gds':
        field_root = root = 'gds-grism-j033236m2748'
    elif key == 'gdn':
        field_root = root = 'gdn-grism-j123656p6215'

    # Optical copy
    for filt in ['f606w', 'f775w', 'f814w', 'f850lp']:
        sci_file = '{0}-050mas-{1}_drz_sci.fits.gz'.format(key, filt)
        if not os.path.exists(sci_file):
            continue
        
        print(sci_file)
         
        #Symlink to individual bands
        os.system('cp -f {0}-050mas-{2}_drz_sci.fits.gz {1}-{2}_drz_sci.fits.gz'.format(key, root, filt))
        os.system('cp -f {0}-050mas-{2}_drz_wht.fits.gz {1}-{2}_drz_wht.fits.gz'.format(key, root, filt))
    
    
    # IR-combined 
    num = None
    for filt in ['f140w', 'f105w', 'f125w', 'f160w']:
        sci_file = '{0}-100mas-{1}_drz_sci.fits.gz'.format(key, filt)
        if not os.path.exists(sci_file):
            continue
        
        print(sci_file)
         
        #Symlink to individual bands
        os.system('cp -f {0}-100mas-{2}_drz_sci.fits.gz {1}-{2}_drz_sci.fits.gz'.format(key, root, filt))
        os.system('cp -f {0}-100mas-{2}_drz_wht.fits.gz {1}-{2}_drz_wht.fits.gz'.format(key, root, filt))
 
        im_i = pyfits.open(sci_file)
        wht_i = pyfits.open(sci_file.replace('_sci','_wht'))
        photflam = im_i[0].header['PHOTFLAM']
        if num is None:
            ref_photflam = photflam*1
            ref_h = im_i[0].header
            num = im_i[0].data*0
            den = num*0
        
        scl = photflam/ref_photflam
        den_i = wht_i[0].data/scl**2
        num += im_i[0].data*scl*den_i
        den += den_i
    
    sci = num/den
    wht = den
    mask = (~np.isfinite(sci)) | (den == 0)
    sci[mask] = 0
    wht[mask] = 0
    
    pyfits.PrimaryHDU(data=sci, header=ref_h).writeto('{0}-{1}_drz_sci.fits'.format(root, 'ir'), overwrite=True, output_verify='fix')
    pyfits.PrimaryHDU(data=wht, header=ref_h).writeto('{0}-{1}_drz_wht.fits'.format(root, 'ir'), overwrite=True, output_verify='fix')
    
    # Fill filters with combination
    for filt in ['f140w', 'f105w', 'f125w', 'f160w']:
        sci_file = '{0}-{1}_drz_sci.fits.gz'.format(root, filt)
        if not os.path.exists(sci_file):
            continue
                 
        im_i = pyfits.open(sci_file, mode='update')
        wht_i = pyfits.open(sci_file.replace('_sci','_wht'))
        photflam = im_i[0].header['PHOTFLAM']
        mask = (wht_i[0].data <= 0) & (wht > 0)
        print(sci_file, mask.sum())
        
        im_i[0].data[mask] = sci[mask]*ref_photflam/photflam
        im_i.flush()
        
    del(sci); del(wht); del(mask); del(im_i); del(wht_i); del(den_i)
    
    # Gzip and sync
    os.system('gzip {0}*fits'.format(root))
    os.system('aws s3 sync --exclude "*" --include "{0}*" ./ s3://grizli-v1/Mosaics/ --acl public-read'.format(root))
     
    ########## Grism prep
    
    #out_root = 'uds-j021732m0512'
    out_root = root.replace('-grism','')
    
    all_visits, all_groups, all_info = np.load('{0}_visits.npy'.format(out_root))
    
    bucket = 'grizli-v1'
    base_path='Exposures'
    
    for v in all_visits:
        # Change aws path
        v['awspath'] = []
        for file in v['files']:
            file_root = os.path.basename(file)
            prog = file_root[:4]
            dataset = file_root[:9]
            
            aws = os.path.join(bucket, base_path, prog, dataset)
            v['awspath'].append(aws)
            
    groups = all_groups#[:28]
    visits = []
    for g in groups:
        visits.append(g['direct'])
        visits.append(g['grism'])

    all_info['keep'] = False
    for v in visits:
        for f in v['files']:
            ix = all_info['FILE'] == f
            all_info['keep'][ix] = True
            
    np.save(root+'_visits.npy', [visits, groups, all_info[all_info['keep']]])
    
    ####### Run model
    grisms = ['g141', 'g102', 'g800l']
    
    if not os.path.exists('{0}-ir_seg.fits'.format(root)):
        os.system('aws s3 sync s3://grizli-v1/GrismMosaics/ ./ --exclude "*" --include "{0}*"'.format(root))
        os.system('aws s3 sync s3://grizli-v1/Mosaics/ ./ --exclude "*" --include "{0}*"'.format(root))
        
        print('gunzip {0}*gz'.format(root))
        os.system('gunzip {0}*gz'.format(root))
    
    visits, groups, info = np.load(root+'_visits.npy')
            
    if not os.path.exists('{0}_phot.fits'.format(root)):
       
        kwargs['multiband_catalog_args']['get_all_filters'] = True
        tab = auto_script.multiband_catalog(field_root=root,
                                        **kwargs['multiband_catalog_args'])
    
        os.system('rm *bkg.fits')
        os.system('gzip {0}*seg.fits'.format(root))
        os.system('aws s3 sync --exclude "*" --include "{0}*seg.fits.gz" --include "{0}*cat.fits" --include "{0}*phot.fits" ./ s3://grizli-v1/Mosaics/ --acl public-read'.format(root))
        os.system('gunzip {0}*seg.fits.gz'.format(root))
    
    grism_files = {}
    for g in grisms:
        grism_files[g] = []
        
    for v in visits:
        has_grism = False
        for g in grisms:
            has_grism |= ('-'+g in v['product'])

        if not has_grism:
            continue
        else:
            print('Add visit {0}'.format(v['product']))
        
        for g in grisms:
            if '-'+g in v['product']:
                grism_files[g].extend(v['files'])        
                 
        for aws, f in zip(v['awspath'], v['files']):
            if os.path.exists(f):
                continue
            else:
                os.system('aws s3 cp s3://{0}/{1} .'.format(aws, f))
            
    for g in grisms:
        os.chdir('../Prep')
        print(g, len(grism_files[g]), len(np.unique(grism_files[g])))
        grism_files[g].sort()
        
        kwargs['grism_prep_args']['files'] = grism_files[g]
        grp = auto_script.grism_prep(field_root=root, **kwargs['grism_prep_args'])
        del(grp)
        
        grism_flt = []
        for file in grism_files[g]:
            grism_flt.extend(glob.glob(file.split('_fl')[0]+'*GrismFLT.fits'))
            
        catalog = glob.glob('{0}-*.cat.fits'.format(root))[0]
        try:
            seg_file = glob.glob('{0}-*_seg.fits'.format(root))[0]
        except:
            seg_file = None
            
        grp = multifit.GroupFLT(grism_files=grism_flt, direct_files=[], ref_file=None, seg_file=seg_file, catalog=catalog, cpu_count=-1, sci_extn=1, pad=256)
        
        # Make drizzle model images
        grp.drizzle_grism_models(root=root, kernel='point', scale=0.15)
        del(grp)
        os.system('aws s3 sync --exclude "*" --include "{0}*_grism*" ./ s3://grizli-v1/GrismMosaics/ --acl public-read'.format(key))
        os.system('aws s3 sync --exclude "*" --include "*GrismFLT*" --include "*wcs.fits" ./ s3://grizli-v1/GrismMosaics/ --acl public-read'.format(key))
    
    
    ################
      
    # Extractions
    pline = auto_script.DITHERED_PLINE
    auto_script.generate_fit_params(field_root=root, prior=None, MW_EBV=0.019, pline=pline, fit_only_beams=True, run_fit=True, poly_order=7, fsps=True, min_sens=0.001, sys_err=0.03, fcontam=0.2, zr=[0.05, 3.4], save_file='fit_args.npy', fit_trace_shift=False, include_photometry=True, use_phot_obj=False)
    
    ids = [27115]
    auto_script.extract(field_root=root, maglim=[13, 24], prior=None, MW_EBV=0.019, ids=ids, pline=pline, fit_only_beams=True, run_fit=False, poly_order=7, oned_R=30, master_files=None, grp=grp, bad_pa_threshold=None, fit_trace_shift=False, size=32, diff=True, min_sens=0.01, fcontam=0.2, min_mask=0.01, sys_err=0.03, skip_complete=True, args_file='fit_args.npy', get_only_beams=False)
    fitting.run_all_parallel(ids[0], verbose=True) 
    
    
def define_tiles(ra=109.3935148, dec=37.74934031, size=(24, 24), tile_size=6, overlap=0.3, field='macs0717', pixscale=0.03, theta=0):
    """
    Make tiles
    """
    from grizli import utils
    import astropy.wcs as pywcs
    from collections import OrderedDict
            
    tile_wcs = OrderedDict()
    
    size_per = tile_size-overlap
    nx = int(np.ceil(size[0]/size_per))
    ny = int(np.ceil(size[1]/size_per))
    
    sx = nx*tile_size-(nx-1)*overlap
    sy = ny*tile_size-(ny-1)*overlap
    px = int(tile_size*60/pixscale)

    # Even number of pixels
    header, parent_wcs = utils.make_wcsheader(ra=ra, dec=dec, size=(sx*60, sy*60), pixscale=pixscale, theta=theta)
    for i in range(nx):
        x0 = int(i*px-i*(overlap/tile_size*px))
        slx = slice(x0, x0+px)
        for j in range(ny):
            y0 = int(j*px-j*(overlap/tile_size*px))
            sly = slice(y0, y0+px)
            
            slice_header = utils.get_wcs_slice_header(parent_wcs, slx, sly)
            slice_wcs = pywcs.WCS(slice_header)
            slice_wcs._header = slice_header
            #tile=(i+1)*10+(j+1)
            tile = '{0:02d}.{1:02d}'.format(i+1, j+1)
            tile_wcs[tile] = slice_wcs
    
    np.save('{0}_{1:02d}mas_tile_wcs.npy'.format(field, int(pixscale*1000)),[tile_wcs])
    
    fpr = open('{0}_{1:02d}mas_tile_wcs.reg'.format(field, int(pixscale*1000)),'w')
    fpr.write('fk5\n')
    for t in tile_wcs:
        fp = tile_wcs[t].calc_footprint()
        pstr = 'polygon('+','.join(['{0:.6f}'.format(i) for i in fp.flatten()])+') # text={{{0}}}\n'.format(t)
        fpr.write(pstr)
    fpr.close()
    
    return tile_wcs
    
####
def drizzle_tiles(visits, tiles, prefix='gdn', filts=['f160w','f140w','f125w','f105w','f814w','f098m','f606w','f475w','f850lp'], pixfrac=0.5, output_bucket=None):   
    """
    mkdir ~/CosmosMosaic
    cd ~/CosmosMosaic
    
    aws s3 cp s3://grizli-preprocess/CosmosMosaic/cosmos_acs_tile_wcs.npy .
    aws s3 cp s3://grizli-preprocess/CosmosMosaic/macs0717_30mas_tile_wcs.npy .
    aws s3 sync --exclude "*" --include "cosmos_visits*" s3://grizli-preprocess/CosmosMosaic/ ./ 

    aws s3 sync --exclude "*" --include "grizli_visits*" s3://grizli-preprocess/CosmosMosaic/ ./ 
    
    """ 
    import numpy as np
    import os
    import copy
    import glob
    
    import astropy.io.fits as pyfits
    import astropy.wcs as pwcs
    from grizli import utils, prep
                
    # By filter
    groups = {}
    for visit in visits:
        if 'footprints' not in visit:
            continue
        
        filt = visit['product'].split('-')[-1]
        if filt not in groups:
            groups[filt] = {'filter':filt, 'files':[], 'footprints':[], 'awspath':[]}
        
        for l in ['files', 'footprints', 'awspath']:
            groups[filt][l].extend(visit[l])
                    
    for t in list(tiles.keys()):
        tile = tiles[t]
        
        psquare = np.sqrt(tile.wcs.cd[0,0]**2+tile.wcs.cd[0,1]**2)
        fine_mas = np.abs(int(np.round(psquare*3600*1000)))
        
        
        h = tile._header
        tile._naxis = [h['NAXIS1'],h['NAXIS2']]
        
        wh = utils.to_header(tile)
        for k in ['NAXIS1','NAXIS2']:
            wh[k] //= 2
        
        wh['CRPIX1'] = (wh['CRPIX1']+0.5)/2.
        wh['CRPIX2'] = (wh['CRPIX2']+0.5)/2.
        
        for k in ['CD1_1','CD1_2', 'CD2_1','CD2_2']:
            if k in wh:
                wh[k] *= 2
        
        #w60 = pywcs.WCS(wh)
        h0 = utils.to_header(tile)
        d0 = np.zeros((h0['NAXIS2'], h0['NAXIS1']), dtype=np.uint8)
        d0[::2,:] = 1
        d0[:,::2] = 1
        pyfits.writeto('fine_grid.fits', data=d0, header=h0, overwrite=True, output_verify='fix')
        
        h1 = wh
        d1 = np.zeros((h1['NAXIS2'], h1['NAXIS1']), dtype=np.uint8)
        d1[::2,:] = 1
        d1[:,::2] = 1
        pyfits.writeto('coarse_grid.fits', data=d1, header=h1, overwrite=True, output_verify='fix')
                
        for filt in filts:
            if filt not in groups:
                continue
            
            visits = [copy.deepcopy(groups[filt])]
            try:
                os.remove('astrodrizzle.log')
            except:
                pass
                      
            version = 'v0.03'
            if filt.startswith('f1') | filt.startswith('f098m') | filt.startswith('g1'):                    
                visits[0]['reference'] = 'coarse_grid.fits'
                visits[0]['product'] = '{0}-{1}-{2:03d}mas-{3}'.format(prefix, t, fine_mas*2, filt)
                
            else:
                visits[0]['reference'] = 'fine_grid.fits'
                visits[0]['product'] = '{0}-{1}-{2:03d}mas-{3}'.format(prefix, t, fine_mas, filt)
            
            if len(glob.glob('{0}*_dr*gz'.format(visits[0]['product']))) > 0:
                continue
                
            old_files = glob.glob('{0}*'.format(visits[0]['product']))
            for file in old_files:
                os.remove(file)
                
            try:
                prep.drizzle_overlaps(visits, parse_visits=False, check_overlaps=True, pixfrac=pixfrac, skysub=False, final_wcs=True, final_wht_type='IVM', static=False, max_files=260, fix_wcs_system=True)
            except:
                os.system('date > {0}.failed'.format(visits[0]['product']))
                continue
                
            #os.system('rm *_fl*fits')

            # Combine split mosaics
            tile_files = glob.glob(visits[0]['product']+'-0*sci.fits')
            if len(tile_files) > 0:
                tile_files.sort()

                im = pyfits.open(visits[0]['reference'])
                img = np.zeros_like(im[0].data, dtype=np.float32)
                wht = np.zeros_like(im[0].data, dtype=np.float32)

                exptime = 0.
                ndrizim = 0.
                ext = 'sci'

                for i, tile_file in enumerate(tile_files):
                    im = pyfits.open(tile_file)
                    wht_i = pyfits.open(tile_file.replace('_sci.f', '_wht.f'))
                    print(i, filt, tile_file, wht_i.filename())

                    exptime += im[0].header['EXPTIME']
                    ndrizim += im[0].header['NDRIZIM']

                    if i == 0:
                        h = im[0].header

                    img += im[0].data*wht_i[0].data
                    wht += wht_i[0].data

                sci = img/wht
                sci[wht == 0] = 0

                h['EXPTIME'] = exptime
                h['NDRIZIM'] = ndrizim

                sci_file = '{0}_drz_sci.fits'.format(visits[0]['product'])
                wht_file = '{0}_drz_wht.fits'.format(visits[0]['product'])
                
                pyfits.writeto(sci_file, data=sci, header=h, overwrite=True)
                pyfits.writeto(wht_file, data=wht, header=h, overwrite=True)
            
            print('gzip {0}_dr*fits'.format(visits[0]['product']))
            os.system('gzip {0}_dr*fits'.format(visits[0]['product']))
            
            if output_bucket is not None:
                os.system('aws s3 sync --exclude "*" --include "{0}_dr*" ./ {1} --acl public-read'.format(visits[0]['product'], output_bucket))
            
            new_files = glob.glob('{0}*fits'.format(visits[0]['product']))
            for file in new_files:
                print('rm '+file)
                os.remove(file)
            
            