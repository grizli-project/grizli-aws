import os
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import glob

from grizli import prep, utils

IR_FILTERS = ['f098m', 'f110w', 'f140w', 'f105w', 'f160w','f125w' ]
OPTICAL_FILTERS = ['f606w', 'f775w', 'f814w', 'f850lp']
GRISMS = ['g102', 'g141', 'g800l']

def test_egs():
    
    size=(80,20)
    theta = 49.699
    dy = -1./60
    
    fields = {}
    fields['egs-rot'] = {'field':'egs-rot', 'ra':214.8288, 'dec':52.8234+dy, 'size':size, 'tile_size':6, 'overlap':0.3, 'pixscale':0.05, 'theta':theta}
    
    key = 'egs-rot'
    tiles = define_tiles(**fields[key])
    os.system('xpaset -p ds9 regions delete all; ds9_reg egs_50mas_tile_wcs_magenta.reg ; ds9_reg egs-rot_50mas_tile_wcs.reg')

def tile_parameters(pixscale=0.05, overlap=0.3, tile_size=6):
    """
    pixscale = pixel scale of optical images.  IR will be double
    
    overlap = overlap of tiles in arcmin
    
    tile_size = tile size in arcmin
    
    """
    fields = {}
    fields['gdn'] = {'field':'gdn', 'ra':189.236, 'dec':62.257, 'size':(24,24), 'tile_size':tile_size, 'overlap':overlap, 'pixscale':pixscale}
    fields['uds'] = {'field':'uds', 'ra':34.361, 'dec':-5.163, 'size':(36,24), 'tile_size':tile_size, 'overlap':overlap, 'pixscale':pixscale}
    fields['egs'] = {'field':'egs', 'ra':214.8288, 'dec':52.8234, 'size':(54, 60), 'tile_size':tile_size, 'overlap':overlap, 'pixscale':pixscale}
    fields['gds'] = {'field':'gds', 'ra':53.1180-0.5/60, 'dec':-27.8216-0.5/60, 'size':(32, 32), 'tile_size':tile_size, 'overlap':overlap, 'pixscale':pixscale}
    fields['cos'] = {'field':'cos', 'ra':150.125-0.041, 'dec':2.25-0.018, 'size':(90, 90), 'tile_size':tile_size, 'overlap':overlap, 'pixscale':pixscale}
    
    # Rotated EGS along strip
    size=(80,20)
    theta = 49.699
    dy = -1./60
    fields['egs'] = {'field':'egs', 'ra':214.8288, 'dec':52.8234+dy, 'size':size, 'tile_size':tile_size, 'overlap':overlap, 'pixscale':pixscale, 'theta':theta}
    
    return fields
#
def define_tiles(ra=109.3935148, dec=37.74934031, size=(24, 24), tile_size=6, overlap=0.3, field='macs0717', pixscale=0.05, theta=0):
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
def drizzle_tiles(visits, tiles, prefix='gdn', filts=['f160w','f140w','f125w','f105w','f814w','f098m','f606w','f475w','f850lp'], pixfrac=0.5, output_bucket=None, clean_intermediate=False, use_keys=None):   
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
      
    if use_keys is None:
        use_keys = list(tiles.keys())
                  
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
        if t not in use_keys:
            continue
        
        tile = tiles[t]
        
        try:
            psquare = np.sqrt(tile.wcs.cd[0,0]**2+tile.wcs.cd[0,1]**2)
        except:
            psquare = np.sqrt(tile.wcs.pc[0,0]**2+tile.wcs.pc[0,1]**2)
            
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
                print('Filter {0} not found in the visit groups'.format(filt))
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
                is_ir=True
                ref_header = h1
                dr = 'drz'
                
            else:
                visits[0]['reference'] = 'fine_grid.fits'
                visits[0]['product'] = '{0}-{1}-{2:03d}mas-{3}'.format(prefix, t, fine_mas, filt)
                is_ir=False
                ref_header = h0
                dr = 'drc'
            
            old_files = glob.glob('{0}*_dr*gz'.format(visits[0]['product']))    
            if len(old_files) > 0:
                print('Files found: {0}'.format(' '.join(old_files)))
                continue
                
            old_files = glob.glob('{0}*'.format(visits[0]['product']))
            for file in old_files:
                os.remove(file)
            
            print('\n\n\n#####\nDrizzle mosaic: {0}\n#####\n\n\n'.format(visits[0]['product']))
            
            try:
                #prep.drizzle_overlaps(visits, parse_visits=False, check_overlaps=True, pixfrac=pixfrac, skysub=False, final_wcs=True, final_wht_type='IVM', static=False, max_files=260, fix_wcs_system=True)
                #
                
                # Use compact drizzler
                ref_header = pyfits.open(visits[0]['reference'])[0].header
                status = utils.drizzle_from_visit(visits[0], ref_header, pixfrac=pixfrac, clean=clean_intermediate, include_saturated=is_ir, kernel='square') 
                outsci, outwht, outh = status
                
                prod = visits[0]['product']
                
                pyfits.writeto('{0}_{1}_sci.fits'.format(prod, dr), 
                               data=outsci, header=outh, overwrite=True, 
                               output_verify='fix')

                pyfits.writeto('{0}_{1}_wht.fits'.format(prod, dr), 
                               data=outwht, header=outh, overwrite=True, 
                               output_verify='fix')
                
                #print(visits[0])
                
            except:
                os.system('date > {0}.failed'.format(visits[0]['product']))
                continue
                
            #os.system('rm *_fl*fits')

            # Combine split mosaics
            # tile_files = glob.glob(visits[0]['product']+'-0*sci.fits')
            # if len(tile_files) > 0:
            #     tile_files.sort()
            # 
            #     im = pyfits.open(visits[0]['reference'])
            #     img = np.zeros_like(im[0].data, dtype=np.float32)
            #     wht = np.zeros_like(im[0].data, dtype=np.float32)
            # 
            #     exptime = 0.
            #     ndrizim = 0.
            #     ext = 'sci'
            # 
            #     for i, tile_file in enumerate(tile_files):
            #         im = pyfits.open(tile_file)
            #         wht_i = pyfits.open(tile_file.replace('_sci.f', '_wht.f'))
            #         print(i, filt, tile_file, wht_i.filename())
            # 
            #         exptime += im[0].header['EXPTIME']
            #         ndrizim += im[0].header['NDRIZIM']
            # 
            #         if i == 0:
            #             h = im[0].header
            # 
            #         img += im[0].data*wht_i[0].data
            #         wht += wht_i[0].data
            # 
            #     sci = img/wht
            #     sci[wht == 0] = 0
            # 
            #     h['EXPTIME'] = exptime
            #     h['NDRIZIM'] = ndrizim
            # 
            #     sci_file = '{0}_drz_sci.fits'.format(visits[0]['product'])
            #     wht_file = '{0}_drz_wht.fits'.format(visits[0]['product'])
            #     
            #     pyfits.writeto(sci_file, data=sci, header=h, overwrite=True)
            #     pyfits.writeto(wht_file, data=wht, header=h, overwrite=True)
            
            print('gzip {0}_dr*fits'.format(visits[0]['product']))
            os.system('gzip {0}_dr*fits'.format(visits[0]['product']))
            
            if output_bucket is not None:
                os.system('aws s3 sync --exclude "*" --include "{0}_dr*" ./ {1} --acl public-read'.format(visits[0]['product'], output_bucket))
            
            new_files = glob.glob('{0}*fits'.format(visits[0]['product']))
            for file in new_files:
                print('rm '+file)
                os.remove(file)
        
def make_candels_tiles(key='gdn', filts=OPTICAL_FILTERS, pixfrac=0.33, bucket='grizli-v1', output_bucket='s3://grizli-v1/Mosaics/', clean_intermediate=False):
    
    import numpy as np
    import matplotlib.pyplot as plt
    import astropy.io.fits as pyfits
    import glob
    from grizli import prep, utils
    import os
    
    from grizli_aws.field_tiles import define_tiles, drizzle_tiles, tile_parameters
    
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
    
    fields = tile_parameters(pixscale=0.05, overlap=0.3, tile_size=6)
    
    tiles = define_tiles(**fields[key])
            
    #all_visits, all_groups, all_info = np.load('goodss-j033236m2748_visits.npy')
    
    #filts = ['f098m', 'f110w', 'f140w', 'f105w', 'f160w','f125w' ]
    
    #filts = ['f606w', 'f775w', 'f814w', 'f850lp']
    
    tiles = np.load('{0}_50mas_tile_wcs.npy'.format(key))[0]
    
    visit_file = glob.glob('{0}-j*visits.npy'.format(key))
    if len(visit_file) == 0:
        os.system('aws s3 sync s3://{0}/Mosaics/ ./ --exclude "*" --include "{1}*npy"'.format(bucket, key))
        visit_file = glob.glob('{0}-j*visits.npy'.format(key))
    
    all_visits, all_groups, info = np.load(visit_file[0])
    
    # Replace aws path for IR
    if bucket == 'grizli-v1':
        for v in all_visits:
            if '_flt' in v['files'][0]:
                v['awspath'] = []
                for file in v['files']:
                    file_root = os.path.basename(file)
                    prog = file_root[:4]
                    dataset = file_root[:9]

                    aws = os.path.join('grizli-v1', 'Exposures', prog, dataset)
                    v['awspath'].append(aws)
            
    if False:
        # test single tile
        t = '01.01'
        one_tile = {}
        one_tile[t] = tiles[t]
        drizzle_tiles(all_visits, one_tile, filts=filts, prefix=key, pixfrac=pixfrac, output_bucket=output_bucket)
    
    drizzle_tiles(all_visits, tiles, filts=filts, prefix=key, pixfrac=pixfrac, output_bucket=output_bucket, clean_intermediate=clean_intermediate)

def combine_tile_filters(key='egs', skip_existing=True, combine_fnu=True, bands=['ir', 'opt']):
    """
    Make combined images weighted by photflam
    """
    tile_groups = []
    for ix in range(20):
        for iy in range(20):
            tile_files = glob.glob('{0}-{1:02d}.{2:02d}-*-f*sci.fits*'.format(key, ix, iy))
            if len(tile_files) == 0:
                continue
            tile_groups.append(tile_files)
            
    if len(tile_groups) == 0:
        tile_groups = [glob.glob('{0}-???mas-f*sci.fits*'.format(key))]
        
    for tile_files in tile_groups:    
        # existing = glob.glob('{0}-{1:02d}.{2:02d}-*-ir*sci.fits*'.format(key, ix, iy))
        # existing += glob.glob('{0}-{1:02d}.{2:02d}-*-opt*sci.fits*'.format(key, ix, iy))
        # 
        # if (len(existing) > 0) & (skip_existing):
        #     continue
            
        #print(tile_files[0])
        tile_files.sort()
        
        ####
        num = {'ir':None, 'opt':None}
        den = {'ir':None, 'opt':None}
        
        if False:
            im = pyfits.open('egs-05.02-050mas-f814w_drc_sci.fits.gz')
            im = pyfits.open('egs-05.02-100mas-f140w_drz_sci.fits.gz')
            
            wh = {}
            for k in ['INSTRUME', 'DETECTOR', 'PHOTFNU', 'PHOTFLAM', 'PHOTBW','PHOTZPT', 'PHOTMODE', 'PHOTPLAM', 'FILTER', 'FILTER1', 'FILTER2']:
                if k in im[0].header:
                    wh[k] = im[0].header[k]
        
        ref_h = {}
        ref_h['opt'] = {'INSTRUME': 'ACS', 'DETECTOR': 'WFC', 'PHOTFLAM': 7.0178627203125e-20, 'PHOTBW': 653.24393453125, 'PHOTZPT': -21.1, 'PHOTMODE': 'ACS WFC1 F814W MJD#56438.5725', 'PHOTPLAM': 8045.415190625002, 'FILTER1': 'CLEAR1L', 'FILTER2': 'F814W'}
        ref_h['ir'] = {'INSTRUME': 'WFC3', 'DETECTOR': 'IR', 'PHOTFNU': 9.5291135e-08, 'PHOTFLAM': 1.4737148e-20, 'PHOTBW': 1132.39, 'PHOTZPT': -21.1, 'PHOTMODE': 'WFC3 IR F140W', 'PHOTPLAM': 13922.907, 'FILTER': 'F140W'}
        
        output_sci = {}
        head = {}
        
        for sci_file in tile_files:
            filt_i = sci_file.split('_dr')[0].split('-')[-1]
            if filt_i in IR_FILTERS:
                band = 'ir'
            else:
                band = 'opt'
            
            if band not in bands:
                continue
                    
            print(sci_file, band)
            output_sci[band] = sci_file.replace(filt_i, band)
            
            im_i = pyfits.open(sci_file)
            wht_i = pyfits.open(sci_file.replace('_sci','_wht'))
            photflam = im_i[0].header['PHOTFLAM']
            ref_photflam = ref_h[band]['PHOTFLAM']

            photplam = im_i[0].header['PHOTPLAM']
            ref_photplam = ref_h[band]['PHOTPLAM']
            
            head[band] = im_i[0].header
            for k in ref_h[band]:
                head[band][k] = ref_h[band][k]
                
            if num[band] is None:
                num[band] = im_i[0].data*0
                den[band] = num[band]*0
                
            scl = photflam/ref_photflam
            if combine_fnu:
                scl *= photplam**2/ref_photplam**2
                
            den_i = wht_i[0].data/scl**2
            num[band] += im_i[0].data*scl*den_i
            den[band] += den_i

        for band in ['opt', 'ir']:
            if num[band] is not None:
                sci = num[band]/den[band]
                wht = den[band]
                
                mask = (~np.isfinite(sci)) | (den == 0)
                sci[mask] = 0
                wht[mask] = 0

                print('Write {0}'.format(output_sci[band]))
                
                pyfits.PrimaryHDU(data=sci, header=head[band]).writeto(output_sci[band], overwrite=True, output_verify='fix')
                pyfits.PrimaryHDU(data=wht, header=head[band]).writeto(output_sci[band].replace('_sci', '_wht'), overwrite=True, output_verify='fix')
            
        ### TBD - full combination of OPT + IR
            
def make_all_tile_catalogs(key='egs', threshold=3, filts=OPTICAL_FILTERS+IR_FILTERS, combine_catalogs=True):

    # Make catalogs
    import os
    import glob
    import astropy.table
    import numpy as np
    
    from grizli import utils, prep
    import sep

    for filt_i in filts:
        files = glob.glob(key+'-[0-9][0-9].[0-9][0-9]*-'+filt_i+'*sci.fits.gz')
        files.sort()
        for i, file in enumerate(files):

            froot = file.split('_dr')[0]
            
            print('\n\n#######\n{0}\n#######\n'.format(froot))
            
            if os.path.exists('{0}.cat.fits'.format(froot)):
                continue
            
            try:
                prep.make_SEP_catalog(froot, threshold=threshold,
                                  rescale_weight=False) 
                os.system('rm {0}_bkg.fits'.format(froot))
                #os.system('rm {0}_seg.fits'.format(froot))
            except:
                print('\n\n !!!! Failed: ', froot)
                pass
        
        if combine_catalogs:
            
            cat_files = glob.glob(key+'-??.??-*{0}*cat.fits'.format(filt_i))
            if len(cat_files) == 0:
                continue

            cat_files.sort()

            cats = [utils.read_catalog(file) for file in cat_files]

            ref = cats[0].copy()

            id_max = [ref['NUMBER'].max()]

            for ic in range(1, len(cats)):
                print(cat_files[ic])
                cat = cats[ic]
                not_edge = (cat['FLAG'] & sep.OBJ_TRUNC) == 0 

                idx, dr = cat.match_to_catalog_sky(ref)     
                new = dr.value > 0.1

                id_max.append(cat['NUMBER'].max())

                cat['NUMBER'] += np.sum(id_max[:-1])

                ref = astropy.table.vstack([ref[new], cat])

            out_root = '{0}-{1}_tiles'.format(key, filt_i)
            ref.write(out_root+'.fits', overwrite=True)
            label = [int(i) for i in ref['NUMBER']]
            prep.table_to_regions(ref, out_root+'.reg', comment=label)
                  
def combine_tile_mosaics(key='gdn', filts=OPTICAL_FILTERS, use_ref='*', extensions = ['sci','wht'], tiles=None, sync=True, bucket='grizli-v1', use_files=None, ref_pixscale=50):
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    import astropy.io.fits as pyfits
    import glob

    from grizli import prep, utils
        
    ##########
    # Combine final mosaics
    
    filt_i ='f814w'
    
    # Compute tiles needed for full mosaic extent
    init = True
    #use_ref = 'all'
    
    if use_ref == 'wfc3ir':
        # wfc3/ir
        files = glob.glob(key+'*-[0-9][0-9].[0-9][0-9]*{0}mas*-*sci.fits.gz'.format(ref_pixscale*2))
    elif use_ref == 'acswfc':
        # acs
        files = glob.glob(key+'*-[0-9][0-9].[0-9][0-9]*{0}mas*-*sci.fits.gz'.format(ref_pixscale))
    elif use_files is not None:
        files = use_files
    else:
        files = glob.glob(key+'*-[0-9][0-9].[0-9][0-9]-{0}sci.fits.gz'.format(use_ref))

    files.sort()

    tile_pos = np.array([int(file.split('-')[-3].replace('.','')) for file in files])
    xp = tile_pos // 100
    yp = tile_pos % 100
    
    left, bot = xp.min(), yp.min()

    nx = (xp.max()-xp.min())+1
    ny = (yp.max()-yp.min())+1
    print('Define mosaic for {0} / {1}: {2} x {3}'.format(key, use_ref, nx, ny))
    
    # WCS defined in lower left corner
    ll = '{0:02d}.{1:02d}'.format(left, bot)
    if tiles is None:
        tiles = np.load('{0}_{1}mas_tile_wcs.npy'.format(key, ref_pixscale))[0]
    
    wcs = tiles[ll]

    ## Loop over filters
    #filts = ['f098m','f110w','f140w','f105w','f160w','f125w']
    for filt_i in filts:
    
        files = glob.glob('*-[0-9][0-9].[0-9][0-9]*-'+filt_i+'*sci.fits.gz')
        files.sort()
        if len(files) == 0:
            continue
            
        tile_pos = np.array([int(file.split('-')[-3].replace('.','')) for file in files])
        xp = tile_pos // 100
        yp = tile_pos % 100
        
        ref = pyfits.open(files[0])[0]
        sh = ref.data.shape
        the_filter = utils.get_hst_filter(ref.header)
    
        is_fine = ('{0}mas'.format(ref_pixscale) in files[0])
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
        coarse_pix = 2*ref_pixscale/1000 #0.1
        overlap_arcmin = 0.3    
        olap = int(overlap_arcmin*60/coarse_pix*(1+is_fine))
    
        for ext in extensions:
            if ext == 'seg':
                data = np.zeros((sh[0]*ny-olap*(ny-1), sh[1]*nx-olap*(nx-1)), dtype=np.int32)
                out_ext = 'seg'
            else:
                data = np.zeros((sh[0]*ny-olap*(ny-1), sh[1]*nx-olap*(nx-1)), dtype=np.float32)
                out_ext = 'drz_'+ext
                
            data_sh = data.shape
            for i, file in enumerate(files):
                x0 = (xp[i]-left)*(sh[1]-olap)
                y0 = (yp[i]-bot)*(sh[0]-olap)
                                    
                slx = slice(x0,x0+sh[1])
                sly = slice(y0,y0+sh[0])

                file_ext = file.replace('drz_sci', out_ext)    
                file_ext = file_ext.replace('drc_sci', out_ext.replace('drz','drc'))    

                if (x0 < 0) | (y0 < 0):
                    print('skip', file_ext)
                    continue
                elif (slx.stop > data_sh[1]) | (sly.stop > data_sh[0]):
                    print('skip', file_ext)
                    continue
                else:
                    print(file_ext)
                
                if not os.path.exists(file_ext):
                    continue
                    
                sci_i = pyfits.open(file_ext)
                data[sly, slx] = sci_i[0].data
        
            mos = pyfits.PrimaryHDU(data=data, header=wh)
            outfile = '{0}-{1:03d}mas-{2}_{3}.fits'.format(key, ref_pixscale*(2-is_fine), filt_i.lower(), out_ext)
            mos.writeto(outfile, overwrite=True)
        
            del(data)
    
    # Gzip and sync
    if sync:
        os.system('gzip --force {0}-???mas*fits'.format(key))
        os.system('aws s3 sync --exclude "*" --include "{0}-???mas*gz" ./ s3://{1}/Mosaics/ --acl public-read'.format(key, bucket))

def full_processing():
    """
    Make mosaics and full-field catalogs
    """
    
    import os
    import glob
    
    import numpy as np
    from grizli_aws import field_tiles
    
    key='egs'
    
    bucket='grizli-v1'
    key='cos'; bucket='grizli-cosmos-v2'
    
    os.system('aws s3 sync s3://{1}/Mosaics/ ./ --exclude "*" --include "{0}-??.??-*" --include "{0}*wcs.*"'.format(key, bucket))
    
    if False:
        os.system('aws s3 sync s3://grizli-v1/Mosaics/ ./ --exclude "*" --include "{0}-??.??-*f814*sci.fits.gz" --include "{0}*wcs.*"'.format(key))
        os.system('aws s3 sync s3://grizli-v1/Mosaics/ ./ --exclude "*" --include "{0}-??.??-*f435*.fits.gz" --include "{0}*wcs.*"'.format(key))
        
        
    fields = field_tiles.tile_parameters(pixscale=0.05, overlap=0.3, tile_size=6)
    tiles = field_tiles.define_tiles(**fields[key])
    
    # Original drizzle
    field_tiles.make_candels_tiles(key=key, filts=field_tiles.IR_FILTERS, pixfrac=0.33, output_bucket='s3://grizli-v1/Mosaics/', bucket='grizli-v1', clean_intermediate=False)

    field_tiles.make_candels_tiles(key=key, filts=['f606w'], pixfrac=0.33, output_bucket='s3://grizli-cosmos-v2/Mosaics/', bucket='grizli-cosmos-v2', clean_intermediate=True)
    field_tiles.make_candels_tiles(key=key, filts=['f350lp', 'f438w', 'f435w', 'f475w', 'f850lp'], pixfrac=0.33, output_bucket='s3://grizli-cosmos-v2/Mosaics/', bucket='grizli-cosmos-v2', clean_intermediate=True)
        
    # Combined band images
    field_tiles.combine_tile_filters(key=key, skip_existing=True)
    
    if False:
        os.system('aws s3 sync s3://{1}/Mosaics/ ./ --exclude "*" --include "{0}-??.??-*" --include "{0}*wcs.*"'.format(key, bucket))
        
        field_tiles.combine_tile_filters(key=key, skip_existing=False)
        field_tiles.combine_tile_filters(key=key, skip_existing=False, bands=['opt'])

        field_tiles.combine_tile_filters(key=key, skip_existing=False, bands=['ir'])
        os.system('aws s3 sync ./ s3://{1}/Mosaics/ --exclude "*" --include "{0}-??.??-100mas-ir*" --acl public-read'.format(key, bucket))

        field_tiles.combine_tile_filters(key=key, skip_existing=False, bands=['opt'])
        
        field_tiles.make_all_tile_catalogs(key='cos', filts=['f350lp', 'f438w', 'f435w', 'f475w', 'f606w','f850lp']+field_tiles.IR_FILTERS)
        field_tiles.make_all_tile_catalogs(key='cos', filts=['f814w'])
        os.system('aws s3 sync ./ s3://{1}/Mosaics/ --exclude "*" --include "{0}*tiles.fits" --acl public-read'.format(key, bucket))
        
        
        use_ref = '*'
        filts = field_tiles.IR_FILTERS
        field_tiles.combine_tile_mosaics(key=key, filts=filts, use_ref=use_ref, extensions = ['sci','wht'], sync=True, bucket='grizli-cosmos-v2')
        
        # Region around CANDELS
        use_files = glob.glob('cos-07.07*sci.fits.gz')
        use_files += glob.glob('cos-07.12*sci.fits.gz')
        use_files += glob.glob('cos-09.07*sci.fits.gz')
        use_files += glob.glob('cos-09.12*sci.fits.gz')
        use_ref = '*'
        
        filts = field_tiles.IR_FILTERS
        filts = ['f814w']
        
        field_tiles.combine_tile_mosaics(key=key, filts=filts, use_ref=use_ref, extensions = ['sci','wht'], sync=True, bucket='grizli-cosmos-v2', use_files=use_files)
        
    # Tile mosaics
    use_ref='wfc3ir'
    field_tiles.combine_tile_mosaics(key=key, filts=field_tiles.IR_FILTERS, use_ref=use_ref, extensions = ['sci','wht'])
    field_tiles.combine_tile_mosaics(key=key, filts=field_tiles.OPTICAL_FILTERS, use_ref=use_ref, extensions = ['sci','wht'])
    field_tiles.combine_tile_mosaics(key=key, filts=['ir'], use_ref=use_ref, extensions = ['sci','wht'])
    field_tiles.combine_tile_mosaics(key=key, filts=['opt'], use_ref=use_ref, extensions = ['sci','wht'])
    
    # Full Catalog
    from grizli.pipeline import auto_script
    kwargs = auto_script.get_yml_parameters()
    
    keys = [key]
    keys = ['cos-cnd']
    
    keys = [f.split('-100mas')[0] for f in glob.glob('*-??.??-100mas-ir*sci.fits.gz')]
    keys.sort()
    
    for ckey in keys:
        kwargs['multiband_catalog_args']['detection_root'] = ckey+'-100mas-ir'
        #kwargs['multiband_catalog_args']['detection_root'] = ckey
        kwargs['multiband_catalog_args']['field_root'] = ckey+'-???mas'
        kwargs['multiband_catalog_args']['output_root'] = ckey+'-mosaic'

        kwargs['multiband_catalog_args']['rescale_weight'] = False
        kwargs['multiband_catalog_args']['det_err_scale'] = 1.

        kwargs['multiband_catalog_args']['filters'] = ['f160w','f140w','f125w','f110w','f105w','f098m','f850lp','f814w','f775w','f606w','f475w','f438w','f435w','f350lp'][::-1]
        auto_script.multiband_catalog(**kwargs['multiband_catalog_args'])
    
    bucket='grizli-cosmos-v2'
    
    os.system('aws s3 sync ./ s3://{0}/Mosaics/ --exclude "*" --include "*mosaic*" --include "*ir.cat.fits" --acl "public-read"'.format(bucket))
    
    os.system('gzip *ir_seg.fits') # *bkg.fits')
    #os.system('aws s3 sync ./ s3://{0}/Mosaics/ --exclude "*" --include "*bkg.fits.gz" --include "*-ir_seg.fits.gz" --acl "public-read"'.format(bucket))
    os.system('aws s3 sync ./ s3://{0}/Mosaics/ --exclude "*" --include "*-ir_seg.fits.gz" --acl "public-read"'.format(bucket))
    
def g800l_prep():
    
    import os
    import glob
    
    import numpy as np
    from grizli_aws import field_tiles
    
    from grizli.pipeline import auto_script
    kwargs = auto_script.get_yml_parameters()
    
    try:
        xxx = key
    except:
        key='egs'
    
    key0 = key
    if key == 'uds':
        field_root = root = 'uds-grism-j021732m0512'
        ref_filt = 'f814w'
    elif key == 'egs':
        field_root = root = 'egs-grism-j141956p5255'
        ref_filt = 'f814w'
    elif key == 'gds':
        field_root = root = 'gds-grism-j033236m2748'
        ref_filt = 'f850lp'
    elif key == 'gdn':
        field_root = root = 'gdn-grism-j123656p6215'
        ref_filt = 'f814w'
    elif key == 'cos':
        field_root = root = 'cos-grism-j100012p0210'
        ref_filt = 'f814w'
        key0 = 'cos-cnd'
        
    new_root = root.replace('-grism', '-g800l')
    
    bucket = 'grizli-v1'
    if key == 'cos':
        bucket = 'grizli-cosmos-v2'
        
    
    if key == 'egs':
        os.system('aws s3 sync s3://grizli-v1/Mosaics/ ./ --exclude "*" --include "egs-0[6-9].0[2-4]*{0}*fits.gz" --include "egs*npy"'.format(ref_filt))
    else:
        os.system('aws s3 sync s3://{2}/Mosaics/ ./ --exclude "*" --include "{0}-050mas-{1}*fits.gz" --include "{0}*npy"'.format(key0, ref_filt, bucket))
        os.system('gunzip {0}-050mas-{1}*fits.gz'.format(key0, ref_filt))
    
    
    os.system('aws s3 cp s3://{1}/Pipeline/{0}/Extractions/{0}-ir.cat.fits ./'.format(field_root, bucket))
    os.system('aws s3 cp s3://{1}/Pipeline/{0}/Extractions/{0}_phot.fits ./'.format(field_root, bucket))
    os.system('aws s3 cp s3://{1}/Pipeline/{0}/Extractions/{0}-ir_seg.fits.gz ./'.format(field_root, bucket))
    os.system('aws s3 cp s3://{1}/Pipeline/{0}/Prep/{0}-ir_seg.fits.gz ./'.format(field_root, bucket))
    os.system('gunzip *seg.fits.gz')
        
    # Tile mosaics
    if key == 'egs':
        fields = field_tiles.tile_parameters(pixscale=0.05, overlap=0.3, tile_size=6)
        tiles = field_tiles.define_tiles(**fields[key])
    
        use_ref='*'
        field_tiles.combine_tile_mosaics(key=key, filts=field_tiles.OPTICAL_FILTERS, use_ref=use_ref, extensions = ['sci','wht'], sync=False)
    
    os.system('cp {0}-050mas-{2}_drz_sci.fits {1}-{2}_drc_sci.fits'.format(key0, new_root, ref_filt))
    os.system('cp {0}-050mas-{2}_drz_wht.fits {1}-{2}_drc_wht.fits'.format(key0, new_root, ref_filt))

    os.system('cp {0}-050mas-{2}_drz_sci.fits {1}-ir_drc_sci.fits'.format(key0, new_root, ref_filt))
    os.system('cp {0}-050mas-{2}_drz_wht.fits {1}-ir_drc_wht.fits'.format(key0, new_root, ref_filt))
    
    if key == 'cos':
        os.system('cp {0}-100mas-ir_seg.fits {1}-ir_seg.fits'.format(key0, new_root))
        os.system('cp {0}-100mas-ir.cat.fits {1}-ir.cat.fits'.format(key0, new_root))
        os.system('cp {0}-mosaic_phot.fits {1}_phot.fits'.format(key0, new_root))
    else:
        os.system('cp {0}-ir_seg.fits {1}-ir_seg.fits'.format(root, new_root))
        os.system('cp {0}-ir.cat.fits {1}-ir.cat.fits'.format(root, new_root))
        os.system('cp {0}_phot.fits {1}_phot.fits'.format(root, new_root))
        
    #kwargs['multiband_catalog_args']['detection_root'] = new_root+'-f814w'
    # kwargs['multiband_catalog_args']['field_root'] = new_root+'-f*'
    # kwargs['multiband_catalog_args']['output_root'] = new_root+'-ir'

    # kwargs['multiband_catalog_args']['rescale_weight'] = True
    # #kwargs['multiband_catalog_args']['det_err_scale'] = 1.
    # 
    # kwargs['multiband_catalog_args']['filters'] = ['f814w']
    # kwargs['multiband_catalog_args']['threshold'] = 2.0
    # 
    # auto_script.multiband_catalog(field_root=new_root, **kwargs['multiband_catalog_args'])
    
    
    #################
    # Visits
    
    all_visits, all_groups, all_info = np.load(glob.glob('{0}-j*_visits.npy'.format(key))[0])
    
    # bucket = 'grizli-v1'
    # base_path='Exposures'
    # 
    # for v in all_visits:
    #     # Change aws path
    #     v['awspath'] = []
    #     for file in v['files']:
    #         file_root = os.path.basename(file)
    #         prog = file_root[:4]
    #         dataset = file_root[:9]
    #         
    #         aws = os.path.join(bucket, base_path, prog, dataset)
    #         
    #         v['awspath'].append(aws)
            
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
            
    np.save(new_root+'_visits.npy', [visits, groups, all_info[all_info['keep']]])
    
    os.system('aws s3 cp {0}_visits.npy s3://grizli-v1/Pipeline/{0}/Extractions/ --acl public-read'.format(new_root))
    
    #########################
    ####### Run model
    grisms = ['g800l']
    
    visits, groups, info = np.load(new_root+'_visits.npy')
                
    from collections import OrderedDict
    grism_files = OrderedDict()
    for g in grisms:
        grism_files[g] = []
        
    for v in visits:
        has_grism = False
        for g in grisms:
            #has_grism |= ('-'+g in v['product'])
            has_grism |= v['product'].endswith('-'+g)
        
        if not has_grism:
            continue
        else:
            print('Add visit {0}'.format(v['product']))
        
        for g in grisms:
            #if '-'+g in v['product']:
                
            if v['product'].endswith('-'+g) :
                grism_files[g].extend(v['files'])        
        
        # Copy file to local dir         
        for aws, f in zip(v['awspath'], v['files']):
            if os.path.exists(f):
                continue
            else:
                os.system('aws s3 cp s3://{0}/{1} .'.format(aws, f))
    
    from grizli import multifit, utils
    
    file_roots = {}
    for g in grisms:
        file_roots[g] = np.array([f[1:4] for f in grism_files[g]])
        for r in np.unique(file_roots[g]):
            print(g, r, (file_roots[g] == r).sum())
    
    from grizli.utils import column_string_operation as colstr
    # Split GDN G800L
    if (key == 'gdn'):
        figs  = colstr(file_roots[g], ['coi'], 'count', 'or')
        barro = colstr(file_roots[g], ['cat'], 'count', 'or')
        pears = colstr(file_roots[g], ['9fa'], 'startswith', 'or')
        hdfn  = colstr(file_roots[g], ['8eb','8n1'], 'startswith', 'or')

        other = ~(figs | barro | pears)
        #grism_files['g800l-figs'] = list(np.array(grism_files[g])[figs])        
        #grism_files['g800l-barro'] = list(np.array(grism_files[g])[barro])
        #grism_files['g800l-pears'] = list(np.array(grism_files[g])[pears])
        grism_files['g800l-hdfn'] = list(np.array(grism_files[g])[hdfn])
        #grism_files['g800l-other'] = list(np.array(grism_files[g])[other])
        grism_files.pop('g800l')
        
    elif (key == 'gds'):
        figs  = colstr(file_roots[g], ['coi'], 'count', 'or')
        vd    = colstr(file_roots[g], ['bhj'], 'count', 'or')
        udf   = colstr(file_roots[g], ['8qq'], 'startswith', 'or')
        pears = colstr(file_roots[g], ['9fa'], 'startswith', 'or')
        
        other = ~(figs | vd | pears)
        #grism_files['g800l-figs'] = list(np.array(grism_files[g])[figs])
        #grism_files['g800l-pears'] = list(np.array(grism_files[g])[pears])
        grism_files['g800l-vd'] = list(np.array(grism_files[g])[vd])
        grism_files['g800l-vd'].sort()
        #grism_files['g800l-vd1'] = grism_files['g800l-vd'][:76]
        grism_files['g800l-vd2'] = grism_files['g800l-vd'][76:]
        grism_files.pop('g800l-vd')
        
        #grism_files['g800l-udf'] = list(np.array(grism_files[g])[udf])
        #grism_files['g800l-other'] = list(np.array(grism_files[g])[other])
        
        grism_files.pop('g800l')
    
    print('\n\n')
            
    for g in grism_files:
        print(g, len(grism_files[g]), len(np.unique(grism_files[g])))
        os.chdir('../Prep')
        grism_files[g].sort()
        
        # Files causing memory errors
        pop_files = ['j8qq32p0q_flc.fits', 'j8qq32p6q_flc.fits', 'j9faj8xxq_flc.fits', 'j8n1t1qqq_flc.fits', 'jbhj33tqq_flc.fits']
        
        for p in pop_files:
            if p in grism_files[g]:
                if p in grism_files[g]:
                    pp  = grism_files[g].pop(grism_files[g].index(p))
                
        kwargs['grism_prep_args']['files'] = grism_files[g]
        kwargs['grism_prep_args']['refine_niter'] = 2
        grp = auto_script.grism_prep(field_root=new_root, **kwargs['grism_prep_args'], cpu_count=1)
        del(grp)
        
        grism_flt = []
        for file in kwargs['grism_prep_args']['files']:
            grism_flt.extend(glob.glob(file.split('_fl')[0]+'*GrismFLT.fits'))
            
        catalog = glob.glob('{0}-*.cat.fits'.format(new_root))[0]
        try:
            seg_file = glob.glob('{0}-*_seg.fits'.format(new_root))[0]
        except:
            seg_file = None
            
        grp = multifit.GroupFLT(grism_files=grism_flt, direct_files=[], ref_file=None, seg_file=seg_file, catalog=catalog, cpu_count=-1, sci_extn=1, pad=256)
        
        # Make drizzle model images
        grp.drizzle_grism_models(root=new_root, kernel='point', scale=0.15)
        del(grp)
        
        os.system('aws s3 sync --exclude "*" --include "{0}*_grism*" ./ s3://{1}/Pipeline/{0}/Extractions/ --acl public-read'.format(new_root, bucket))
        os.system('aws s3 sync --exclude "*" --include "*GrismFLT*" --include "*wcs.fits" ./ s3://{1}/Pipeline/{0}/Extractions/ --acl public-read'.format(new_root, bucket))

    ########
    os.system('aws s3 sync --exclude "*" --include "{0}*ir.cat.fits*" --include "{0}*phot.fits*" ./ s3://grizli-v1/Pipeline/{0}/Extractions/ --acl public-read'.format(new_root))
        
    MW_EBV = {'gdn':0.011, 'gds':0.007, 'egs': 0.007, 'uds':0.019, 'cos':0.0148}

    pline = auto_script.DITHERED_PLINE.copy()
    pline['pixscale'] = 0.05
    pline['size'] = 4
    
    auto_script.generate_fit_params(field_root=new_root, prior=None, MW_EBV=MW_EBV[key], pline=pline, fit_only_beams=True, run_fit=True, poly_order=7, fsps=True, min_sens=0.00, sys_err=0.03, fcontam=0.2, zr=[0.05, 3.4], save_file='fit_args.npy', fit_trace_shift=True, include_photometry=True, use_phot_obj=False, scale_photometry=False)
    args = np.load('fit_args.npy')[0]
    args['sed_args'] = {'bin': 4, 'xlim': [0.3, 2.5]}
    np.save('fit_args.npy', [args])
    
    os.system('aws s3 sync --exclude "*" --include "{0}*ir.cat.fits*" --include "{0}*phot.fits*" --include "fit_args.npy*" ./ s3://{1}/Pipeline/{0}/Extractions/ --acl public-read'.format(new_root, bucket))
    
    from grizli.aws import lambda_handler
    lambda_handler.extract_beams_from_flt(new_root, bucket, id, clean=False)
     
    os.system('aws s3 sync --exclude "*" --include "{0}_{1:05d}*" ./ s3://grizli-v1/Pipeline/{0}/Extractions/ --acl public-read'.format(new_root, id))
    
    # New Extractions
    os.system('aws s3 sync s3://grizli-v1/Pipeline/{0}/Extractions/ ./ --exclude "*" --include "*phot.fits" --include "*wcs.fits"'.format(new_root))
    
    from grizli import utils
    import numpy as np
    import glob
    
    phot = utils.read_catalog('{0}_phot.fits'.format(new_root))
    phot['has_grism'] = 0
    for f in glob.glob('*wcs.fits'):
        w = utils.WCSFootprint(f, ext=0)
        phot['has_grism'] += w.path.contains_points(np.array([phot['ra'], phot['dec']]).T)
        
    imag = 23.9-2.5*np.log10(phot['f814w_flux_aper_1']*phot['flux_auto']/phot['flux_aper_1']*phot['f814w_tot_corr'])
    sel = np.isfinite(imag) & (imag > 20) & (imag < 26) & (phot['has_grism'] > 8)
    from grizli.aws import fit_redshift_lambda
    fit_redshift_lambda.fit_lambda(root=new_root, beams=[], ids=phot['number'][sel], newfunc=False, bucket_name='grizli-v1', skip_existing=True, sleep=False, skip_started=True, quasar_fit=False, output_path=None, show_event=False, zr=[0.01,2])
    
def grism_prep():
    
    import numpy as np
    from grizli import utils, multifit
    from grizli.pipeline import auto_script
    import astropy.io.fits as pyfits
    import os
    import glob
    
    from grizli.pipeline import auto_script
    kwargs = auto_script.get_yml_parameters()
    
    if False:
        # Test cdatalog parameters
        kwargs['multiband_catalog_args']['detection_root'] = 'gds-03.04-100mas-ir'
        kwargs['multiband_catalog_args']['field_root'] = 'gds-03.04-*mas'
        kwargs['multiband_catalog_args']['output_root'] = 'gds-03.04'
        kwargs['multiband_catalog_args']['filters'] = ['f160w','f140w','f125w','f110w','f105w','f098m','f850lp','f814w','f775w','f606w'][::-1]
        auto_script.multiband_catalog(**kwargs['multiband_catalog_args'])

        # Full-field multiband catalog
        key = 'gds'
        
        kwargs['multiband_catalog_args']['detection_root'] = key+'-100mas-ir'
        kwargs['multiband_catalog_args']['field_root'] = key+'-???mas'
        kwargs['multiband_catalog_args']['output_root'] = key+'-mosaic'

        kwargs['multiband_catalog_args']['rescale_weight'] = False
        kwargs['multiband_catalog_args']['det_err_scale'] = 1.

        kwargs['multiband_catalog_args']['filters'] = ['f160w','f140w','f125w','f110w','f105w','f098m','f850lp','f814w','f775w','f606w'][::-1]
        auto_script.multiband_catalog(**kwargs['multiband_catalog_args'])

        os.system('aws s3 sync ./ s3://grizli-v1/Mosaics/ --exclude "*" --include "*mosaic*" --include "*ir.cat.fits"')
        
        os.system('gzip *ir_seg.fits') # *bkg.fits')
        os.system('aws s3 sync ./ s3://grizli-v1/Mosaics/ --exclude "*" --include "*bkg.fits.gz" --include "*-ir_seg.fits.gz"')
    
    key0 = key    
    if key == 'uds':
        field_root = root = 'uds-grism-j021732m0512'
    elif key == 'egs':
        field_root = root = 'egs-grism-j141956p5255'
    elif key == 'gds':
        field_root = root = 'gds-grism-j033236m2748'
    elif key == 'gdn':
        field_root = root = 'gdn-grism-j123656p6215'
    elif key == 'cos':
        field_root = root = 'cos-grism-j100012p0210'
        key0 = 'cos-cnd'
    
    bucket = 'grizli-v1'
    if key == 'cos':
        bucket = 'grizli-cosmos-v2'
    
    os.system('aws s3 sync s3://{0}/Mosaics/ ./ --exclude "*" --include "{1}*"'.format(bucket, key0))
    os.system('aws s3 sync s3://{0}/Mosaics/ ./ --exclude "*" --include "{1}*"'.format(bucket, field_root))
    
    os.system('cp -f {0}-100mas-ir_seg.fits.gz {1}-ir_seg.fits.gz'.format(key0, root, filt))
    os.system('cp -f {0}-100mas-ir_drz_sci.fits.gz {1}-ir_drz_sci.fits.gz'.format(key0, root, filt))
    os.system('cp -f {0}-100mas-ir_drz_wht.fits.gz {1}-ir_drz_wht.fits.gz'.format(key0, root, filt))
    os.system('cp -f {0}-100mas-ir_seg.fits.gz {1}-ir_seg.fits.gz'.format(key0, root, filt))
    os.system('cp -f {0}-mosaic_phot.fits {1}_phot.fits'.format(key0, root, filt))

    os.system('cp -f {0}-100mas-ir.cat.fits {1}-ir.cat.fits'.format(key0, root, filt))

    # Optical copy
    filts = ['f606w', 'f775w', 'f814w', 'f850lp']
    
    for filt in filts:
        sci_file = '{0}-050mas-{1}_drz_sci.fits.gz'.format(key0, filt)
        if not os.path.exists(sci_file):
            continue
        
        print(sci_file)
         
        #Symlink to individual bands
        os.system('cp -f {0}-050mas-{2}_drz_sci.fits.gz {1}-{2}_drc_sci.fits.gz'.format(key0, root, filt))
        os.system('cp -f {0}-050mas-{2}_drz_wht.fits.gz {1}-{2}_drc_wht.fits.gz'.format(key0, root, filt))
    
    
    # IR-combined 
    num = None
    for filt in ['f140w', 'f105w', 'f125w', 'f160w']:
        sci_file = '{0}-100mas-{1}_drz_sci.fits.gz'.format(key0, filt)
        if not os.path.exists(sci_file):
            continue
        
        print(sci_file)
         
        #Symlink to individual bands
        os.system('cp -f {0}-100mas-{2}_drz_sci.fits.gz {1}-{2}_drz_sci.fits.gz'.format(key0, root, filt))
        os.system('cp -f {0}-100mas-{2}_drz_wht.fits.gz {1}-{2}_drz_wht.fits.gz'.format(key0, root, filt))
 
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
    if key == 'cos':
        im = pyfits.open('cos-grism-j100012p0210-ir_drz_sci.fits')
        wht = pyfits.open('cos-grism-j100012p0210-ir_drz_wht.fits')[0].data
        sci = im[0].data
        ref_photflam = im[0].header['PHOTFLAM']
        
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
    
    # Reference WCS
    ref_im = pyfits.open('{0}-ir_drz_sci.fits'.format(root))
    ref_wcs =  utils.WCSFootprint(ref_im, ext=0)
    
    if bucket == 'grizli-v1':
        #bucket = 'grizli-v1'
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
    products = [g['direct']['product'] for g in groups]
    so = np.argsort(products)
    visits = []
    #for g in groups:
    for i in so:
        g = groups[i]
        if g['grism']['footprints'][0].intersects(ref_wcs.polygon):
            print(g['direct']['product'], g['grism']['product'])
            visits.append(g['direct'])
            visits.append(g['grism'])

    all_info['keep'] = False
    for v in visits:
        for f in v['files']:
            ix = all_info['FILE'] == f
            all_info['keep'][ix] = True
            
    np.save(root+'_visits.npy', [visits, groups, all_info[all_info['keep']]])

def grism_model():
    import numpy as np
    from grizli import utils, multifit
    from grizli.pipeline import auto_script
    import astropy.io.fits as pyfits
    import os
    import glob
    
    from grizli.pipeline import auto_script
    kwargs = auto_script.get_yml_parameters()
    
    key0 = key    
    if key == 'uds':
        field_root = root = 'uds-grism-j021732m0512'
    elif key == 'egs':
        field_root = root = 'egs-grism-j141956p5255'
    elif key == 'gds':
        field_root = root = 'gds-grism-j033236m2748'
    elif key == 'gdn':
        field_root = root = 'gdn-grism-j123656p6215'
    elif key == 'cos':
        field_root = root = 'cos-grism-j100012p0210'
        key0 = 'cos-cnd'
    
    bucket = 'grizli-v1'
    if key == 'cos':
        bucket = 'grizli-cosmos-v2'
        
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
            #has_grism |= ('-'+g in v['product'])
            has_grism |= v['product'].endswith('-'+g)

        if not has_grism:
            continue
        else:
            print('Add visit {0}'.format(v['product']))
        
        for g in grisms:
            #if '-'+g in v['product']:
            if v['product'].endswith('-'+g):
                grism_files[g].extend(v['files'])        
        
        # Copy file to local dir         
        for aws, f in zip(v['awspath'], v['files']):
            if os.path.exists(f):
                continue
            else:
                os.system('aws s3 cp s3://{0}/{1} .'.format(aws, f))
            
    for g in grisms:
        print(g, len(grism_files[g]), len(np.unique(grism_files[g])))
        os.chdir('../Prep')
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
        
        # Report
        try:
            auto_script.make_report(root, make_rgb=False, gzipped_links=False)
            os.system('cp -f ../Prep/*html ../Extractions/')
        except:
            pass
            
        os.chdir('../Extractions')
        
        os.system('aws s3 sync --exclude "*" --include "{0}*_grism*" --include "{0}*" ./ s3://{1}/Pipeline/{0}/Extractions/ --acl public-read'.format(root, bucket))
        os.system('aws s3 sync --exclude "*" --include "*GrismFLT*" --include "*wcs.fits" ./ s3://{1}/Pipeline/{0}/Extractions/ --acl public-read'.format(root, bucket))
    
    if False:
        # Copy gds g141 files to correct directory
        for i, file in enumerate(grism_files['g141']):
            print('# {0} {1}'.format(i, file))
            os.system('aws s3 sync --exclude "*" --include "{0}*" s3://grizli-v1/GrismMosaics/ s3://grizli-v1/Pipeline/{1}/Extractions/'.format(file.split('_fl')[0], root))
        
        # Wrong directory gdn >> gds
        old = 's3://grizli-v1/Pipeline/gdn-grism-j123656p6215/Extractions'
        new = 's3://grizli-v1/Pipeline/gds-grism-j033236m2748/Extractions'
        
        for g in grisms:
            for i, file in enumerate(grism_files[g]):
                filer = file.split('_fl')[0]
                files = glob.glob(filer+'*')
                if len(files) > 0:
                    print('# {0} {1}'.format(i, file))
                    for file_i in files:
                        print('  '+file_i)
                        os.system('mv {0} ../../gds-grism-j033236m2748/Prep/'.format(file_i))
                        os.system('aws s3 mv {0}/{2} {1}/{2} --acl public-read'.format(old, new, file_i))
        
        # Check missing
        os.system('aws s3 ls s3://grizli-v1/Pipeline/{0}/Extractions/ |grep GrismFLT.fits | sed "s/.01.GrismFLT.fits//" | awk \'{{{{print $4}}}}\' > {1}.files'.format(root, key))
        remote = np.loadtxt('{0}.files'.format(key), dtype=str)
        all_files = []
        for g in grism_files:
            print('---'+g)
            missing = []
            for f in grism_files[g]:
                if f.split('_fl')[0] not in remote:
                    missing.append(f)
            
            print('    Missing: '+' '.join(missing))
            all_files.extend(grism_files[g])
            
        extra = []
        for f in remote:
            if f+'_flt.fits' not in all_files:
                print('Extra: '+f)
                
    ################
      
    # Extractions
    MW_EBV = {'gdn':0.011, 'gds':0.007, 'egs': 0.007, 'uds':0.019, 'cos': 0.0148}
    
    pline = auto_script.DITHERED_PLINE
    auto_script.generate_fit_params(field_root=root, prior=None, MW_EBV=MW_EBV[key], pline=pline, fit_only_beams=True, run_fit=True, poly_order=7, fsps=True, min_sens=0.001, sys_err=0.03, fcontam=0.2, zr=[0.05, 3.4], save_file='fit_args.npy', fit_trace_shift=False, include_photometry=True, use_phot_obj=False)
    
    os.system('aws s3 sync --exclude "*" --include "{0}*ir.cat.fits*" --include "{0}*phot.fits*" --include "fit_args.npy*" ./ s3://{1}/Pipeline/{0}/Extractions/ --acl public-read'.format(root, bucket))
    
    #####
    
    ids = [27115]
    auto_script.extract(field_root=root, maglim=[13, 24], prior=None, MW_EBV=0.019, ids=ids, pline=pline, fit_only_beams=True, run_fit=False, poly_order=7, oned_R=30, master_files=None, grp=grp, bad_pa_threshold=None, fit_trace_shift=False, size=32, diff=True, min_sens=0.01, fcontam=0.2, min_mask=0.01, sys_err=0.03, skip_complete=True, args_file='fit_args.npy', get_only_beams=False)
    fitting.run_all_parallel(ids[0], verbose=True) 
    
def setup_figs():
    import numpy as np
    import os
    
    root = 'gdn-grism-j123656p6215'
    #root = 'gds-grism-j033236m2748'
    
    new_root = root.replace('-grism', '-figs')
    
    for ext in ['-ir.cat.fits', '_phot.fits', '_visits.npy', '-ir_seg.fits.gz']:
        os.system('aws s3 cp s3://grizli-v1/Pipeline/{0}/Extractions/{0}{2} s3://grizli-v1/Pipeline/{1}/Extractions/{1}{2}'.format(root, new_root, ext))

    for ext in ['-ir.cat.fits', '_phot.fits', '_visits.npy', '-ir_seg.fits', '-ir_seg.fits.gz']:
        os.system('aws s3 cp s3://grizli-v1/Pipeline/{0}/Prep/{0}{2} s3://grizli-v1/Pipeline/{1}/Prep/{1}{2}'.format(root, new_root, ext))
    
    # Fit args
    os.system('aws s3 cp s3://grizli-v1/Pipeline/{0}/Extractions/fit_args.npy .'.format(root))
    args = np.load('fit_args.npy')[0]
    args['root'] = args['group_name'] = new_root
    
    args['use_phot_obj'] = False
    args['scale_photometry'] = False
    args['fit_trace_shift'] = True
    np.save('fit_args.npy', [args])

    os.system('aws s3 cp fit_args.npy s3://grizli-v1/Pipeline/{0}/Extractions/fit_args.npy --acl public-read'.format(new_root))

    os.system('aws s3 sync s3://grizli-v1/Pipeline/{0}/Extractions/FIGS/ s3://grizli-v1/Pipeline/{1}/Extractions/ --acl public-read'.format(root, new_root))

    os.system('aws s3 rm --recursive s3://grizli-v1/Pipeline/{0}/Extractions/FIGS/'.format(root))
    
    
def grism_regions():
    from grizli import utils
    import glob
    
    wfc3 = True
    
    if wfc3:
        reg_file = 'candels_grism_wfc3.reg'
        fp = open(reg_file,'w')
        files = glob.glob('i*wcs.fits')
    else:
        reg_file = 'candels_grism_acs.reg'
        fp = open(reg_file,'w')
        files = glob.glob('j*wcs.fits')

    files.sort()

    fp.write('fk5\n')
    
    foots = OrderedDict()
    
    for file in files:
        print(file)
        foot = utils.WCSFootprint(file, ext=0)
        fp.write(foot.region+'\n')
        foots[file] = foot
        
    fp.close()
    np.save(reg_file.replace('.reg', '.npy'), [foots])
    
    ## Photometric objects in grism exposures
    cat_file = 'gdn-mosaic_phot_apcorr.fits'
    cat_file = 'gds-mosaic_phot_apcorr.fits'
    #cat_file = 'uds-mosaic_phot_apcorr.fits'
    cat_file = 'egs-mosaic_phot_apcorr.fits'

    phot = utils.read_catalog(cat_file)
    for band in ['acs', 'wfc3']:
        print(band)
        phot['has_{0}_grism'.format(band)] = 0
        coords = np.array([phot['ra'], phot['dec']]).T
        foots = np.load('candels_grism_{0}.npy'.format(band))[0]
        for ik, k in enumerate(foots):
            print(ik)
            phot['has_{0}_grism'.format(band)] += foots[k].path.contains_points(coords)
    
    has_grism = (phot['has_acs_grism'] > 1)*1 + (phot['has_wfc3_grism'] > 1)*2
    phot['has_grism'] = has_grism
    
    phot.write(cat_file.replace('.fits', '.grism.fits'), overwrite=True)
    
    # Make plot
    from matplotlib.colors import LogNorm
    from mastquery import overlaps
    
    hmag = 23.9-2.5*np.log10(phot['f160w_tot_2'])
    hmag = phot['mag_auto']
    
    sel = (hmag < 24)
    
    fig = plt.figure(figsize=[10,5])
    
    for i, col in enumerate(['has_acs_grism', 'has_wfc3_grism']):
        ax = fig.add_subplot(121+i)
        gg = phot[col][sel] > 0
        
        ax.scatter(phot['ra'][sel][~gg], phot['dec'][sel][~gg], alpha=0.1, marker='+', color='k')
        ax.text(0.95, 0.95, col.split('_')[1].upper(), ha='right', va='top', transform=ax.transAxes)
        
        ax.scatter(phot['ra'][sel][gg], phot['dec'][sel][gg], c=phot[col][sel][gg], alpha=0.1, marker='s', s=20, norm=LogNorm(vmin=1, vmax=200))
        ax.set_xlim(ax.get_xlim()[::-1])
        overlaps.draw_axis_labels(ax=ax, nlabel=3, format='latex')
        ax.set_aspect(1/np.cos(phot['dec'][0]/180*np.pi))
        ax.grid()
    
    fig.axes[1].set_yticklabels([])
    fig.axes[0].text(0.05, 0.05, cat_file.split('-')[0].upper(), ha='left', va='bottom', transform=fig.axes[0].transAxes, fontsize=16)
    
    fig.tight_layout(w_pad=1, pad=0.1)
    fig.tight_layout(w_pad=1, pad=0.1)
    fig.savefig(cat_file.replace('.fits', '.grism.png'))
    
def finkelstein():
    
    from grizli.aws import lambda_handler, fit_redshift_lambda
    ids=[37344] # z=7.5
    
    ids=[37344,37308]
    
    for root in ['gdn-grism-j123656p6215', 'gdn-figs-j123656p6215', 'gdn-g800l-j123656p6215'][1:]:
        events = fit_redshift_lambda.fit_lambda(root=root, ids=ids, bucket_name='grizli-v1', show_event=2, zr=[6,8.5], skip_started=False, scale_photometry=False, use_phot_obj=False, verbose=True, run_fit='False', clean=False)
        for event in events:
            lambda_handler.redshift_handler(event, {})

def highz_matches():
         
    ### Matches in high-z catalogs
    import numpy as np
    import numpy as np
    key = 'gdn'
    
    if key == 'uds':
        field_root = root = 'uds-grism-j021732m0512'
        ref_filt = 'f814w'
    elif key == 'egs':
        field_root = root = 'egs-grism-j141956p5255'
        ref_filt = 'f814w'
    elif key == 'gds':
        field_root = root = 'gds-grism-j033236m2748'
        ref_filt = 'f850lp'
    elif key == 'gdn':
        field_root = root = 'gdn-grism-j123656p6215'
    
    phot = utils.read_catalog(key+'-mosaic_phot_apcorr.grism.fits')
    
    newcols = []
    hasmat = phot['id']*0
    
    for f in [('finkelstein_2015.fits','f15',['zphot','_1500Mag']),
              ('bouwens2015.fits', 'b15', ['zset', 'zphot', 'F160W']),
              ('harikane_2016.fits', 'h16', ['ID'])]:
                  
        ext = utils.read_catalog('3DHST_Catalogs/'+f[0])
        #idx, dr = phot.match_to_catalog_sky(ext, other_radec=['_RAJ2000', '_DEJ2000']) 

        idx, dr = ext.match_to_catalog_sky(phot, ['_RAJ2000', '_DEJ2000']) 
        mat = dr.value < 0.4
        hasmat += mat*1
        
        phot['dr_'+f[1]] = dr.value
        phot['dr_'+f[1]].format = '.2f'
        
        newcols += ['dr_'+f[1]]
        for c in f[2]:
            print(f[1], c)
            newcol = '{0}_{1}'.format(f[1], c)
            newcols.append(newcol)
            
            phot[newcol] = np.ma.array(ext[c][idx], mask=~mat).filled()
    
    phot['ra'].format = phot['dec'].format = '.5f'
    phot['mag_auto'].format = phot['flux_radius'].format = '.2f'
    
    full_cols = ['id', 'ra', 'dec', 'mag_auto', 'flux_radius', 'has_acs_grism', 'has_wfc3_grism']+newcols
    filter_columns = full_cols[3:]
    
    # PNGs
    png_ext = 'stack'
    for g in ['g800l', 'figs', 'grism']:
        if g == 'figs':
            stack_list = ['stack', 'R30']
        else:
            stack_list = ['stack', 'full']
        
        for png_ext in stack_list:
            png = ['https://s3.amazonaws.com/grizli-v1/Pipeline/{0}/Extractions/{0}_{1:05d}.{2}.png'.format(field_root.replace('grism', g), id, png_ext) for id in phot['id']]
            png_col = '{0}_{1}'.format(g, png_ext)
            phot[png_col] = ['<a href="{0}"><img src="{0}" height=220px></a>'.format(p) for p in png]
            full_cols += [png_col]
        
    phot[full_cols][hasmat > 0].write_sortable_html('{0}-highz.html'.format(key),
                 replace_braces=True, localhost=False, 
                 max_lines=len(phot)+10, table_id=None, 
                 table_class='display compact', css=None, 
                 filter_columns=filter_columns, use_json=True)
    