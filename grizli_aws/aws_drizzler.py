#!/bin/env python

def group_by_filter():
    """
    aws s3 sync --exclude "*" --include "cosmos_visits*" s3://grizli-preprocess/CosmosMosaic/ ./ 
    
    """
    from grizli import prep, utils
    import numpy as np
    
    master='cosmos'
    master='grizli-jan2019'
    
    tab = utils.read_catalog('{0}_visits.fits'.format(master))
    all_visits = np.load('{0}_visits.npy'.format(master))[0]
    
    # By filter
    
    # Exclude DASH
    dash = utils.column_string_operation(tab['product'], 'icxe', 'startswith')
    
    # Don't exclude DASH
    dash = utils.column_string_operation(tab['product'], 'xxxx', 'startswith')
    
    groups = {}
    fpstr = {}
    
    for filt in np.unique(tab['filter']):
        mat = (tab['filter'] == filt) & (~dash)
        groups[filt] = {'filter':filt, 'files':[], 'awspath':[], 'footprints':[]}
        fpstr[filt] = 'fk5\n'
        
        for ix in np.where(mat)[0]:
            fp = all_visits[ix]['footprint']
            if hasattr(fp, '__len__'):
                fps = fp
            else:
                fps = [fp]
            for fp in fps: 
                xy = fp.boundary.xy
                pstr = 'polygon('+','.join(['{0:.6f}'.format(i) for i in np.array([xy[0].tolist(), xy[1].tolist()]).T.flatten()])+') # text={{{0}}}\n'.format(all_visits[ix]['product'])
            
                fpstr[filt] += pstr
            
            for k in ['files', 'awspath','footprints']:
                groups[filt][k].extend(all_visits[ix][k])
        
        fp = open('{0}-pointings-{1}.reg'.format(master, filt),'w')
        fp.write(fpstr[filt])
        fp.close()
    
        print('{0} {1:>3d} {2:>4d}'.format(filt, mat.sum(), len(groups[filt]['files'])))
    
    np.save('{0}_filter_groups.npy'.format(master), [groups])
    
DEFAULT_RGB = {'output_dpi': 75, 'add_labels':False, 'output_format':'png', 'show_ir':False}

def drizzle_images(label='cp561356', ra=150.208875, dec=1.850241667, pixscale=0.06, size=40, pixfrac=0.8, theta=0, half_optical_pixscale=False, filts=['f160w','f814w', 'f140w','f125w','f105w','f110w','f098m','f850lp', 'f775w', 'f606w','f475w'], remove=True, rgb_params=DEFAULT_RGB, master='grizli-jan2019', aws_bucket='s3://grizli/CutoutProducts/'):
    """
    label='cp561356'; ra=150.208875; dec=1.850241667; size=40; filts=['f160w','f814w', 'f140w','f125w','f105w','f606w','f475w']
    
    
    """
    from grizli import prep, utils
    from grizli.pipeline import auto_script
    
    import numpy as np
    import astropy.io.fits as pyfits
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    
    import glob
    
    import copy
    import os
    
    if isinstance(ra, str):
        coo = SkyCoord('{0} {1}'.format(ra, dec), unit=(u.hour, u.deg))
        ra, dec = coo.ra.value, coo.dec.value
    
    if label is None:
        try:
            import mastquery.utils
            label = mastquery.utils.radec_to_targname(ra=ra, dec=dec, round_arcsec=(1/15, 1), targstr='j{rah}{ram}{ras}{sign}{ded}{dem}{des}')
        except:
            label = 'grizli-cutout'
            
    #master = 'cosmos'
    #master = 'grizli-jan2019'
    
    if master == 'grizli-jan2019':
        parent = 's3://grizli/MosaicTools/'
    elif master == 'cosmos':
        parent = 's3://grizli-preprocess/CosmosMosaic/'
    else:
        parent = None
        
    for ext in ['_visits.fits', '_visits.npy', '_filter_groups.npy']:
        if (not os.path.exists('{0}{1}'.format(master, ext))) & (parent is not None):
            os.system('aws s3 cp {0}{1}{2} ./'.format(parent, master, ext))
            
    tab = utils.read_catalog('{0}_visits.fits'.format(master))
    all_visits = np.load('{0}_visits.npy'.format(master))[0]
    groups = np.load('{0}_filter_groups.npy'.format(master))[0]
        
    #filts = ['f160w','f814w', 'f110w', 'f098m', 'f140w','f125w','f105w','f606w', 'f475w']
    
    has_filts = []
    
    for filt in filts:
        if filt not in groups:
            continue
        
        visits = [copy.deepcopy(groups[filt])]
        #visits[0]['reference'] = 'CarlosGG/ak03_j1000p0228/Prep/ak03_j1000p0228-f160w_drz_sci.fits'

        hdu = utils.make_wcsheader(ra=ra, dec=dec, size=size, pixscale=pixscale, get_hdu=True, theta=theta)
        visits[0]['product'] = label+'-'+filt

        h = hdu.header

        if (filt[:2] in ['f0', 'f1', 'g1']) | (not half_optical_pixscale):
            data = hdu.data  
        else:
            for k in ['NAXIS1','NAXIS2','CRPIX1','CRPIX2']:
                h[k] *= 2

            h['CRPIX1'] -= 0.5
            h['CRPIX2'] -= 0.5

            for k in ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']:
                h[k] /= 2

            data = np.zeros((h['NAXIS2'], h['NAXIS1']), dtype=np.int16)
                        
        pyfits.PrimaryHDU(header=h, data=data).writeto('ref.fits', overwrite=True, output_verify='fix')
        visits[0]['reference'] = 'ref.fits'
                    
        status = prep.drizzle_overlaps(visits, parse_visits=False, check_overlaps=True, pixfrac=pixfrac, skysub=False, final_wcs=True, final_wht_type='IVM', static=True, max_files=260, fix_wcs_system=True)
        
        if len(glob.glob('{0}-{1}*sci.fits'.format(label, filt))):
            has_filts.append(filt)
            
        if remove:
            os.system('rm *_fl*fits')
    
        # Combine split mosaics
        tile_files = glob.glob(visits[0]['product']+'-0*sci.fits')
        if len(tile_files) > 0:
            tile_files.sort()
        
            im = pyfits.open(visits[0]['reference'])
            img = np.zeros_like(im[0].data)
            wht = np.zeros_like(im[0].data)
    
            exptime = 0.
            ndrizim = 0.
            ext = 'sci'
                
            if len(tile_files) == 0:
              continue

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
        
            pyfits.writeto('{0}_drz_sci.fits'.format(visits[0]['product']), data=sci, header=h, overwrite=True)
            pyfits.writeto('{0}_drz_wht.fits'.format(visits[0]['product']), data=wht, header=h, overwrite=True)
     
    if len(has_filts) == 0:
        return []
           
    if rgb_params:
        auto_script.field_rgb(root=label, HOME_PATH=None, filters=has_filts, **rgb_params)
     
    if aws_bucket:   
        #aws_bucket = 's3://grizli-cosmos/CutoutProducts/'
        #aws_bucket = 's3://grizli/CutoutProducts/'
    
        os.system('aws s3 sync --exclude "*" --include "{0}*" ./ {1} --acl public-read'.format(label, aws_bucket))
    
        os.system("""echo "<pre>" > index.html; aws s3 ls AWSBUCKETX --human-readable | sort -k 1 -k 2 | grep -v index | awk '{printf("%s %s",$1, $2); printf(" %6s %s ", $3, $4); print "<a href="$5">"$5"</a>"}'>> index.html; aws s3 cp index.html AWSBUCKETX --acl public-read""".replace('AWSBUCKETX', aws_bucket))
    
    return has_filts
    
if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 5:        
        print('Usage: aws_drizzler.py cp561356 150.208875 1.850241667 40 ')
        print(sys.argv)
        exit()
    
    #print('xxx')    
    drizzle_images(label=sys.argv[1], ra=float(sys.argv[2]), dec=float(sys.argv[3]), size=float(sys.argv[4]))
    
def go():
    """
    grep -v "\#" QGs.txt | awk '{print "./aws_drizzler.py",$1,$2,$3,"60"}' > run.sh
    grep -v "\#" gomez.txt | awk '{print "./aws_drizzler.py",$1,$2,$3,"60"}' >> run.sh
    
    """
    pass