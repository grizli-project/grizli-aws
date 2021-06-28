
# pip install git+https://github.com/gbrammer/fitsmap

import os
import glob
import astropy.io.fits as pyfits
import astropy.wcs as pywcs

import numpy as np

import matplotlib.pyplot as plt
plt.ioff()

def all_fields():
    from grizli.aws import db
    engine = db.get_db_engine()
    ch = db.from_sql("select * from charge_fields where log like 'Finish%%' and field_ra > 0 and filters like '%%F1%%' and nfilt > 1 and field_root like 'j%%'", engine)
    
    so = np.argsort(ch['field_root'])
    with open('fields.txt','w') as fp:
        for r in ch['field_root'][so]:
            fp.write(f'{r}\n')
            
    bash = """
    roots=`cat fields.txt | shuf`
    for root in ${roots}; do
        if [ ! -e "${root}" ]; then
            date > ./log/${root}.log
            rm -rf ${root}/output/catalog_assets ${root}/output/js
            python make_fitsmap.py ${root}
            rm ${root}/*fits*
        fi
    done
    
    """
    
    import json
    
    js = json.load(open('CHArGE-Mar2020.json.bkup'))
    count = 0
    for row in js['data']:
        root=row[0]
        map_href = f"<a href=https://s3.amazonaws.com/grizli-v1/Pipeline/{root}/Map/index.html>{root}</a>"
        if os.path.exists(root):
            row[0] = map_href
            count += 1
            
    with open('CHArGE-Mar2020.json','w') as fp:
        json.dump(js, fp)
    
    os.system('aws s3 cp CHArGE-Mar2020.json s3://grizli-v1/Master/ --acl public-read')
    os.system('aws s3 cp CHArGE-Mar2020.json.bkup s3://grizli-v1/Master/ --acl public-read')
    
def make_seg(segfile, outfile='seg.png', xs=8, seed=1):
    """
    Make colorful segmentation image
    """
    segim = pyfits.open(segfile)[0].data
    
    ids = np.unique(segim)
    np.random.seed(seed)
    rnd = np.random.rand(segim.max()+1)
    
    rnd_seg = rnd[segim.flatten()].reshape(segim.shape)
    
    sh = segim.shape
    
    fig, ax = plt.subplots(1,1,figsize=(xs, xs*sh[0]/sh[1]))

    dpi = sh[1]/xs
    norm = None
                       
    rnd_seg[segim == 0] = np.nan
    ax.imshow(rnd_seg, cmap='terrain_r', vmin=0, vmax=1, origin='lower')
    
    ax.axis('off')
    fig.tight_layout(pad=0)
    fig.savefig(outfile, dpi=dpi)
    return fig
    
    
def run_root(root='j002532m1223', min_zoom=2, get_grism=True):
    """
    Prepare images for fitsmap.convert
    """        
    from grizli.pipeline import auto_script
    from grizli import utils
    import eazy.utils
    
    from fitsmap import convert

    print('sync')
    
    os.system(f'aws s3 sync s3://grizli-v1/Pipeline/{root}/Prep/ {root}/ --exclude "*" --include "*sci.fits.gz" --include "*phot.fits" --include "*seg.fits.gz"')
    os.system(f'aws s3 sync s3://grizli-v1/Pipeline/{root}/IRAC/ {root}/ --exclude "*" --include "*sci.fits*" --include "*model.fits"')
    os.system(f'aws s3 sync s3://grizli-v1/Pipeline/{root}/Map/ {root}/ --exclude "*" --include "{root}.*png"')
    
    os.chdir(root)
    
    if not os.path.exists(f'{root}.rgb.png'):
        _ = auto_script.field_rgb(root=root, xsize=6, full_dimensions=True, HOME_PATH=None, gzext='*', suffix='.rgb', output_format='png')
    
    # IR
    files = glob.glob(f'{root}-[if][r01]*sci.fits*')
    files.sort()
    filts = [file.split(f'{root}-')[1].split('_')[0] for file in files]
    for filt in filts:
        if os.path.exists(f'{root}.{filt}.png'):
            continue
            
        _ = auto_script.field_rgb(root=root, xsize=6, full_dimensions=True, HOME_PATH=None, gzext='*', filters=[filt], suffix=f'.{filt}', output_format='png', invert=True, scl=2)
    
    # Optical, 2X pix
    files = glob.glob(f'{root}-[f][2-8]*sci.fits*')
    files.sort()
    filts = [file.split(f'{root}-')[1].split('_')[0] for file in files]
    for filt in filts:
        if os.path.exists(f'{root}.{filt}.png'):
            continue
            
        _ = auto_script.field_rgb(root=root, xsize=6, full_dimensions=2, HOME_PATH=None, gzext='*', filters=[filt], suffix=f'.{filt}', output_format='png', invert=True, scl=2)

    # Spitzer
    if glob.glob(f'{root}-ch*fits.gz'):
        import reproject
        out_img = pyfits.open(f'{root}-ir_drz_sci.fits.gz')
        repr_hdu = out_img[0]
        # repr_hdu = utils.make_maximal_wcs([out_wcs], pixel_scale=0.2, 
        #                                   verbose=False, pad=0, poly_buffer=0)
        repr_wcs = pywcs.WCS(repr_hdu.header)
        
        mosaics = glob.glob(f'{root}-ch[12]*sci.fits.gz')
        mosaics.sort()
        for mos in mosaics:
            ch = mos.split(f'{root}-')[1].split('_')[0]
            if os.path.exists(f'{root}.{ch}.png'):
                continue
            
            print(f'Reproject {ch}')
            in_img = pyfits.open(f'{root}-{ch}_drz_sci.fits.gz')
            in_wcs = pywcs.WCS(in_img[0].header)
            
            reproj = utils.blot_nearest_exact(in_img[0].data, 
                                              in_wcs, repr_wcs, 
                                              scale_by_pixel_area=False)
                                  
            pyfits.writeto(f'{root}-{ch}_drz_sci.fits', data=reproj, 
                       header=repr_hdu.header, overwrite=True)
                    
            ext = [ch]
            
            if os.path.exists(f'{root}-{ch}_model.fits'):
                # resid
                print(f' {ch} model')
                in_img = pyfits.open(f'{root}-{ch}_model.fits')
                reproj = utils.blot_nearest_exact(in_img[1].data, 
                                                  in_wcs, repr_wcs, 
                                                  scale_by_pixel_area=False)
                pyfits.writeto(f'{root}-{ch}m_drz_sci.fits', data=reproj, 
                      header=repr_hdu.header, overwrite=True)
                ext.append(ch+'m')
                                
            for filt in ext:
                _ = auto_script.field_rgb(root=root, xsize=6, full_dimensions=True, HOME_PATH=None, gzext='', filters=[filt], suffix=f'.{filt}', output_format='png', invert=True, scl=2)
    
    if not os.path.exists(f'{root}.seg.png'):
        sfig = make_seg(f'{root}-ir_seg.fits.gz', outfile=f'{root}.seg.png')
    
    filelist = []
    for q in ['*f[2-8]', '*f[01]*', '*ir*', '*ch[12]','*seg', '*rgb']:
        l_i = glob.glob(q+'*png')
        l_i.sort()
        filelist.extend(l_i)
    
    ph = utils.read_catalog(f'{root}_phot.fits')
    ph['id'] = ph['number']
    ph['ra'].format = '.6f'
    ph['dec'].format = '.6f'
    ph['mag'] = ph['mag_auto']
    ph['mag'].format = '.2f'
    
    ph['query'] = [eazy.utils.query_html(r, d).split(') ')[1]
                   for r, d in zip(ph['ra'], ph['dec'])]
    
    ph['id','ra','dec','query','mag'].write('phot.cat', format='ascii.csv', overwrite=True)
    
    filelist += ['phot.cat']
    
    if get_grism:
        from grizli.aws import db
        engine = db.get_db_engine()
        gr = db.from_sql(f"select root, id, ra, dec, z_map from redshift_fit where root='{root}'", engine)
        
        print(f'grism.cat: {len(gr)} sources')
        
        if len(gr) > 0:
            gr['query'] = [eazy.utils.query_html(r, d).split(') ')[1]
                           for r, d in zip(gr['ra'], gr['dec'])]
                       
            gr['stack'] = [f'<img src="https://s3.amazonaws.com/grizli-v1/Pipeline/{root}/Extractions/{root}_{id:05d}.stack.png"  height="100px"/>' for id in gr['id']]
            gr['full'] = [f'<img src="https://s3.amazonaws.com/grizli-v1/Pipeline/{root}/Extractions/{root}_{id:05d}.full.png"  height="100px"/>' for id in gr['id']]
            gr['line'] = [f'<img src="https://s3.amazonaws.com/grizli-v1/Pipeline/{root}/Extractions/{root}_{id:05d}.line.png" height="80px"/>' for id in gr['id']]
            
            gr['ra'].format = '.6f'
            gr['dec'].format = '.6f'
            gr['z_map'].format = '.4f'
            
            gr['id','ra','dec','query','z_map', 'stack','full','line'].write('grism.cat', format='grism.csv', overwrite=True)

            filelist += ['grism.cat']
    
    
    convert.MPL_CMAP = 'gray_r'
    convert.cartographer.MARKER_HTML_WIDTH = '650px'
    convert.cartographer.MARKER_HTML_HEIGHT = '440px'
    convert.POPUP_CSS = [
        "span { text-decoration:underline; font-weight:bold; line-height:12pt; }",
        "tr { line-height: 7pt; }",
        "table { width: 100%; }",
        "img { height: 100px; width: auto; }",
    ]
    
    convert.dir_to_map(
        "./",
        filelist=filelist,
        out_dir="output",
        cat_wcs_fits_file=f"{root}-ir_drz_sci.fits.gz",
        catalog_delim=',',
        min_zoom=min_zoom,
        task_procs=False,
        image_engine='MPL'
    )
    plt.close('all')
    
    if os.path.exists('output/index.html'):
        os.system(f'aws s3 sync output/ s3://grizli-v1/Pipeline/{root}/Map/ --acl public-read --quiet')
        os.system(f'aws s3 sync ./ s3://grizli-v1/Pipeline/{root}/Map/ --exclude "*" --include "{root}.*png" --acl public-read')
        
        
if __name__ == '__main__':
    import sys
    import os
    pwd = os.getcwd()
    root = sys.argv[1]
    print(root)
    run_root(root=root)
    
def testing():
    
    from fitsmap import convert, cartographer
    from importlib import reload
    reload(cartographer); reload(convert)
    
    root = 'j004404m2034'
    
    norm_kwargs = dict(stretch='log', min_cut=-0.02,
                           max_cut=5)
    #
    norm_kwargs = dict(stretch='linear', min_percent=0.1,
                           max_percent=99.9)
    
    # norm_kwargs = None
                           
    os.system('rm -rf output/j004404m2034_f160w_drz_sci')
    
    convert.MPL_CMAP = 'gray_r'
    convert.dir_to_map(
        "./",
        out_dir="output",
        cat_wcs_fits_file=f"{root}-f160w_drz_sci.fits",
        catalog_delim=',',
        min_zoom=min_zoom,
        task_procs=False,
        image_engine='MPL',
        norm_kwargs=norm_kwargs
    )
    