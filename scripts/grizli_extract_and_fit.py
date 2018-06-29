#!/usr/bin/env python
def auto_extract(root='j023507-040202', maglim=[16.5,26]):
    
    from grizli import utils
    from grizli.pipeline import auto_script
    utils.set_warnings()
    
    tab = utils.GTable.gread('{0}_footprint.fits'.format(root))
    
    DITHERED_PLINE = {'kernel': 'point', 'pixfrac': 0.2, 'pixscale': 0.1, 'size': 8, 'wcs': None}
    
    auto_script.extract(field_root=root, maglim=maglim, ids=None, run_fit=False, MW_EBV=tab.meta['MW_EBV'], pline=DITHERED_PLINE)
    
    # Test
    #auto_script.extract(field_root=root, maglim=[16.5,26], ids=[403], run_fit=False, MW_EBV=tab.meta['MW_EBV'], pline=DITHERED_PLINE)
    
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
    
    if sys.argv[2] == 'init':
        
        #generate_scripts(root=root)
        
        print("""
aws s3 sync s3://aws-grivam/Pipeline/{0}/Extractions/ ./
aws s3 cp s3://aws-grivam/Pipeline/{0}_footprint.fits ./

grizli_extract_and_fit.py {0} run

mpiexec -n 16 python -m mpi4py.futures $GRIZLICODE/grizli/pipeline/run_MPI.py

grizli_extract_and_fit.py {0} summary

aws s3 sync --exclude "*" --include "{0}_*fits" ./ s3://aws-grivam/Pipeline/{0}/Extractions/
aws s3 sync --exclude "*" --include "{0}_*png" --include "*html" --acl public-read ./ s3://aws-grivam/Pipeline/{0}/Extractions/

""".format(root))

    elif sys.argv[2] == 'summary':
        auto_script.summary_catalog(field_root=root, dzbin=0.01, use_localhost=False, filter_bandpasses=None)
    else: # "run"
        if len(sys.argv) > 3:
            maglim = np.cast[float](sys.argv[3].split(','))
        else:
            maglim = [16.5, 26]
        
        print(root, maglim)
            
        auto_extract(root=root, maglim=maglim) 