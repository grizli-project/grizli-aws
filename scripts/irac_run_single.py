#!/usr/bin/env python
    
def auto_run(root='j023507-040202', args=[]):
    import os
    import yaml
    import time
    
    import matplotlib.pyplot as plt
    plt.ioff()
    
    from grizli import utils
    utils.set_warnings()
    
    import golfir.pipeline
    
    # Run query again
    
    kwargs = {'home':'/GrizliImaging/', 'pixfrac':0.2, 'kernel':'square', 'initial_pix':1.0, 'final_pix':0.5, 'pulldown_mag':15.2, 'sync_xbcd':True, 'skip_fetch':False, 'radec':None, 'mosaic_pad':2.5, 'drizzle_ref_file':''}
        
    master_radec = '{0}/{1}_master.radec'.format(os.getcwd(), root)
    parent_radec = '{0}/{1}_parent.radec'.format(os.getcwd(), root)
    if os.path.exists(master_radec):
        kwargs['radec'] = master_radec
    elif os.path.exists(parent_radec):
        kwargs['radec'] = parent_radec
        
    # Command line arguments
    if args:
        for arg in args:
            if arg.startswith('--'):
                if arg in ['--grism', '--sync', '--noclean', '--lambda_verbose']:
                    continue
                    
                pspl = arg.strip('--').split('=')[0]
                val = arg.split('=')[1]
                
                if pspl in ['radec']:
                    continue
                                        
                # Split nested dictionaries by '.'
                if '.' in pspl:
                    valid = False
                    ps = pspl.split('.')
                    d = kwargs
                    for p in ps:
                        if p in d:
                            valid = True
                            if isinstance(d[p], dict):
                                d = d[p]
                else:
                    d = kwargs
                    valid = pspl in kwargs
                    p = pspl
                                    
                if valid:
                    if val.strip() in ['None', 'null']:
                        d[p] = None
                        print('Runtime argument: {0} = {1}/[None]'.format(p, val))
                    elif isinstance(d[p], list):
                        lval = val.replace('[','').replace(']','').split(',')
                        
                        # Item shoud be a list
                        if (len(lval) < len(d[p])) & ('filter' not in arg):
                            msg = 'Parameter {0} should be a list like {1}'
                            raise(ValueError(msg.format(arg, d[p])))
                        
                        try:
                            lval = list(np.cast[float](lval))
                        except:
                            pass
                            
                        d[p] = lval
                        print('Runtime argument: {0} = {1}'.format(p, lval))
                    else:
                        if isinstance(d[p], bool):
                            if val.isdigit():
                                d[p] = int(val)
                            else:
                                d[p] = val.lower() == 'true'
                        else:
                            try:
                                d[p] = d[p].__class__(val)
                            except:
                                try:
                                    d[p] = float(val)
                                except:
                                    d[p] = val
                                    
                        print('Runtime argument: {0} = {1}'.format(p, d[p]))
    
    golfir.pipeline.irac_mosaics(root=root, **kwargs)                        
        
if __name__ == "__main__":
    import sys
    import numpy as np
    from grizli import utils
    utils.set_warnings()
    
    # if len(sys.argv) < 3:
    #     print('Usage: aws.py {field} {init/run/summary} [min_mag,max_mag]')
    #     exit 
        
    root = sys.argv[1]
    auto_run(root=root, args=sys.argv[2:])