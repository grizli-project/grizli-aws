"""
Check drizzled images from prep and log checks to a file
"""

def run_check(check_root='GrizliPrep', query='./j*/Prep/*phot.fits'):
    
    import glob
    import os
    import time
    
    qfiles = glob.glob(query)
    
    try:
        checked = [line.strip() for line in open(check_root+'Passed.log').readlines()]
    except:
        checked = []
    
    try:
        failed = [line.strip() for line in open(check_root+'Failed.log').readlines()]
    except:
        failed = []
        
    for file in qfiles:
        #root=os.path.basename(file)[:14]
        root=os.path.basename(file).split('_phot')[0]#[:14]
        
        if root in checked:
            print('Already checked: {0}'.format(root))
            continue
            
        path = os.path.dirname(file)
        os.system('ds9 {0}/{1}*sci.fits {0}/../Extractions/{1}*grism*.fits & '.format(path, root))
        
        time.sleep(2)
        
        os.system('scale_threedhst; ds9_match wcs')
        os.system('open {0}/*fine.png'.format(path))
        
        x = input('{0}, OK? [y/n] '.format(root))
        
        if x == 'n':
            if root not in failed:
                failed.append(root)
                
                fp = open(check_root+'Failed.log','w')
                fp.writelines([item+'\n' for item in failed])
                fp.close()
        else:
            checked.append(root)
            
            fp = open(check_root+'Passed.log','w')
            fp.writelines([item+'\n' for item in checked])
            fp.close()

def sync_new():
    
    checked = [line.strip() for line in open(check_root+'Passed.log').readlines()]
    
    pass
    