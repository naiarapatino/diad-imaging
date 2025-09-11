############# Define your own Paths here!###############
#local_path_to_au    = '/projectnb/bu-disks/nzpatino/ALMA/analysis_scripts/'
local_path_to_utils = '/projectnb/bu-disks/nzpatino/ALMA/Helper-Scripts/'

#######################################################


import os
import sys
import numpy as np
from casatasks import tclean, exportfits
# sys.path.append(local_path_to_au)
# import analysisUtils as au

# utilsdir = local_path_to_utils
# execfile(utilsdir+'reduction_utils.py')

vis = sys.argv[1]
cellsize = sys.argv[2] # '0.04arcsec'
imsize = int(sys.argv[3]) # 512 - must be integer for tclean
mask_pa = float(sys.argv[4]) #-16 #position angle of mask in degrees
mask_maj = float(sys.argv[5]) #4. #semimajor axis of mask in arcsec
mask_min = float(sys.argv[6]) #1  #semiminor axis of mask in arcsec
mask_ra = str(sys.argv[7]) #'04h20m21.464s'
mask_dec = str(sys.argv[8]) #'+28.13.48.829'
LB_scales = eval(sys.argv[9]) #[0, 5, 10, 20, 30]
robust = float(sys.argv[10]) #0.5 - must be float for tclean
threshold = str(sys.argv[11]) #'0.08mJy'

print(f"Imaging {vis} with cellsize={cellsize}, imsize={imsize}, mask_pa={mask_pa}, mask_maj={mask_maj}, mask_min={mask_min}, mask_ra={mask_ra}, mask_dec={mask_dec}, LB_scales={LB_scales}, robust={robust}, threshold={threshold}")

common_mask = 'ellipse[[%s, %s], [%.1farcsec, %.1farcsec], %.1fdeg]' % (mask_ra, mask_dec, mask_maj, mask_min, mask_pa)
imname = f'{vis[:-3]}.284ghz.R{robust}'

os.system('rm -rf '+imname+'*')

tclean(vis=vis,
       imagename=imname,
       specmode='mfs',
       spw='',
       deconvolver = 'multiscale',
       scales = LB_scales,
       weighting = 'briggs',
       robust = robust,
       restoringbeam = 'common',
       cell=cellsize,
       imsize=imsize,
       smallscalebias = 0.6,
       niter = 50000,
       threshold = threshold,
       perchanweightdensity=False,
       mask = common_mask,
       interactive=False)

exportfits(f'{imname}.image', f'{imname}.fits')