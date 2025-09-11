import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="casaconfig")

# load the visibility array
import numpy as np
import sys
from casatools import table
from casatasks import split
import os

vis_mod = np.load(sys.argv[1])

os.system('rm -rf '+sys.argv[3])
# duplicate the MS structure and DATA column into a new MS
split(vis=sys.argv[2], outputvis=str(sys.argv[3]), datacolumn='data')

# modify the output MS and insert the new 'vis' NumPy array into the measurement set

tb = table()
tb.open(sys.argv[3], nomodify=False)
data0 = tb.getcol('DATA')
npol, nchan, nvis = data0.shape
vis_model = np.broadcast_to(vis_mod[np.newaxis, np.newaxis, :], (npol, nchan, nvis))
vis_model = vis_model.astype(data0.dtype, copy=False)
tb.putcol('DATA', vis_model)
tb.close()