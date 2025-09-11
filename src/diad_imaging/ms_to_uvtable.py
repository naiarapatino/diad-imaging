import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="casaconfig")

from uvplot import export_uvtable
import sys
import casatools

# Create CASA table tool
tb = casatools.table()

# Arguments: sys.argv[1] = ms_file, sys.argv[2] = uvtable_file
ms_file = str(sys.argv[1])
uvtable_file = str(sys.argv[2])

export_uvtable(uvtable_file, tb, vis=ms_file, datacolumn='DATA')

print(f"UV table saved to {uvtable_file}")