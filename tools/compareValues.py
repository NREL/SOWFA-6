#!/usr/bin/env python
"""
Helper function for verifying that two different sets of data (e.g.,
coordinates) are (nearly) equal
"""
import sys
import numpy as np
import pandas as pd

fname1 = sys.argv[1]
fname2 = sys.argv[2]
coords1 = pd.read_csv(fname1)
coords2 = pd.read_csv(fname2)
if not np.allclose(coords1, coords2):
    sys.exit('WARNING: '+fname1+' and '+fname2+' are significantly different!')
