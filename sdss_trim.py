# SDSS catalog for DESCQA2

from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import sys
import fitsio
from astropy.table import Table
from os.path import expanduser
home = expanduser("~")+'/'
sys.path.append(home+'git/Python/user_modules/')

mag_col = ['modelMag_u', 'modelMag_g', 'modelMag_r', 'modelMag_i', 'modelMag_z']
ext_col = ['extinction_u', 'extinction_g', 'extinction_r', 'extinction_i', 'extinction_z']

cat = Table.read('/Users/roz18/Documents/Data/DESCQA-Color-Test/SpecPhoto20161107.fit')
print(len(cat))

# # Keep objects with valid photometry
# mask = np.ones(len(cat), dtype=bool)
# for index in range(len(mag_col)):
#     col = mag_col[index]
#     mask = mask & cat[col]>0
# cat = cat[mask]
# print(len(cat))

# Select main sample galaxies from the Legacy Survey
mask = (cat['legacy_target1']&64>0)
cat = cat[mask]
print(len(cat))

# Extinction correction
mask = np.ones(len(cat), dtype=bool)
for index in range(len(mag_col)):
    col = mag_col[index]
    col1 = ext_col[index]
    mask = cat[col]>0
    cat[col][mask] = cat[col][mask] - cat[col1][mask]

cat.remove_columns(['extinction_u', 'extinction_g', 'extinction_r', 'extinction_i', 'extinction_z', 'class', 'subClass', 'flags'])

# Fill masked value
cat = cat.filled(-99)

cat.write('/Users/roz18/Documents/Data/DESCQA-Color-Test/SpecPhoto_sdss_mgs_extinction_corrected.fits')