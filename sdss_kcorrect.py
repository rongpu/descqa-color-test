# This script is used to create the k-corrected SDSS magnitudes in DESCQA v1.

from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import kcorrect
from astropy.cosmology import Planck13

# k-correct to this redshift
kcorrect_z = 0.06

####### Load and preprocess SDSS catalog #######

# Load SDSS catalog in FITS format
cat = Table.read('SpecPhoto_sdss_extinction_corrected_trimmed.fits')
print(len(cat))

# Make sure all magnitudes are valid
mask = (cat['modelMag_u']>0) & (cat['modelMag_g']>0) & (cat['modelMag_r']>0) & (cat['modelMag_i']>0) & (cat['modelMag_z']>0)
cat = cat[mask]
print(len(cat))

# Redshift cut
mask = (cat['z']>0.06) & (cat['z']<0.09)
cat = cat[mask]
print(len(cat))

####### Convert the catalog to kcorrect format #######

redshifts = np.array(cat['z'])
mags = np.array([cat['modelMag_u'], cat['modelMag_g'], cat['modelMag_r'], cat['modelMag_i'], cat['modelMag_z']]).transpose()
magerrs = np.array([cat['modelMagErr_u'], cat['modelMagErr_g'], cat['modelMagErr_r'], cat['modelMagErr_i'], cat['modelMagErr_z']]).transpose()

maglo = mags - magerrs
maghi = mags + magerrs

# convert magnitudes and errors to maggies
maggies = 10**(-0.4*mags)
errmaggies = 0.5*(10**(-0.4*maglo)-10**(-0.4*maghi))

# Save catalog in kcorrect format
data = np.concatenate((redshifts[:, None], maggies, errmaggies), axis=1)
data = Table(data)
data.write('sdss_maggies_z_0.06_0.09.dat', format='ascii.no_header')

################## Run kcorrect #####################

# Load default kcorrect templates and filters
kcorrect.load_templates()
kcorrect.load_filters()

# Fit the SED templates to observed photometry
print('Fitting SED templates to photometry')
kcorrect.fit_coeffs_from_file('sdss_maggies_z_0.06_0.09.dat', outfile='sdss_output_coeffs_z_0.06_0.09.dat')

# Apply kcorrection and save results
print('Computing k-corrected photometry')
kcorrect.reconstruct_maggies_from_file('sdss_output_coeffs_z_0.06_0.09.dat', redshift=kcorrect_z, outfile='sdss_k_corrected_maggies_z_0.06_0.09_kcorrect_{:.2}.dat'.format(kcorrect_z))

####### Convert kcorrected maggies to magnitudes #######
cat1 = Table.read('sdss_k_corrected_maggies_z_0.06_0.09_kcorrect_{:.2}.dat'.format(kcorrect_z), format='ascii.no_header', names=('redshift', 'maggies_u', 'maggies_g', 'maggies_r', 'maggies_i', 'maggies_z'))

cat_kcorrect = Table()
cat_kcorrect['redshift'] = redshifts

# k-corrected apparent magnitudes
cat_kcorrect['u'] = -2.5*np.log10(cat1['maggies_u'])
cat_kcorrect['g'] = -2.5*np.log10(cat1['maggies_g'])
cat_kcorrect['r'] = -2.5*np.log10(cat1['maggies_r'])
cat_kcorrect['i'] = -2.5*np.log10(cat1['maggies_i'])
cat_kcorrect['z'] = -2.5*np.log10(cat1['maggies_z'])

# compute distance modulus with Planck 2013 cosmology;
# beware that some references use h=1 for absolute magnitude
dm = np.array(Planck13.distmod(cat_kcorrect['redshift']))
# compute absolute magnitude
cat_kcorrect['M_u'] = cat_kcorrect['u'] - dm
cat_kcorrect['M_g'] = cat_kcorrect['g'] - dm
cat_kcorrect['M_r'] = cat_kcorrect['r'] - dm
cat_kcorrect['M_i'] = cat_kcorrect['i'] - dm
cat_kcorrect['M_z'] = cat_kcorrect['z'] - dm

cat_kcorrect.write('sdss_k_corrected_magnitudes_z_0.06_0.09_kcorrect_{:.2}.fits'.format(kcorrect_z))
