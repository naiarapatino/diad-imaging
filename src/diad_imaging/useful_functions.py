##
import galario
import numpy as np
import astropy.io.fits as pyfits
import astropy.io.ascii as ascii
import sys
import shutil
import pandas as pd
import uvplot
from uvplot import UVTable, arcsec
from uvplot import COLUMNS_V0  
from galario.double import sampleImage
from galario import deg, arcsec
from astropy.io import fits
from emcee import EnsembleSampler
import corner

def read_mod(mod_file):
    
#Opens DIAD model image and makes it 
#compatible with galario

    hdu    = fits.open(mod_file)[0]
    img0   = hdu.data.squeeze()
    image  = img0.astype(np.float64, copy=False)
    # image dimension is 1025, have to change it to even number 1024
    nx, ny = image.shape
    nx_even = nx - (nx % 2)
    ny_even = ny - (ny % 2)
    image   = image[:nx_even, :ny_even]

    image = np.ascontiguousarray(image)
    dxy    = np.deg2rad(abs(hdu.header['CDELT1']))
    
    
    return image, dxy

def resample_mod(mod, obs, dxy, wle, PA, inc, dRA=0.0, dDEC=0.0):

''' Takes DIAD model image and resamples it to match 
the visibilities of an observation.
INPUTS:
    mod: image file obtained from read_mod function
    obs: observed visibilities as a uvtable with .txt extension
    dxy: 
    wle: wavelength in meters
    PA: PA in degrees
    inc: inclination of the disk in degrees
    
'''

    u, v, Re, Im, w = np.require(np.loadtxt("HD163296_uvtable.txt", unpack=True), requirements='C')
    u /= wle
    v /= wle


    PA_rad = np.deg2rad(PA)
    inc_rad = np.deg2rad(inc)

    vis_mod = sampleImage(
        image, dxy, u, v,
        dRA=dRa, dDec=dDec,
        PA=PA_rad,
        check=False
    )

    uv = UVTable(uvtable=[u*wle, v*wle, Re, Im, w], wle=wle, columns=COLUMNS_V0)
    uv.deproject(inc_rad, PA_rad)

    uv_mod = UVTable(uvtable=[u*wle, v*wle, vis_mod.real, vis_mod.imag, w], wle=wle, columns=COLUMNS_V0)
    uv_mod.deproject(inc_rad, PA_rad)

    return uv, uv_mod
