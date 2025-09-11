##
import numpy as np
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
import os
from pathlib import Path

package_root = Path(__file__).parent.parent.parent
scripts_root = str(package_root / 'src' / 'diad_imaging')

def ms_to_uvtable(ms_file, uvtable_file, verbose=False):
    ''' Converts a measurement set to a uvtable text file.
    INPUTS:
    ms_file: path to the measurement set
    uvtable_file: path to the output uvtable text file
    '''

    ms_to_uvtable_script = str(scripts_root + '/ms_to_uvtable.py')

    # Clean up CASA lock file
    casa_data_lock = f"{package_root}/.pixi/envs/default/share/casa/data/data_update.lock"
    os.system(f"rm -f '{casa_data_lock}'")
    
    # Run CASA script using pixi's Python environment with CASA libraries
    # Suppress pkg_resources warning
    casa_cmd = f"cd {package_root} && PYTHONWARNINGS=ignore::UserWarning pixi run python '{ms_to_uvtable_script}' '{ms_file}' '{uvtable_file}'"


    if verbose:
        print(f"Running command: {casa_cmd}")
    os.system(casa_cmd)
    os.system(f"cd {package_root};rm -f *.log")

    return uvtable_file

def mod_to_ms(mod_file, obs_ms, output_ms, verbose=False):
    ''' Converts a resampled model image to a measurement set.
    INPUTS:
    mod_file: path to the model image file
    obs_ms: path to the observed measurement set (to copy metadata)
    output_ms: path to the output measurement set
    '''
    mod_to_ms_script = str(scripts_root + '/mod_to_ms.py')

    # Clean up CASA lock file and any existing output
    casa_data_lock = f"{package_root}/.pixi/envs/default/share/casa/data/data_update.lock"
    os.system(f"rm -f '{casa_data_lock}'")
    os.system(f"rm -rf '{output_ms}'")
    
    # Suppress pkg_resources warning
    casa_cmd = f"cd {package_root} && PYTHONWARNINGS=ignore::UserWarning pixi run python '{mod_to_ms_script}' '{mod_file}' '{obs_ms}' '{output_ms}'"
    if verbose:
        print(f"Running command: {casa_cmd}")
    os.system(casa_cmd)
    os.system(f"cd {package_root};rm -f *.log")

    # os.system(f"arch -x86_64 casa --nologger --nogui -c \"execfile('{mod_to_ms_script}')\" {mod_file} {obs_ms} {output_ms}")
    # os.system(f"rm *.log")

    return output_ms

def inputs_for_imaging(object_name, freq):
    ''' Reads the input parameters for imaging from a CSV file.
    INPUTS:
    object_name: name of the object to be imaged (should match the entry in the CSV file)
    freq: frequency of the observation
    '''

    imaging_inputs_csv = str(f'{scripts_root}/tables/imaging_inputs.csv')

    df = pd.read_csv(imaging_inputs_csv)
    row = df[df['obs'] == f'{object_name}_{freq}'].iloc[0]

    cellsize = row['cellsize']
    imsize = row['imsize']
    threshold = row['threshold']
    mask_ra = row['mask_ra']
    mask_dec = row['mask_dec']
    mask_maj = row['mask_maj']
    mask_min = row['mask_min']
    mask_pa = row['mask_pa']
    robust = row['robust']
    LB_scales = eval(row['LB_scales'])  # Convert string representation of list to actual list

    return cellsize, imsize, threshold, mask_ra, mask_dec, mask_maj, mask_min, mask_pa, robust, LB_scales

def image_mod(mod_file, cellsize, imsize, threshold, mask_ra, mask_dec, mask_maj, mask_min, mask_pa, robust=0.5, LB_scales=[0, 5, 10, 20, 30],verbose=False):

    ''' Creates a CLEANed image from a model measurement set.
    INPUTS:
    mod_file: path to the model image file
    cellsize: cell size for the image
    imsize: image size
    threshold: threshold for cleaning
    mask_ra: right ascension of the mask
    mask_dec: declination of the mask
    mask_maj: major axis of the mask
    mask_min: minor axis of the mask
    mask_pa: position angle of the mask
    robust: robust parameter for cleaning
    LB_scales: list of scales for multi-scale cleaning
    '''

    # Clean up CASA lock file
    casa_data_lock = f"{package_root}/.pixi/envs/default/share/casa/data/data_update.lock"
    os.system(f"rm -f '{casa_data_lock}'")

    image_mod_script = str(scripts_root + '/image_mod_script.py')
    casa_cmd = f"cd {package_root} && PYTHONWARNINGS=ignore::UserWarning pixi run python '{image_mod_script}' '{mod_file}' '{cellsize}' '{imsize}' '{mask_pa}' '{mask_maj}' '{mask_min}' '{mask_ra}' '{mask_dec}' '{LB_scales}' '{robust}' '{threshold}'"
    if verbose:
        print(f"Running command: {casa_cmd}")
    os.system(casa_cmd)
    os.system(f"cd {package_root};rm -f *.log")

    return


def read_mod(mod_file):

    ''' Opens DIAD model image and makes it compatible with galario '''

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

def resample_mod(savename,mod, obs, dxy, wle, PA, inc, dRA=0.0, dDEC=0.0):

    ''' Takes DIAD model image and resamples it to match 
    the visibilities of an observation.
    INPUTS:
    savename: path/name to save the resampled model image (without extension)
    mod: image file obtained from read_mod function
    obs: observed visibilities as a uvtable with .txt extension
    dxy: pixel scale in radians
    wle: wavelength in meters
    PA: PA in degrees
    inc: inclination of the disk in degrees
    
    '''

    u, v, Re, Im, w = np.require(np.loadtxt(obs, unpack=True), requirements='C')
    u /= wle
    v /= wle


    PA_rad = np.deg2rad(PA)
    inc_rad = np.deg2rad(inc)

    vis_mod = sampleImage(
        mod, dxy, u, v,
        dRA=dRA, dDec=dDEC,
        PA=PA_rad,
        check=False
    )

    uv = UVTable(uvtable=[u*wle, v*wle, Re, Im, w], wle=wle, columns=COLUMNS_V0)
    uv.deproject(inc_rad, PA_rad)

    uv_mod = UVTable(uvtable=[u*wle, v*wle, vis_mod.real, vis_mod.imag, w], wle=wle, columns=COLUMNS_V0)
    uv_mod.deproject(inc_rad, PA_rad)

    np.save(f'{savename}.npy', vis_mod)

    return uv, uv_mod
