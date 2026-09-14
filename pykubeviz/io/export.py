import os
from astropy.io import fits
import numpy as np

def export_results_to_fits(filepath: str, results_cube: np.ndarray, errors_cube: np.ndarray, header: fits.Header = None):
    """
    Exports the fitting results to a multi-extension FITS file, matching the IDL kubeviz format.
    
    Extension 0: Primary (Empty or basic info)
    Extension 1: Results Data Cube
    Extension 2: Errors Data Cube
    """
    primary_hdu = fits.PrimaryHDU(header=header)
    primary_hdu.header['CREATOR'] = 'Kubeviz Python Port'
    
    results_hdu = fits.ImageHDU(results_cube, name='RESULTS')
    errors_hdu = fits.ImageHDU(errors_cube, name='ERRORS')
    
    hdul = fits.HDUList([primary_hdu, results_hdu, errors_hdu])
    
    # Overwrite if exists
    if os.path.exists(filepath):
        os.remove(filepath)
        
    hdul.writeto(filepath)
    print(f"Results successfully exported to {filepath}")
