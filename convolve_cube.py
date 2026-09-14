import numpy as np
from astropy.io import fits
import sys
from scipy.interpolate import interp1d

def convolve_muse_cube(cube_path, filter_path, output_path):
    """
    Convolves a MUSE datacube with a filter transmission curve.
    
    Parameters:
    - cube_path: Path to the MUSE datacube FITS file.
    - filter_path: Path to the filter transmission file (text file with 2 columns: wavelength, transmission).
    - output_path: Path to save the resulting 2D image FITS file.
    """
    
    print(f"Reading MUSE datacube: {cube_path}")
    with fits.open(cube_path) as hdul:
        # MUSE data is typically in extension 1 (DATA)
        # We also need the header to get the wavelength grid and WCS
        data = hdul[1].data
        header = hdul[1].header
        
        # Get wavelength grid from header
        # Check standard keywords for MUSE
        try:
            crval3 = header['CRVAL3']
            cd3_3 = header.get('CD3_3', header.get('CDELT3'))
            crpix3 = header['CRPIX3']
            naxis3 = header['NAXIS3']
        except KeyError as e:
            raise ValueError(f"Missing WCS keyword in header: {e}")
            
        # Create wavelength array (1D)
        # Pixel indices are 1-based in FITS (CRPIX), so we subtract 1 for 0-based Python
        wave_cube = crval3 + (np.arange(naxis3) - (crpix3 - 1)) * cd3_3
        
        # Read filter transmission curve
        print(f"Reading filter transmission curve: {filter_path}")
        filter_data = np.loadtxt(filter_path)
        wave_filter = filter_data[:, 0]
        trans_filter = filter_data[:, 1]
        
        # Interpolate the transmission curve onto the datacube's wavelength grid
        print("Interpolating filter transmission...")
        interp_func = interp1d(wave_filter, trans_filter, bounds_error=False, fill_value=0.0)
        trans_interp = interp_func(wave_cube)
        
        # Normalize the transmission curve so that the integral is 1 (optional, depends on use case, but standard for preserving flux units)
        trans_interp /= np.trapz(trans_interp, wave_cube)
        
        # Convolve datacube with filter
        # Reshape transmission curve for broadcasting (lambda, 1, 1)
        trans_reshaped = trans_interp[:, np.newaxis, np.newaxis]
        
        print("Convolving datacube...")
        # Replace NaNs with 0 to prevent propagation if desired, but here we just multiply
        # Assuming data units are flux per wavelength unit (e.g. erg/s/cm2/A)
        # The convolution is an integral over wavelength: Integral(Flux * Trans * dlambda)
        # Since we normalized Trans such that Integral(Trans dlambda) = 1, we can just do sum over lambda of (Flux * Trans_reshaped) * dlambda
        
        dlambda = cd3_3
        image = np.nansum(data * trans_reshaped, axis=0) * dlambda
        
        # Create output header by copying WCS info from cube header, but dropping 3rd axis
        out_header = header.copy()
        
        # Remove 3rd axis keywords
        keys_to_remove = ['NAXIS3', 'CRVAL3', 'CRPIX3', 'CD3_3', 'CDELT3', 'CTYPE3', 'CUNIT3', 'CD1_3', 'CD2_3', 'CD3_1', 'CD3_2']
        for k in keys_to_remove:
            if k in out_header:
                del out_header[k]
                
        out_header['NAXIS'] = 2
        
        # Save output FITS
        print(f"Saving image to: {output_path}")
        out_hdu = fits.PrimaryHDU(image, header=out_header)
        out_hdu.writeto(output_path, overwrite=True)
        print("Done.")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python convolve_cube.py <cube_fits> <filter_txt> <output_fits>")
        print("Example: python convolve_cube.py MUSE_datacube.fits filter_V.txt output_image.fits")
    else:
        convolve_muse_cube(sys.argv[1], sys.argv[2], sys.argv[3])
