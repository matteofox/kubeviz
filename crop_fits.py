import astropy.io.fits as fits

input_file = "ADP.2016-06-17T14_40_46.135.fits"
output_file = "ADP.2016-06-17T14_40_46.135_cropped.fits"

with fits.open(input_file) as hdul:
    for hdu in hdul:
        if hdu.is_image and hdu.data is not None and hdu.data.ndim == 3:
            z, y, x = hdu.data.shape
            
            # Central 100x100
            y_start = max(0, y // 2 - 50)
            y_end = min(y, y // 2 + 50)
            x_start = max(0, x // 2 - 50)
            x_end = min(x, x // 2 + 50)
            
            hdu.data = hdu.data[:, y_start:y_end, x_start:x_end]
            
            # Update header CRPIX
            if 'CRPIX1' in hdu.header:
                hdu.header['CRPIX1'] -= x_start
            if 'CRPIX2' in hdu.header:
                hdu.header['CRPIX2'] -= y_start
                
    hdul.writeto(output_file, overwrite=True)
    
print(f"Successfully cropped {input_file} and saved to {output_file}")
