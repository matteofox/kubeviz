import os
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from pykubeviz.core.state import KubevizState

def load_datacube(state: KubevizState, filepath: str) -> bool:
    """
    Loads a FITS datacube and populates the KubevizState.
    Equivalent to the IDL fxread_cube and related setup logic.
    """
    if not os.path.exists(filepath):
        print(f"Error: File {filepath} not found.")
        return False
        
    state.filename = os.path.basename(filepath)
    state.indir = os.path.dirname(filepath) or "."
    
    try:
        with fits.open(filepath) as hdul:
            # Typically, primary HDU has no data or a 2D image, and ext 1 has the cube
            # We'll try to find the first 3D datacube
            cube_hdu = None
            pri_hdu = hdul[0]
            
            for hdu in hdul:
                if hdu.is_image and hdu.data is not None and hdu.data.ndim == 3:
                    cube_hdu = hdu
                    break
                    
            if cube_hdu is None:
                print("Error: No 3D datacube found in the FITS file.")
                return False
                
            state.indatacube = cube_hdu.data
            state.inprihead = pri_hdu.header
            state.indatahead = cube_hdu.header
            
            # WCS parsing (simplistic for now)
            header = cube_hdu.header
            state.nwpix, state.nrow, state.ncol = state.indatacube.shape
            
            # Wavelength solution from WCS standard CDELT3/CRVAL3/CRPIX3
            state.dlambda = header.get('CDELT3', header.get('CD3_3', 1.0))
            state.lambda0 = header.get('CRVAL3', 0.0)
            state.pix0 = header.get('CRPIX3', 1.0)
            
            # Create wavelength array
            pixels = np.arange(state.nwpix)
            state.wave = state.lambda0 + (pixels + 1 - state.pix0) * state.dlambda
            
            # WCS Object
            try:
                state.wcs = WCS(header)
            except Exception:
                state.wcs = None
            
            # --- Instrument Specific Parsing ---
            prihead = pri_hdu.header
            instrume = str(prihead.get('INSTRUME', '')).strip().lower()
            state.instr = instrume
            
            noise_hdu = None
            
            if state.instr == 'muse':
                print("[KUBEVIZ] Detected instrument: MUSE")
                state.noiseisvar = True
                state.band = str(prihead.get('INS MODE', prihead.get('HIERARCH ESO INS MODE', ''))).strip()
                state.fluxfac = 1.0
                
                # Instrumental resolution polynomial (fallback to manual)
                state.instrres_extpoly = np.array([499.446, -0.201608, 0.000130290, -7.87479e-09, 0.0])
                
                # Look for STAT extension
                for hdu in hdul:
                    if hdu.name == 'STAT' and hdu.data is not None and hdu.data.ndim == 3:
                        noise_hdu = hdu
                        break
                        
            elif state.instr in ['kmos', 'kmos3d']:
                print(f"[KUBEVIZ] Detected instrument: {state.instr.upper()}")
                state.band = str(prihead.get('INS FILT1 NAME', prihead.get('HIERARCH ESO INS FILT1 NAME', ''))).strip()
                state.fluxfac = 0.1
                
                # Try to detect IFU
                for j in range(1, 25):
                    key = f'OCS ARM{j} NAME'
                    hier_key = f'HIERARCH ESO OCS ARM{j} NAME'
                    val = prihead.get(key, prihead.get(hier_key, ''))
                    if str(val).strip() > '0':
                        state.ifu = j
                        break
                        
                # KMOS usually has noise in ext 2 if datacube is ext 1
                if len(hdul) > 2 and hdul[2].is_image and hdul[2].data is not None and hdul[2].data.ndim == 3:
                    noise_hdu = hdul[2]

            # Initialize smoothed cubes
            state.datacube = state.indatacube.copy() * state.fluxfac
            
            if noise_hdu is not None:
                state.innoisehead = noise_hdu.header
                state.innoisecube = noise_hdu.data
                
                if state.noiseisvar:
                    # Avoid sqrt of negative variance
                    var = np.clip(state.innoisecube, 0, None)
                    state.noisecube = np.sqrt(var) * state.fluxfac
                else:
                    state.noisecube = state.innoisecube.copy() * state.fluxfac
            else:
                state.noisecube = None
            
            print(f"Loaded {state.filename}: Shape {state.datacube.shape}")
            return True
            
    except Exception as e:
        print(f"Error loading datacube: {e}")
        return False

def save_fit_results(filepath: str, rescube: np.ndarray, z_base: float, header: fits.Header = None):
    """
    Saves the fit results as a 3D FITS datacube.
    rescube is expected to have shape (n_params, ny, nx).
    The layers correspond to:
    0: Flux_Ha
    1: dFlux_Ha
    2: Flux_NII6583
    3: dFlux_NII6583
    4: Velocity
    5: dVelocity
    6: Sigma
    7: dSigma
    """
    if header is None:
        header = fits.Header()
        
    from pykubeviz.core.lines import LINES_DB
        
    header['DESC'] = 'KubeViz Fit Results'
    header['ZBASE'] = (z_base, 'Base redshift used for fitting')
    
    for i, line in enumerate(LINES_DB):
        header[f'LAYER{i*2}'] = f"Flux_{line['name']}"
        header[f'LAYER{i*2+1}'] = f"dFlux_{line['name']}"
        
    header['LAYER34'] = 'Velocity_kms'
    header['LAYER35'] = 'dVelocity_kms'
    header['LAYER36'] = 'Sigma_kms'
    header['LAYER37'] = 'dSigma_kms'
    
    for i in range(1, 10):
        header[f'LAYER{38+(i-1)}'] = f'Continuum_Set{i}'
        header[f'LAYER{47+(i-1)}'] = f'dContinuum_Set{i}'
    
    hdu = fits.PrimaryHDU(data=rescube, header=header)
    hdul = fits.HDUList([hdu])
    
    try:
        hdul.writeto(filepath, overwrite=True)
        return True
    except Exception as e:
        print(f"Error saving fit results: {e}")
        return False

def load_fit_results(filepath: str):
    """
    Loads fit results from a FITS file.
    Returns (rescube, z_base)
    """
    from astropy.io import fits
    with fits.open(filepath) as hdul:
        rescube = hdul[0].data
        z_base = hdul[0].header.get('ZBASE', 0.0)
        return rescube, z_base
