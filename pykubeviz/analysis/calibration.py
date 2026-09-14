import numpy as np

def compute_instrumental_resolution(wave: np.ndarray, instrument: str, mode: str = 'polynomial', coeffs: np.ndarray = None) -> np.ndarray:
    """
    Computes the instrumental resolution as a function of wavelength.
    Matches the IDL implementation for handling spectral resolution of different instruments (e.g. KMOS, SINFONI, MUSE).
    
    Returns:
        np.ndarray: The instrumental resolution (sigma in km/s or Angstroms depending on calibration)
    """
    if mode == 'polynomial' and coeffs is not None:
        # Polynomial fit
        res = np.polyval(coeffs[::-1], wave) # Note: IDL uses c[0] + c[1]*x, np.polyval uses c[0]*x^n + ... + c[n]
        return res
    elif mode == 'constant':
        # Default constant resolution
        return np.ones_like(wave) * (coeffs[0] if coeffs is not None else 1.0)
    else:
        # Fallback
        return np.zeros_like(wave)

def vacuum_to_air(wave_vac: np.ndarray) -> np.ndarray:
    """
    Converts vacuum wavelengths to air wavelengths using the Morton (1991) formula.
    """
    s = 1e4 / wave_vac
    n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s**2) + 0.0001599740894897 / (38.92568793293 - s**2)
    return wave_vac / n

def air_to_vacuum(wave_air: np.ndarray) -> np.ndarray:
    """
    Converts air wavelengths to vacuum wavelengths using the Ciddor (1996) formula.
    """
    s = 1e4 / wave_air
    n = 1 + 0.0000834254 + 0.02406147 / (130 - s**2) + 0.00015998 / (38.9 - s**2)
    return wave_air * n
