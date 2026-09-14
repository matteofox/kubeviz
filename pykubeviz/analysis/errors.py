import numpy as np
from typing import Callable, Dict, List

def bootstrap_error_analysis(flux: np.ndarray, wave: np.ndarray, fit_function: Callable, n_iterations: int = 100) -> Dict[str, np.ndarray]:
    """
    Performs bootstrap resampling error analysis on the spectrum.
    """
    results = {
        'amplitude': [],
        'center': [],
        'sigma': []
    }
    
    n_pixels = len(flux)
    for _ in range(n_iterations):
        # Resample with replacement
        indices = np.random.randint(0, n_pixels, n_pixels)
        resampled_wave = wave[indices]
        resampled_flux = flux[indices]
        
        # Sort back by wavelength so fitting works smoothly
        sort_idx = np.argsort(resampled_wave)
        resampled_wave = resampled_wave[sort_idx]
        resampled_flux = resampled_flux[sort_idx]
        
        try:
            _, params = fit_function(resampled_wave, resampled_flux)
            results['amplitude'].append(params['amplitude'])
            results['center'].append(params['center'])
            results['sigma'].append(params['sigma'])
        except Exception:
            pass
            
    return {k: np.array(v) for k, v in results.items()}

def montecarlo_error_analysis(flux: np.ndarray, wave: np.ndarray, noise: np.ndarray, fit_function: Callable, n_iterations: int = 100) -> Dict[str, np.ndarray]:
    """
    Performs Monte Carlo error analysis by injecting Gaussian noise into the spectrum.
    Matches IDL Kubeviz MC1/MC2 concepts.
    """
    results = {
        'amplitude': [],
        'center': [],
        'sigma': []
    }
    
    for _ in range(n_iterations):
        # Inject noise based on the noise cube/array
        noisy_flux = flux + np.random.normal(0, noise)
        
        try:
            _, params = fit_function(wave, noisy_flux)
            results['amplitude'].append(params['amplitude'])
            results['center'].append(params['center'])
            results['sigma'].append(params['sigma'])
        except Exception:
            pass
            
    return {k: np.array(v) for k, v in results.items()}
