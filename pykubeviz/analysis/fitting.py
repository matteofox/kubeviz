import numpy as np
from lmfit import Model, Parameters
from typing import Optional, Tuple, Dict
from multiprocessing import Pool

def gaussian(x: np.ndarray, amplitude: float, center: float, sigma: float) -> np.ndarray:
    """ 1D Gaussian function. """
    return amplitude * np.exp(-((x - center) ** 2) / (2 * sigma ** 2))

def fit_emission_line(wave: np.ndarray, flux: np.ndarray, noise: Optional[np.ndarray] = None, 
                      init_amp: float = 1.0, init_center: float = 5000.0, init_sigma: float = 2.0) -> Tuple[np.ndarray, Dict]:
    """
    Fits a single Gaussian emission line to the provided spectrum using lmfit.
    Replaces a simplified IDL mpfit invocation.
    """
    gmodel = Model(gaussian)
    params = gmodel.make_params(amplitude=init_amp, center=init_center, sigma=init_sigma)
    
    # Optional constraints based on IDL parameter tying:
    params['sigma'].min = 0.0  # Sigma must be positive
    params['amplitude'].min = 0.0 # Emission line
    
    weights = None
    if noise is not None:
        valid = (noise > 0)
        weights = np.zeros_like(noise)
        weights[valid] = 1.0 / noise[valid]
        
    result = gmodel.fit(flux, params, x=wave, weights=weights)
    
    # Return best fit flux array and a dict of parameters
    best_fit = result.best_fit
    fit_params = {
        'amplitude': result.params['amplitude'].value,
        'center': result.params['center'].value,
        'sigma': result.params['sigma'].value,
        'amplitude_err': result.params['amplitude'].stderr,
        'center_err': result.params['center'].stderr,
        'sigma_err': result.params['sigma'].stderr,
        'chi2': result.chisqr,
        'redchi': result.redchi
    }
    
    return best_fit, fit_params

def fit_tied_emission_lines(wave: np.ndarray, flux: np.ndarray, noise: Optional[np.ndarray] = None,
                            main_center: float = 5007.0, main_amp: float = 1.0, init_sigma: float = 2.0,
                            tie_ratio: float = 2.98, secondary_offset: float = -48.0) -> Tuple[np.ndarray, Dict]:
    """
    Fits two emission lines where the second line is tied to the first line by a fixed ratio and fixed wavelength offset.
    e.g. [OIII] 5007 and 4959.
    """
    gmodel1 = Model(gaussian, prefix='g1_')
    gmodel2 = Model(gaussian, prefix='g2_')
    
    model = gmodel1 + gmodel2
    params = model.make_params()
    
    # Line 1 (Main)
    params['g1_amplitude'].set(value=main_amp, min=0.0)
    params['g1_center'].set(value=main_center)
    params['g1_sigma'].set(value=init_sigma, min=0.0)
    
    # Line 2 (Tied)
    # Amplitude is tied
    params['g2_amplitude'].set(expr=f'g1_amplitude / {tie_ratio}')
    # Center is tied (offset)
    params['g2_center'].set(expr=f'g1_center + {secondary_offset}')
    # Sigma is tied
    params['g2_sigma'].set(expr='g1_sigma')
    
    weights = None
    if noise is not None:
        valid = (noise > 0)
        weights = np.zeros_like(noise)
        weights[valid] = 1.0 / noise[valid]
        
    result = model.fit(flux, params, x=wave, weights=weights)
    
    fit_params = {
        'amplitude': result.params['g1_amplitude'].value,
        'center': result.params['g1_center'].value,
        'sigma': result.params['g1_sigma'].value,
    }
    
    return result.best_fit, fit_params

def compute_moments(wave: np.ndarray, flux: np.ndarray, window_mask: np.ndarray) -> Dict[str, float]:
    """
    Computes 0th, 1st, and 2nd moments of an emission line.
    Equivalent to the IDL moments computation.
    """
    w = wave[window_mask]
    f = flux[window_mask]
    
    if len(w) == 0 or np.sum(f) <= 0:
        return {'flux': 0.0, 'vel': 0.0, 'disp': 0.0}
        
    # 0th moment (Total Flux)
    m0 = np.sum(f)
    
    # 1st moment (Velocity/Center)
    m1 = np.sum(f * w) / m0
    
    # 2nd moment (Velocity Dispersion / Sigma)
    # Variance = sum(f * (w - m1)^2) / sum(f)
    m2_var = np.sum(f * (w - m1)**2) / m0
    m2 = np.sqrt(m2_var) if m2_var > 0 else 0.0
    
    return {
        'flux': m0,
        'vel': m1,
        'disp': m2
    }

def estimate_continuum_sidebands(wave: np.ndarray, flux: np.ndarray, center_wave: float, 
                                 minoff: float = 200.0, maxoff: float = 500.0, 
                                 minperc: float = 40.0, maxperc: float = 60.0) -> Tuple[float, float]:
    """
    Estimates continuum using the SDSS method:
    Selects side bands [center - maxoff, center - minoff] and [center + minoff, center + maxoff].
    Calculates the mean and std of fluxes between minperc and maxperc percentiles in these bands.
    """
    left_mask = (wave >= center_wave - maxoff) & (wave <= center_wave - minoff)
    right_mask = (wave >= center_wave + minoff) & (wave <= center_wave + maxoff)
    mask = left_mask | right_mask
    
    valid = np.isfinite(flux)
    band_flux = flux[mask & valid]
    
    if len(band_flux) < 2:
        return 0.0, 0.0
        
    p_low = np.percentile(band_flux, minperc)
    p_high = np.percentile(band_flux, maxperc)
    
    perc_mask = (band_flux >= p_low) & (band_flux <= p_high)
    filtered_flux = band_flux[perc_mask]
    
    if len(filtered_flux) == 0:
        return np.mean(band_flux), np.std(band_flux)
        
    return np.mean(filtered_flux), np.std(filtered_flux)


from pykubeviz.core.lines import LINES_DB
from collections import defaultdict

def _get_in_range_sets(wave_min, wave_max, z_init):
    in_range = defaultdict(list)
    for line in LINES_DB:
        obs_wave = line['rest_wave'] * (1.0 + z_init)
        if wave_min <= obs_wave <= wave_max:
            in_range[line['set']].append(line)
    return in_range

def _build_and_fit_lineset(f, w, weights, lineset_lines, z_init, cont_mode, scale_covar, 
                           ckms, instrres_extpoly, vmin, vmax, smin, smax,
                           fixed_vel=None, fixed_sigma_kms=None, fit_window=None,
                           lines_config=None):
    from lmfit.models import GaussianModel, ConstantModel
    import numpy as np
    
    ref_line = lineset_lines[0]
    ref_prefix = ref_line['name'] + '_'
    ref_rest = ref_line['rest_wave']
    lambda_expected = ref_rest * (1 + z_init)
    
    mask = np.ones_like(w, dtype=bool)
    if fit_window is not None:
        mask = (w >= lambda_expected - fit_window) & (w <= lambda_expected + fit_window)
        if not np.any(mask):
            return None, ref_prefix, lambda_expected, mask
        w = w[mask]
        f = f[mask]
        if weights is not None:
            weights = weights[mask]
            
    models = []
    for line in lineset_lines:
        models.append(GaussianModel(prefix=line['name'] + '_'))
        
    model = models[0]
    for m in models[1:]:
        model = model + m
        
    if cont_mode == 1:
        model = model + ConstantModel(prefix='cont_')
        
    params = model.make_params()
    
    if cont_mode == 1:
        params['cont_c'].set(value=np.nanmedian(f))
        
    if fixed_vel is not None:
        params[ref_prefix + 'center'].set(value=lambda_expected * (1 + fixed_vel/ckms), vary=False)
    else:
        params[ref_prefix + 'center'].set(value=lambda_expected)
        if vmin is not None: params[ref_prefix + 'center'].set(min=lambda_expected * (1 + vmin/ckms))
        if vmax is not None: params[ref_prefix + 'center'].set(max=lambda_expected * (1 + vmax/ckms))
        
    if fixed_sigma_kms is not None:
        params.add(ref_prefix + 'sigma_kms', value=fixed_sigma_kms, vary=False)
    else:
        s_init = 2.0 / lambda_expected * ckms
        params.add(ref_prefix + 'sigma_kms', value=np.nanmax([s_init, 10.0]), min=smin if smin is not None else 0.0)
        if smax is not None:
            params[ref_prefix + 'sigma_kms'].set(max=smax)
            
    params.add(ref_prefix + 'sigma', expr=f"{ref_prefix}sigma_kms * {ref_prefix}center / {ckms}")
    
    fmin, fmax = 0.0, None
    if lines_config:
        lc = lines_config.get(ref_line['name'], {})
        if lc.get('min') is not None: fmin = lc['min']
        if lc.get('max') is not None: fmax = lc['max']
        
    params[ref_prefix + 'amplitude'].set(value=np.nanmax(f)*5.0, min=fmin, max=fmax)
    
    for line in lineset_lines[1:]:
        prefix = line['name'] + '_'
        l_exp = line['rest_wave'] * (1 + z_init)
        
        params[prefix + 'center'].set(expr=f"{ref_prefix}center * {line['rest_wave']} / {ref_rest}")
        params.add(prefix + 'sigma', expr=f"{ref_prefix}sigma_kms * {prefix}center / {ckms}")
            
        lfmin, lfmax = 0.0, None
        if lines_config:
            lc = lines_config.get(line['name'], {})
            if lc.get('min') is not None: lfmin = lc['min']
            if lc.get('max') is not None: lfmax = lc['max']
            
        params[prefix + 'amplitude'].set(value=np.nanmax(f)*1.0, min=lfmin, max=lfmax)
        
    line_names = [l['name'] for l in lineset_lines]
    if 'n2_b' in line_names and 'n2_r' in line_names:
        params['n2_b_amplitude'].set(expr='n2_r_amplitude / 3.071')
    if 'o3_b' in line_names and 'o3_r' in line_names:
        params['o3_b_amplitude'].set(expr='o3_r_amplitude / 2.98')
    if 'o1_b' in line_names and 'o1_r' in line_names:
        params['o1_r_amplitude'].set(expr='o1_b_amplitude / 3.0')
    try:
        result = model.fit(f, params, x=w, weights=weights, scale_covar=scale_covar)
        return result, ref_prefix, lambda_expected, mask
    except Exception:
        return None, None, None, mask

def fit_all_lines(wave: np.ndarray, flux: np.ndarray, noise: Optional[np.ndarray] = None,
                  z_init: Optional[float] = None, init_sigma: float = 2.0,
                  init_amp_ha: Optional[float] = None,
                  cont_mode: int = 0, cont_minoff: float = 200.0, cont_maxoff: float = 500.0,
                  cont_minperc: float = 40.0, cont_maxperc: float = 60.0,
                  lines_config: Optional[dict] = None, instrres_extpoly: Optional[np.ndarray] = None,
                  error_method: str = 'none', fit_window: Optional[float] = None) -> Tuple[np.ndarray, dict]:
    import numpy as np
    valid = np.isfinite(flux)
    if not np.any(valid):
        return np.full_like(flux, np.nan), {}

    w = wave[valid]
    f = flux[valid]
    

    weights = None
    if noise is not None:
        v_n = noise[valid]
        w_valid = (v_n > 0)
        weights = np.zeros_like(v_n)
        weights[w_valid] = 1.0 / v_n[w_valid]
    
    ckms = 299792.458
    
    if lines_config:
        vmin = lines_config.get('vel', {}).get('min', -1500.0)
        vmax = lines_config.get('vel', {}).get('max', 1500.0)
        smin = lines_config.get('sigma', {}).get('min', 10.0)
        smax = lines_config.get('sigma', {}).get('max', 500.0)
    else:
        vmin, vmax, smin, smax = -1500.0, 1500.0, 10.0, 500.0

    z = z_init if z_init is not None else 0.0
    in_range_sets = _get_in_range_sets(w.min(), w.max(), z)
    
    if not in_range_sets:
        return np.full_like(flux, np.nan), {}
        
    ordered_sets = sorted(in_range_sets.keys())
    
    primary_set_idx = None
    
    best_fit_total = np.full_like(f, np.nan)
    fit_params = {}
    
    for s_idx in ordered_sets:
        lineset_lines = in_range_sets[s_idx]
        
        if lines_config:
            lineset_lines = [l for l in lineset_lines if lines_config.get(l['name'], {}).get('enabled', True)]
            
        if not lineset_lines:
            continue
            
        fixed_vel = None
        fixed_sigma_kms = None
        if primary_set_idx is not None:
            fixed_vel = fit_params.get('vel')
            fixed_sigma_kms = fit_params.get('sigma_kms')
            
        lam_exp_main = lineset_lines[0]['rest_wave'] * (1 + z)
        
        if cont_mode == 0:
            cont_val, _ = estimate_continuum_sidebands(
                wave, flux, lam_exp_main,
                minoff=cont_minoff, maxoff=cont_maxoff,
                minperc=cont_minperc, maxperc=cont_maxperc
            )
            f_local = f - cont_val
        else:
            cont_val = 0.0
            f_local = f
            
        res, ref_prefix, lam_exp, mask = _build_and_fit_lineset(
            f_local, w, weights, lineset_lines, z, cont_mode, error_method == 'noise_scaling',
            ckms, instrres_extpoly, vmin, vmax, smin, smax,
            fixed_vel=fixed_vel, fixed_sigma_kms=fixed_sigma_kms,
            fit_window=fit_window, lines_config=lines_config
        )
        
        if res is not None:
            # Initialize with 0 in the masked region if it was nan, so we can overwrite or add
            if np.isnan(best_fit_total[mask]).all():
                best_fit_total[mask] = 0.0
                
            if cont_mode == 0:
                # Add continuum back to the fitted lines
                best_fit_total[mask] = res.best_fit + cont_val
            else:
                # ConstantModel is already included in res.best_fit
                best_fit_total[mask] = res.best_fit
                c_val = res.params['cont_c'].value
                from lmfit.models import GaussianModel
                for line in lineset_lines:
                    prefix = line['name'] + '_'
                    best_fit_total[mask] += GaussianModel().eval(x=w[mask], amplitude=res.params[prefix+'amplitude'].value, 
                                                           center=res.params[prefix+'center'].value, 
                                                           sigma=res.params[prefix+'sigma'].value)
            
            if cont_mode == 1:
                c_val = res.params['cont_c'].value
                c_err = res.params['cont_c'].stderr or 0.0
            else:
                c_val = cont_val
                c_err = 0.0
                
            fit_params[f'continuum_set{s_idx}'] = c_val
            fit_params[f'continuum_set{s_idx}_err'] = c_err
            
            if primary_set_idx is None:
                primary_set_idx = s_idx
                primary_result = res
                vel = (res.params[ref_prefix + 'center'].value / lam_exp - 1.0) * ckms
                vel_err = ((res.params[ref_prefix + 'center'].stderr or 0.0) / lam_exp) * ckms
                
                sigma_kms_total = res.params[ref_prefix + 'sigma_kms'].value
                sigma_kms_err = res.params[ref_prefix + 'sigma_kms'].stderr or 0.0
                
                if instrres_extpoly is not None:
                    R = np.polyval(instrres_extpoly[::-1], res.params[ref_prefix + 'center'].value)
                    sigma_instr_kms = ckms / (R * 2.35482)
                    
                    val2 = sigma_kms_total**2 - sigma_instr_kms**2
                    if val2 < 0:
                        sigma_kms = 0.0
                    else:
                        sigma_kms = np.sqrt(val2)
                else:
                    sigma_kms = sigma_kms_total
                
                fit_params['vel'] = vel
                fit_params['vel_err'] = vel_err
                fit_params['sigma_kms'] = sigma_kms
                fit_params['sigma_kms_err'] = sigma_kms_err
                
            for line in lineset_lines:
                prefix = line['name'] + '_'
                amp = res.params[prefix + 'amplitude'].value
                damp = res.params[prefix + 'amplitude'].stderr or 0.0
                
                line_flux = amp
                line_flux_err = damp
                    
                fit_params[f"flux_{line['name']}"] = line_flux
                fit_params[f"flux_{line['name']}_err"] = line_flux_err
                    
            if cont_mode == 1 and s_idx == primary_set_idx:
                best_fit_total += c_val
                
    best_fit = np.full_like(flux, np.nan)
    best_fit[valid] = best_fit_total
    
    return best_fit, fit_params

def reconstruct_fit_spectrum(wave: np.ndarray, rescube_params: np.ndarray, z_base: float, fit_window: float = 80.0, instrres_extpoly: Optional[np.ndarray] = None) -> np.ndarray:
    import numpy as np
    from lmfit.models import GaussianModel
    from pykubeviz.core.lines import LINES_DB
    
    if not np.isfinite(rescube_params[0]) and not np.isfinite(rescube_params[34]):
        return np.full_like(wave, np.nan)
        
    vel = rescube_params[34]
    sigma_intrinsic_kms = rescube_params[36]
    
    if np.isnan(vel) or np.isnan(sigma_intrinsic_kms):
        return np.full_like(wave, np.nan)
        
    ckms = 299792.458
    tot_spec = np.full_like(wave, np.nan)
    
    # Group lines by set
    from collections import defaultdict
    sets = defaultdict(list)
    for i, line in enumerate(LINES_DB):
        sets[line['set']].append((i, line))
        
    for s_idx, lines in sets.items():
        cont = rescube_params[38 + (s_idx - 1)]
        if not np.isfinite(cont):
            cont = 0.0
            
        # Find expected lambda of first line in set to define window
        ref_rest = lines[0][1]['rest_wave']
        lambda_expected = ref_rest * (1.0 + z_base) * (1.0 + vel/ckms)
        
        mask = np.zeros_like(wave, dtype=bool)
        if fit_window is not None:
            # We use a broad mask covering all lines in the set +/- fit_window
            for i, line in lines:
                le = line['rest_wave'] * (1.0 + z_base) * (1.0 + vel/ckms)
                mask |= (wave >= le - fit_window) & (wave <= le + fit_window)
        else:
            mask = np.ones_like(wave, dtype=bool)
            
        if np.isnan(tot_spec[mask]).all():
            tot_spec[mask] = cont
        else:
            tot_spec[mask] = np.where(np.isnan(tot_spec[mask]), cont, tot_spec[mask])
            
        for i, line in lines:
            flux = rescube_params[i*2]
            if np.isnan(flux) or flux <= 0: continue
            
            rest = line['rest_wave']
            center = rest * (1.0 + z_base) * (1.0 + vel/ckms)
            
            if instrres_extpoly is not None:
                R = np.polyval(instrres_extpoly[::-1], center)
                sigma_instr_kms = ckms / (R * 2.35482)
                sigma_total_kms = np.sqrt(sigma_intrinsic_kms**2 + sigma_instr_kms**2)
            else:
                sigma_total_kms = sigma_intrinsic_kms
                
            sigma = sigma_total_kms * center / ckms
            
            amp = flux
            tot_spec[mask] += GaussianModel().eval(x=wave[mask], amplitude=amp, center=center, sigma=sigma)
            
    return tot_spec

def _fit_single_spaxel(args):
    x, y, wave, flux, noise, kwargs = args
    _, params = fit_all_lines(wave, flux, noise, **kwargs)
    return x, y, params

def fit_map(datacube: np.ndarray, wave: np.ndarray, noisecube: Optional[np.ndarray] = None, 
            z_init: Optional[float] = None, n_jobs: int = 4, progress_callback=None,
            cont_mode: int = 0, cont_minoff: float = 200.0, cont_maxoff: float = 500.0,
            cont_minperc: float = 40.0, cont_maxperc: float = 60.0,
            lines_config: Optional[dict] = None, instrres_extpoly: Optional[np.ndarray] = None,
            error_method: str = 'none', fit_window: Optional[float] = None) -> np.ndarray:
    import numpy as np
    nz, ny, nx = datacube.shape
    
    n_params = 56
    rescube = np.full((n_params, ny, nx), np.nan)
    
    kwargs = {
        'z_init': z_init,
        'cont_mode': cont_mode,
        'cont_minoff': cont_minoff,
        'cont_maxoff': cont_maxoff,
        'cont_minperc': cont_minperc,
        'cont_maxperc': cont_maxperc,
        'lines_config': lines_config,
        'instrres_extpoly': instrres_extpoly,
        'error_method': error_method,
        'fit_window': fit_window
    }
    
    args_list = []
    for y in range(ny):
        for x in range(nx):
            flux = datacube[:, y, x]
            if np.any(np.isfinite(flux)):
                noise = noisecube[:, y, x] if noisecube is not None else None
                args_list.append((x, y, wave, flux, noise, kwargs))
                
    if not args_list:
        return rescube
        
    total_spaxels = len(args_list)
    completed = 0
    
    from multiprocessing import Pool
    with Pool(processes=n_jobs) as pool:
        for x, y, params in pool.imap_unordered(_fit_single_spaxel, args_list):
            completed += 1
            
            if params:
                for i, line in enumerate(LINES_DB):
                    rescube[i*2, y, x] = params.get(f"flux_{line['name']}", np.nan)
                    rescube[i*2+1, y, x] = params.get(f"flux_{line['name']}_err", np.nan)
                    
                rescube[34, y, x] = params.get('vel', np.nan)
                rescube[35, y, x] = params.get('vel_err', np.nan)
                rescube[36, y, x] = params.get('sigma_kms', np.nan)
                rescube[37, y, x] = params.get('sigma_kms_err', np.nan)
                
                for set_id in range(1, 10):
                    rescube[38 + (set_id-1), y, x] = params.get(f'continuum_set{set_id}', np.nan)
                    rescube[47 + (set_id-1), y, x] = params.get(f'continuum_set{set_id}_err', np.nan)
                    
            if progress_callback:
                try:
                    progress_callback(completed, total_spaxels, rescube)
                except Exception as e:
                    if str(e) == "Cancelled":
                        break
                    else:
                        raise
                        
    return rescube

def guess_from_adjacent_spaxels(rescube: np.ndarray, x: int, y: int, sn_thresh: float, maxvelerr: float, z_base: float) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    import numpy as np
    ny, nx = rescube.shape[1], rescube.shape[2]
    dx = [-1, 0, 1, 1, 1, 0, -1, -1]
    dy = [1, 1, 1, 0, -1, -1, -1, 0]
    
    vel_vals = []
    vel_errs = []
    sig_vals = []
    sig_errs = []
    amp_vals = []
    amp_errs = []
    
    ckms = 299792.458
    sqrt_2pi = np.sqrt(2 * np.pi)
    
    for i in range(8):
        nx_i, ny_i = x + dx[i], y + dy[i]
        if 0 <= nx_i < nx and 0 <= ny_i < ny:
            if check_autoflag(rescube, nx_i, ny_i, sn_thresh, maxvelerr):
                vel = rescube[34, ny_i, nx_i]
                vel_err = rescube[35, ny_i, nx_i]
                
                if vel_err > 0 and not np.isnan(vel):
                    vel_vals.append(vel)
                    vel_errs.append(vel_err)
                
                sig = rescube[36, ny_i, nx_i]
                sig_err = rescube[37, ny_i, nx_i]
                
                if sig_err > 0 and not np.isnan(sig):
                    sig_vals.append(sig)
                    sig_errs.append(sig_err)
                    
                flux_ha = rescube[0, ny_i, nx_i]
                flux_ha_err = rescube[1, ny_i, nx_i]
                
                if not np.isnan(flux_ha) and not np.isnan(sig):
                    center_ha = 6562.819 * (1.0 + z_base) * (1.0 + vel/ckms)
                    sigma_ha = sig * center_ha / ckms
                    amp = flux_ha / (sigma_ha * sqrt_2pi) if sigma_ha > 0 else 0
                    
                    amp_err = amp * (flux_ha_err / flux_ha) if flux_ha > 0 else 0
                    
                    if amp > 0 and amp_err > 0:
                        amp_vals.append(amp)
                        amp_errs.append(amp_err)
                
    if not vel_vals:
        return None, None, None
        
    vel_vals = np.array(vel_vals)
    vel_errs = np.array(vel_errs)
    sig_vals = np.array(sig_vals)
    sig_errs = np.array(sig_errs)
    
    w_vel = 1.0 / (vel_errs ** 2)
    vel_guess = np.sum(vel_vals * w_vel) / np.sum(w_vel)
    
    w_sig = 1.0 / (sig_errs ** 2)
    sig_guess = np.sum(sig_vals * w_sig) / np.sum(w_sig)
    
    z_guess = (1.0 + z_base) * (1.0 + vel_guess / ckms) - 1.0
    
    amp_guess = None
    if amp_vals:
        amp_vals = np.array(amp_vals)
        amp_errs = np.array(amp_errs)
        w_amp = 1.0 / (amp_errs ** 2)
        amp_guess = np.sum(amp_vals * w_amp) / np.sum(w_amp)
    
    return z_guess, sig_guess, amp_guess

def fit_adj_all(datacube: np.ndarray, wave: np.ndarray, noisecube: Optional[np.ndarray], rescube: np.ndarray,
                kwargs: dict, sn_thresh: float = 3.0, maxvelerr: float = 50.0,
                progress_callback=None, max_iter: int = 100000, z_base: float = 0.0):
    import numpy as np
    from typing import Tuple, Optional
    nz, ny, nx = datacube.shape
    new_rescube = rescube.copy()
    
    hope = np.ones((ny, nx), dtype=bool)
    total_fixed = 0
    
    for it in range(max_iter):
        fixed_this_iter = 0
        bad_spaxels_by_neighbors = {i: [] for i in range(1, 9)}
        
        for y in range(ny):
            for x in range(nx):
                if hope[y, x] and not check_autoflag(new_rescube, x, y, sn_thresh, maxvelerr):
                    dx = [-1, 0, 1, 1, 1, 0, -1, -1]
                    dy = [1, 1, 1, 0, -1, -1, -1, 0]
                    ok_count = 0
                    for i in range(8):
                        nx_i, ny_i = x + dx[i], y + dy[i]
                        if 0 <= nx_i < nx and 0 <= ny_i < ny:
                            if check_autoflag(new_rescube, nx_i, ny_i, sn_thresh, maxvelerr):
                                ok_count += 1
                                
                    if ok_count > 0:
                        bad_spaxels_by_neighbors[ok_count].append((x, y))
        
        bad_spaxels = []
        for i in range(8, 0, -1):
            if len(bad_spaxels_by_neighbors[i]) > 0:
                bad_spaxels = bad_spaxels_by_neighbors[i]
                break
                
        if progress_callback:
            total_bad = sum(len(v) for v in bad_spaxels_by_neighbors.values())
            progress_callback(it, total_fixed, total_bad, new_rescube)
            
        if not bad_spaxels:
            break
            
        for x, y in bad_spaxels:
            z_guess, sig_guess, amp_guess = guess_from_adjacent_spaxels(new_rescube, x, y, sn_thresh, maxvelerr, z_base)
            if z_guess is not None:
                flux = datacube[:, y, x]
                noise = noisecube[:, y, x] if noisecube is not None else None
                
                kw = kwargs.copy()
                kw['z_init'] = z_guess
                kw['init_sigma'] = sig_guess * 6562.819 / 299792.458 if sig_guess else 2.0
                kw['init_amp_ha'] = amp_guess
                
                _, params = fit_all_lines(wave, flux, noise, **kw)
                if params:
                    for i, line in enumerate(LINES_DB):
                        new_rescube[i*2, y, x] = params.get(f"flux_{line['name']}", np.nan)
                        new_rescube[i*2+1, y, x] = params.get(f"flux_{line['name']}_err", np.nan)
                        
                    new_rescube[34, y, x] = params.get('vel', np.nan)
                    new_rescube[35, y, x] = params.get('vel_err', np.nan)
                    new_rescube[36, y, x] = params.get('sigma_kms', np.nan)
                    new_rescube[37, y, x] = params.get('sigma_kms_err', np.nan)
                    
                    for set_id in range(1, 10):
                        new_rescube[38 + (set_id-1), y, x] = params.get(f'continuum_set{set_id}', np.nan)
                        new_rescube[47 + (set_id-1), y, x] = params.get(f'continuum_set{set_id}_err', np.nan)
                        
                    if check_autoflag(new_rescube, x, y, sn_thresh, maxvelerr):
                        fixed_this_iter += 1
                        total_fixed += 1
                    else:
                        hope[y, x] = False
                else:
                    hope[y, x] = False
                    
        if fixed_this_iter == 0:
            break
            
    if progress_callback:
        progress_callback(it, total_fixed, 0, new_rescube)
        
    return new_rescube, total_fixed

def check_autoflag(rescube: np.ndarray, x: int, y: int, sn_thresh: float, maxvelerr: float) -> bool:
    import numpy as np
    flux_ha = rescube[0, y, x]
    flux_ha_err = rescube[1, y, x]
    vel_err = rescube[35, y, x]
    
    if np.isnan(flux_ha) or flux_ha <= 0: return False
    if np.isnan(flux_ha_err) or flux_ha_err <= 0: return False
    
    sn = flux_ha / flux_ha_err
    if sn < sn_thresh: return False
    
    if vel_err > maxvelerr: return False
    
    return True
