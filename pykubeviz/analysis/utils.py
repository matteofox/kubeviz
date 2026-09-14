import numpy as np

def kubeviz_weighted_median(array: np.ndarray, weight: np.ndarray) -> float:
    """
    Computes the weighted median of a 1D array.
    Equivalent to the IDL kubeviz_weighted_median.
    """
    array = np.asarray(array)
    weight = np.asarray(weight)
    
    # Sort array and weights
    ind = np.argsort(array)
    array_sort = array[ind]
    weight_sort = weight[ind]
    
    tot_weight = np.sum(weight_sort)
    if tot_weight == 0:
        return np.median(array)
        
    cum_weight = np.cumsum(weight_sort)
    median_idx = np.searchsorted(cum_weight, tot_weight / 2.0)
    
    return array_sort[median_idx]

def kubeviz_percentile(array: np.ndarray, perc: float) -> np.ndarray:
    """
    Computes percentiles of an array.
    Using np.percentile natively, which is an exact analog for IDL's custom percentiles.
    """
    return np.percentile(array, perc)

def kubeviz_sigma_clip(array: np.ndarray, nsig: float = 3.0, nIter: int = 3) -> np.ndarray:
    """
    Iterative sigma clipping.
    Removes outliers beyond nsig standard deviations, iterating nIter times.
    """
    valid_mask = np.ones(array.shape, dtype=bool)
    
    for _ in range(nIter):
        valid_data = array[valid_mask]
        if len(valid_data) == 0:
            break
            
        m = np.median(valid_data)
        s = np.std(valid_data)
        
        # Keep elements within nsig * s
        valid_mask &= np.abs(array - m) < (nsig * s)
        
    return array[valid_mask]

def kubeviz_remove_badvalues(array: np.ndarray, repval: float = 0.0) -> np.ndarray:
    """
    Replaces NaNs and Infs with a replacement value.
    """
    out = array.copy()
    out[~np.isfinite(out)] = repval
    return out

def kubeviz_dataclip(image: np.ndarray, percentage: float = 95.0):
    """
    Clips an image within the requested percentile range.
    e.g. percentage=95 returns range=[2.5%ile to 97.5%ile]
    """
    img = image.flatten()
    valid = img[np.isfinite(img)]
    
    if len(valid) > 0:
        lower = 50.0 - 0.5 * percentage
        upper = 50.0 + 0.5 * percentage
        return np.percentile(valid, [lower, upper])
    return [0.0, 0.0]

def create_mask(shape: tuple, center: tuple, radius: float, mode: str = 'circle') -> np.ndarray:
    """
    Creates a boolean 2D mask.
    mode: 'circle' or 'square'
    """
    ny, nx = shape
    y, x = np.ogrid[:ny, :nx]
    cy, cx = center
    
    if mode == 'circle':
        mask = (x - cx)**2 + (y - cy)**2 <= radius**2
    elif mode == 'square':
        mask = (np.abs(x - cx) <= radius) & (np.abs(y - cy) <= radius)
    else:
        mask = np.zeros(shape, dtype=bool)
        
    return mask
