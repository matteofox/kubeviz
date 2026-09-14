from dataclasses import dataclass, field
from typing import Optional, List, Any
import numpy as np

@dataclass
class KubevizState:
    """
    Central state manager for the Kubeviz application.
    Corresponds to the `kubeviz_state` COMMON block from the original IDL code.
    """
    version: str = 'K3.0-Python'
    debug: bool = False
    ckms: float = 299792.458  # speed of light in km/s
    percent_step: int = 10
    filename: str = ""
    indir: str = ""
    cwdir: str = ""
    outdir: str = ""
    
    # Data arrays
    indatacube: Optional[np.ndarray] = None
    innoisecube: Optional[np.ndarray] = None
    inprihead: Optional[Any] = None  # Header
    indatahead: Optional[Any] = None
    innoisehead: Optional[Any] = None
    datacube: Optional[np.ndarray] = None
    noisecube: Optional[np.ndarray] = None
    noise: Optional[np.ndarray] = None
    badpixelmask: Optional[np.ndarray] = None
    badpixelimg: Optional[np.ndarray] = None
    
    fluxfac: float = 1.0
    noiseisvar: bool = False
    
    # Instrument Specific
    instr: str = ""
    band: str = ""
    instrres_extpoly: Optional[np.ndarray] = None
    ifu: int = -1
    
    # Wavelength calibration
    lambda0: float = 0.0
    dlambda: float = 0.0
    pix0: float = 0.0
    wave: Optional[np.ndarray] = None
    wcs: Optional[Any] = None
    
    # Continuum fitting parameters
    cont_mode: int = 0  # 0: SDSS method (side bands), 1: MPFIT (constant fit)
    cont_minoff: float = 200.0
    cont_maxoff: float = 500.0
    cont_minperc: float = 40.0
    cont_maxperc: float = 60.0
    
    # Autoflag thresholds
    mask_sn_thresh: float = 3.0
    mask_maxvelerr: float = 50.0
    
    # Fit window size around mainline (in Angstroms)
    fit_window: float = 80.0
    
    # Error method
    error_method: str = 'none'
    
    # Fit Lines config
    fit_lines_config: dict = field(default_factory=lambda: {
        'vel': {'min': -500.0, 'max': 500.0},
        'sigma': {'min': None, 'max': None},
        'ha': {'enabled': True, 'min': None, 'max': None},
        'n2r': {'enabled': True, 'min': None, 'max': None},
        'n2b': {'enabled': True, 'min': None, 'max': None}
    })
    
    # Spectral extraction
    medspec: Optional[np.ndarray] = None
    nmedspec: Optional[np.ndarray] = None
    
    # Error estimation modes
    domontecarlo: int = 0
    plotMonteCarlodistrib: bool = False
    saveMonteCarlodistrib: bool = False
    useMonteCarlonoise: bool = False
    scaleNoiseerrors: bool = False
    
    # Display settings
    zmin_ima: float = 0.0
    zmax_ima: float = 0.4
    zmin_spec: float = 0.0
    zmax_spec: float = 0.4
    
    # Instrument and IFU
    instr: str = ''
    band: str = ''
    ifu: int = 0
    vacuum: bool = False
    pixscale: int = 0
    
    # Navigation
    cubesel: int = 0
    col: int = 0
    row: int = 0
    wpix: int = 0
    startcol: int = 0
    startrow: int = 0
    startwpix: int = 0
    ncol: int = 0
    nrow: int = 0
    nwpix: int = 0
    
    # UI flags
    zoommap: bool = False
    linefitmap: bool = False
    scroll: bool = False
    marker: int = -1
    scale: int = -1
    zoomrange: int = 0
    zcuts: int = 0
    
    # Masking and selections
    spaxselect: Optional[np.ndarray] = None
    cursormode: int = 1
    specmode: int = 0
    imgmode: int = 0
    
    # Images for display
    curr_ima: Optional[np.ndarray] = None
    
    # Interactive variables
    wavsel: int = 1
    wavrange1: List[int] = field(default_factory=lambda: [0, 0])
    wavrange2: List[int] = field(default_factory=lambda: [0, 0])
    spaxrange: List[int] = field(default_factory=lambda: [0, 0, 0, 0])
    
    # Zooming
    zoomfac: int = 0
    zoomspax: int = 1
    
    # More state variables will be added here as needed for parity
