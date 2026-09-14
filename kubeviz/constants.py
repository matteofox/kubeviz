"""Physical constants and sentinel values shared by the whole package.

Sentinels follow the IDL code exactly so that result files stay interchangeable.
"""

CKMS = 299792.458          # speed of light [km/s]           (state.ckms)
SIGTOFWHM = 2.35482        # sigma -> FWHM for a Gaussian

# Sentinel values used inside the results containers (IDL convention)
NOT_FIT = -999.0           # parameter or error never fitted / reset
NO_ERRORS = -998.0         # value available but no error (moments without MC)

# Flag bits stored in plane 0 of the narrow/broad/continuum result cubes.
# Flag is the sum of the bits that apply; 0 means OK.
FLAG_OK = 0
FLAG_LOW_SN = 1            # S/N of the (main) line below mask_sn_thresh
FLAG_VELERR = 2            # velocity or dispersion error above mask_maxvelerr
FLAG_INVALID_VEL = 4       # velocity / dispersion errors invalid (0 or NOT_FIT)
FLAG_SIGCLIP = 8           # velocity / dispersion outside 4 sigma of the map (fitadj only)
FLAG_MANUAL = 16           # manual flag / bad spaxel / reset spaxel

# Error methods (state.domontecarlo)
ERR_NOISE = 0
ERR_BOOTSTRAP = 1
ERR_MC1 = 2
ERR_MC2 = 3
ERR_MC3 = 4
ERR_METHOD_NAMES = {
    ERR_NOISE: "Noise Cube",
    ERR_BOOTSTRAP: "Bootstrap",
    ERR_MC1: "MonteCarlo 1",
    ERR_MC2: "MonteCarlo 2",
    ERR_MC3: "MonteCarlo 3",
}

# Fit types (state.linefit_type) and modes (state.linefit_mode)
FIT_GAUSS = 0
FIT_MOMENTS = 1
MODE_SPAXEL = 0
MODE_MASK = 1

# Spectrum extraction modes (state.specmode)
SPEC_SLICE = 0
SPEC_SUM = 1
SPEC_MEDIAN = 2
SPEC_WAVG = 3
SPEC_MEDSUB = 4
SPEC_OPTIMAL = 5

# Image modes (state.imgmode)
IMG_SLICE = 0
IMG_SUM1 = 1
IMG_MED1 = 2
IMG_WAVG1 = 3
IMG_WMED1 = 4
IMG_MEDSUB1 = 5
IMG_SUM2 = 6
IMG_MED2 = 7
IMG_WAVG2 = 8
IMG_WMED2 = 9
IMG_MEDSUB2 = 10
IMG_MED2_MINUS_MED1 = 11
IMG_MED1_MINUS_MED2 = 12

# Cube selection in the spaxel viewer (state.cubesel)
CUBE_DATA = 0
CUBE_NOISE = 1
CUBE_BADPIX = 2
CUBE_SN = 3
CUBE_LINEFIT = 4
CUBE_LINEFIT_ERR = 5
CUBE_LINEFIT_SN = 6

# Percentiles kept for the Monte Carlo error distributions (IDL montecarlo_percs).
# The first two MUST be +1 sigma and -1 sigma.
_MED = 50.0
_SIG1 = 50.0 * 0.6826894850
_SIG2 = 50.0 * 0.9544997241
MONTECARLO_PERCS = (
    _MED + _SIG1, _MED - _SIG1,
    _MED + _SIG2, _MED - _SIG2,
    _MED, 10.0, 90.0, 25.0, 75.0,
)
N_MONTECARLO_PERCS = len(MONTECARLO_PERCS)

# Instrumental resolution modes (state.instrres_mode)
INSTRRES_VARPOLY = 0       # polynomial fitted to skylines in the noise cube
INSTRRES_EXTPOLY = 1       # external polynomial (header / templates / hardcoded)
INSTRRES_TEMPLATE = 2      # sigma from a template line profile
MAX_POLYCOEFF_INSTRRES = 5

# Maximum number of lines supported by the fixed-size parameter arrays
MAXLINES = 25
