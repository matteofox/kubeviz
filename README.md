# kubeviz

Interactive visualisation and emission-line fitting of IFU datacubes.

This repository contains two implementations:

* **`kubeviz/` (Python 3, in development)** - a full port of the IDL program to Python
  (numpy/scipy/astropy, Qt GUI). The headless engine and the batch mode are complete;
  the GUI is under construction. See `PORTING.md` for the routine-by-routine status.
* **`kubeviz.pro` (IDL, V2.2)** - the original program, unchanged, with its `addons/`,
  `templates/` and `doc/` directories. The IDL instructions follow below.

## Python quick start

```bash
pip install -e .            # engine only (numpy, scipy, astropy)
pip install -e ".[gui]"     # plus PyQt6 and pyqtgraph for the GUI (phase 3)
python -m pytest            # run the test-suite on synthetic cubes
```

Batch mode (fit every spaxel of a cube, write the results FITS and a session file):

```bash
kubeviz cube.fits --redshift 0.85 --lineset 1 --batch --fit-all-lines --outdir results/
```

Every keyword of the IDL procedure is available as a command line flag (`kubeviz --help`).
Results files use the same plane layout and header keywords as the IDL version, so they can be
exchanged between the two implementations. Python sessions are saved as `.kvz` files.

From Python:

```python
from kubeviz.session import start_session
from kubeviz.fitting.linefit import dofit
from kubeviz.fitting.fitall import fitall
from kubeviz.fitting.flags import autoflag
from kubeviz.io.results import saveres

state = start_session("cube.fits", redshift=0.85, lineset=1)
state.pndofit[:state.Nlines] = 1      # fit the narrow component of every line in the set
state.pcdofit[:state.Nlines] = 1      # and the continuum
fitall(state)                          # or dofit(state) for the current spaxel (state.col/row)
autoflag(state)
saveres(state, "cube_res.fits")
res = state.res                        # ResultSet: res.n[row, col] = [flag, dv, sigma, fluxes...]
```

---

# IDL version (V2.2) - original README

Acknowledgements:

If you use KUBEVIZ in your work please add the following sentence in the ackowledgement section of any resulting publication:

This work made use of the KUBEVIZ software which is publicly available at https://github.com/matteofox/kubeviz/.
The developement of the KUBEVIZ code code was supported by the Deutsche Forschungsgemeinschaft via Project IDs: 
WI3871/1-1 and WI3871/1-2

And please cite Fossati et al. (2016) http://adsabs.harvard.edu/abs/2016MNRAS.455.2028F.

*****************

Setting up kubeviz (V2.0)

1. Setup appropriate directory structure to store kubeviz versions
   outside IDL_PATH 

   mkdir kubeviz

2. Unpack the gzipped tarfile in the directory you have created:

   cd kubeviz
   tar xvfz kubeviz_v2.0.tar.gz

3. This creates a subdirectory with the specific version number 
   e.g. kubeviz/v2.0 Now link the current version inside IDL_PATH (assumed
   IDL root directory $IDL_ROOT):

   ln -s kubeviz/v2.0 $IDL_ROOT/kubeviz_current
   
   Make sure your IDL_PATH runs recusrsively into the subdirectories 
   (addons/ must be in the IDL_PATH).
   
3. Run kubeviz!

   idl
   > kubeviz, 'cube.fits'

4. Documentation
   
   The doc/ directory of the package includes an instruction file
   and a detailed changelog. We suggest you to have a quick look 
   at the instructions.txt file before using the code. This describes
   the calling sequence and the features of the interactive environment.
   
   From the IDL prompt run:
   > doc_library, 'kubeviz'
   to see the header of the code and the available keywords.

   In the GUI, select 
      'help' -> 'Instructions' 
      'help' -> 'What's new' 
   to see the documentation.   
   

!! IMPORTANT !!
  
   Kubeviz requires the following external libraries:
   * Astrolib Library (updated after 2012, or late 2014 if you load MUSE cubes)
         http://idlastro.gsfc.nasa.gov/contents.html
   * Some routines in the astrolib make use of programs in the Coyote Library. 
     Although they have recently bundled into the astrolib package it is 
     recommended to have a full and up-to-date  installation of the Coyote 
     library which must be downloaded separately.
   * Mpfit suite of code (V1.70 or greater)
         http://cow.physics.wisc.edu/~craigm/idl/fitting.html

   If kubeviz does not work properly and stops when calling one of those 
   third-party libraries, the first thing to do is to update the library. 
   Then restart IDL and try again.

(BASIC) TROUBLESHOOTING

1) If you see this message:
    
    [WARNING] ##       This is a 32 bit IDL version.        ##
    [WARNING] ##   Some features might now work correctly!  ##
   
   There is no reason to worry unless you are trying to load a very large 
   datacube. Almost certainly MUSE datacubes require a 64bit IDL version 
   for the load to be successful.

2) If you see errors like:
      
      % Error occurred at: VALID_NUM   line 71     
      % AIRTOVAC: Incorrect number of arguments.
   
   It means you have not updated the Astrolib library since 2012. 
   It is important to keep Astrolib updated frequently as bugs 
   are fixed and the library is gradually being converted to handle 
   64bit variables whenever possible. 

3) If you get an error when you try to load data from an instrument 
   which is not in the list of supported instrument please report it 
   to mfossati at mpe.mpg.de. We will try to solve the issue.
