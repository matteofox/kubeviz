;+
;#############################################################################
;
; Copyright (C) 2014-2019, David J. Wilman, Matteo Fossati, Joris Gerssen
; E-mail: matteo.fossati_at_durham.ac.uk
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;#############################################################################
;
;
; NAME: kubeviz
; 
; VERSION: 2.2
; 
; AUTHORS: David J. Wilman, Matteo Fossati, & Joris Gerssen
; 
; PURPOSE: Interactively examine a FITS datacube (x,y,lambda).
;          kubeviz also provides interactive functionality to fit emission lines 
;          using mpfit, & select sub-sections of the cube for line fitting and 
;          further arithmetic. 
; 
; ORIGIN: Based on cubeviz by Joris Gerssen which hosts the interactive visualization 
;         tool with some of the functionality and no fitting.
; 
; CALLING SEQUENCE:
;         kubeviz [,datafile] [,noisefile=noise_file] [,ext=ext] [,noise_ext=noise_ext] 
;         	  [,trim=trim] [,fluxfac=fluxfac] [,smooth=spax_smooth] [,specsmooth=spec_smooth] 
;         	  [,traspose=transpose] [,help=help] [,version=version] [,debug=debug]
;         	  [,vacuum=vacuum] [,logarithmic=logarithmic] [,waveunit=waveunit] 
;         	  [,redshift=redshift] [,do_mc_errors=do_mc_errors] [,bootstrap=bootstrap_file] 
;         	  [,nmontecarlo=nmontecarlo]  [,plot_mc_pdf=plot_mc_pdf] [,save_mc_pdf=save_mc_pdf]   
;         	  [,use_mc_noise=use_mc_noise] [,scroll=scroll] [,lineset=lineset]  
;         	  [,mask_sn_thresh=mask_sn_thresh] [,mask_maxvelerr=mask_maxvelerr]  
;         	  [,mom_thresh=mom_thresh] [,spmask=spmask] [,instr=instr] [,band=band]
;         	  [,batch=batch] [,fit_all_lines=fit_all_lines] [fix_ratios=fix_ratios]
;         	  [,resfile=results_file] [,fittype=fittype] [,fitpars=fitpars] 
;                 [,momwindowfile=momwindowfile] [,gaussinitfile=gaussinitfile]
;
; INPUT:    Datacube (x, y, lambda) in FITS file format or kubeviz saved session
;
; OUTPUT:   None
; 
; KEYWORDS: datafile:    name of the datacube FITS file.
;                        If no path is given, the file is assumed to be in the current directory.
;                        If the file has extension '.sav' it is assumed to be a previously saved
;                        session (all other keywords are disregarded)
;           noisefile:   (optional) Name of the separate noise file. 
;           ext:         (optional) FITS extension of data cube (default=1, if noisefile set default=0)
;           noise_ext:   (optional) FITS extension of noise cube (default=2, if noisefile set default=0,1)
;           trim:        (optional) formatted string to define which portion of the datacube needs to be
;                        loaded. The string must have the following syntax '[x1:x2,y1:y2,z1:z2]'. The
;                        wildcard * can be used to select all element along a dimension ('[x1:x2,*,*]').
;                        The indices are 1 based as in IRAF and ds9 and both extremes are included.  
;           fluxfac:     (optional) factor to multiply flux-scale (default instrument dependent)
;           smooth:      (optional) spatial smoothing kernel in spaxels (default 1)
;           transpose:   (optional) transpose data? 
;           help:        (optional) Show instructions and exit
;           version:     (optional) Show version and exit
;           logarithmic: (optional) Use if wavelength solution in header is logarithmic
;           vacuum:      (optional) Use if wavelength solution in header is in vacuum. For known instruments 
;                                   it is automatically set.
;           waveunit:    (optional) Unit of wavelength if not in CUNIT3 keyword (overrides). Converted to Angstroms.
;           redshift:    (foptional) redshift of target (default 0)
;           do_mc_errors:(optional) choice of error type for linefitting. Default 0 (noise cube)
;                         0: Single fit with errors from noise cube
;                         1: Bootstrap cubes (requires multi-extension bootstrap cubes-frame as input)
;                         2: Monte Carlo 1 : generates monte-carlo cubes via adding random gaussian-distributed 
;                                            errors with the noise cube as amplitude and THEN smooth.
;                         3: Monte Carlo 2 : As Monte Carlo 1 but smooths before generating monte-carlo cubes.
;                         4: Monte Carlo 3 : As Monte Carlo 2 but with residuals from a simple fit instead of the noise. 
;                                            The residuals are randomly shifted along the wavelenght axis
;           bootstrap:   (optional) bootstrap fits file = multi-extension fits datacube
;                        Primary header should contain NEXT parameter containing number of extensions (realisations). 
;                        Then each extension should contain one realisation of the datacube. 
;                        The file is searched in the datacube directory if no path is given.  
;           nmontecarlo: (optional) number of monte-carlo realisations where required (default 100)
;           plot_mc_pdf: (optional) save diagnostic plots of distribution of results with multiple realisations
;           save_mc_pdf: (optional) save the distribution of montecarlo results into a multiextension FITS table.
;           use_mc_noise:(optional) use the noise cube computed from the monte-carlo cubes
;           scroll:      (optional) force display of scrollbar in linefit window
;           lineset:     (optional) select only one set of lines to fit. Default=0 (selects all lines).
;                        1: Ha+[NII]
;                        2: [OIII]
;                        3: [OII]
;                        4: Hbeta
;                        5: [SII]
;                        6: [OI]
;   	    	    	 7: [SIII]
;   	    	    	 8: [HeI]
;                        9: Lya
;           debug:       (optional) verbose
;           mask_sn_thresh: (optional) S/N (amplitude) threshold for autoflag acceptance of a fit.  Default=3.
;           mask_maxvelerr: (optional) velocity (or dispersion) error threshold for autoflag accepatance of a fit. Default=50 km/s.
;           mom_thresh:  (optional) flux error threshold for 1st and higher moments computation. 
;                        Not used for the 0th moment (total flux). Default 0.5
;           spmask:      (optional) name of the spaxel mask FITS file to be loaded.
;           instr:       (optional) instrument. If datafile header contains INSTRUME keyword it overrides.
;           band:        (optional) band (for info on instrumental resolution and wavelength range). 
;                        For ESO instruments taken from header.
;           batch:       (optional) run kubeviz without the GUI. Run fit of all spaxels, save the results file 
;                        and the session file. If spmask is given the fit loops over all the masks, saves the 
;                        results and the session file.     
;           fit_all_lines: (optional, works only with /batch) fit all the lines in the selected lineset. 
;                        If not set fit only the mainline.
;           fix_ratios:  (optional) fix the intensity ratios for NII, OIII and OI lines. 
;                        Ratios from Storey & Zeippen (2000)
;           resfile:     (optional) load a kubeviz formatted results file. It is searched in the current 
;                        working dir if no path is given.
;           fittype:     (optional) choice of type of fitting. Can be 'gauss' (Default) or 'moments'. 
;           outdir:      (optional) the directory where output files are written
;           logunit:     (optional) logical unit number where the terminal output is redirect. The unit must be 
;                                   already open. If not set the output is sent to the standard output.
;           fitpars:     (optional) a vector of 6 or 8 elements specifying (1,2 optional) fitrange; 3,4 cont 
;                                   min/max offset; 5,6 cont min/max percentiles; 7 continuum polynomial order; 
;                                   8 continuum mode (0 = use offset windows, 1= use MPFIT in fitrange)  
;           momwindowfile: (optional) a file containing a datacube with 2 planes. The first defines the mean velocity and 
;                                     the second the range of velocity for the moments window used to integrate the 
;                                     flux (zeroth order moment)
;           gaussinitfile: (optional) a file containing a datacube for which each plane corresponds to the initial 
;                                     estimates for gaussian parameters
;
; COMMON_BLOCKS:
;           kubeviz_state:     structure to pass variables/data between 
;                              program modules
;           kubeviz_winsize:   structure listing the physical sizes in pixels
;                              of the drawing windows
;           kubeviz_widgetids: structure to store/pass the IDs of the graphical widgets
;           kubeviz_fit:       structure to store/pass fit parameters between program
;                              modules
;           kubeviz_linesdb    structure to store the rest frame wavelengths of 
;                              emission lines and their names
;
; RESTRICTIONS:
;           - requires IDL 6.1 or higher to dynamically resize the specrum windows
;           - modifies colour table
;
; EXAMPLES:  
;
; kubeviz, 'combined_gs3_19791.fits', redshift=2.225, lineset=1, smooth=2, bootstrap='bootstrap_gs3_19791.fits'
; 
; EXTERNAL CODE:
; 
;           Astrolib Library: http://idlastro.gsfc.nasa.gov/contents.html 
;                   (Requires a version released later than 2011)
;           Some routines in the astrolib make use of programs in the Coyote Library, 
;                   which must be downloaded separately.
;           Mpfit suite of code:  http://cow.physics.wisc.edu/~craigm/idl/fitting.html 
;                   (Requires Version 1.70 or greater)
;           xpdmenu.pro, zscale_range.pro, fxpar_sp.pro, fxread_cube.pro, plotimage.pro 
;                   (Kubeviz specific versions included in addons)
;           
; ENVIRONMENT:
;           Requires about 1Gb of RAM when working on a 20x20x2048 cube. Can require up to 3Gb of RAM 
;           if all error methods are used. Included calibrations: KMOS polynomial fit to spectral 
;           resolution (M.Fossati, E. Wisnioski), Night sky spectrum for the OH lines identification.
;
; MODIFICATION HISTORY:  
; See the full changelog in Help -> What's new
;   
; 10-02-2014 : Release of V1.0
; 28-02-2014 : Release of V1.1
;              Add the batch mode and the possibility to load a session or a result file
;              Add montecarlo 3 method for error estimate. Add display header option in file menu
;              Add a button to fit bad spaxels using initial guess from adjacent spaxels
;              Add a button to set an automatic mask based on S/N in the linefit range 
;              Countless bug fixes including redshift change, save results and fit parameters
; 26-03-2014 : Release of V1.2
;              Add the moments as an alternative way of measuring emission lines parameters
;              The noise can be correctly computed from the mc cubes if the appropriate option is on. 
;              Several options can now be changed interactively (rather than only from the command line)
;              Various performance and memory improvements. 
; 31-10-2014 : Release of V1.3
;              Works natively with SINFONI and MUSE data. Performance improvements to handle large IFUs.
;              Flux ratios of OI, OIII, NII lines can be fixed when fitting.
;              Fully redesigned routine to evaluate/read instrumental resolution values. 
;              The code now uses vacuum wavelengths for cryo instruments.
; 24-12-2014 : Release of V1.3.1
;              More robust method to derive the line initial guess.
;              Memory improvements
; 26-06-2015 : Release of V1.3.2
;              Spaxel masks are now processed as 2D fits images rather than 1D vectors
;              GUI improvements, improvements in the moments algorithm.
;              Several minor bug fixes
; 22-10-2015 : Release of V1.3.3
;              Addons are now kubeviz specific to avoid conflicts
;              Several minor bug fixes and stability improvements
; 24-11-2015 : Release of V1.4
;              Improved spectrum viewer title bar, minor bug fixes.
;              Test that 3rd party libraries are present and up-to-date.
;              Save current image as FITS now available.
; 16-12-2015 : Release of V1.4.1
;              Critical bug fix preventing the results to be saved if lineset=0
;              Improved intial guess on the line position.
; 12-06-2016 : Release of V2.0!
;              Zooming on the spaxel viewer is now possible.
;              Circular and square spaxel masks can be now conventiently defined
;              A second gauss comp. can be forced to be at higher z or larger width in the fit
;              The range of spaxels for fitall and fitadjall can be modified by the user. Fit runs
;              can now be interrupted.
; 06-11-2017 : Release of V2.1
;              Several bug fixes. Significantly improved algorithm for the iterative fitadjall.
;              Two component fits have a "smart" option to decide whether two offset kinematical 
;              components are present in a given spaxel (see documentation for more details)
;              (BETA) User initial guess maps (gaussian) and window parameter
;              maps (moments) can be loaded.
; 11-06-2019 : Release of V2.2
;              Several bug fixes. Added O2 to the database of lines. Lineset hierarchy has been
;              re-defined.
;
; BUGS         Needs to reset the results container but NOT the widgets
;              unless we are definitely running in interactive mode!!
; FUTURE WORK: improve consistency of initial guess maps (switch on/off, load from GUI)
;              improve handling of initial guess maps for gaussian fits
;-------------------------------------------------------------------------
;-

FUNCTION kubeviz_str, num, format=format
if n_elements(format) gt 0 then return, strcompress(string(num,format=format),/REMOVE_ALL) else return, strcompress(string(num),/REMOVE_ALL)
end

FUNCTION kubeviz_SetUnion, a, b
IF a[0] LT 0 THEN RETURN, b    ;A union NULL = a
IF b[0] LT 0 THEN RETURN, a    ;B union NULL = b
RETURN, Where(Histogram([a,b], OMin = omin)) + omin ; Return combined set
END

FUNCTION kubeviz_SetIntersection, a, b
minab = Min(a, Max=maxa) > Min(b, Max=maxb) ;Only need intersection of ranges
maxab = maxa < maxb

   ; If either set is empty, or their ranges don't intersect: result = NULL.

IF maxab LT minab OR maxab LT 0 THEN RETURN, -1
r = Where((Histogram(a, Min=minab, Max=maxab) NE 0) AND  $
          (Histogram(b, Min=minab, Max=maxab) NE 0), count)

IF count EQ 0 THEN RETURN, -1 ELSE RETURN, r + minab
END

FUNCTION kubeviz_SetDifference, a, b

   ; = a and (not b) = elements in A but not in B

mina = Min(a, Max=maxa)
minb = Min(b, Max=maxb)
IF (minb GT maxa) OR (maxb LT mina) THEN RETURN, a ;No intersection...
r = Where((Histogram(a, Min=mina, Max=maxa) NE 0) AND $
          (Histogram(b, Min=mina, Max=maxa) EQ 0), count)
IF count eq 0 THEN RETURN, -1 ELSE RETURN, r + mina
END

FUNCTION kubeviz_weighted_median, array, weight

tot_weight = total(weight)

ind = sort(array)
array_sort = array[ind]
weight_sort = weight[ind]
sum = 0
n=0

for i=0,size(array, /N_elements)-1 do begin
  sum += weight_sort[i]
  n += 1
  if sum gt tot_weight/2. then break
endfor

return, array_sort[n-1]

end

FUNCTION kubeviz_remove_badvalues, array, repval=repval

if n_elements(repval) eq 0 then repval=0.
badvalues = where(finite(array) eq 0, Nnotvalid)

if Nnotvalid gt 0 then begin
  array[badvalues] = repval
  return, array
endif else return, array

end

FUNCTION kubeviz_percentile, X, perc

Nperc = N_Elements(perc)
N = N_Elements(X)
s = sort(X)
Xs = X(s)
pos = .01*perc*(N-1) ; 0,1,2...
pos1 = fix(pos,type=3)
pos2 = pos1+1
frac1 = pos2-pos
frac2 = pos-pos1
okpos1 = where(pos1 ge 0, Nokpos1)
okpos2 = where(pos2 le N-1, Nokpos2)
contr1 = replicate(0.,Nperc)
contr2 = replicate(0.,Nperc)
if Nokpos1 gt 0 then contr1(okpos1) = frac1(okpos1)*Xs([pos1(okpos1)])
if Nokpos2 gt 0 then contr2(okpos2) = frac2(okpos2)*Xs([pos2(okpos2)])
Xperc = contr1+contr2

return, Xperc
end

FUNCTION kubeviz_closest, array, value

indexvec = 0
if (n_elements(value) eq 0) then indexvec=[indexvec,-1] else begin
  if (n_elements(array) eq 0) then indexvec=[indexvec,-1] else begin
    for i=0,n_elements(value)-1 do begin
 	abdiff = abs(array-value[i])	 ;form absolute difference
 	mindiff = min(abdiff,index)	 ;find smallest difference
	indexvec=[indexvec,index]
    endfor	
  endelse
endelse
return,indexvec[1:n_elements(indexvec)-1]
end

pro kubeviz_sigma_clip,array,nsig=nsig,nIter=nIter,index=index

if n_elements(nsig) eq 0 then nsig=3
if n_elements(nIter) eq 0 then nIter=3

wif=n_elements(array)
index=lindgen(wif)

for i=0,nIter-1 do begin
    m=median(array[index])
    s=stddev(array[index])
    clip=nsig*s
    w=where(abs(array[index]-m) lt clip,wif)
    index=index(w)
endfor
array=array[index]
end

pro kubeviz_cumdistrplot, x, xr=xr, yr=yr, oplot=oplot, col=col, xtit=xtit, ytit=ytit, tit=tit, thick=thick, linestyle=linestyle, ylog=ylog, ynorm=ynorm, renorm=renorm, charsize=charsize, charthick=charthick, reverse=reverse, yoff=yoff, title=title, xout=xout, yout=yout

Nx = N_Elements(x)
if N_Elements(xr) eq 0 then xr=[min(x),max(x)]
if N_Elements(renorm) eq 0 then renorm=0
if N_Elements(yr) eq 0 then if renorm eq 0 then yr=[0.,1.] else yr=[0.,1.*total(x ge xr[0] and x le xr[1])/Nx]
    
if N_Elements(oplot) eq 0 then oplot = 0
if N_Elements(col) eq 0 then col=2
if N_Elements(xtit) eq 0 then xtit=''
if N_Elements(ytit) eq 0 then ytit=''
if N_Elements(title) eq 0 then title=''
if N_Elements(thick) eq 0 then thick=1.
if N_Elements(linestyle) eq 0 then linestyle=0
if N_Elements(ylog) eq 0 then ylog=0
if N_Elements(ynorm) eq 0 then ynorm=1.
if N_Elements(yoff) eq 0 then yoff=0.
if N_elements(charsize) eq 0 then charsize=1.
if N_elements(charthick) eq 0 then charthick=1.
if N_elements(reverse) eq 0 then reverse=0

y = yoff+(ynorm*findgen(Nx)/((1.0-yoff)*Nx))
if reverse eq 1 then sx = x(reverse(sort(x))) else sx = x(sort(x))

!x.style = 1
!y.style = 1

if oplot eq 0 then begin
    plot, sx, y, col=col, xr=xr,yr=yr, xtit=xtit, ytit=ytit, title=title, xthick=thick, ythick=thick, thick=thick, linestyle=linestyle, ylog=ylog, charsize=charsize, charthick=charthick
endif else begin
    oplot, sx, y, col=col, thick=thick, linestyle=linestyle
endelse

xout=sx
yout=y

end


function kubeviz_dataclip, image, percentage=pclip

; clip image within the percentile range requested
; e.g. percentage=95 requires range=[2.5%ile to 97.5%ile]

; default 95: 
if n_elements(pclip) eq 0 then pclip=95.

; Remove nans
img = reform(image, (size(image))[1]*(size(image))[2])
valid = where( finite(img) eq 1, Nfinite)

if Nfinite gt 0 then begin

img = img[valid]
range = kubeviz_percentile(img,50.-[1.,-1.]*.5*pclip)
return, range

endif else return, [0,0]

end

;-------------------------------------------------------------------
;  THE REAL CODE STARTS HERE
;-------------------------------------------------------------------

pro kubeviz_common

; initialize the state vector
common kubeviz_state, state
common kubeviz_winsize, winsize
common kubeviz_widgetids, widgetids
common kubeviz_fit, fit
common kubeviz_linesdb, linesdb
common kubeviz_pars, pars


state = {  $
          version: 'K2.2',                 $  ; version number 
	  log_lun: -1,                     $  ; log file unit, if -1 print to the terminal
	  debug: 0L,                       $  ; if set print additional information
          ckms: 299792.458D,               $  ; speed of light in km/s
          percent_step: 10,                $  ; step for progress stats for long calculations
	  filename: "",                    $  ; name of input file
          indir: "",                       $  ; input directory
	  cwdir: "",                       $  ; working directory
          outdir: "",                      $  ; output directory
	  indatacube: ptr_new(),           $  ; pointer to the (unsmoothed) input datacube array
          innoisecube: ptr_new(),          $  ; pointer to the (unsmoothed) input noisecube array
          inprihead: ptr_new(),            $  ; pointer to the primary header (if exists)
	  indatahead: ptr_new(),           $  ; pointer to the header of the input datacube array
          innoisehead: ptr_new(),          $  ; pointer to the header of the input datacube array
	  datacube: ptr_new(),             $  ; pointer to the smoothed datacube array
          noisecube: ptr_new(),            $  ; pointer to the smoothed noisecube array
          noise : ptr_new(),               $  ; pointer to the current noisecube (can be from data or from montecarlo)
          badpixelmask: ptr_new(),         $  ; bad pixel mask (3d cube)
          badpixelimg:  ptr_new(),         $  ; pointer to the badpixel image (2d image)
          fluxfac: 1.,                     $  ; multiplicative flux factor to be applied to cubes
          noiseisvar: 0B,                  $  ; is the original noise extension the variance or the STD?
	  lambda0: 0.D,                    $  ; wavelength calibration parameter
          dlambda: 0.D,                    $  ;  ditto
          pix0: 0.D,                       $  ;  ditto
          wave: ptr_new(),                 $  ; pointer to the (wcs) wavelength array
          medspec: ptr_new(),              $  ; pointer to the median spectrum array (can also be summed spec)
          nmedspec: ptr_new(),             $  ; pointer to the median spectrum noise array
          domontecarlo: 0L,                $  ; do montecarlo error estimate?
          montecarlocubes: replicate(ptr_new(),1000),    $ ; array of pointers to cubes for montecarlo errors , current max 1000
          medspec_montecarlo: replicate(ptr_new(),1000), $ ; array of pointers to median/summed spec, current max 1000
          Nmontecarlo : 0L,                $  ; Number of cubes of the selected MonteCarlo method
          Nmontecarlo_percs: 9L,           $  ; Number of percentiles to keep results for MonteCarlo cubes
          montecarlo_percs: ptr_new(),     $  ; pointer to array containing the percentiles to keep results for MonteCarlo cubes
          plotMonteCarlodistrib: 0B,       $  ; boolean: plot the distribution of bootstrap fit values?
          saveMonteCarlodistrib: 0B,       $  ; boolean: plot the distribution of bootstrap fit values?
	  useMonteCarlonoise: 0B,          $  ; boolean: use the noise obtained from the montecarlo cubes?
	  scaleNoiseerrors:0B,             $  ; boolean: scale mpfit errors assuming the chi-sq is 1? Only for noisecube errors
	  zmin_ima: 0.0,                   $  ; minimum intensity to use in the spaxel viewer
          zmax_ima: 0.4,                   $  ; maximum intensity to use in the spaxel viewer
          zmin_spec: 0.0,                  $  ; minimum intensity to use in the spectrum viewer
          zmax_spec: 0.4,                  $  ; maximum intensity to use in the spectrum viewer
	  instr: '',                       $  ; instrument name
          band: '',                        $  ; band
          ifu: 0L,                         $  ; IFU number for multi ifu instruments
          vacuum: 0B,                      $  ; is the wavelength solution computed in vacuum?
	  pixscale: 0L,                    $  ; pixel scale for sinfoni cubes
	  cubesel: 0L,                     $  ; which cube? (0=data/1=noise/..)
          col: 0L,                         $  ; column number of current spaxel (image coordinates)
          row: 0L,                         $  ; row number of current spaxel    (image coordinates)
          wpix: 0L,                        $  ; current wavelength pixel number
	  Startcol: 0L,                    $  ; first column extracted from the original image (phys coordinates)
	  Startrow: 0L,                    $  ; first row    extracted from the original image (phys coordinates)
          Startwpix: 0L,                   $  ; first wpixel extracted from the original image (phys coordinates)
	  Ncol: 0L,                        $  ; number of columns
          Nrow: 0L,                        $  ; number of row
          Nwpix: 0L,                       $  ; number of wavelength pixels
          zoommap:0B,                      $  ; zoom panel mapped (0/1)
          linefitmap:0B,                   $  ; linefit window mapped (0/1)
	  scroll: 0B,                      $  ; scroll linefit window? 
          marker: 0L,                      $  ; over plot marker lines (-1/1)
          scale: 0L,                       $  ; scale spectrum window to zmin,zmax (-1/1)
          zoomrange: 0L,                   $  ; half-range of spectrum zoom window 
          zcuts: 0L,                       $  ; display cuts (User/MinMax/Zscale/HistEq)
          spaxselect: ptr_new(),           $  ; array to enumerate user selected spaxels
          cursormode: 1L,                  $  ; cursor mode: crosshair or spaxel select
          specmode: 0L,                    $  ; plot single, median, sum or single - median or optimal
          imgmode: 0L,                     $  ; plot slice, median, sum or slice - median
          img1: ptr_new(),                 $  ; pointer to the median/summed image
          img2: ptr_new(),                 $  ; pointer to the median/summed image
          nimg1: ptr_new(),                $  ; pointer to the median/combined noise image 
          nimg2: ptr_new(),                $  ; pointer to the median/combined noise image 
          bimg1: ptr_new(),                $  ; pointer to the median/combined badpixel image 
          bimg2: ptr_new(),                $  ; pointer to the median/combined badpixel image 
          curr_ima: ptr_new(),             $  ; pointer to the currently displayed image
	  curr_byteima: ptr_new(),         $  ; pointer to the currently displayed image
	  curr_byteima_grid: ptr_new(),    $  ; pointer to the currently displayed and regridded image
	  lineresimg: ptr_new(),           $  ; pointer to the linefit masked results image
          lineerrresimg: ptr_new(),        $  ; pointer to the linefit masked error results image
          unm_lineresimg: ptr_new(),       $  ; pointer to the linefit unmasked results image
          unm_lineerrresimg: ptr_new(),    $  ; pointer to the linefit unmasked error results image
	  flagmode : 1B,                   $  ; flag or unflag the linefit result image
          wavsel: 1L,                      $  ; selected wavelength range indicator (Default = 1)
          wavrange1:intarr(2),             $  ; range of selected wavelength
          wavrange2:intarr(2),             $  ; range of selected wavelength
          spaxrange:intarr(4),             $  ; range of currently displayed spaxels
	  drag_spax: 0L,                   $  ; flag set when mouse button hold down in spax
          drag_spec: 0L,                   $  ; flag set when mouse button hold down in spec
          zoomfac: 0L,                     $  ; Zoom factor for the spaxel viewer
	  zoomspax: 1L,                    $  ; Zoom factor for interactive rescaling of spaxel maps (can be included in the previous?)
	  maxcol: 0L,                      $  ; number of color levels
          red: bytarr(256),                $  ; red colour table indices
          green: bytarr(256),              $  ; green colour table indices
          blue: bytarr(256),               $  ; blue colour table indices
          invert:-1L,                      $  ; switch to invert colour tables
          ctab: 0L,                        $  ; store the (selected) colour table
          smooth: 1L,                      $  ; spatial smoothing size in pixels
          specsmooth: 1L,                  $  ; spectral smoothing size in pixels
          transpose: 0L,                   $  ; transpose frames? (0/1)
          npress: 0L,                      $  ; counter for the number of 's' key presses 
          Nmask: 1L,                       $  ; number of masks (initialize to one)
          maxNmask: 50000L,                $  ; maximum number of masks
          imask: 1L,                       $  ; index to keep track of the mask to use (base 1)
	  maskradius :1L,                  $  ; radius of the mask (for circular) or size (for square)
          maskmode: 0L,                    $  ; mask mode 0: single spaxel, 1: circular, 2:square
	  nrescube: ptr_new(),             $  ; pointer to the mask linefit results data cube for narrow lines
          brescube: ptr_new(),             $  ; pointer to the mask linefit results data cube for broad lines
          crescube: ptr_new(),             $  ; pointer to the mask linefit results data cube for underlying continuum
          mrescube: ptr_new(),             $  ; pointer to the mask moment results
          nerrrescube: ptr_new(),          $  ; pointer to the errors on mask linefit results data cube for narrow lines
          berrrescube: ptr_new(),          $  ; pointer to the errors on mask linefit results data cube for broad lines
          cerrrescube: ptr_new(),          $  ; pointer to the errors on mask linefit results data cube for underlying continuum
          merrrescube: ptr_new(),          $  ; pointer to the errors on mask moment results
          sp_nrescube: ptr_new(),          $  ; pointer to the spaxel linefit results data cube for narrow lines
          sp_brescube: ptr_new(),          $  ; pointer to the spaxel linefit results data cube for broad lines
          sp_crescube: ptr_new(),          $  ; pointer to the spaxel linefit results data cube for underlying continuum
          sp_mrescube: ptr_new(),          $  ; pointer to the spaxel moment results
          sp_nerrrescube: ptr_new(),       $  ; pointer to the errors on spaxel linefit results data cube for narrow lines
          sp_berrrescube: ptr_new(),       $  ; pointer to the errors on spaxel linefit results data cube for broad lines
          sp_cerrrescube: ptr_new(),       $  ; pointer to the errors on spaxel linefit results data cube for underlying continuum
          sp_merrrescube: ptr_new(),       $  ; pointer to the errors on spaxel moment results
          continuumfit_mode: 0D,           $  ; Select continuum fitting mode. 0: SDSS method, 1: fit internally using MPFIT
	  continuumfit_minoff: 200.D,      $  ; min offset in wavelength units to define continuum region
          continuumfit_maxoff: 500.D,      $  ; max offset in wavelength units to define continuum region
          continuumfit_minperc: 40.D,      $  ; min %ile of spectrum in continuum region to estimate continuum level
          continuumfit_maxperc: 60.D,      $  ; max %ile of spectrum in continuum region to estimate continuum level	    
          continuumfit_order: 0L,          $  ; order of fit to continuum in continuum region
          maxwoffb: 80.D,                  $  ; max wavelength offset from lineset to fit bluewards (Angstroms)
          maxwoffr: 80.D,                  $  ; max wavelength offset from lineset to fit redwards (Angstroms)
          secondcomp_mode: 0L,             $  ; Select second component mode. 0: free, 1: larger vel, 2: larger sigma
	  secondcomp_smart: 0B,            $  ;Use advanced rejection method for second component?
	  instrres_mode: 0L,               $  ; mode for instrumental resolution = 0 (fit to variance), 1 (use polynomial) 2 (use spec. templates)
          instrres: ptr_new(),             $  ; instrumental resolution per mask
          sp_instrres: ptr_new(),          $  ; instrumental resolution per spaxel
          instrres_tplsig: 0D,             $  ; sigma from templates for instrumental resolution
	  instrres_extpoly: replicate(0.D,7), $  ; array of polynomial coefficients for instrumental resolution from Emily's fit
          instrres_varpoly: replicate(0.D,7), $  ; array of polynomial coefficients for isntrumental resolution from noise spectra
	  max_polycoeff_instrres: 5L,      $  ; max. polynomial order to describe instrumental resolution as f(wavelength)
          fitallrange:replicate(0,4),      $  ; range of spaxels where fitall should loop
	  Nloopadj:0L,                     $  ; number of loops for fit adj (0= until no improvement)
	  gnfit: replicate(-999.,25),      $  ; array storing narrow line linefit user-supplied parameters
          gbfit: replicate(-999.,25),      $  ; array storing broad line linefit user-supplied parameters
          pnfix: bytarr(25),               $  ; narrow line linefit parameters fitted(0) or kept fixed (1)
          pbfix: bytarr(25),               $  ; broad line linefit parameters fitted(0) or kept fixed (1)
          gnlims: replicate(-999.,25,2),   $  ; array storing narrow line linefit param limits
          gblims: replicate(-999.,25,2),   $  ; array storing broad line linefit param limits
          pndofit: bytarr(25),             $  ; array storing  narrow line linefit selection of lines to fit
          pbdofit: bytarr(25),             $  ; array storing broad line linefit selection of lines to fit
          pcdofit: bytarr(25),             $  ; array storing selection for which lines to fit a continuum 
          nshow: replicate(0B,25),         $  ; show narrow linefit in speczoom window?
          bshow: replicate(0B,25),         $  ; show broad linefit in speczoom window?
          cshow: replicate(0B,25),         $  ; show continuum fit in speczoom window?
          par_imagebutton: '',             $  ; Encoded parameter of selected image button
          fitconstr: 0B,                   $  ; Choice to fit with constraints or without.
          fitfixratios: 0B,                $  ; Choice to lock line ratios of (some) lines.
	  Nlines_all: 0L,                  $  ; Number of emission lines in the spectrum
          Nlines: 0L,                      $  ; Number of emission lines in the selected lineset
          lines_all: fltarr(25),           $  ; Estimated central wavelengths of all lines
          lines: fltarr(25),               $  ; Estimated central wavelengths of all lines in the current session
	  linesets: ptr_new(),             $  ; Pointer to array associating each line in the session to a lineset
          linenames: ptr_new(),            $  ; Pointer to array of line names.
          linefancynames: ptr_new(),       $  ; Pointer to array of line names to be used in linefit window.
          line_bisectors: fltarr([25,2]),  $  ; Wavelength bisectors of lines used to define the wavelength range used for each line 
	  lineset_ind: replicate(ptr_new(),100), $  ; Pointer to arrays of lines to be fit togheter
	  lineset_max: 0L,                 $  ; The number of linesets in the current database of lines
	  redshift: 0D,                    $  ; user defined redshift of the galaxy
          selected_lineset: 0L,            $  ; use a subset only of all [DEFAULT] available lines
          chisq: ptr_new(),                $  ; pointer to the mask chi^2 values 
          sp_chisq: ptr_new(),             $  ; pointer to the spaxel chi^2 values 
          linefit_mode: 0B,                $  ; 0=spaxels, 1=masks
          linefit_type: 0B,                $  ; 0=gaussian fit, 1=moments
          mom_windowmap: ptr_new(),        $  ; pointer to datacube with 2 planes (images) epf v,dv, defined to set the centre and width of the moments window
          gauss_initmap: ptr_new(),        $  ; pointer to datacube with each plane used to define the starting estimate of each parameter.
	  mom_windowmap_bootstrap: replicate(ptr_new(),1000), $  ; pointer to bootstrap cubes, one for each bootstrap equivalent to mom_windowmap 
          mask_sn_thresh: 3.0D,            $  ; max S/N of line flux to mask out with AUTOFLAG
          mask_maxvelerr: 50.0D,           $  ; max velocity (or dispersion) error to allow without masking out using AUTOFLAG
	  mom_thresh: 0.5D,                $  ; min S/N of flux values to be used in the computation of moments
          Nbootstrap: 0L,                  $  ; Number of bootstrap cubes
          inbootstrapcubes: replicate(ptr_new(),1000), $ ; array of pointers to unsmoothed bootstrap cubes
          bootstrapcubes: replicate(ptr_new(),1000),   $ ; array of pointers to cubes, current max 1000
          bootstrapnoise: ptr_new(),       $  ; pointer to the noise cube from bootstrap cubes
          Nmc1: 0L,                        $  ; number of iterations for MonteCarlo1 method 
          inmc1cubes: replicate(ptr_new(),1000), $ ; array of pointers to mc1cubes, current max 1000
          mc1cubes: replicate(ptr_new(),1000),   $ ; array of pointers to mc1cubes, current max 1000
          mc1noise: ptr_new(),             $  ; pointer to the noise cube from mc1cubes
          Nmc2: 0L,                        $  ; number of iterations for MonteCarlo2 method
          inmc2cubes: replicate(ptr_new(),1000), $ ; array of pointers to mc2cubes, current max 1000
          mc2cubes: replicate(ptr_new(),1000),   $ ; array of pointers to mc2cubes, current max 1000
          mc2noise: ptr_new(),             $  ; pointer to the noise cube from mc2cubes
          Nmc3: 0L,                        $  ; number of iterations for MonteCarlo3 method 
          inmc3cubes: replicate(ptr_new(),1000), $ ; array of pointers to mc3cubes, current max 1000
          mc3cubes: replicate(ptr_new(),1000),   $ ; array of pointers to mc3cubes, current max 1000
          mc3noise: ptr_new(),             $  ; pointer to the noise cube from mc3cubes
          noi_sp_nrescube: ptr_new(),      $ ; pointer to the spaxel linefit results data cube for narrow lines
          noi_sp_brescube: ptr_new(),      $ ; pointer to the spaxel linefit results data cube for broad lines
          noi_sp_crescube: ptr_new(),      $ ; pointer to the spaxel linefit results data cube for underlying continuum
          noi_sp_mrescube: ptr_new(),      $ ; pointer to the spaxel linefit results data cube for moments
          noi_sp_nerrrescube: ptr_new(),   $ ; pointer to the errors on spaxel linefit results data cube for narrow lines
          noi_sp_berrrescube: ptr_new(),   $ ; pointer to the errors on spaxel linefit results data cube for broad lines
          noi_sp_cerrrescube: ptr_new(),   $ ; pointer to the errors on spaxel linefit results data cube for underlying continuum
          noi_sp_merrrescube: ptr_new(),   $ ; pointer to the errors on spaxel linefit results data cube for moments
          noi_nrescube: ptr_new(),         $ ; pointer to the mask linefit results data cube for narrow lines
          noi_brescube: ptr_new(),         $ ; pointer to the mask linefit results data cube for broad lines
          noi_crescube: ptr_new(),         $ ; pointer to the mask linefit results data cube for underlying continuum
          noi_mrescube: ptr_new(),         $ ; pointer to the mask linefit results data cube for moments
          noi_nerrrescube: ptr_new(),      $ ; pointer to the errors on mask linefit results data cube for narrow lines
          noi_berrrescube: ptr_new(),      $ ; pointer to the errors on mask linefit results data cube for broad lines
          noi_cerrrescube: ptr_new(),      $ ; pointer to the errors on mask linefit results data cube for underlying continuum
          noi_merrrescube: ptr_new(),      $ ; pointer to the errors on mask linefit results data cube for moments
          boot_sp_nrescube: ptr_new(),     $ ; pointer to the spaxel linefit results data cube for narrow lines
          boot_sp_brescube: ptr_new(),     $ ; pointer to the spaxel linefit results data cube for broad lines
          boot_sp_crescube: ptr_new(),     $ ; pointer to the spaxel linefit results data cube for underlying continuum
          boot_sp_mrescube: ptr_new(),     $ ; pointer to the spaxel linefit results data cube for moments
          boot_sp_nerrrescube: ptr_new(),  $ ; pointer to the errors on spaxel linefit results data cube for narrow lines
          boot_sp_berrrescube: ptr_new(),  $ ; pointer to the errors on spaxel linefit results data cube for broad lines
          boot_sp_cerrrescube: ptr_new(),  $ ; pointer to the errors on spaxel linefit results data cube for underlying continuum
          boot_sp_merrrescube: ptr_new(),  $ ; pointer to the errors on spaxel linefit results data cube for moments
          boot_nrescube: ptr_new(),        $ ; pointer to the mask linefit results data cube for narrow lines
          boot_brescube: ptr_new(),        $ ; pointer to the mask linefit results data cube for broad lines
          boot_crescube: ptr_new(),        $ ; pointer to the mask linefit results data cube for underlying continuum
          boot_mrescube: ptr_new(),        $ ; pointer to the mask linefit results data cube for moments
          boot_nerrrescube: ptr_new(),     $ ; pointer to the errors on mask linefit results data cube for narrow lines
          boot_berrrescube: ptr_new(),     $ ; pointer to the errors on mask linefit results data cube for broad lines
          boot_cerrrescube: ptr_new(),     $ ; pointer to the errors on mask linefit results data cube for underlying continuum
          boot_merrrescube: ptr_new(),     $ ; pointer to the errors on mask linefit results data cube for moments
          mc1_sp_nrescube: ptr_new(),      $ ; pointer to the spaxel linefit results data cube for narrow lines
          mc1_sp_brescube: ptr_new(),      $ ; pointer to the spaxel linefit results data cube for broad lines
          mc1_sp_crescube: ptr_new(),      $ ; pointer to the spaxel linefit results data cube for underlying continuum
          mc1_sp_mrescube: ptr_new(),      $ ; pointer to the spaxel linefit results data cube for moments
          mc1_sp_nerrrescube: ptr_new(),   $ ; pointer to the errors on spaxel linefit results data cube for narrow lines
          mc1_sp_berrrescube: ptr_new(),   $ ; pointer to the errors on spaxel linefit results data cube for broad lines
          mc1_sp_cerrrescube: ptr_new(),   $ ; pointer to the errors on spaxel linefit results data cube for underlying continuum
          mc1_sp_merrrescube: ptr_new(),   $ ; pointer to the errors on spaxel linefit results data cube for moments
          mc1_nrescube: ptr_new(),         $ ; pointer to the mask linefit results data cube for narrow lines
          mc1_brescube: ptr_new(),         $ ; pointer to the mask linefit results data cube for broad lines
          mc1_crescube: ptr_new(),         $ ; pointer to the mask linefit results data cube for underlying continuum
          mc1_mrescube: ptr_new(),         $ ; pointer to the mask linefit results data cube for moments
          mc1_nerrrescube: ptr_new(),      $ ; pointer to the errors on mask linefit results data cube for narrow lines
          mc1_berrrescube: ptr_new(),      $ ; pointer to the errors on mask linefit results data cube for broad lines
          mc1_cerrrescube: ptr_new(),      $ ; pointer to the errors on mask linefit results data cube for underlying continuum
          mc1_merrrescube: ptr_new(),      $ ; pointer to the errors on mask linefit results data cube for moments
          mc2_sp_nrescube: ptr_new(),      $ ; pointer to the spaxel linefit results data cube for narrow lines
          mc2_sp_brescube: ptr_new(),      $ ; pointer to the spaxel linefit results data cube for broad lines
          mc2_sp_crescube: ptr_new(),      $ ; pointer to the spaxel linefit results data cube for underlying continuum
          mc2_sp_mrescube: ptr_new(),      $ ; pointer to the spaxel linefit results data cube for moments
          mc2_sp_nerrrescube: ptr_new(),   $ ; pointer to the errors on spaxel linefit results data cube for narrow lines
          mc2_sp_berrrescube: ptr_new(),   $ ; pointer to the errors on spaxel linefit results data cube for broad lines
          mc2_sp_cerrrescube: ptr_new(),   $ ; pointer to the errors on spaxel linefit results data cube for underlying continuum
          mc2_sp_merrrescube: ptr_new(),   $ ; pointer to the errors on spaxel linefit results data cube for moments
          mc2_nrescube: ptr_new(),         $ ; pointer to the mask linefit results data cube for narrow lines
          mc2_brescube: ptr_new(),         $ ; pointer to the mask linefit results data cube for broad lines
          mc2_crescube: ptr_new(),         $ ; pointer to the mask linefit results data cube for underlying continuum
          mc2_mrescube: ptr_new(),         $ ; pointer to the mask linefit results data cube for moments
          mc2_nerrrescube: ptr_new(),      $ ; pointer to the errors on mask linefit results data cube for narrow lines
          mc2_berrrescube: ptr_new(),      $ ; pointer to the errors on mask linefit results data cube for broad lines
          mc2_cerrrescube: ptr_new(),      $ ; pointer to the errors on mask linefit results data cube for underlying continuum
          mc2_merrrescube: ptr_new(),      $ ; pointer to the errors on mask linefit results data cube for moments
          mc3_sp_nrescube: ptr_new(),      $ ; pointer to the spaxel linefit results data cube for narrow lines
          mc3_sp_brescube: ptr_new(),      $ ; pointer to the spaxel linefit results data cube for broad lines
          mc3_sp_crescube: ptr_new(),      $ ; pointer to the spaxel linefit results data cube for underlying continuum
          mc3_sp_mrescube: ptr_new(),      $ ; pointer to the spaxel linefit results data cube for moments
          mc3_sp_nerrrescube: ptr_new(),   $ ; pointer to the errors on spaxel linefit results data cube for narrow lines
          mc3_sp_berrrescube: ptr_new(),   $ ; pointer to the errors on spaxel linefit results data cube for broad lines
          mc3_sp_cerrrescube: ptr_new(),   $ ; pointer to the errors on spaxel linefit results data cube for underlying continuum
          mc3_sp_merrrescube: ptr_new(),   $ ; pointer to the errors on spaxel linefit results data cube for moments
          mc3_nrescube: ptr_new(),         $ ; pointer to the mask linefit results data cube for narrow lines
          mc3_brescube: ptr_new(),         $ ; pointer to the mask linefit results data cube for broad lines
          mc3_crescube: ptr_new(),         $ ; pointer to the mask linefit results data cube for underlying continuum
          mc3_mrescube: ptr_new(),         $ ; pointer to the mask linefit results data cube for moemnts
          mc3_nerrrescube: ptr_new(),      $ ; pointer to the errors on mask linefit results data cube for narrow lines
          mc3_berrrescube: ptr_new(),      $ ; pointer to the errors on mask linefit results data cube for broad lines
          mc3_cerrrescube: ptr_new(),      $ ; pointer to the errors on mask linefit results data cube for underlying continuum
          mc3_merrrescube: ptr_new()       $ ; pointer to the errors on mask linefit results data cube for moments
       
        }


winsize = { $
            xwin1: 0L,    $  ; xsize of spaxel window  
            ywin1: 0L,    $  ; ysize of spaxel window  
            xwin2: 500,   $  ; xsize of spectrum window       - hard coded but can be resized
            ywin2: 250,   $  ; ysize of spectrum window       - hard coded but can be resized
            xwin3: 500,   $  ; xsize of spectrum zoom window  - hard coded but can be resized
            ywin3: 250    $  ; ysize of spectrum zoom window  - hard coded but can be resized
            ;xwin5: 500,   $  ; xsize of linefit window       - Not required because it does 
            ;ywin5: 500    $  ; ysize of linefit window       - not contain draw_widget calls
          }

widgetids = { $
	      base1: 0L,                     $ ; ID of spaxel base
              base2: 0L,                     $ ; ID of spectrum base
              base3: 0L,                     $ ; ID of spectrum zoom base
              base5: 0L,                     $ ; ID of linefit control base
	      wid1: 0L, 		     $ ; ID of spaxel window 
              wid2: 0L, 		     $ ; ID of spectrum window
              wid3: 0L, 		     $ ; ID of spectrum zoom draw window
              wid4: 0L, 		     $ ; ID of the colourbar widget
              wid5: 0L, 		     $ ; ID of the scale (zcut) parameters draw widget
              draw2: 0L,                     $ ; ID of spectrum draw window (required for dyn resize)
	      draw3: 0L,                     $ ; ID of spectrum zoom draw window (required for dyn resize)
              colmin_id: 0L,		     $ ; ID of min color label
              colmax_id: 0L,		     $ ; ID of max color label
              slid_id: 0L,		     $ ; ID of the wavelength slider
              pixima_id: 0L,		     $ ; ID of image pixel label
	      pixphys_id: 0L,		     $ ; ID of physical pixel label
              wcs_id: 0L,		     $ ; ID of wcs label
              pixval_id: 0L,		     $ ; ID of pixel value label
              smooth_id: 0L,                 $ ; ID of the current smoothing label
	      cube_id: 0L,		     $ ; ID of cube view
              imask_id: 0L,		     $ ; ID of imask label
              imask_lf_id: 0L,  	     $ ; ID of linefit imask label
              stepmask_id: 0L,  	     $ ; ID of the prev/next mask buttons parent base
              zminsp_id: 0L,		     $ ; ID of zmin label in spec window
              zmaxsp_id: 0L,		     $ ; ID of zmax label in spec window
              zminim_id: 0L,		     $ ; ID of zmin label in image user params
              zmaxim_id: 0L,		     $ ; ID of zmax label in image user params
	      spatzoom_id: 0L,               $ ; ID of spatial zoom label
	      zoom_id: 0L,                   $ ; ID of zoom-range label
              imgmode_id: 0L,                $ ; ID of image mode label
              continuummode_button: 0L,      $ ; ID of button to select internal continuum mode
	      continuumfit_minoff_id: 0L,    $ ; ID for min cont region offset 
              continuumfit_maxoff_id:  0L,   $ ; ID for max cont region offset
              continuumfit_minperc_id: 0L,   $ ; ID for min cont region %ile for fitting
              continuumfit_maxperc_id: 0L,   $ ; ID for max cont region  %ile for fitting
              continuumfit_order_id:   0L,   $ ; ID for order of fit to continuum
              maxwoffb_id: 0L,               $ ; ID of widget for maxwoff blue
              maxwoffr_id: 0L,               $ ; ID of widget for maxwoff red
              instrres_mode_id: 0L,          $ ; ID of instrumental resolutiion mode
              instrres_id: 0L,               $ ; ID of instrumental resolution label
              nlamb_id: lonarr(25),	     $ ; IDs of narrow line centroids
              blamb_id: lonarr(25),	     $ ; IDs of broad line centroids
              gnfit_id: lonarr(25),	     $ ; IDs of narrow line linefit user-supplied param.
              gbfit_id: lonarr(25),	     $ ; IDs of broad line linefit user-supplied param.
              pnfit_id: lonarr(25),	     $ ; IDs of narrow line linefit fitted param.
              pbfit_id: lonarr(25),	     $ ; IDs of broad line linefit fitted param.
              pcfit_id: lonarr(25),	     $ ; IDs of continuum linefit fitted param.
              ncomp_id: lonarr(25),	     $ ; IDs of narrow line component label
              bcomp_id: lonarr(25),	     $ ; IDs of broad line  component label
              ccomp_id: lonarr(25),	     $ ; IDs of continuum  component label
              gnlims_id: lonarr(25,2),       $ ; IDs of narrow line linefit param limits
              gblims_id: lonarr(25,2),       $ ; IDs of broad line linefit param limits.
	      nfixbutton: lonarr(25),	     $ ; IDs of narrow line linefit fix buttons
              bfixbutton: lonarr(25),	     $ ; IDs of broad line linefit fix buttons
              nresetbutton: lonarr(25),      $ ; IDs of narrow line reset buttons
              bresetbutton: lonarr(25),      $ ; IDs of broad line reset buttons
              nimagebutton: lonarr(25),      $ ; IDs of narrow line image buttons
              bimagebutton: lonarr(25),      $ ; IDs of broad line image buttons
              cimagebutton: lonarr(25),      $ ; IDs of continuum image buttons
              pndofitbutton: lonarr(25),     $ ; IDs of narrow line linefit dofit buttons
              pbdofitbutton: lonarr(25),     $ ; IDs of broad line linefit dofit buttons
              pcdofitbutton: lonarr(25),     $ ; IDs of buttons for continuum fit selection
              nshowbutton: lonarr(25),       $ ; IDs of narrow line show buttons
              bshowbutton: lonarr(25),       $ ; IDs of broad line show buttons
              cshowbutton: lonarr(25),       $ ; IDs of continuum show buttons
              flag_button: 0L,               $ ; ID of good/bad fit flag button
              flag_imagebutton: 0L,          $ ; ID of flag image button
              chisq_imagebutton: 0L,         $ ; ID of chisq image button
	      fitconstrbutton: 0L,	     $ ; ID of button for fit with constraints
              fitconstrlabel: 0L,	     $ ; ID of label for fit with constraints
              fixratiosbutton: 0L,	     $ ; ID of button for tie line ratios
              fixratioslabel: 0L,	     $ ; ID of label for tie line ratios
              secondcomplabel: 0L,           $ ; ID of label defining second component
	      secondcompbuttons: lonarr(3),  $ ; ID of buttons defining second component
	      secondcompsmartlabel : 0L,     $ ; ID of label defining the smart rejection os 2nd component
              secondcompsmartbutton : 0L,    $ ; ID of button defining the smart rejection os 2nd component
	      columns_id : lonarr(25),       $ ; IDs of column descriptions
              fitoption_id: lonarr(25),      $ ; IDs of the buttons in the last row of the linefit window
              redshift_id: 0L,               $ ; ID of redshift label
              chisq_id: 0L,                  $ ; ID of chi-sq label
              linefit_mode_id: 0L,	     $ ; ID of linefit mode label
              linefit_type_id: 0L,	     $ ; ID of linefit type label
              errormethod_id: 0L,	     $ ; ID of errormethod label.
              montecarloplot_id: 0L,	     $ ; ID of montecarlo plots label.
              montecarlosave_id: 0L,	     $ ; ID of montecarlo save label.
	      montecarlonoise_id: 0L,	     $ ; ID of montecarlo noise label.
	      scalenoise_id: 0L,	     $ ; ID of scale noise label.
	      inputmap_flag: 0L,    	     $ ; ID of input map label.
	      masksnthresh_id: 0L,	     $ ; ID of sn threshold for autoflag
              maskmaxvelerr_id: 0L,	     $ ; ID of max vel error for autoflag
              mom_thresh_id: 0L,             $ ; ID of min SN for moments
	      cancelfit_id: 0L,              $ ; ??
	      polycoeffres_id: lonarr(7)     $ ; IDs for polynomial coefficients for instrumental resolution

            }

fit = { $
      spec: ptr_new(),      	   $  ; Pointer to spectrum
      noise: ptr_new(),     	   $  ; Pointer to noise spectrum
      wave: ptr_new(),      	   $  ; Pointer to wavelength axis
      weights: ptr_new(),   	   $  ; Pointer to fit weights
      spec_bootstrap: ptr_new(),   $  ; Pointer to bootstrap spectra
      mfitpars: intarr(25), 	   $  ; binary array of lines to evaluate moments (1=fit)
      nfitpars: intarr(25), 	   $  ; binary array of narrow line parameters to fit (1=fit)
      bfitpars: intarr(25), 	   $  ; binary array of broad line parameters to fit (1=fit)
      cfitpars: 0B,         	   $  ; binary value for continuum to fit?
      cpar: 0,              	   $  ; index of continuum pars array
      nparstart: -1,        	   $  ; index of first narrow line in pars array
      nparend: -1,          	   $  ; index of last narrow line in pars array
      bparstart:  -1,       	   $  ; index of first broad line in pars array
      bparend:  -1,         	   $  ; index of last broad line in pars array
      contfitregion: ptr_new(),    $  ; array of indices for continuum fitting region
      okcontfit: ptr_new(), 	   $  ; array of wavelength indices to fit continuum
      contpars: ptr_new(),  	   $  ; continuum fit parameters (results)
      contpars_bootstrap: ptr_new(),  $  ; Pointer to continuum fit parameter arrays based on bootstrap cubes
      sigcontpars: ptr_new(),      $  ; errors on contpars (formal errors from the fit)
      sigcont: ptr_new(),          $  ; error on the continuum value (from standard deviation of residuals)
      pars: ptr_new(),             $  ; results array with all line-fitting parameters
      sigpars: ptr_new(),          $  ; errors on result parameters
      pars_bootstrap: ptr_new(),   $  ; Pointer to parameter arrays based on bootstrap cubes
      gpars: ptr_new(),    	   $  ; initial guess of all parameters
      parinfo: ptr_new(),  	   $  ; parinfo structure to be passed to mpcurvefit
      chisq: 0D,           	   $  ; best chi^2 fit result
      maxwoffb: 0D,        	   $  ; max wavelength offset from lineset to fit (Angstroms)
      maxwoffr: 0D,        	   $  ; max wavelength offset from lineset to fit (Angstroms)
      okk: ptr_new(),      	   $  ; array of wavelength indices to fit line
      specfit: ptr_new(),  	   $  ; fit to spectrum
      status: 0L,          	   $  ; fit status parameter
      errmsg: "",          	   $  ; error message parameter
      iter: 0L,             	   $  ; number of iterations for fit
      firstset: 0B                 $  ; is this the first lineset with useful lines when we fit multiple sets?
      }

linesdb = { $
          lyalpha: 1215.670D,   $         ;ALL THOSE LINES ARE AIR WAVE
          halpha : 6562.819D,   $
          n2_b   : 6548.050D,   $
          n2_r   : 6583.450D,   $
          hbeta  : 4861.333D,   $
          o3_b   : 4958.911D,   $
          o3_r   : 5006.843D,   $
          s2_b   : 6716.440D,   $
          s2_r   : 6730.810D,   $
          o1_b   : 6300.304D,   $          ; table 13.2 Osterbrock
          o1_r   : 6363.776D,   $          ; table 13.2 Osterbrock
          he1_b  : 5875.621D,   $
	  he1_r  : 6678.152D,   $
          s3_b   : 9068.600D,   $
	  s3_r   : 9530.600D,   $
	  o2_b   : 3726.030D,   $
	  o2_r   : 3728.820D,   $
	  linesets   : fltarr(25), $    
	  lines_rest : fltarr(25), $      ; array of all lines
          linenames  : strarr(25), $      ; array of line names
          linefancynames : strarr(25) $
          }

linesdb.linesets = [1,1,1,2,2,3,3,4,5,5,6,6,7,7,8,8,9] ;Associate to each line its corresponding lineset
linesdb.lines_rest = [linesdb.halpha,linesdb.n2_b,linesdb.n2_r,linesdb.o3_b,linesdb.o3_r, $
                      linesdb.o2_b,linesdb.o2_r,linesdb.hbeta, $
                      linesdb.s2_b,linesdb.s2_r,linesdb.o1_b,linesdb.o1_r, linesdb.s3_b, linesdb.s3_r,  $
		      linesdb.he1_b, linesdb.he1_r, linesdb.lyalpha] ;array of all lines
		      
linesdb.linenames = ['Ha', 'n2_b', 'n2_r', 'o3_b', 'o3_r','o2_b', 'o2_r','Hb','s2_b', 's2_r','o1_b', 'o1_r','s3_b','s3_r','He1_b','He1_r','Lya']

linesdb.linefancynames = ['Ha', '[NII] blue', '[NII] red', '[OIII] blue', '[OIII] red', '[OII] blue', '[OII] red','Hb','[SII] blue', $
                         '[SII] red','[OI] blue', '[OI] red', '[SIII] blue', '[SIII] red', 'HeI blue', 'HeI red', 'Lya']

pars = { $
         gotomask  : 0L,      $
	 spatsmooth: 0L,      $ 
	 specsmooth: 0L,      $
	 zmin_ima  : 0D,      $ 
	 zmax_ima  : 0D,      $
	 selhead   : -1L,     $
	 maskradius: 0L,      $
	 maskmode  : 0L,      $
	 fitallrange : replicate(0,4), $
	 Nloopadj  : 0L       $
       }


end

; --------------------------------------------------------------------------------

pro kubeviz_linefit_define_flags
common flags, flag_ok, flag_bad
  
flag_bad =      [                                               $
                [255B, 255B, 255B, 255B],                       $
                [255B, 255B, 255B, 255B],                       $
                [015B, 060B, 014B, 252B],                       $
                [015B, 028B, 012B, 248B],                       $
                [207B, 201B, 201B, 249B],                       $
                [207B, 201B, 201B, 249B],                       $
                [207B, 201B, 201B, 249B],                       $
                [015B, 008B, 200B, 249B],                       $
                [015B, 008B, 200B, 249B],                       $
                [207B, 201B, 201B, 249B],                       $
                [207B, 201B, 201B, 249B],                       $
                [207B, 201B, 201B, 249B],                       $
                [015B, 200B, 009B, 248B],                       $
                [015B, 204B, 009B, 252B],                       $
                [255B, 255B, 255B, 255B],                       $
                [255B, 255B, 255B, 255B]                        $
                ]

flag_ok =       [                                               $
                [000B, 000B, 000B, 000B],                       $
                [192B, 143B, 049B, 006B],                       $
                [224B, 159B, 049B, 006B],                       $
                [096B, 152B, 049B, 006B],                       $
                [096B, 152B, 057B, 006B],                       $
                [096B, 152B, 031B, 006B],                       $
                [096B, 152B, 007B, 006B],                       $
                [096B, 152B, 007B, 006B],                       $
                [096B, 152B, 015B, 006B],                       $
                [096B, 152B, 029B, 006B],                       $
                [096B, 152B, 057B, 006B],                       $
                [096B, 152B, 049B, 000B],                       $
                [096B, 152B, 049B, 000B],                       $
                [224B, 159B, 049B, 006B],                       $
                [192B, 143B, 049B, 006B],                       $
                [000B, 000B, 000B, 000B]                        $
                ]

end

function kubeviz_excludewave, wave
; evaluate indices for which wave vector lies in hardcoded "bad"
; regions of the spectrum which are to be excluded.

Nwave = N_Elements(wave)
bad = replicate(0B,Nwave)

; Atmospheric O2 feature, leaves strong residuals: 12680-12710 AA 
bad_start = [12680.]
bad_end = [12710.]
Nbad = N_Elements(bad_start)
if N_Elements(bad_end) ne Nbad then begin
    print, 'ERROR: MISMATCHED NUMBER OF ELEMENTS IN BAD_START AND BAD_END ARRAYS WHEN DEFINING BAD SPECTRAL REGIONS!'
    stop
endif

for i=0, Nbad-1 do begin
    good = (1-bad)*(wave lt bad_start[i] or wave gt bad_end[i])
    bad = (1-good)
endfor

return, bad

end

;---------------------------------------------------------
function kubeviz_resultsstring, value, error, format=format

if N_Elements(format) eq 0 then begin
   format1='(F9.2)'
   format2='(F8.2)'
endif else begin
   format1=format
   format2=format
endelse   

if value gt 1E6 then begin
format1='(F9.1)'
format2='(F8.1)'
endif

case error[0] of 
 -999 : resultstr = 'Not Fit' 
 -998 : resultstr =  kubeviz_str(value,format=format1)+' +/- No Errors' 
 else : begin
        case n_elements(error) of
        1: resultstr = kubeviz_str(value,format=format1)+'+/-'+kubeviz_str(error,format=format2) 
        2: resultstr = kubeviz_str(value,format=format1)+'+'+kubeviz_str(error[0],format=format2)+'/'+kubeviz_str(error[1],format=format2)
        else: resultstr =  'BUG: Wrong number of elements in error for results string!'
       endcase
       end   
endcase

return, resultstr
end

;-----------------------------------------------------------
function kubeviz_userparstring, value, format=format
if N_Elements(format) eq 0 then format='(F7.2)'

if value eq -999 then parstr = 'Not Set' else parstr = kubeviz_str(value,format=format)

return, parstr
end

;-----------------------------------------------------------
function kubeviz_getsnimage, im, nim, abs=abs

sn_image = (N_elements(abs) eq 0) ? im/nim : abs(im/nim)
sn_image = kubeviz_remove_badvalues(sn_image)

return, sn_image
end

;-----------------------------------------------------------
function kubeviz_decodetrimstr, trimstr, fcbaxis
common kubeviz_state

;First process the trim string to make it understood by IDL
x_size = fcbaxis[0]
y_size = fcbaxis[1]
z_size = fcbaxis[2]
oktrim = 1

xytrim = strsplit(strmid(trimstr, 1, strlen(trimstr)-2), ',', /extract)
if N_elements(xytrim) eq 3 then begin
  ; Start with x
  case 1 of 
    strlen(xytrim[0]) gt 1 and n_elements(strsplit(xytrim[0], ':', /extract)) eq 2 : begin
  	temp = fix(strsplit(xytrim[0], ':', /extract))
  	if temp[0] ge 1 and temp[1] le x_size then  xtrim = [(temp[0]-1), (temp[1]-1)] else oktrim = 0
    end 
    strlen(xytrim[0]) eq 1 and xytrim[0] eq '*' : xtrim = [0,x_size-1]
    else: oktrim = 0
  endcase
  ; Then process y
  case 1 of 
    strlen(xytrim[1]) gt 1 and n_elements(strsplit(xytrim[1], ':', /extract)) eq 2 : begin
  	temp = fix(strsplit(xytrim[1], ':', /extract))
  	if temp[0] ge 1 and temp[1] le y_size then ytrim = [(temp[0]-1), (temp[1]-1)] else oktrim = 0
    end 
    strlen(xytrim[1]) eq 1 and xytrim[1] eq '*' : ytrim = [0,y_size-1]
    else: oktrim = 0
  endcase
  ; And then z
  case 1 of 
    strlen(xytrim[2]) gt 1 and n_elements(strsplit(xytrim[2], ':', /extract)) eq 2 : begin
  	temp = fix(strsplit(xytrim[2], ':', /extract))
  	if temp[0] ge 1 and temp[1] le z_size then ztrim = [(temp[0]-1), (temp[1]-1)] else oktrim = 0
    end 
    strlen(xytrim[2]) eq 1 and xytrim[2] eq '*' : ztrim = [0,z_size-1]
    else: oktrim = 0
  endcase
endif else oktrim = 0

if oktrim eq 1 then begin 
  state.Startcol  = xtrim[0]
  state.Startrow  = ytrim[0]
  state.Startwpix = ztrim[0]
  return, [xtrim,ytrim,ztrim]
endif else begin
  printf, state.log_lun, '[WARNING] Trim keyword syntax not understood. No trim will take place'
  return, replicate(-1, 6) 
endelse

end

;-----------------------------------------------------------
pro kubeviz_statusline, step_count, abs_step, t_start, close=close
common kubeviz_state

if keyword_set(close) then begin
 kubeviz_statuslinecore, /close
 printf, state.log_lun, ''
 return
endif

step_count += 1

if step_count mod abs_step eq 0 then begin
   curr_percent = state.percent_step*step_count/abs_step
   kubeviz_statuslinecore,  /clear
   t_here = systime(/seconds)
   if curr_percent lt 100 then begin
    eta = ((t_here-t_start)/(curr_percent)*(100-curr_percent))
    format='(I02)'
    case 1 of
      eta le 60 : etastr=kubeviz_str(fix(eta))+'s'
      eta gt 60 and eta le 3600 : etastr=kubeviz_str(fix(eta/60))+':'+kubeviz_str(fix(eta mod 60), format=format)+' mm:ss'
      eta gt 3600 : etastr=kubeviz_str(fix(eta/3600))+':'+kubeviz_str(fix((eta mod 3600)/60), format=format)+':'+kubeviz_str(fix((eta mod 3600) mod 60), format=format)+' hh:mm:ss'  
    endcase
    kubeviz_statuslinecore, '[PROGRES] '+kubeviz_str(curr_percent)+'% processed - ETA: '+etastr,0 
   endif else begin
    eta = ((t_here-t_start))
    format='(I02)'
    case 1 of
      eta le 60 : etastr=kubeviz_str(fix(eta))+'s'
      eta gt 60 and eta le 3600 : etastr=kubeviz_str(fix(eta/60))+':'+kubeviz_str(fix(eta mod 60), format=format)+' mm:ss'
      eta gt 3600 : etastr=kubeviz_str(fix(eta/3600))+':'+kubeviz_str(fix((eta mod 3600)/60), format=format)+':'+kubeviz_str(fix((eta mod 3600) mod 60), format=format)+' hh:mm:ss'  
    endcase
    kubeviz_statuslinecore, '[PROGRES] 100% processed - Time elapsed: '+etastr,0  
   endelse 
endif

end

;-----------------------------------------------------------
function kubeviz_timediff, datetime1, datetime2, ret_unit = ret_unit

if N_elements(ret_unit) eq 0 then ret_unit = 's'

y1 = fix(strmid(datetime1, 0, 4))
m1 = fix(strmid(datetime1, 5, 2))
d1 = fix(strmid(datetime1, 8, 2))
hh1 = fix(strmid(datetime1, 11, 2))
mm1 = fix(strmid(datetime1, 14, 2))
ss1 = fix(strmid(datetime1, 17, 2))

y2 = fix(strmid(datetime2, 0, 4))
m2 = fix(strmid(datetime2, 5, 2))
d2 = fix(strmid(datetime2, 8, 2))
hh2 = fix(strmid(datetime2, 11, 2))
mm2 = fix(strmid(datetime2, 14, 2))
ss2 = fix(strmid(datetime2, 17, 2))

juldate, [y1, m1, d1, hh1, mm1, ss1], jd1
juldate, [y2, m2, d2, hh2, mm2, ss2], jd2

diff_s = (jd2-jd1)

case ret_unit of
  's': timediff = (jd2-jd1)*3600*24
  'm': timediff = (jd2-jd1)*60*24
  'h': timediff = (jd2-jd1)*24
  else: begin
   printf, state.log_lun, '[WARNING] Unrecognized time format. Returning the time difference in seconds instead.
   timediff = (jd2-jd1)*3600*24
  end
endcase

return, timediff

end


;---------------------------------------------------------
function kubeviz_getmainline, redshift=redshift, vacuum=vacuum, linesarr=linesarr
common kubeviz_state, state
common kubeviz_linesdb

; default to selected lineset
if N_Elements(redshift) eq 0 then redshift = state.redshift
if N_elements(vacuum)   eq 0 then vacuum   = state.vacuum 
if N_elements(linesarr) eq 0 and      ptr_valid(state.linenames) then linesarr = (*state.linenames)
if N_elements(linesarr) eq 0 and  not ptr_valid(state.linenames) then return, linesdb.halpha 

if state.Nlines eq 0 then mainline = linesdb.halpha else begin
     case 1 of 
       total(strmatch(linesarr, 'Ha'))    eq 1 : mainline = linesdb.halpha 
       total(strmatch(linesarr, 'o3_r'))  eq 1 : mainline = linesdb.o3_r
       total(strmatch(linesarr, 'o2_r'))  eq 1 : mainline = linesdb.o2_r
       total(strmatch(linesarr, 'Hb'))    eq 1 : mainline = linesdb.hbeta 
       total(strmatch(linesarr, 's2_r'))  eq 1 : mainline = linesdb.s2_r 
       total(strmatch(linesarr, 'o1_b'))  eq 1 : mainline = linesdb.o1_b
       total(strmatch(linesarr, 'He1_b')) eq 1 : mainline = linesdb.he1_b
       total(strmatch(linesarr, 's3_b'))  eq 1 : mainline = linesdb.s3_b
       total(strmatch(linesarr, 'Lya'))   eq 1 : mainline = linesdb.lyalpha
       else: mainline = linesdb.halpha
     endcase
endelse  

if vacuum eq 1 then begin
  airtovac, mainline, mainline_vac
  mainline = mainline_vac
endif 

return, mainline*(1.+redshift)

end

;----------------------------------------------------------------------
function kubeviz_moment_reshape, momentcube, error=error
;Reshapes the moment results cube such that it can be used in the linefit_image_update routine
common kubeviz_state

if N_elements(error) eq 0 then error = 0 else error = 1
main_ind = where(abs(state.lines - kubeviz_getmainline()) lt 1)

if error eq 0 then begin
  rescube = replicate(0.D, state.Ncol, state.Nrow , 3+state.Nlines)
  rescube[*,*,0] = momentcube[*,*,main_ind*6+5]
  rescube[*,*,1] = momentcube[*,*,main_ind*6+1]
  rescube[*,*,2] = momentcube[*,*,main_ind*6+2]
  for line= 0 , state.Nlines-1 do rescube[*,*,line+3] = momentcube[*,*,line*6]
  return, rescube
endif else begin
  errrescube = replicate(0.D, state.Ncol, state.Nrow , 3+state.Nlines, state.Nmontecarlo_percs)
  errrescube[*,*,0,*] = momentcube[*,*,main_ind*6+5,*]
  errrescube[*,*,1,*] = momentcube[*,*,main_ind*6+1,*]
  errrescube[*,*,2,*] = momentcube[*,*,main_ind*6+2,*]
  for line= 0 , state.Nlines-1 do errrescube[*,*,line+3,*] = momentcube[*,*,line*6,*]
  return, errrescube
endelse

end

;----------------------------------------------
function kubeviz_getinstrres, lambda=lambda
common kubeviz_state
; get the instrumental resolution at lambda:
; resolution in "R" units, d(lambda)/lambda

if N_elements(lambda) eq 0 then lambda = kubeviz_getmainline()
instrres = 0.

if state.Nlines gt 0 then begin
  case state.instrres_mode of
     0: for i=0, state.max_polycoeff_instrres do instrres += state.instrres_varpoly[i]*lambda^i ;Variance fit
     1: for i=0, state.max_polycoeff_instrres do instrres += state.instrres_extpoly[i]*lambda^i ;Polynomial fit
     2: instrres = lambda / (2.35482*state.instrres_tplsig)  ;Template fit
  endcase  
endif 

return, instrres
end


;------------------------------------------
function kubeviz_linefit_estimatekinematics, gpars, lambda0
common kubeviz_state
common kubeviz_fit

;MF. HERE WE USE THE MPFITPEAK METHOD FOR THE STARTING VALUES 
;WHICH PROVES TO BE VERY ROBUST
;Note: no need to do a continuum subtraction in this routine as 
;MPFIT takes care of it 

glambda = lambda0*(1.+gpars[fit.nparstart]/state.ckms)
checkrange = 40/(1.+state.redshift)   ; angstroms to avoid very large bins
nearpeak = where( abs((*fit.wave)-glambda) le checkrange, nx)

if nx lt 5 then goto, GUESSFAIL

x = (*fit.wave)[nearpeak]
y = (*fit.spec)[nearpeak]

is = sort(x)
xs = x[is] & ys = y[is]
maxx = max(xs, min=minx) & maxy = max(ys, min=miny, nan=nan)
dx = 0.5 * [xs[1]-xs[0], xs[2:*] - xs, xs[nx-1] - xs[nx-2]]
totarea = total(dx*ys, nan=nan)       ;; Total area under curve
av = totarea/(maxx - minx)  ;; Average height

;; Degenerate case: all flat with no noise
if miny EQ maxy then begin
    est = ys[0]*0.0 + [0,xs[nx/2],(xs[nx-1]-xs[0])/2, ys[0]]
    guess = 1
    return, est
endif

;; Compute the spread in values above and below average... we
;; take the narrowest one as the one with the peak
wh1 = where(y GE av, ct1)
wh2 = where(y LE av, ct2)
if ct1 EQ 0 OR ct2 EQ 0 then begin
    print, 'ERROR: average Y value should fall within the range of Y data values but does not'
    return, !values.d_nan
endif
sd1 = total(x[wh1]^2)/ct1 - (total(x[wh1])/ct1)^2
sd2 = total(x[wh2]^2)/ct2 - (total(x[wh2])/ct2)^2
    
;; Compute area above/below average
;; Search for a positive peak

cent  = x[where(y EQ maxy)] & cent = cent[0]
peak  = maxy - av
peakarea = totarea - total(dx*(ys<av), nan=nan)
if peak EQ 0 then peak = 0.5*peakarea
width = peakarea / (2*abs(peak))
if width EQ 0 OR finite(width) EQ 0 then width = median(dx)

est = [cent, width]
return, est

GUESSFAIL:
 cent  = -1
 width = -1
 if state.debug eq 1 then begin
   print, '[DEBUG] GUESS FOR SIGMA FAILED.'
   print, '[DEBUG] SIGMA AND CENTROID FIXED '
 endif
 return, [cent, width]

end

;--------------------------------------------------------------------------------
function kubeviz_get_residual_spec, col, row
common kubeviz_state
common kubeviz_fit

sigtofwhm = 2.35482

if N_elements(col) eq 0 then col = state.col
if N_elements(row) eq 0 then row = state.row

spec = (*state.datacube)[col,row,*]
residual=spec

; Is there a spectrum to fit?
;okspec = where(spec ne 0., Nokspec)
okspec = where(finite(spec) eq 1, Nokspec)
if Nokspec gt 0 then begin

     kubeviz_linefit_dofit, col=col, row=row, /nosave

     Nnfitpars = total(fit.nfitpars)
     Nbfitpars = total(fit.bfitpars) 
     Nnfitlines = max([Nnfitpars-2,0])
     Nbfitlines = max([Nbfitpars-2,0])

     if Nnfitlines gt 0 then begin

        ;subtract continuum at mainline
        lambda_line = state.lines[0]*(1.+(*fit.pars)[1]/state.ckms)
        residual -= kubeviz_getcontatlambda((*fit.contpars),lambda_line)
        xx = (*state.wave)[0:state.Nwpix-1]

        for i=0,Nnfitlines-1 do begin
             dv = (*fit.pars)[1]
             if dv eq -999. then dv = 0.
             pos = state.lines[i]*(1.+dv/state.ckms)
             sigv = (*fit.pars)[2]
             if sigv eq -999. then sigv = 0.
             width = state.lines[i]*(sigv/state.ckms)
             instrres_A = pos/(sigtofwhm*kubeviz_getinstrres(lambda=pos)) ; convert from R to Angstroms
             totwidthsq = width^2 + instrres_A^2
             norm = (*fit.pars)[i+3]
             gn = (norm/sqrt(2.*!DPI*totwidthsq))*exp(-((xx-pos)^2/totwidthsq)/2.)
             residual -= gn
        endfor
     endif 

     if Nbfitlines gt 0 then begin
         for i=0,Nbfitlines-1 do begin
             dv = (*fit.pars)[Nnfitlines+3]
             if dv eq -999. then dv = 0.
             pos = state.lines[i]*(1.+dv/state.ckms)
             sigv = (*fit.pars)[Nnfitlines+4]
             if sigv eq -999. then sigv = 0.
             width = state.lines[i]*(sigv/state.ckms)
             instrres_A = pos/(sigtofwhm*kubeviz_getinstrres(lambda=pos)) ; convert from R to Angstroms
             totwidthsq = width^2 + instrres_A^2
             norm = (*fit.pars)[Nnfitlines+i+5]
             gn = (norm/sqrt(2.*!DPI*totwidthsq))*exp(-((xx-pos)^2/totwidthsq)/2.)
             residual -= gn
        endfor
     endif
     ;Free the fit pointers 
     heap_free, fit
endif

return, residual
end

;--------------------------------------------------------------------------------
function kubeviz_linefit_getredshift
common kubeviz_state
common kubeviz_fit

sigtofwhm = 2.35482

; Save state of some fitting variables
montecarlo = state.domontecarlo
pndofit = state.pndofit
pcdofit = state.pcdofit
maxoffwb = state.maxwoffb
maxoffwr = state.maxwoffr
continuumfit_minoff = state.continuumfit_minoff 
continuumfit_maxoff = state.continuumfit_maxoff  
continuumfit_minperc = state.continuumfit_minperc   
continuumfit_maxperc = state.continuumfit_maxperc   
continuumfit_order  = state.continuumfit_order   

; Set the default parameters for this automatic fit
for iline=0, state.Nlines-1 do begin
  state.pndofit[iline] = 1
  state.pbdofit[iline] = 0
  state.pcdofit[iline] = 1
endfor

state.domontecarlo = 0
kubeviz_set_continuumfit_defaults
kubeviz_set_linefitrange_defaults


kubeviz_medianspec, /sum, /redsh
spec = *state.medspec
noise = *state.nmedspec
wave = *state.wave
spec_bootstrap = replicate(-1.,1,state.Nwpix)
okspec = where(finite(spec) eq 1, Nokspec)
firstset=0

; save status of the fix button and start parameters for wavelength
; parameters (to be reset after fitting all lines):
pnfix = state.pnfix
gnfit = state.gnfit
pbfix = state.pbfix
gbfit = state.gbfit

if Nokspec gt 0 then begin
   for set=0, state.lineset_max-1 do kubeviz_linefit_fitset, wave(okspec), spec(okspec), noise(okspec), (*state.lineset_ind[set]),  spec_bootstrap[*,okspec], 0, firstset
endif

; reset the fix button and start parameters:
state.pnfix = pnfix
state.gnfit = gnfit
state.pbfix = pbfix
state.gbfit = gbfit

main_ind = (where(abs(state.lines - kubeviz_getmainline()) lt 1))[0]

fitres   =  state.lines[main_ind]*(1.+(*fit.pars)[1]/state.ckms)
if ( (*fit.pars)[3]/(*fit.sigpars)[3] gt state.mask_sn_thresh and (*fit.pars)[3] gt 0) $ 
   then redshift =  (fitres / kubeviz_getmainline(redshift=0)) - 1 else redshift = -1

;Free the fit pointers 
heap_free, fit

;Set original values of the parameters
state.domontecarlo = montecarlo 
state.pndofit =  pndofit 
state.pcdofit =  pcdofit 
state.maxwoffb = maxoffwb
state.maxwoffr = maxoffwr
state.continuumfit_minoff = continuumfit_minoff 
state.continuumfit_maxoff = continuumfit_maxoff 
state.continuumfit_minperc = continuumfit_minperc
state.continuumfit_maxperc = continuumfit_maxperc
state.continuumfit_order = continuumfit_order  

return, redshift
end

;---------------------------------------------------------------------------
function kubeviz_linefit_startminmaxconsistencycheck, guess, par, linetype, silent=silent
common kubeviz_state
common kubeviz_fit

if n_elements(silent) eq 0 then silent=0 else silent=1

setguess=guess
linefit_update=0

; consistency check that the start value is within the [min,max] range
; defined by the user. ONLY WHERE SET CONSTRAINTS BUTTON IS PRESSED.
if state.fitconstr eq 1 then begin
    case linetype of
        'N': begin
            if guess lt state.gnlims[par,0] and state.gnlims[par,0] ne -999. then begin
                if silent eq 0 then begin
		  printf, state.log_lun, '[WARNING] Narrow Line Parameter '+kubeviz_str(par,format='(I2)')+' initial guess value = '+kubeviz_str(guess,format='(f7.2)')+' is below user-defined min value of '+kubeviz_str(state.gnlims[par,0],format='(f7.2)')
                  printf, state.log_lun, '[WARNING] Initial guess set to '+kubeviz_str(state.gnlims[par,0],format='(f7.2)')
                endif
		setguess = state.gnlims[par,0]
                state.gnfit[par] = setguess
                if not silent then linefit_update = 1
            endif
            if guess gt state.gnlims[par,1] and state.gnlims[par,1] ne -999. then begin
                if silent eq 0 then begin
		  printf, state.log_lun, '[WARNING] Narrow Line Parameter '+kubeviz_str(par,format='(I2)')+' initial guess value = '+kubeviz_str(guess,format='(f7.2)')+' is above user-defined max value of '+kubeviz_str(state.gnlims[par,1],format='(f7.2)')
                  printf, state.log_lun, '[WARNING] Initial guess set to '+kubeviz_str(state.gnlims[par,1],format='(f7.2)')
                endif
		setguess = state.gnlims[par,1]
                state.gnfit[par] = setguess
                if not silent then linefit_update = 1
            endif
        end
        'B': begin
            if guess lt state.gblims[par,0] and state.gblims[par,0] ne -999. then begin
                if silent eq 0 then begin
		  printf, state.log_lun, '[WARNING] Broad Line Parameter '+kubeviz_str(par,format='(I2)')+' initial guess value = '+kubeviz_str(guess,format='(f7.2)')+' is below user-defined min value of '+string(state.gblims[par,0],format='(f7.2)')
                  printf, state.log_lun, '[WARNING] Initial guess set to '+string(state.gblims[par,0],format='(f7.2)')
                endif
		setguess = state.gblims[par,0]
                state.gbfit[par] = setguess
                if not silent then linefit_update = 1
            endif
            if guess gt state.gblims[par,1] and state.gblims[par,1] ne -999. then begin
                if silent eq 0 then begin
		  printf, state.log_lun, '[WARNING] Broad Parameter '+kubeviz_str(par,format='(I2)')+' initial guess value = '+kubeviz_str(guess,format='(f7.2)')+' is above user-defined max value of '+string(state.gblims[par,1],format='(f7.2)')
                  printf, state.log_lun, '[WARNING] Initial guess set to '+string(state.gblims[par,1],format='(f7.2)')
                endif
		setguess = state.gblims[par,1]
                state.gbfit[par] = setguess
                if not silent then linefit_update = 1
            endif
        end
    endcase
endif

if linefit_update then  kubeviz_linefit_update, /update_userpars

return, setguess
end

;------------------------------------------------------------------------
function kubeviz_linefit_Nokadjacentspax, x_fit, y_fit, flag 
common kubeviz_state

; get surrounding spaxel indices:
dx = [-1,0,1,1,1,0,-1,-1]
dy = [1,1,1,0,-1,-1,-1,0]
x_neighbours = x_fit + dx
y_neighbours = y_fit + dy

neighbours_ok = where(x_neighbours ge 0 and x_neighbours le state.Ncol-1 and y_neighbours ge 0 and y_neighbours le state.Nrow-1 and flag[x_neighbours,y_neighbours] eq 0, Nneighbours_ok)

return, Nneighbours_ok

end

;------------------------------------------------------------------------
function kubeviz_linefit_guessfromadjacentspax, x_fit, y_fit, flag, nrescube, brescube, nerrrescube, berrrescube
common kubeviz_state
; Take a guess at the parameters in a given spaxel given the fit
; parameters in adjacent spaxels with flag=ok
; Take the inverse-variance weighted mean of each parameter value

; We simply do this by populating the gnfit, gbfit arrays in the
; common state structure. These are then propogated into the fit in
; the usual way.

; get surrounding spaxel indices:
dx = [-1,0,1,1,1,0,-1,-1]
dy = [1,1,1,0,-1,-1,-1,0]
x_neighbours = x_fit + dx
y_neighbours = y_fit + dy

Npar = 3+state.Nlines

neighbours_ok = where(x_neighbours ge 0 and x_neighbours le state.Ncol-1 and y_neighbours ge 0 and y_neighbours le state.Nrow-1 and flag[x_neighbours,y_neighbours] eq 0, Nneighbours_ok)

if Nneighbours_ok gt 0 then begin
   x_neighbours_ok = x_neighbours[neighbours_ok]
   y_neighbours_ok = y_neighbours[neighbours_ok]
   for par=1,Npar-1 do begin ; par=0 is the flag
      if total(nerrrescube[x_neighbours_ok,y_neighbours_ok,replicate(par,Nneighbours_ok),replicate(0,Nneighbours_ok)] eq -999.) eq 0 then begin 
         values_ok = nrescube[x_neighbours_ok,y_neighbours_ok,replicate(par,Nneighbours_ok)]
	 err_ok    = 0.5*(abs(nerrrescube[x_neighbours_ok,y_neighbours_ok,replicate(par,Nneighbours_ok),replicate(0,Nneighbours_ok)])+abs(nerrrescube[x_neighbours_ok,y_neighbours_ok,replicate(par,Nneighbours_ok),replicate(1, Nneighbours_ok)]))
	 guess_adj = total(values_ok/(err_ok)^2, /nan)/total(1./(err_ok)^2, /nan)
	 if finite(guess_adj) $
	    then state.gnfit[par] = kubeviz_linefit_startminmaxconsistencycheck(guess_adj,par, 'N', /silent) $
	    else state.gnfit[par] = kubeviz_linefit_startminmaxconsistencycheck(0.,       par, 'N', /silent)
      endif  
      if total(berrrescube[x_neighbours_ok,y_neighbours_ok,replicate(par,Nneighbours_ok),replicate(0,Nneighbours_ok)] eq -999.) eq 0 then begin 
         values_ok = brescube[x_neighbours_ok,y_neighbours_ok,replicate(par,Nneighbours_ok)]
	 err_ok    = 0.5*(abs(berrrescube[x_neighbours_ok,y_neighbours_ok,replicate(par,Nneighbours_ok),replicate(0,Nneighbours_ok)])+abs(berrrescube[x_neighbours_ok,y_neighbours_ok,replicate(par,Nneighbours_ok),replicate(1, Nneighbours_ok)]))
	 guess_adj = total(values_ok/(err_ok)^2, /nan)/total(1./(err_ok)^2, /nan)
	 if finite(guess_adj) $
	     then state.gbfit[par] = kubeviz_linefit_startminmaxconsistencycheck(guess_adj,par, 'B', /silent) $
	     else state.gbfit[par] = kubeviz_linefit_startminmaxconsistencycheck(0.,       par, 'B', /silent)
      endif
   endfor
   dofit = 1
   
endif else dofit=0

return, dofit

end

;-----------------------------------------------------------
pro kubeviz_linefit_skylines
common kubeviz_state

printf, state.log_lun, '[KUBEVIZ] Computing the polynomial coefficients for the instrumental resolution.'
printf, state.log_lun, '[KUBEVIZ] This might take a few minutes, but will only happen once..'

findpro, 'kubeviz', dirlist=kubeviz_dir, /noprint 

ohspec_file = kubeviz_dir[0]+'/templates/ohspec/kmos_oh_spec.fits'

ohspec=readfits(ohspec_file,hdroh, /silent)
crvaloh = sxpar(hdroh,'CRVAL1') * 1.e4
cdeltoh = sxpar(hdroh,'CDELT1') * 1.e4
crpixoh = sxpar(hdroh,'CRPIX1')
naxisoh = sxpar(hdroh,'NAXIS1')
indxoh  = indgen(naxisoh)
lamboh = ((dindgen(naxisoh)+1-crpixoh)*cdeltoh+crvaloh)

;First check that at least half of the spectral range of the data is mapped in the ohspec
;if not estimate from variance spectrum, if yes cut the ohspec to the useful part

Noverlap = N_elements(kubeviz_SetIntersection((*state.wave),lamboh))
if Noverlap gt fix(state.Nwpix/2.) then begin
  ohspec = ohspec[where(lamboh gt (*state.wave)[0] and lamboh lt (*state.wave)[state.Nwpix-1])]
  lamboh = lamboh[where(lamboh gt (*state.wave)[0] and lamboh lt (*state.wave)[state.Nwpix-1])]
  
  nbins = 5
  zpos2 = 0
  binslim = long(1.*N_elements(lamboh)/nbins) 

  ;Find the lines above a local threshold (currently assumes maxval/50 and 5 bins)
  for b=0,nbins-1 do begin  
    zpos = 0
    zspec=ohspec[b*binslim:((b+1)*binslim)-1]
    maxval = max(zspec[where(finite(zspec))])
    zline=zspec
    zline[where(zline lt maxval/50)] = 0
    zdiff = zline[0:n_elements(zline)-2] - zline[1:n_elements(zline)-1]
    for i=10L,n_elements(zdiff)-10 do if (zdiff[i] gt 0 and zdiff[i-1] lt 0) then zpos=[zpos,i]
    
    ; Reject those lines with an asymmetric profile down to 1/50 of the peak value
    for i=1,n_elements(zpos)-1 do begin
      left  = zspec[0:zpos[i]-1]
      right = zspec[zpos[i]+1:n_elements(zspec)-1]
      leftlim  = where(reverse(left) lt zspec[zpos[i]]/50)
      rightlim = where(right lt zspec[zpos[i]]/50)
      if abs(leftlim[0]-rightlim[0]) lt 2 then zpos2=[zpos2,zpos[i]+b*binslim]
    endfor
  endfor
  zpos2 = zpos2[1:n_elements(zpos2)-1] ;Clip 1st element because it is 0
  skywave_all = lamboh[zpos2]          ;Convert into basolute wavelength
  if state.debug eq 1 then printf, state.log_lun, '[KUBEVIZ] The number of identified OH lines is: '+kubeviz_str(N_elements(skywave_all))

endif else begin
  
  ;Now calculate how much of the spectrum is above 6000A, there are no lines below this value for a meaningful fit
  Nwave = n_elements(*state.wave)
  dummy = where((*state.wave) gt 6000, Nabove6000)
  
  if Nabove6000 lt Nwave/2. then begin
    printf, state.log_lun, '[KUBEVIZ] The spectrum is too blue for a fit of the instr. resolution from skylines.'
    printf, state.log_lun, '[KUBEVIZ] Input the instrumental resolution manually.'
    return
  endif
  
  printf, state.log_lun, '[WARNING] The template OH spec does not overlap enough with the data spectral range.'
  printf, state.log_lun, '[KUBEVIZ] The position of the skylines is guessed automatically.'
  
  kubeviz_medianspec, /sum, /all
  wave = *state.wave
  var = reform(*state.nmedspec)^2 ; to access skylines

  ; First truncate spectrum to useful part:
  var = var[where(var gt 0.)]
  wave = wave[where(var gt 0.)]
  Nw = n_elements(wave)
  
  ; using only values between the 30 & 60 %iles,
  ; smooth it fairly heavily to get a "sky continuum"
  smth_size = 20
  var_smth = smooth(var,smth_size,/EDGE_TRUNCATE)
  var_sub = var-var_smth

  ; histogram to find the peak and rms of most values in var_sub
  ; initial estimate using most of range:
  perc_10 = (kubeviz_percentile(var_sub,10))[0]
  perc_80 = (kubeviz_percentile(var_sub,80))[0]
  Nbins = (Nw/15.)
  hh = histogram(var_sub,Nbins=Nbins,min=perc_10,max=perc_80)
  bb = perc_10+(perc_80-perc_10)*indgen(Nbins)/(1.*(Nbins-1))

  ; repeat with improved binsize - iterate.
  for it=1,5 do begin
      hh_peakval = max(hh,hh_peak)
      i=hh_peak
      while hh[i] ge .5*hh_peakval and i gt 0 do i--
      hh_hm_min = i
      i=hh_peak
      while hh[i] ge .5*hh_peakval do i++
      hh_hm_max = i
      hh_hm_fwhm_estimate = .5*(bb[hh_hm_max]+bb[hh_hm_max-1]) - .5*(bb[hh_hm_min]+bb[hh_hm_min+1])
      binsize = hh_hm_fwhm_estimate/10.

      Nbins = (perc_80-perc_10)/binsize
      if finite(Nbins) eq 0 then return
      hh = histogram(var_sub,Nbins=Nbins,min=perc_10,max=perc_80)
      bb = perc_10+(perc_80-perc_10)*indgen(Nbins)/(1.*(Nbins-1))

      maxhh = max(hh,maxbb)
      peak = bb[maxbb]
      gg = gaussfit(bb,hh,A,Nterms=6)
      mode = A[1]
      rms = A[2]
  endfor

  ; focus on the positive clipped values, divide into a set of "lines",
  ; excluding other features which are obviously too broad.
  skylines_start = [-1]
  skylines_end = [-1]
  currentline = -1
  inline = 0
  skylinewidth_max = 15 ; pixels above 3-sigma
  ; exclude 10 pixels at either end:
  for i=10, Nw-11 do begin
      if var_sub[i]-mode ge 10.*rms then begin
  	  if inline eq 0 then begin
  	      ; new line
  	      inline=1
  	      currentline++
  	      if skylines_start[0] eq -1 then skylines_start[0]=i else skylines_start = [skylines_start,i]
  	  endif
      endif else begin
  	  if inline eq 1 then begin
  	      ; end of line
  	      inline=0
  	      if skylines_end[0] eq -1 then skylines_end[0]=i else skylines_end = [skylines_end,i]
  	      skyline_width = skylines_end[currentline] - skylines_start[currentline]
  	      ; remove if too broad to be real line:
  	      if skyline_width gt skylinewidth_max then begin
  		  currentline--
  		  if currentline eq -1 then begin
  		      skylines_start = [-1] 
  		      skylines_end = [-1]
  		  endif else begin
  		      skylines_start = skylines_start[indgen(currentline+1)]
  		      skylines_end = skylines_end[indgen(currentline+1)]
  		  endelse
  	      endif
  	  endif
      endelse
  endfor
  Nskylines = currentline+1
  skywave_all = (wave[skylines_start]+wave[skylines_end])/2
  if n_elements(skywave_all) eq 0 then begin
    printf, state.log_lun, '[WARNING] The number of identified lines in the variance spec is zero.'
    printf, state.log_lun, '[KUBEVIZ] Input the instrumental resolution manually.'
    return
  endif
endelse

; now loop through these lines and fit them:
skylinepos = 0
skylineres = 0

;Step is required to avoid hours of CPU time for large IFUs (SINFONI)
;Still loops over all the spaxels for KMOS. 25 takes dithering into account
step = max([1,round(sqrt(state.Ncol*state.Nrow)/25.)])
for col=0,state.Ncol-1, step do begin
  for row=0, state.Nrow-1, step do begin
   
   ; get spectrum and wavelength vector
   wave  = *state.wave
   var = reform((*state.noise)[col,row,*])^2 ; to access skylines
   
   ; First truncate spectrum to useful part:
   keep_trunc = where(var gt 0., Nokk)
   if Nokk gt 0 then begin

       var = var(keep_trunc)
       wave = wave(keep_trunc)
       
       ;now clip the skylines outside the valid data range and convert to pixel values
       if n_elements(wave) lt 10 then Nskylines =0 else $
         oklines = where(skywave_all gt wave[fix(n_elements(wave)/100.*5.)] and skywave_all lt wave[fix(n_elements(wave)/100.*95.)], Nskylines)
       if Nskylines gt 0 then begin
         skywave = skywave_all[oklines]
         skyindx = kubeviz_closest(wave, skywave)
       
         for i=0, Nskylines-1 do begin
   	   ; need extended wavelength range to fit:
   	   wrfit = [skyindx[i]-7,skyindx[i]+7]
   	   gskyline = mpfitpeak(wave[wrfit[0]:wrfit[1]], var[wrfit[0]:wrfit[1]], A, Nterms=4, chisq=chisq, sigma=sigma, /positive)
   	   resid  = total(abs((var[wrfit[0]:wrfit[1]] - A[3]- A[0]*exp(-0.5*((wave[wrfit[0]:wrfit[1]] - A[1])/A[2])^2)))/total(var[wrfit[0]:wrfit[1]])*100)
	   ; Keep only the results from good fit. Line roughly centered and small residuals
	   if resid lt 10 and abs(skywave[i]-A[1]) lt 2*state.dlambda then begin
	     skylinepos = [skylinepos,A[1]]
	     ; convert to a resolution R value:
             skylineres = [skylineres,A[1]/(A[2]*2.35482)]
           endif
         endfor
       endif 
   endif  
endfor
endfor

skylinepos_bin = 0
skylineres_bin = 0
skylineres_sig = 0

if state.debug eq 1 then begin
    !x.style = 1
    !y.style = 1
    set_plot, 'ps'
    filename = state.cwdir+'skynoise_res.eps'
    device, filename = filename, /encapsulate, /color
    if state.ifu gt 0 then title = '!6Instr: '+strupcase(state.instr)+' IFU: '+strtrim(state.ifu,2)+' Band: '+strupcase(state.band) else $
                           title = '!6Instr: '+strupcase(state.instr)+' Band: '+strupcase(state.band)
    plotsym, 0, 0.5
    plot, [(*state.wave)[0], (*state.wave)[state.Nwpix-1]], [0,7000], title=title, xtitle='!6Wavelength (!6!sA!r!u!9 %!6!n)', $
          ytitle='!6Resolution', xthick=3, ythick=3, thick=3, charthick=3, /nodata
endif

for i=0,n_elements(skywave_all)-1 do begin
  thisline = where(abs(skylinepos-skywave_all[i]) lt 2*state.dlambda, Nok)
  if Nok gt 5 then begin
    ;For each line with more than 5 good measurements iterative median sigma clipping
    thislineres = skylineres[thisline]
    kubeviz_sigma_clip, thislineres, nsig=2
    plotsym, 0, 0.5
    oplot, replicate(skywave_all[i], n_elements(thislineres)), thislineres, psym=8, col=248
    ; Save position, median resolution and sigma of the mean (median treated as the mean)
    skylinepos_bin = [skylinepos_bin,skywave_all[i]]
    skylineres_bin = [skylineres_bin,median(thislineres, /even)]
    skylineres_sig = [skylineres_sig,(.5*((kubeviz_percentile(thislineres,84))[0]-(kubeviz_percentile(thislineres,16)))[0])/sqrt(n_elements(thislineres))]
  endif
endfor

if n_elements(skylinepos_bin) gt 1 then begin 
 skylinepos_bin = skylinepos_bin[1:n_elements(skylinepos_bin)-1]
 skylineres_bin = skylineres_bin[1:n_elements(skylineres_bin)-1]
 skylineres_sig = skylineres_sig[1:n_elements(skylineres_sig)-1]
 ; fit a poly function and save the coefficients as those from variance
 polypars = poly_fit(skylinepos_bin, skylineres_bin, 4, measure_errors = skylineres_sig, /double)
 state.instrres_varpoly[0:4] = polypars
endif 

if state.debug eq 1 then begin
    oploterror, skylinepos_bin, skylineres_bin, skylineres_sig, psym=4, col=253, errcolor=253, thick=3
    oplot, wave,  state.instrres_extpoly[0]+state.instrres_extpoly[1]*wave+state.instrres_extpoly[2]*wave^2+$
    		  state.instrres_extpoly[3]*wave^3+state.instrres_extpoly[4]*wave^4, col=254, thick=5
    oplot, wave,  polypars[0]+polypars[1]*wave+polypars[2]*wave^2+$
    		  polypars[3]*wave^3+polypars[4]*wave^4, col=249, thick=5
    device,/close
    set_plot, 'X'
    !p.multi = 0	   
endif

end

;------------------------------------------------
pro kubeviz_ngaussmodel, x, a, fx
common kubeviz_state
common kubeviz_fit

sigtofwhm = 2.35482

; Use the proper definition of a Gaussian: 

; gauss = norm / ( sigma * sqrt(2.0 * !DPI) ) * exp(-0.5 *(f)^2)
;     where f = (x - xcen) / width
;           xcen = lambda0*(1.+dv/c)
;           width = lambda0*(sig(v)/c)

; Fit a single lineset together.
; Allows narrow and broad components.

Nx = (size(x))[1]

; start with continuum:
fx = replicate(a[fit.cpar],Nx)

; number of parameters:
Nnfitpars = total(fit.nfitpars)
Nbfitpars = total(fit.bfitpars) 
Npars = 1 + Nnfitpars + Nbfitpars

; Lines to fit: 
if Nnfitpars gt 0 then begin
    nfitlines = fit.nfitpars[2:n_elements(fit.nfitpars)-1]
    nfit = where(nfitlines eq 1, Nnfitlines)
endif

if Nbfitpars gt 0 then begin
    bfitlines = fit.bfitpars[2:n_elements(fit.bfitpars)-1]
    bfit = where(bfitlines eq 1, Nbfitlines)
endif

; add narrow components:
if Nnfitpars gt 0 then begin
    dv = a[fit.nparstart]
    sigv = a[fit.nparstart+1]
    for i=0, Nnfitlines-1 do begin
        lambda0 = state.lines[nfit[i]]
        xcen = lambda0*(1.+dv/state.ckms)
        width = lambda0*(sigv/state.ckms)
        instrres_A = xcen/(sigtofwhm*kubeviz_getinstrres(lambda=xcen)) ; convert from R to Angstroms
	totwidthsq = width^2 + instrres_A^2
        norm = a[fit.nparstart+2+i]
        gn = (norm/sqrt(2.*!DPI*totwidthsq))*exp(-((x-xcen)^2/totwidthsq)/2.)
        fx += gn
    endfor
endif

; add broad components:
if Nbfitpars gt 0 then begin
    dv = a[fit.bparstart]
    sigv = a[fit.bparstart+1]
    for i=0, Nbfitlines-1 do begin
        lambda0 = state.lines[bfit[i]]
        xcen = lambda0*(1.+dv/state.ckms)
        width = lambda0*(sigv/state.ckms)
        instrres_A = xcen/(sigtofwhm*kubeviz_getinstrres(lambda=xcen)) ; convert from R to Angstroms
        totwidthsq = width^2 + instrres_A^2
        norm = a[fit.bparstart+2+i]
        gb = (norm/sqrt(2.*!DPI*totwidthsq))*exp(-((x-xcen)^2/totwidthsq)/2.)
        fx += gb
    endfor
endif

end


;------------------------------------------------
pro kubeviz_linefit_getstartvals, col=col, row=row, silent=silent
common kubeviz_state
common kubeviz_fit

; Set the start "guess" values for fitting.

if N_elements(silent) eq 0 then silent=0

; fit.pars array contains all and only the parameters to fit.
; described by the fit.nfitpars & fit.bfitpars & fit.cfitpars arrays and
; and indices described by the cfit.par, fit.nparstart, fit.nparend,
; fit.bparstart & fit.bparend arrays.

Nnfitpars = total(fit.nfitpars)
Nbfitpars = total(fit.bfitpars)
Npars = 1 + Nnfitpars + Nbfitpars
Nnfitlines = max([Nnfitpars-2,0])
Nbfitlines = max([Nbfitpars-2,0])

; Lines to fit: 
if Nnfitlines gt 0 then begin
    nfitlines = fit.nfitpars[2:n_elements(fit.nfitpars)-1]
    nfit = where(nfitlines eq 1)
    nmainline = kubeviz_getmainline(linesarr=(*state.linenames)[nfit])
    nmainlineind = where(abs(state.lines[nfit] - nmainline) lt 2)
endif

if Nbfitlines gt 0 then begin
    bfitlines = fit.bfitpars[2:n_elements(fit.bfitpars)-1]
    bfit = where(bfitlines eq 1)
    bmainline = kubeviz_getmainline(linesarr=(*state.linenames)[bfit])
    bmainlineind = where(abs(state.lines[bfit] - bmainline) lt 2)
endif

; starting values:
gpars = replicate(0.D,Npars)

; Continuum (depends on the mode)
cpar = fit.cpar
; set to zero if SDSS method is in use.
if state.continuumfit_mode eq 0 then gpars[cpar] = 0.
; else use the SDSS method result as an estimate
if state.continuumfit_mode eq 1 then gpars[cpar] = (*fit.contpars)[0]

; For a given spaxel [col,row] Establish whether a user-defined central 
; velocity and window has been defined 
if n_elements(col) gt 0 and n_elements(row) gt 0 and ptr_valid(state.gauss_initmap) then begin
  
  nestimates = (size((*state.gauss_initmap)))[3]
  offsetinkms    = (*state.gauss_initmap)[col,row,0]
  nwidthinkms    = (*state.gauss_initmap)[col,row,1]
  ; Flux estimates for NARROW LINES ONLY (CAN BE SOMETHING TO UPDATE)
  fluxestimates0 = reform((*state.gauss_initmap)[col,row,2+indgen(nestimates-2)])
  validfluxestimates0 = where(fluxestimates0 ne -999, Nvalidfluxestimates0)
  
  if state.debug eq 1 then begin 
    print, '[DEBUG] Starting values from the external map:'  
    print, reform((*state.gauss_initmap)[col,row,*])
  endif
  
endif else begin
  ; estimate from spectrum:
  kinguess = kubeviz_linefit_estimatekinematics(gpars,nmainline)
  if kinguess[0] eq -1 then begin ;if estimate fails use those rough values
     offsetinkms = 0.
     nwidthinkms = 25. 
  endif else begin
     offsetinkms = (kinguess[0]-nmainline)/nmainline*state.ckms
     nwidthinkms = (kinguess[1]/nmainline)*state.ckms
     if offsetinkms lt 1E-3 then offsetinkms = 0.
     if nwidthinkms lt 1E-3 then nwidthinkms = 0.
  endelse
  Nvalidfluxestimates0=0
endelse

; Line offset in velocity units:
par = 1
npar = fit.nparstart
bpar = fit.bparstart
;Separate the lines if 2nd component is set to larger offset
if Nbfitlines gt 0 and state.secondcomp_mode le 1 then begin
  noffsetinkms = offsetinkms -100.
  boffsetinkms = offsetinkms +100.
endif else begin
  noffsetinkms = offsetinkms
  boffsetinkms = offsetinkms
endelse

; If initial guess is user defined in GUI then use it
if state.gnfit[par] ne -999. then noffsetinkms = state.gnfit[par] 
if state.gbfit[par] ne -999. then boffsetinkms = state.gbfit[par] 
; Only if the set constraints button is ON, we then check whether the
; starting value is within the range of the limits:
if Nnfitlines gt 0 then gpars[npar] = kubeviz_linefit_startminmaxconsistencycheck(noffsetinkms,par, 'N',silent=silent)
if Nbfitlines gt 0 then gpars[bpar] = kubeviz_linefit_startminmaxconsistencycheck(boffsetinkms,par, 'B',silent=silent)

; Line width in velocity units
par = 2
npar = fit.nparstart+1
bpar = fit.bparstart+1
;Separate the lines if 2nd component is set to larger width
if Nbfitlines gt 0 and state.secondcomp_mode eq 2 then begin
  bwidthinkms = 3.*nwidthinkms
endif else begin
  bwidthinkms = nwidthinkms
endelse

; If initial guess is user defined in GUI then use it
if state.gnfit[par] ne -999. then nwidthinkms = state.gnfit[par] 
if state.gbfit[par] ne -999. then bwidthinkms = state.gbfit[par]
; Only if the set constraints button is ON, we then check whether the
; starting value is within the range of the limits:
if Nnfitlines gt 0 then gpars[npar] = kubeviz_linefit_startminmaxconsistencycheck(nwidthinkms,par, 'N',silent=silent)
if Nbfitlines gt 0 then gpars[bpar] = kubeviz_linefit_startminmaxconsistencycheck(bwidthinkms,par, 'B',silent=silent) 

; Line fluxes:

; 1st component (Narrow):
if Nnfitlines gt 0 then begin
    for i=0,Nnfitlines-1 do begin
        line = state.lines[nfit[i]]
        par = nfit[i]+3
        npar = fit.nparstart+2+i
        if state.gnfit[par] ne -999. then begin
	  gpars[npar] = kubeviz_linefit_startminmaxconsistencycheck( state.gnfit[par], par , 'N',silent=silent) 
	endif else begin
	    makeestimate=1
	    fluxestimate=0
	    if Nvalidfluxestimates0 gt 0 then begin
	      if fluxestimates0[i] ne -999 then begin
	        gpars[npar]  = fluxestimates0[i] 
	        makeestimate = 0
	      endif	
	    endif 
	    if makeestimate eq 1 then begin
	      if i eq nmainlineind then begin
                goff_A  = (gpars[fit.nparstart]/state.ckms)*line
                gwidth_A = (gpars[fit.nparstart+1]/state.ckms)*nmainline
	        wr_line = where((*fit.wave) ge line+goff_A-10. and (*fit.wave) le line+goff_A+10., Nwrline)
                if Nwrline gt 0 then fluxestimate = max((*fit.spec)[wr_line])*sqrt(2.*!DPI)*gwidth_A
	      endif else begin
	        fluxestimate = 0.33 * gpars[fit.nparstart+2+nmainlineind]
	      endelse
	     gpars[npar] = kubeviz_linefit_startminmaxconsistencycheck( fluxestimate, par , 'N' ,silent=silent)
	    endif
        endelse
    endfor
endif

; 2nd component (Broad):
if Nbfitlines gt 0 then begin
    for i=0,Nbfitlines-1 do begin
        line = state.lines[bfit[i]]
        par = bfit[i]+3
        bpar = fit.bparstart+2+i
        if state.gbfit[par] ne -999. then begin
	  gpars[bpar] = kubeviz_linefit_startminmaxconsistencycheck( state.gbfit[par], par, 'B' ,silent=silent)
        endif else begin
            makeestimate=1
	    fluxestimate=0
	    if Nvalidfluxestimates0 gt 0 then begin
	      if fluxestimates0[i] ne -999 then begin
	        gpars[npar]  = fluxestimates0[i] 
	        makeestimate = 0
	      endif	
	    endif
	    if makeestimate eq 1 then begin
               if i eq bmainlineind then begin
	         goff_A  = (gpars[fit.nparstart]/state.ckms)*line
	         gwidth_A = (gpars[fit.bparstart+1]/state.ckms)*bmainline
	         wr_line = where((*fit.wave) ge line+goff_A-10. and (*fit.wave) le line+goff_A+10., Nwrline)
                 if Nwrline gt 0 then fluxestimate = max((*fit.spec)[wr_line])*sqrt(2.*!DPI)*gwidth_A
               endif else begin
	        fluxestimate= 0.33 * gpars[fit.bparstart+2+bmainlineind]
	       endelse
	     gpars[bpar] = kubeviz_linefit_startminmaxconsistencycheck( fluxestimate, par, 'B' ,silent=silent)
	    endif
        endelse
    endfor
endif

fit.gpars = ptr_new(gpars)

if state.debug eq 1 then begin 
  print, '[DEBUG] Starting values:'  
  print, gpars
endif

end

;----------------------------------------------------------------
pro kubeviz_linefit_setparinfo
common kubeviz_state
common kubeviz_fit

; setup the parinfo array for fitting which sets constraints on the parameters

; fit.pars array contains all and only the parameters to fit.
; described by the fit.nfitpars, fit.bfitpars & fit.cfitpars arrays and
; indices described by the cfit.par, fit.nparstart, fit.nparend,
; fit.bparstart & fit.bparend arrays.

Nnfitpars = total(fit.nfitpars)
Nbfitpars = total(fit.bfitpars) 
Npars = 1 + Nnfitpars + Nbfitpars
Nnfitlines = max([Nnfitpars-2,0])
Nbfitlines = max([Nbfitpars-2,0])

; PARAMETERS to fit: 
if Nnfitlines gt 0 then nfitp = where(fit.nfitpars eq 1)
if Nbfitlines gt 0 then bfitp = where(fit.bfitpars eq 1)

parinfo = replicate({fixed:0, limited:[0,0], limits:[-9.e9,9.e9], tied:''}, Npars)

; DJW 14-01-10 continuum is now FIXED regardless of constraints option, as we have already subtracted it.
; MF 14-12-04 now depends on the continuum method in use
parinfo[fit.cpar].fixed = 1-state.continuumfit_mode

; choice: fit with or without user-constraints. 
; If we are fitting multiple linesets we need to go through this procedure
; to fix kinematics to the strongest lines
if state.fitconstr eq 1 or state.selected_lineset eq 0 then begin
    ; narrow lines:
    if Nnfitpars gt 0 then begin
        for i=0,Nnfitpars-1 do begin
            par = nfitp[i]+1
            npar=fit.nparstart+i
            if state.pnfix[par] eq 1 then parinfo[npar].fixed = 1
            if state.gnlims[par,0] ne -999. then begin
                parinfo[npar].limited[0] = 1
                parinfo[npar].limits[0] = state.gnlims[par,0]
            endif else begin
            ; min flux and width is always zero:
                if par gt 1 then begin
                    parinfo[npar].limited[0] = 1
                    parinfo[npar].limits[0] = 0.
                endif
            endelse
            if state.gnlims[par,1] ne -999. then begin
                parinfo[npar].limited[1] = 1
                parinfo[npar].limits[1] = state.gnlims[par,1]
            endif
        endfor
    endif
    ; broad lines:
    if Nbfitpars gt 0 then begin
        for i=0,Nbfitpars-1 do begin
            par = bfitp[i]+1
            bpar=fit.bparstart+i
            if state.pbfix[par] eq 1 then parinfo[bpar].fixed = 1
            if state.gblims[par,0] ne -999. then begin
                parinfo[bpar].limited[0] = 1
                parinfo[bpar].limits[0] = state.gblims[par,0]
            endif else begin
            ; min flux and width is always zero:
                if par gt 1 then begin
                    parinfo[bpar].limited[0] = 1
                    parinfo[bpar].limits[0] = 0.
                endif
            endelse
            if state.gblims[par,1] ne -999. then begin
                parinfo[bpar].limited[1] = 1
                parinfo[bpar].limits[1] = state.gblims[par,1]
            endif      
        endfor
    endif
endif else begin
    ; Flux and width >= 0 even without constraints option set.
    if Nnfitpars gt 2 then begin
        for i=1,Nnfitpars-1 do begin
            npar=fit.nparstart+i
            parinfo[npar].limited[0] = 1
            parinfo[npar].limits[0] = 0.
        endfor
    endif
    if Nbfitpars gt 2 then begin
        for i=1,Nbfitpars-1 do begin
            bpar=fit.bparstart+i
            parinfo[bpar].limited[0] = 1
            parinfo[bpar].limits[0] = 0.
        endfor
    endif  
endelse

; choice: fix line ratios
if state.fitfixratios eq 1 then begin
    ; narrow lines:
    if Nnfitpars gt 0 then begin
        linepos = replicate(-1,6) ; Respectively NII, OIII,OI blue and red
        linenames = (*state.linenames)[where(fit.nfitpars[2:*] eq 1)]
        for i=2,Nnfitpars-1 do begin
          npar=fit.nparstart+i
	  case linenames[i-2] of
	   'n2_b' : linepos[0] = npar
	   'n2_r' : linepos[1] = npar
	   'o3_b' : linepos[2] = npar
	   'o3_r' : linepos[3] = npar
	   'o1_b' : linepos[4] = npar
	   'o1_r' : linepos[5] = npar 
	   else: 
	  endcase
	endfor
	if linepos[0] ge 0 and linepos[1] ge 0 then parinfo[linepos[0]].tied='P['+kubeviz_str(linepos[1])+']/3.071' ; NII  from Storey & Zeippen 2000
	if linepos[2] ge 0 and linepos[3] ge 0 then parinfo[linepos[2]].tied='P['+kubeviz_str(linepos[3])+']/3.013' ; OIII from Storey & Zeippen 2000
	if linepos[4] ge 0 and linepos[5] ge 0 then parinfo[linepos[5]].tied='P['+kubeviz_str(linepos[4])+']/2.997' ; OI   from Storey & Zeippen 2000
    endif
    
    ; broad lines:
    if Nbfitpars gt 0 then begin
        linepos = replicate(-1,8) ; Respectively NII, OIII, SII, OI blue and red
	linenames = (*state.linenames)[where(fit.bfitpars[2:*] eq 1)]
	for i=2,Nbfitpars-1 do begin
          bpar=fit.bparstart+i
          case linenames[i-2] of
	   'n2_b' : linepos[0] = bpar
	   'n2_r' : linepos[1] = bpar
	   'o3_b' : linepos[2] = bpar
	   'o3_r' : linepos[3] = bpar
	   'o1_b' : linepos[4] = bpar
	   'o1_r' : linepos[5] = bpar 
	   else: 
	  endcase
	endfor
	if linepos[0] ge 0 and linepos[1] ge 0 then parinfo[linepos[0]].tied='P['+kubeviz_str(linepos[1])+']/3.071' ; NII  from Storey & Zeippen 2000
	if linepos[2] ge 0 and linepos[3] ge 0 then parinfo[linepos[2]].tied='P['+kubeviz_str(linepos[3])+']/3.013' ; OIII from Storey & Zeippen 2000
	if linepos[4] ge 0 and linepos[5] ge 0 then parinfo[linepos[4]].tied='P['+kubeviz_str(linepos[5])+']/2.997' ; OI   from Storey & Zeippen 2000
     endif	
endif

fit.parinfo = ptr_new(parinfo)

end

;--------------------------------------------------------------------------
pro kubeviz_linefit_fit, wave, spec, weights, gpars, pars, savepars=savepars

common kubeviz_state
common kubeviz_fit

if n_elements(savepars) eq 0 then savepars = 0 ; default OFF

itmax=200
quiet = 1-state.debug

openu, lun, state.cwdir+'fit_errmsg_'+strmid(state.filename,0,strlen(state.filename)-5)+'.txt', /get_lun, /append

; initial estimates:
pars = gpars
specfit = mpcurvefit( wave, spec, weights, pars, sigpars, chisq=chisq, dof=dof, $
                  itmax=itmax, function_name='kubeviz_ngaussmodel', iter=iter, parinfo=(*fit.parinfo), $
                  status=status, /autoderivative, errmsg=errmsg, quiet=quiet)

if strlen(errmsg) gt 0 then begin ;If an error occurred
    case state.linefit_mode of
    0: begin
       printf, lun, 'Spaxel: ', state.col, state.row
       printf, lun, 'Fitting Error Message: ',errmsg
       end
    1: begin
       printf, lun, 'Mask: ', state.imask-1
       printf, lun, 'Fitting Error Message: ',errmsg
       end 
    endcase     
    if state.debug eq 1 then begin
        print, 'Fitting Error Message: ',errmsg
        print, 'Pause to debug: .c to continue'
        stop
    endif
    Npars = 1 + total(fit.nfitpars) + total(fit.bfitpars)
    if savepars eq 1 then begin
        fit.pars = ptr_new(replicate(-999.D,Npars))
        fit.sigpars = ptr_new(replicate(-999.D,Npars))
    endif
endif else begin ;If everything was OK
    if savepars eq 1 then begin
        fit.iter = iter
        fit.specfit = ptr_new(specfit)
        fit.pars = ptr_new(pars)
        fit.sigpars = ptr_new(sigpars)
        fit.chisq = chisq/dof
        fit.status = status
        fit.errmsg = errmsg
    endif
endelse

if state.debug eq 1 then begin
    print, 'Starting pars: ',(*fit.gpars)
    print, 'Pars: ',pars
;    print, 'Par Errors: ',sigpars
    print, 'chisq: ',chisq
    print, 'dof: ', dof
    print, 'status: ',status
    print, 'Number of iterations: ', iter
    print, 'Error Message: ', errmsg
endif

free_lun, lun

;stop

end

;-------------------------------------------------------------------------
pro kubeviz_cut_momoutliers, moments, Nok
  
if (n_elements(moments) gt 0) then begin
  if (n_elements(moments) gt 2) then begin
    do_cut = 1
    i = 1
    while ((do_cut eq 1) and (n_elements(moments) gt 2)) do begin
      ; cut a single value at the left?
      if (moments[i]-1 ne moments[i-1]) then moments = moments[i:n_elements(moments)-1] else begin
	; cut a pair of single values at the left?
	if ((do_cut eq 1) and (moments[i+1]-2 ne moments[i-1])) then moments = moments[i+1:n_elements(moments)-1] else do_cut = 0
      endelse
    endwhile

    do_cut = 1
    while ((do_cut eq 1) and (n_elements(moments) gt 2)) do begin
      i = n_elements(moments)-1
      ; cut a single value at the right?
      if (moments[i]-1 ne moments[i-1]) then moments = moments[0:i-1] else begin
	; cut a pair of values at the right?
	if (moments[i]-2 ne moments[i-2]) then moments = moments[0:i-2] else do_cut = 0
      endelse
    endwhile
    if (n_elements(moments) lt 3) then moments = -1
  endif else moments = -1
endif else moments = -1

if moments[0] eq -1 then Nok = 0 else Nok = N_elements(moments)
  
end 

;--------------------------------------------------------------------------
pro kubeviz_linefit_moments, wave, spec, weights, pars, savepars=savepars, col=col, row=row, boot=boot

common kubeviz_state
common kubeviz_fit
sigtofwhm = 2.35482

if n_elements(boot) eq 0 then boot=-1 ; default is NOT a bootstrap cube

if n_elements(savepars) eq 0 then savepars = 0 ; default OFF

; For a given spaxel [col,row] Establish whether a user-defined central velocity and window has been defined 
if n_elements(col) gt 0 and n_elements(row) gt 0 and ptr_valid(state.mom_windowmap) then begin
  fixedwin = 1
; If this is from a bootstrap cube, and there is a valid bootstrap version of the user-defined central velocity and window, use that:
  if boot ge 0 and ptr_valid(state.mom_windowmap_bootstrap[boot]) then begin
    vcen = (*state.mom_windowmap_bootstrap[boot])[col,row,0]
    vwin = (*state.mom_windowmap_bootstrap[boot])[col,row,1]
  endif else begin  
; No valid bootstrap version so use same user-defined central velocity and window for all bootstraps:
    vcen = (*state.mom_windowmap)[col,row,0]
    vwin = (*state.mom_windowmap)[col,row,1]
  endelse
  
endif else begin
  fixedwin = 0
  vcen = 0
  vwin = 300
endelse

lines_mom = where(fit.mfitpars eq 1, Nlinesmom)
pars = replicate(0.D, 5*Nlinesmom)
sigpars = replicate(-998.D, 5*Nlinesmom)

for line = 0, Nlinesmom-1 do begin
; locally keep the window in wavelength as defined using vcen,vwin
   lambdaz = state.lines[lines_mom[line]]
   lambda0 = state.lines[lines_mom[line]]*(1+vcen/state.ckms)
   dlambdawin = lambda0*(vwin/state.ckms) ; half-window size
   bisectors = state.line_bisectors[lines_mom[line],*]*(1+vcen/state.ckms)
   
   if ptr_valid(fit.sigcont) then sigcont = (*fit.sigcont)[0] else sigcont = 0
   
   good = where(wave gt bisectors[0] and wave lt bisectors[1] $ 
           and spec gt state.mom_thresh * sigcont, Ngood )
   kubeviz_cut_momoutliers, good, Ngood	   
   
   if Ngood gt 0 then begin
      mom1   = (total(spec[good]*weights[good]*wave[good])/total(spec[good]*weights[good]))                        ;1st order moment MEAN
      mom2   = (total(spec[good]*weights[good]*(wave[good]-mom1)^2)/total(spec[good]*weights[good]))^0.5           ;2nd order moment STD DEV
      mom3   = (total(spec[good]*weights[good]*((wave[good]-mom1)/mom2)^3)/total(spec[good]*weights[good]))        ;3rd order moment SKEWNESS
      mom4   = (total(spec[good]*weights[good]*((((wave[good]-mom1)/mom2)^4)-3))/total(spec[good]*weights[good]))  ;4th order moment KURTOSIS
      pars[line*5+1] = state.ckms*((mom1/lambdaz)-1) 								   ; Convert into velocity offset
      pars[line*5+2] = (state.ckms/lambdaz)*sqrt(max([(mom2)^2-(mom1/(sigtofwhm*kubeviz_getinstrres(lambda=mom1)))^2,0],/nan))  ; Convert into line width
      pars[line*5+3] = mom3
      pars[line*5+4] = mom4
   endif  
   
   okk_max = where(wave gt bisectors[0] and wave lt bisectors[1], Nokkmax)
   okk_min = where(wave gt lambda0-dlambdawin and wave lt lambda0+dlambdawin, Nokkmin)
   
; if the user-defined central velocity and window has been defined, stick with it:
   if fixedwin eq 1 then okk_final = okk_min else begin
   
; use iterative method, defining a window limited by okk_min and okk_max:
    okk_diff = kubeviz_setdifference(okk_max, okk_min)
    ind_below = where(wave[okk_diff] lt lambda0, Nokkbelow)
    ind_above = where(wave[okk_diff] gt lambda0, Nokkabove)
   
    if Nokkbelow gt 0 then begin
      okk_below = okk_diff[ind_below]
      limit_below = okk_below[Nokkbelow-1]
      while spec[limit_below] gt 0. and limit_below ge okk_below[0] and limit_below ge 1 do limit_below--
      limit_below++
    endif else limit_below = min(okk_min)
    if Nokkabove gt 0 then begin
      okk_above = okk_diff[ind_above]
      limit_above = okk_above[0]
      while spec[limit_above] gt 0. and limit_above le okk_above[Nokkabove-1] and limit_above le n_elements(spec)-2 do limit_above++
      limit_above--
    endif else limit_above = max(okk_min)
   
    okk_final = indgen(limit_above-limit_below+1)+limit_below
   endelse
   
   Nokk_final = n_elements(okk_final)
   
   if Nokk_final  gt 0 then begin
      ;mom0   = (total(spec[okk]*weights[okk])/total(weights[okk])) * float(Nokk)                                         ;0th order moment FLUX
      mom0   = total(spec[okk_final], /nan)                                                                               ;0th order moment FLUX
   endif else mom0 = -99
  
   pars[line*5+0] = mom0 * state.dlambda          ;max([mom0 * state.dlambda,0]) 				          ; Convert into flux 
 
endfor

if savepars eq 1 then begin
  fit.pars=ptr_new(pars)
  fit.sigpars=ptr_new(sigpars)
endif
end

;---------------------------------------------------------
pro kubeviz_linefit_setfixline, linetypes
common kubeviz_state
common kubeviz_fit

; Fix the line position and width.
; Also set the start value to the current best fit value.

Nlinetypes = n_elements(linetypes)
for ilinetype = 0, Nlinetypes-1 do begin
    case linetypes[ilinetype] of
        'N': begin
            state.pnfix[1] = 1
            state.pnfix[2] = 1
            state.gnfit[1] = (*fit.pars)[fit.nparstart+0]  ;nrescube[1]
            state.gnfit[2] = (*fit.pars)[fit.nparstart+1]  ;nrescube[2]
        end
        'B': begin
            state.pbfix[1] = 1
            state.pbfix[2] = 1
            state.gbfit[1] = (*fit.pars)[fit.bparstart+0]  ;brescube[1]
            state.gbfit[2] = (*fit.pars)[fit.bparstart+1]  ;brescube[2]
        end
    endcase
endfor

end

;---------------------------------------------------
pro kubeviz_linefit_saveMonteCarlodistrib, bootpars, cbootpars
common kubeviz_state
common kubeviz_fit

case state.linefit_type of
    0: typestr = 'gau'
    1: typestr = 'mom'
endcase

case state.linefit_mode of
    0: begin
       locstr = 'spax'
       extstr = 'SPX_'+kubeviz_str(state.col+1)+'_'+kubeviz_str(state.row+1)
       nexten = state.col * state.row
    end   
    1: begin
       locstr = 'mask'
       extstr = 'MASK_'+kubeviz_str(state.imask)
       nexten = state.Nmask
    end   
endcase

case state.domontecarlo of
    1: montedir = 'boot'
    2: montedir = 'MC1'
    3: montedir = 'MC2'
    4: montedir = 'MC3'
    else: montedir=''
endcase

savefilename = state.outdir+strmid(state.filename,0,strlen(state.filename)-5)+'_'+typestr+'_'+locstr+'_PDF.fits'

;Transpose to make it a vertical image. The first par in bootpars is empty (flag) in the gaussian fitting mode. 
;Then we append the continuum values.
case state.linefit_type of
    0: pdf_values = [transpose(bootpars[*,1:*]), transpose(cbootpars)]
    1: pdf_values = [transpose(bootpars), transpose(cbootpars)]
endcase

;Prepare the header
mkhdr, hdr, pdf_values
sxdelpar, hdr, 'COMMENT'
sxdelpar, hdr, 'SIMPLE'
sxaddpar, hdr, 'XTENSION', 'IMAGE', before='BITPIX'
sxaddpar, hdr, 'EXTNAME' , kubeviz_str(extstr) 
sxaddpar, hdr, 'DATACUBE', state.filename           	    , 'Name of the original datacube'
case state.linefit_mode of
    0: begin
    sxaddpar, hdr, 'COL', state.col+1          	   	    , 'Spaxel column number'
    sxaddpar, hdr, 'ROW', state.row+1          	   	    , 'Spaxel row number'
    end
    1: begin
    sxaddpar, hdr, 'MASK', state.imask           	    , 'Number of the current mask'
    end
endcase
sxaddpar, hdr, 'FITTYPE' , fix(state.linefit_type)  	    , 'Type of fit: 0: gauss 1: moments'
sxaddpar, hdr, 'RESTYPE' , fix(state.linefit_mode)  	    , 'Type of result: 0: spaxels 1: masks'
sxaddpar, hdr, 'ERMETHOD', state.domontecarlo       	    , 'Err method: 0:Noise, 1:Bootstrap 2,3,4:MC1,2,3'
sxaddpar, hdr, 'MCNOISE' , fix(state.useMonteCarlonoise)    , 'MonteCarlo noise: 0: Off 1: On'
sxaddpar, hdr, 'SPATSMTH', state.smooth             	    , 'Spatial smoothing applied'
sxaddpar, hdr, 'SPECSMTH', state.specsmooth         	    , 'Spectral smoothing applied'
sxaddpar, hdr, 'LINESET' , state.selected_lineset   	    , 'Lineset used'
case state.linefit_type of
  0: begin ;Gaussian fit
     ind=1
     if fit.nfitpars[0] eq 1 then begin
       sxaddpar, hdr, 'COL'+kubeviz_str(ind) , 'narrow line dv'
       ind += 1
     endif  
     if fit.nfitpars[1] eq 1 then begin
       sxaddpar, hdr, 'COL'+kubeviz_str(ind) , 'narrow line sig(v)'
       ind += 1
     endif  
     for i=2, N_elements(fit.nfitpars)-1 do begin
       if fit.nfitpars[i] eq 1 then begin 
        sxaddpar, hdr, 'COL'+kubeviz_str(ind) , 'narrow line '+(*state.linenames)[i-2]+' flux'
        ind += 1
       endif
     endfor
     if fit.bfitpars[0] eq 1 then begin
       sxaddpar, hdr, 'COL'+kubeviz_str(ind) , 'broad line dv'
       ind += 1
     endif
     if fit.bfitpars[1] eq 1 then begin
       sxaddpar, hdr, 'COL'+kubeviz_str(ind) , 'broad line sig(v)'
       ind += 1
     endif
     for i=2, N_elements(fit.bfitpars)-1 do begin
       if fit.bfitpars[i] eq 1 then begin
         sxaddpar, hdr, 'COL'+kubeviz_str(ind) , 'broad line '+(*state.linenames)[i-2]+' flux'
	 ind += 1
       endif
     endfor
     for i=2, N_elements(fit.nfitpars)-1 do begin
       if fit.nfitpars[i] eq 1 then begin 
        sxaddpar, hdr, 'COL'+kubeviz_str(ind) , 'narrow line '+(*state.linenames)[i-2]+' continuum flux'
        ind += 1
       endif
     endfor
     for i=2, N_elements(fit.bfitpars)-1 do begin
       if fit.bfitpars[i] eq 1 then begin
         sxaddpar, hdr, 'COL'+kubeviz_str(ind) , 'broad line '+(*state.linenames)[i-2]+' continuum flux'
	 ind += 1
       endif
     endfor
 
  end
  1: begin ;Moment fit
     ind=1
     for i=0, N_elements(fit.mfitpars)-1 do begin
       if fit.mfitpars[i] eq 1 then begin 
        sxaddpar, hdr, 'COL'+kubeviz_str(ind)   , 'moment '+(*state.linenames)[i]+' 0th (flux)'
	sxaddpar, hdr, 'COL'+kubeviz_str(ind+1) , 'moment '+(*state.linenames)[i]+' 1st (dv)'
	sxaddpar, hdr, 'COL'+kubeviz_str(ind+2) , 'moment '+(*state.linenames)[i]+' 2nd (sigv)'
	sxaddpar, hdr, 'COL'+kubeviz_str(ind+3) , 'moment '+(*state.linenames)[i]+' 3rd (norm Skewness)'
	sxaddpar, hdr, 'COL'+kubeviz_str(ind+4) , 'moment '+(*state.linenames)[i]+' 4th (norm Kurtosis)'
        ind += 5
       endif
     endfor
     for i=0, N_elements(fit.mfitpars)-1 do begin
       if fit.mfitpars[i] eq 1 then begin 
        sxaddpar, hdr, 'COL'+kubeviz_str(ind)   , 'moment '+(*state.linenames)[i]+' continuum'
        ind += 1
       endif
     endfor
  end
endcase

if file_test(savefilename) eq 0 then begin
   
   mkhdr, prihdr, ''
   sxdelpar, prihdr, 'COMMENT' 
   sxaddpar, prihdr, 'NEXT ' , nexten
   sxaddpar, prihdr, 'NMC  ' , state.Nmontecarlo
   
   mwrfits, 0, savefilename, prihdr, /create     
      
endif 

;Add a new extension or replace an existing one ??
fits_info, savefilename, extname=extname, /silent
case total(strmatch(kubeviz_str(extname),extstr, /fold_case)) of 
 0: writefits, savefilename, pdf_values, hdr, /append 
 1: modfits, savefilename, pdf_values, hdr, extname=extstr
 2: stop
endcase 

end

;---------------------------------------------------
pro kubeviz_linefit_plotMonteCarlodistrib, res_bootstrap, percs_bootstrap, val, restype, par
common kubeviz_state
common kubeviz_fit

; make plots of the distributions of line fit params for the bootstrap cubes

case state.linefit_mode of
    0: begin ; spaxels:
        locstr = 'spax_'+kubeviz_str(state.col,format='(I2)')+'_'+kubeviz_str(state.row,format='(I2)')
    end
    1: begin                    ; masks:
        locstr = 'mask_'+kubeviz_str(state.imask-1,format='(I3)')
    end
endcase

case par of
    0: strpar='flag'
    1: if (restype eq 'cont') then strpar='flux@'+(*state.linenames)[par-1] else strpar='dv'
    2: if (restype eq 'cont') then strpar='flux@'+(*state.linenames)[par-1] else strpar='sig(v)'
    else: if (restype eq 'cont') then strpar='flux@'+(*state.linenames)[par-1] else strpar=(*state.linenames)[par-3]+'_flux'
endcase

case state.domontecarlo of
    0: montedir = ''
    1: montedir = 'bootstrap'
    2: montedir = 'montecarlo1'
    3: montedir = 'montecarlo2'
    4: montedir = 'montecarlo3'
    else: print, '[ERROR] Montecarlo method not implemented yet'
endcase

file_mkdir, state.cwdir+montedir
plotfilename = state.cwdir+montedir+'/'+locstr+'_'+restype+'_'+strpar+'_distr.eps'
!x.style = 1
!y.style = 1
!p.multi = 0
set_plot, 'ps'
device,filename = plotfilename, /encapsulate, /color

kubeviz_cumdistrplot, res_bootstrap, col=254, xtit=restype+' '+strpar, ytit='f('+restype+' '+strpar+'>x)', $
              title='Error method: '+montedir, charthick=3., thick=3., linestyle=0
xr = [min(res_bootstrap),max(res_bootstrap)]
oplot, [val,val], [0.,1.], col=254, linestyle=0, thick=5
for iperc=0, state.Nmontecarlo_percs-1 do begin
    case iperc of
        0: begin ; +1sig
            linestyle=1 
            col=254 ; black
        end 
        1: begin ; -1sig
            linestyle=1 
            col=254 ; black
        end
        2: begin ; +2sig
            linestyle=1
            col=252 ; blue
        end
        3: begin ; -2sig
            linestyle=1 
            col=252 ; blue
        end
        4: begin ; median
            linestyle=0
            col=253 ; red
        end
        5: begin ; 10
            linestyle=2
            col=249 ; green
        end
        6: begin ; 90
            linestyle=2
            col=249 ; green
        end
        7: begin ; 25
            linestyle=2
            col=250 ; grey
        end
        8: begin ; 75
            linestyle=2
            col=250 ; grey
        end
        else: begin
            linestyle=2
            col=252
        end
    endcase
    oplot, replicate(percs_bootstrap[iperc],2), [0.,1.], col=col, linestyle=linestyle, thick=2.
    oplot, xr, replicate(.01*(*state.montecarlo_percs)[iperc],2), col=col, linestyle=linestyle, thick=2.
endfor
device,/close
set_plot, 'X'

end

;--------------------------------------------------------------------
function kubeviz_getcontatlambda, contpars, lambda
common kubeviz_state

cont = replicate(0.D,n_elements(lambda))
; estimate continuum level at lambda:
for i=0, state.continuumfit_order do cont += contpars[i]*(lambda^i)

return, cont

end

;-------------------------------------------------------------------
pro kubeviz_linefit_swapnarrowbroad
common kubeviz_fit

;Save narrow pars to buffer
temp     = (*fit.pars)[fit.nparstart:fit.bparstart-1]
temp_err = (*fit.sigpars)[fit.nparstart:fit.bparstart-1]
Nel = n_elements(temp)

;Copy broad onto narrow
(*fit.pars)[fit.nparstart:fit.bparstart-1] = (*fit.pars)[fit.bparstart:fit.bparstart+Nel-1]
(*fit.sigpars)[fit.nparstart:fit.bparstart-1] = (*fit.sigpars)[fit.bparstart:fit.bparstart+Nel-1]

;Copy buffer over broad
(*fit.pars)[fit.bparstart:fit.bparstart+Nel-1] = temp
(*fit.sigpars)[fit.bparstart:fit.bparstart+Nel-1] = temp_err

end


;-------------------------------------------------------------------
pro kubeviz_linefit_keepfit
common kubeviz_state
common kubeviz_fit
; copy the fit results into the results cube

;do we need to scale the errors?
if state.scaleNoiseerrors gt 0 and state.domontecarlo eq 0 then (*fit.sigpars) *= sqrt(fit.chisq)

Nnfitpars = total(fit.nfitpars)
Nbfitpars = total(fit.bfitpars) 
Ncfitpars = total(fit.cfitpars)
Npars = Ncfitpars + Nnfitpars + Nbfitpars
Nnfitlines = max([Nnfitpars-2,0])
Nbfitlines = max([Nbfitpars-2,0])

;if a second component is fit check if results satisfy the condition
;in case of free fit the fainter component is the second in the results
;FIXME can the condition on the equal number of lines be lifted?
if Nbfitlines gt 0 and Nnfitlines eq Nbfitlines then begin
  case state.secondcomp_mode of
    0: if (*fit.pars)[fit.nparstart+2] lt (*fit.pars)[fit.bparstart+2] then kubeviz_linefit_swapnarrowbroad 
    1: if (*fit.pars)[fit.nparstart+0] gt (*fit.pars)[fit.bparstart+0] then kubeviz_linefit_swapnarrowbroad 
    2: if (*fit.pars)[fit.nparstart+1] gt (*fit.pars)[fit.bparstart+1] then kubeviz_linefit_swapnarrowbroad 
  endcase
endif

; pointers to the parameters which have been fit for narrow and broad line components:
if Nnfitlines gt 0 then nfitp = where(fit.nfitpars eq 1)
if Nbfitlines gt 0 then bfitp = where(fit.bfitpars eq 1)

; current values:
case state.linefit_mode of
    0: begin ; spaxels:
        nrescube = reform((*state.sp_nrescube)[state.col,state.row,*])
        nerrrescube = reform((*state.sp_nerrrescube)[state.col,state.row,*,*])
        brescube = reform((*state.sp_brescube)[state.col,state.row,*])
        berrrescube = reform((*state.sp_berrrescube)[state.col,state.row,*,*])
        crescube = reform((*state.sp_crescube)[state.col,state.row,*])
        cerrrescube = reform((*state.sp_cerrrescube)[state.col,state.row,*,*])
    end
    1: begin ; masks:
        nrescube = reform((*state.nrescube)[state.imask-1,*])
        nerrrescube = reform((*state.nerrrescube)[state.imask-1,*,*])
        brescube = reform((*state.brescube)[state.imask-1,*])
        berrrescube = reform((*state.berrrescube)[state.imask-1,*,*])
        crescube = reform((*state.crescube)[state.imask-1,*])
        cerrrescube = reform((*state.cerrrescube)[state.imask-1,*,*])
    end
endcase

if state.domontecarlo gt 0  then begin
    cont_lambda_bootstrap = replicate(0.D, state.Nmontecarlo, Nnfitlines+Nbfitlines)
    cont_boot_par = 0
endif

; NOTE on continuum:
; Continuum is fit separately, but provides a value of the best fit
; continuum at the position of each line. This should be saved at the
; position of the lines which have been fit. Therefore this is done
; below when saving the line parameters.
; NOTE on continuum 2: MF 2014-12-07
; If continuum mode is not 0 (SDSS) the above is not true. We save 
; the MPFIT result.

if Nnfitpars gt 0 then begin
    dv = 0. ; initialize
    for i=0,Nnfitpars-1 do begin
        par = nfitp[i]+1 ; +1 for the flag
        npar=fit.nparstart+i
        if (*fit.pars)[npar] ne -999. and (*fit.sigpars)[npar] ne -999. then begin
	    nrescube[par] = (*fit.pars)[npar]
          ; record dv (1st parameter), useful later for the continuum
	    if i eq 0 then dv=(*fit.pars)[npar] 
          ; If we fit multiple linesets, save the kinematic errors only if this is the first lineset (otherwise errors are 0 due to fix)   
	    if ~(par le 2 and state.selected_lineset eq 0 and fit.firstset eq 0) then begin 
	     if state.domontecarlo gt 0  then begin                
                 percs_bootstrap = kubeviz_percentile(reform((*fit.pars_bootstrap)[*,npar]),*state.montecarlo_percs)
                 nerrrescube[par,*] = percs_bootstrap-(*fit.pars)[npar]
                 if state.plotMonteCarlodistrib eq 1 then kubeviz_linefit_plotMonteCarlodistrib, reform((*fit.pars_bootstrap)[*,npar]), percs_bootstrap, (*fit.pars)[npar], 'narrow', par
             endif else begin
                 ; 1-sig errors in the first 2 elements of the error array for each parameter:                
                 nerrrescube[par,0] = (*fit.sigpars)[npar]
                 nerrrescube[par,1] = -(*fit.sigpars)[npar]
             endelse
	    endif 
        endif
; value of continuum = continuum fit evaluated at the wavelength of the line:
        if par gt 2 and Ncfitpars gt 0 then begin
            cpar = par - 2  
            line = cpar - 1
	    if state.continuumfit_mode eq 0 then begin
	      lambda_line = state.lines[line]*(1.+dv/state.ckms)
              crescube[cpar] = kubeviz_getcontatlambda((*fit.contpars),lambda_line)
            endif else crescube[cpar] = (*fit.pars)[fit.cpar]
	    if state.domontecarlo gt 0  then begin
		for boot=0,state.Nmontecarlo-1 do begin
		   if state.continuumfit_mode eq 0 then cont_lambda_bootstrap[boot,cont_boot_par] = kubeviz_getcontatlambda((*fit.contpars_bootstrap)[boot,*],lambda_line) $
		                                   else cont_lambda_bootstrap[boot,cont_boot_par] = (*fit.pars_bootstrap)[boot,fit.cpar]
		endfor
		percs_bootstrap = kubeviz_percentile(cont_lambda_bootstrap[*,cont_boot_par],*state.montecarlo_percs)
                cerrrescube[cpar,*] = percs_bootstrap-crescube[cpar]
                if state.plotMonteCarlodistrib eq 1 then kubeviz_linefit_plotMonteCarlodistrib, cont_lambda_bootstrap[*,cont_boot_par], percs_bootstrap, crescube[cpar], 'cont', cpar
                cont_boot_par += 1
            endif else begin

	    ; Assuming independent errors to polynomial fit, it is then simple to apply the polynomial errors @ lambda
	    ; cerrrescube[cpar,0] = abs(kubeviz_getcontatlambda((*fit.sigcontpars),lambda_line))
	    ; MF, replaced with standard deviation of residual values used to evaluate the continuum level.
            if state.continuumfit_mode eq 0 then cerrrescube[cpar,0] = (*fit.sigcont) else cerrrescube[cpar,0] = (*fit.sigpars)[fit.cpar]
            cerrrescube[cpar,1] = -cerrrescube[cpar,0]
	    	
            endelse
          endif
    endfor
endif

if Nbfitpars gt 0 then begin
    dv = 0. ; initialize
    for i=0,Nbfitpars-1 do begin
        par = bfitp[i]+1
        bpar = fit.bparstart+i
        if (*fit.pars)[bpar] ne -999. and (*fit.sigpars)[bpar] ne -999. then begin
            brescube[par] = (*fit.pars)[bpar]
            if i eq 0 then dv=(*fit.pars)[bpar] ; record dv (1st parameter)
            if ~(par le 2 and state.selected_lineset eq 0 and state.pbfix[par] eq 1 and berrrescube[par,0] gt 0) then begin 
	     if state.domontecarlo gt 0  then begin                
                 percs_bootstrap = kubeviz_percentile(reform((*fit.pars_bootstrap)[*,bpar]),*state.montecarlo_percs)
                 berrrescube[par,*] = percs_bootstrap-(*fit.pars)[bpar]
                 if state.plotMonteCarlodistrib eq 1 then kubeviz_linefit_plotMonteCarlodistrib, reform((*fit.pars_bootstrap)[*,bpar]), percs_bootstrap, (*fit.pars)[bpar], 'broad', par
             endif else begin
                ; 1-sig errors in the first 2 elements of the error array for each parameter:                
                berrrescube[par,0] = (*fit.sigpars)[bpar]
                berrrescube[par,1] = -(*fit.sigpars)[bpar]
             endelse
	    endif 
        endif
; value of continuum = continuum fit evaluated at the wavelength of the line:
        if par gt 2 and Ncfitpars gt 0 then begin
            cpar = par - 2
            line = cpar - 1
            if state.continuumfit_mode eq 0 then begin
	      lambda_line = state.lines[line]*(1.+dv/state.ckms)
              crescube[cpar] = kubeviz_getcontatlambda((*fit.contpars),lambda_line)
            endif else crescube[cpar] = (*fit.pars)[fit.cpar]
	    if state.domontecarlo gt 0 then begin
                for boot=0,state.Nmontecarlo-1 do begin
		   if state.continuumfit_mode eq 0 then cont_lambda_bootstrap[boot,cont_boot_par] = kubeviz_getcontatlambda((*fit.contpars_bootstrap)[boot,*],lambda_line) $
		                                   else cont_lambda_bootstrap[boot,cont_boot_par] = (*fit.pars_bootstrap)[boot,fit.cpar]
                endfor
		percs_bootstrap = kubeviz_percentile(cont_lambda_bootstrap[*,cont_boot_par],*state.montecarlo_percs)
                cerrrescube[cpar,*] = percs_bootstrap-crescube[cpar]
                if state.plotMonteCarlodistrib eq 1 then kubeviz_linefit_plotMonteCarlodistrib, cont_lambda_bootstrap[*,cont_boot_par], percs_bootstrap, crescube[cpar], 'cont', cpar
		cont_boot_par += 1
            endif else begin
            if state.continuumfit_mode eq 0 then cerrrescube[cpar,0] = (*fit.sigcont) else cerrrescube[cpar,0] = (*fit.sigpars)[fit.cpar]
            cerrrescube[cpar,1] = -cerrrescube[cpar,0]
            endelse
        endif
    endfor
endif

; record in state variable:
case state.linefit_mode of
    0: begin ; spaxels:
        
        if Nnfitpars+Nbfitpars gt 0 then (*state.sp_chisq)[state.col,state.row] = fit.chisq
        
        (*state.sp_nrescube)[state.col,state.row,*] = nrescube
        (*state.sp_brescube)[state.col,state.row,*] = brescube
        (*state.sp_crescube)[state.col,state.row,*] = crescube 
        (*state.sp_nerrrescube)[state.col,state.row,*,*] = nerrrescube
        (*state.sp_berrrescube)[state.col,state.row,*,*] = berrrescube
        (*state.sp_cerrrescube)[state.col,state.row,*,*] = cerrrescube
                
    end
    1: begin ; masks:
         
         if Nnfitpars+Nbfitpars gt 0 then (*state.chisq)[state.imask-1] = fit.chisq
         
        (*state.nrescube)[state.imask-1,*] = nrescube
        (*state.brescube)[state.imask-1,*] = brescube
        (*state.crescube)[state.imask-1,*] = crescube 
        (*state.nerrrescube)[state.imask-1,*,*] = nerrrescube
        (*state.berrrescube)[state.imask-1,*,*] = berrrescube
        (*state.cerrrescube)[state.imask-1,*,*] = cerrrescube
        
    end
endcase
kubeviz_linefit_pointerswitch, /save

;Save the montecarlo distributions. 
if state.saveMonteCarlodistrib eq 1 and state.domontecarlo gt 0 then kubeviz_linefit_saveMonteCarlodistrib, (*fit.pars_bootstrap), cont_lambda_bootstrap

;MF free the pointers to the current fit.
;Needed to avoid thousands of unreferenced heap variables
heap_free, fit

end

;--------------------------------------
pro kubeviz_linefit_keepmom
common kubeviz_state
common kubeviz_fit

; copy the fit results into the results cube

mfitp = where(fit.mfitpars eq 1, Nmfitpars)
cfitp = where(fit.cfitpars eq 1, Ncfitpars)

kubeviz_linefit_pointerswitch, /load

; current values:
case state.linefit_mode of
    0: begin ; spaxels:
        mrescube = reform((*state.sp_mrescube)[state.col,state.row,*])
        merrrescube = reform((*state.sp_merrrescube)[state.col,state.row,*,*])
        crescube = reform((*state.sp_crescube)[state.col,state.row,*])
        cerrrescube = reform((*state.sp_cerrrescube)[state.col,state.row,*,*])
    end
    1: begin ; masks:
        mrescube = reform((*state.mrescube)[state.imask-1,*])
        merrrescube = reform((*state.merrrescube)[state.imask-1,*,*])
        crescube = reform((*state.crescube)[state.imask-1,*])
        cerrrescube = reform((*state.cerrrescube)[state.imask-1,*,*])
    end
endcase

if state.domontecarlo gt 0  then begin
    cont_lambda_bootstrap = replicate(0.D, state.Nmontecarlo, Nmfitpars)
    cont_boot_par = 0
endif

if Nmfitpars gt 0 then begin
    for i=0,Nmfitpars-1 do begin
        for j=0,4 do begin
          par = 6*mfitp[i]+j 
          npar= 5*i+j 
          mrescube[par] = (*fit.pars)[npar]
          if state.domontecarlo gt 0  then begin               
              percs_bootstrap = kubeviz_percentile(reform((*fit.pars_bootstrap)[*,npar]),*state.montecarlo_percs)
              merrrescube[par,*] = percs_bootstrap-(*fit.pars)[npar]
          endif else merrrescube[par,*] = (*fit.sigpars)[npar]
        endfor
        if Ncfitpars gt 0 then begin
            cpar = i+1 
            ; value of continuum = continuum fit evaluated at the wavelength of the line:
            lambda_line = state.lines[i]*(1.+mrescube[6*i+1]/state.ckms)
            crescube[cpar] = kubeviz_getcontatlambda((*fit.contpars),lambda_line)
            if state.domontecarlo gt 0  then begin
                for boot=0,state.Nmontecarlo-1 do cont_lambda_bootstrap[boot,cont_boot_par] = kubeviz_getcontatlambda((*fit.contpars_bootstrap)[boot,*],lambda_line)
		percs_bootstrap = kubeviz_percentile(cont_lambda_bootstrap[*,cont_boot_par],*state.montecarlo_percs)
                cerrrescube[cpar,*] = percs_bootstrap-crescube[cpar]
                if state.plotMonteCarlodistrib eq 1 then kubeviz_linefit_plotMonteCarlodistrib, cont_lambda_bootstrap[*,cont_boot_par], percs_bootstrap, crescube[cpar], 'cont', cpar
                cont_boot_par += 1
            endif else begin
                cerrrescube[cpar,0] = (*fit.sigcont)
                cerrrescube[cpar,1] = -cerrrescube[cpar,0]
            endelse
         endif
    endfor
endif

; record in state variable:
case state.linefit_mode of
    0: begin ; spaxels:
        (*state.sp_mrescube)[state.col,state.row,*] = mrescube
        (*state.sp_merrrescube)[state.col,state.row,*,*] = merrrescube
        (*state.sp_crescube)[state.col,state.row,*] = crescube
        (*state.sp_cerrrescube)[state.col,state.row,*,*] = cerrrescube
        
    end
    1: begin ; masks:
        (*state.mrescube)[state.imask-1,*] = mrescube
        (*state.merrrescube)[state.imask-1,*,*] = merrrescube
        (*state.crescube)[state.imask-1,*] = crescube
        (*state.cerrrescube)[state.imask-1,*,*] = cerrrescube
    end
endcase

kubeviz_linefit_pointerswitch, /save

;Save the montecarlo distributions 
if state.saveMonteCarlodistrib eq 1 and state.domontecarlo gt 0 then kubeviz_linefit_saveMonteCarlodistrib, (*fit.pars_bootstrap), cont_lambda_bootstrap

;MF free the pointers to the current fit.
;Needed to avoid thousands of unreferenced heap variables
heap_free, fit

end

;------------------------------------------------------------------------
pro kubeviz_fitcontinuum, contfitregion, spec, noise, wave, wave_cen, okcontfit, contpars, sigcontpars, sigcont
common kubeviz_state
common kubeviz_fit

; Select the appropriate region of the spectrum:
spec_contfitregion  = spec[contfitregion]
wave_contfitregion  = wave[contfitregion]
noise_contfitregion = noise[contfitregion]


if state.continuumfit_minperc eq 0 and state.continuumfit_maxperc eq 100 then begin ;inverse variance method
    
    weights = 1./noise_contfitregion^2    
; DJW 13-10-2016
; FIX: Set weights to 0 for |noise - median(noise)|  > 90%ile - 10%ile
; Deal with very noisy regions of data which, even with low weights, can
; influence continuum fit.
; ALSO deal with spikes in spectrum, which, for the sake of continuum
; fitting, can be dealt with in exactly the same way.
    okval = where(finite(weights) eq 1 and abs(noise_contfitregion - median(noise_contfitregion)) le (kubeviz_percentile(noise_contfitregion,[90.]))[0]-(kubeviz_percentile(noise_contfitregion,[10.]))[0] and abs(spec_contfitregion - median(spec_contfitregion)) le (kubeviz_percentile(spec_contfitregion,[90.]))[0]-(kubeviz_percentile(spec_contfitregion,[10.]))[0], Nok)
    
    if Nok gt 0 then begin
    	; fit it with a polynomial of the selected order
    	estimates = replicate(0.,state.continuumfit_order+1)
    	estimates[0] = mean(spec_contfitregion[okval])
    	contpars  = svdfit(wave_contfitregion[okval], spec_contfitregion[okval] , state.continuumfit_order+1, A=estimates, chisq=chisq, sigma=sigcontpars, weights=weights[okval], status=status)
    	polyspec  = kubeviz_getcontatlambda(contpars, wave_contfitregion[okval])
    	sigcont   =  0.5*(kubeviz_percentile(spec_contfitregion[okval]-polyspec,[84])+abs(kubeviz_percentile(spec_contfitregion[okval]-polyspec,[16]))) 
    endif else begin
        contpars	= replicate(0.,state.continuumfit_order+1)
        sigcontpars = replicate(0.,state.continuumfit_order+1)
        sigcont	= 0.
    endelse

endif else begin ;SDSS method

   noisecut_perc = 20

   new_range = (state.continuumfit_maxperc-state.continuumfit_minperc)/(0.01*(100-noisecut_perc))
   mid_perc = .5*(state.continuumfit_minperc+state.continuumfit_maxperc)
   min_perc = mid_perc-(0.5*new_range)
   max_perc = mid_perc+(0.5*new_range)

   ; BLUE SIDE

   blue = where(wave_contfitregion lt wave_cen, Nblue)
   if Nblue gt 0 then begin 
    ; First reject the highest 20 percent of noise channels (sky lines) 
    noisepercrange = kubeviz_percentile(noise_contfitregion[blue],100-noisecut_perc)
    oknoiseb = where(noise_contfitregion[blue] lt noisepercrange[0], Nokcontfitb)
    if Nokcontfitb gt 0 then begin
     spec_b = (spec_contfitregion[blue])[oknoiseb]
     wave_b = (wave_contfitregion[blue])[oknoiseb]
     ; Select the percentile range
     specpercrange = kubeviz_percentile(spec_b,[min_perc,max_perc])
     okpercb = where(spec_b ge specpercrange[0] and spec_b le specpercrange[1], Nokcontfitb)
    endif
    if Nokcontfitb gt 0 then begin
      spec_b = spec_b[okpercb]
      wave_b = wave_b[okpercb]
    endif
   endif else begin
    Nokcontfitb = 0
   endelse

   ; RED SIDE

   red = where(wave_contfitregion gt wave_cen, Nred)
   if Nred gt 0 then begin 
    ; First reject the highest 20 percent of noise channels (sky lines) 
    noisepercrange = kubeviz_percentile(noise_contfitregion[red],100-noisecut_perc)
    oknoiser = where(noise_contfitregion[red] lt noisepercrange[0], Nokcontfitr)
    if Nokcontfitr gt 0 then begin
     spec_r = (spec_contfitregion[red])[oknoiser]
     wave_r = (wave_contfitregion[red])[oknoiser]
     ; Select the percentile range
     specpercrange = kubeviz_percentile(spec_r,[min_perc,max_perc])
     okpercr = where(spec_r ge specpercrange[0] and spec_r le specpercrange[1], Nokcontfitr)
    endif
    if Nokcontfitr gt 0 then begin
      spec_r = spec_r[okpercr]
      wave_r = wave_r[okpercr]
    endif 
   endif else begin
    Nokcontfitr = 0
   endelse

   Nokcontfit = Nokcontfitb + Nokcontfitr

   if Nokcontfit gt 0 then begin
    
    case 1 of
     Nokcontfitb gt 0 and Nokcontfitr eq 0 : begin
       spec_fit = spec_b
       wave_fit = wave_b
     end
     Nokcontfitb eq 0 and Nokcontfitr gt 0 : begin
       spec_fit = spec_r
       wave_fit = wave_r
     end
     Nokcontfitb gt 0 and Nokcontfitr gt 0 : begin
       spec_fit = [spec_b,spec_r]
       wave_fit = [wave_b,wave_r]
     end
    endcase 

    ; fit it with a polynomial of the selected order
    estimates = replicate(0.,state.continuumfit_order+1)
    estimates[0] = mean(spec_fit)
    contpars  = svdfit(wave_fit, spec_fit, state.continuumfit_order+1, A=estimates, chisq=chisq, sigma=sigcontpars, status=status)
    polyspec  = kubeviz_getcontatlambda(contpars, wave_fit)
    sigcont   =  0.5*(kubeviz_percentile(spec_fit-polyspec,[84])+abs(kubeviz_percentile(spec_fit-polyspec,[16]))) 

   endif else begin 
    contpars	= replicate(0.,state.continuumfit_order+1)
    sigcontpars = replicate(0.,state.continuumfit_order+1)
    sigcont	= 0.
   endelse

endelse ;end SDSS method

end

;--------------------------------------------------------------------------
pro kubeviz_definecontinuumfitregion, wave, linesinset, wave_cen, Ncont
common kubeviz_state
common kubeviz_fit

; centre continuum fitting region on the main line of the linesinset:
contfitregion_cen = wave_cen

; default offsets from the globally set parameters from the linefit
; window:
contfitregion_minoff = state.continuumfit_minoff
contfitregion_maxoff = state.continuumfit_maxoff

; Define 2 regions: one blue- and one red-ward of the linesinset:
contfitregion_lims = [ [(contfitregion_cen - contfitregion_maxoff), (contfitregion_cen - contfitregion_minoff)], [(contfitregion_cen + contfitregion_minoff), (contfitregion_cen + contfitregion_maxoff)] ]

; array of pixels which are within these limits:
contfitregion_nocuts = where( ( wave ge contfitregion_lims[0,0] and wave le contfitregion_lims[1,0] ) or ( wave ge contfitregion_lims[0,1] and wave le contfitregion_lims[1,1] ), Ncontfitregion_nocuts )

; Define regions which contain any line; which are too close to the
; ends of the spectrum; or which are forcibly excluded from any fit. 
;; Djw 18-10-16: Added Kubeviz_Excludewave Routine / Call
; Start with the whole allowed fitting range. We then avoid all the regions where a line is expected (also those lines not in the current lineset) 
contfitregion_suitable = where (wave ge (*state.wave)[0]+10 and wave le (*state.wave)[state.Nwpix-1]-10 and (kubeviz_excludewave(wave) eq 0) )
for i=0, state.Nlines_all-1 do if N_elements(contfitregion_suitable) gt 1 then $
    contfitregion_suitable = kubeviz_SetDifference(contfitregion_suitable, where (wave ge state.lines_all[i]-state.maxwoffb and wave le state.lines_all[i]+state.maxwoffr))

; Exclude these from the continuum fitting region:
contfitregion = kubeviz_SetIntersection(contfitregion_nocuts, contfitregion_suitable) 
dummy = where(contfitregion ne -1, Ncont)

; save to the fit structure:
fit.contfitregion = ptr_new(contfitregion)

end

;------------------------------------------------------------------------------
pro kubeviz_linefit_fitset, wave, spec, noise, linesinset, spec_bootstrap, keep, firstset, col=col, row=row, silent=silent
common kubeviz_state
common kubeviz_fit

if n_elements(silent) eq 0 then silent = 0

; which lines to fit?
if linesinset[0] ne -1 then begin
    Nlinesinset = n_elements(linesinset)

; max. range to fit around main line (0th in linesinset), or centered
; between two lines for linesets with an even number of lines
    wmainline = median([state.lines[linesinset]],/even)    
    fit.maxwoffb = state.maxwoffb
    fit.maxwoffr = state.maxwoffr
    wmin = wmainline-fit.maxwoffb
    wmax = wmainline+fit.maxwoffr
; determine pixels in wavelength axis to fit. 
    useline = where(state.pndofit[linesinset]+state.pbdofit[linesinset] ge 1, Nuseline)
    if Nuseline gt 0 then begin
        bisectorrange = [ min(state.line_bisectors[linesinset(useline),0]) , max(state.line_bisectors[linesinset(useline),1]) ]
        wmin = max([bisectorrange[0],wmin])
        wmax = min([bisectorrange[1],wmax])
    endif

    ;; DJW 18-10-16: forcibly exclude hardcoded regions of spectrum
    ;; using kubeviz_excludewave call
    okk = where(wave ge wmin and wave le wmax and (kubeviz_excludewave(wave) eq 0), Nokk)
    
    if Nokk gt 0 then begin
       ; build array of narrow and broad line pars to fit (set to 1):
       nfitpars = replicate(0B,2+state.Nlines)
       bfitpars = replicate(0B,2+state.Nlines)
       cfitpars = 0B

       for iline=0, Nlinesinset-1 do begin
           line = linesinset[iline]
           ; flux parameter in narrow / broad-line components:
           npar = line+2
           bpar = line+2
           if state.pndofit[line] eq 1 then nfitpars[npar]=1
           if state.pbdofit[line] eq 1 then bfitpars[bpar]=1
       endfor
       Nnfitlines = total(nfitpars)
       Nbfitlines = total(bfitpars) 
      
       ; First check if we need to bother fitting any lines in the set:
       if Nnfitlines+Nbfitlines gt 0 then begin

           ; if we are fitting any lines then we need to fit the line offset, width
           if Nnfitlines gt 0 then nfitpars[indgen(2)] = 1
           if Nbfitlines gt 0 then bfitpars[indgen(2)] = 1

           ; and continuum if requested:
           cfitpars = (state.pcdofit[linesinset[0]] eq 1)

           ; record parameters to fit:
           fit.nfitpars = nfitpars
           fit.bfitpars = bfitpars
           fit.cfitpars = cfitpars

           fit.okk = ptr_new(okk)
           fit.wave = ptr_new(wave[okk])
           fit.spec = ptr_new(spec[okk])
           fit.noise = ptr_new(noise[okk])
           fit.spec_bootstrap = ptr_new(spec_bootstrap[*,okk])

           ; Arrays containing all and only the parameters to fit:     
           ; e.g. for >= 1 narrow line: 			       
           ; 0: continuum					       
           ; 1: mean_narrow					       
           ; 2: width_narrow					       
           ; 3 - (2+Nnfitlines): amplitude_narrow		       
           ; 3+Nnfitlines: mean_broad				       
           ; 4+Nnfitlines: width_broad  			       
           ; 5+Nnfitlines - (4+Nnfitlines+Nbfitlines): amplitude_broad 
           
	   fit.cpar = 0
           if Nnfitlines gt 0 then begin
               fit.nparstart = 1
               fit.nparend = fit.nparstart+Nnfitlines+1
           endif else begin
               fit.nparstart = -1
               fit.nparend = -1
           endelse
           if Nbfitlines gt 0 then begin
               if fit.nparend eq -1 then fit.bparstart = 1 else fit.bparstart = fit.nparend+1
               fit.bparend = fit.bparstart+Nbfitlines+1
           endif else begin
               fit.bparstart = -1
               fit.bparend = -1
           endelse 
           
            ; define continuum regions (exclude spec ends and other lines) & save in fit struct
           if cfitpars eq 1 then begin
               kubeviz_definecontinuumfitregion, wave, linesinset, wmainline, Ncont
	       if Ncont gt 0 then begin
            ; then fit it
                 kubeviz_fitcontinuum, (*fit.contfitregion), spec, noise, wave, wmainline, okcontfit, contpars, sigcontpars, sigcont
            ; continuum spectrum:
                 contsub = kubeviz_getcontatlambda(contpars,wave[okk])
            ; save to fit structure:
                 fit.okcontfit   = ptr_new(okcontfit)
                 fit.contpars    = ptr_new(contpars)
                 fit.sigcontpars = ptr_new(sigcontpars)
                 fit.sigcont     = ptr_new(sigcont)
		 
            ; if necessary, repeat for bootstrap spectra:
                 if state.domontecarlo gt 0 then begin
        	     contpars_bootstrap = replicate(-999.D,state.Nmontecarlo,state.continuumfit_order+1)
        	     contsub_bootstrap = replicate(0.D,state.Nmontecarlo,Nokk)
        	     for boot=0,state.Nmontecarlo-1 do begin
        	         spec_boot = reform(spec_bootstrap[boot,*])
			 kubeviz_fitcontinuum, (*fit.contfitregion), spec_boot, noise, wave, wmainline, okcontfit_bootstrap, contpars_boot, sigcontpars_boot, sigcont_boot
			 contpars_bootstrap[boot,*] = contpars_boot
        	         contsub_bootstrap[boot,*]  = kubeviz_getcontatlambda(contpars_boot,wave[okk])
          	     endfor
        	     fit.contpars_bootstrap = ptr_new(contpars_bootstrap)
                 endif
                endif else begin ;Ncont = 0
		   if state.continuumfit_mode eq 0 then begin ;If SDSS mode raise a warning
		     printf, state.log_lun, '[WARNING] No suitable continuum region found. No fitting takes place.'
		     return
		   endif else begin ;In internal continuum mode continue
		     fit.contpars = ptr_new(replicate(0.,5))
		     contsub      = replicate(0.,Nokk)
		   endelse  
		endelse 
              endif else begin 
              contsub = replicate(0.,Nokk) ; no continuum subtracted. If necessary return it also for montecarlo spectra
              if state.domontecarlo gt 0 then contsub_bootstrap = replicate(0.D,state.Nmontecarlo,Nokk)
           endelse
           
	   ; inverse variance weighting 
	   ; MF How can we add the continuum subtraction uncertainty?
	   
	   weights = 1.D/((*fit.noise)^2) 
           
	   ; Deal with zero noise pixels (e.g. bad pixels): set weights to zero.
           nonoise = where((*fit.noise) eq 0, Nnonoise)
           if Nnonoise gt 0 then weights[nonoise] = 0.
           fit.weights = ptr_new(weights)

          ; here we split between the two fit types (gaussian fit or moments)
	  ; Gaussian:	  
           if state.linefit_type eq 0 then begin 
          
	  ; choose the starting values:
             kubeviz_linefit_getstartvals, col=col, row=row, silent=silent

          ; choose the parameter limits:
             kubeviz_linefit_setparinfo
  
          ; fit these parameters:
             if state.continuumfit_mode eq 0 then kubeviz_linefit_fit, (*fit.wave), (*fit.spec)-contsub, (*fit.weights), (*fit.gpars), pars, /savepars $
             else                                 kubeviz_linefit_fit, (*fit.wave), (*fit.spec),         (*fit.weights), (*fit.gpars), pars, /savepars
	     
	     
	     ;secondcomp_smart is an advanced rejection routine to identify if two component fits are actually necessary. It should be 
	     ;used only to identify two truly separated kinematical components. Do not set to find broad undelying lines.
	     ;Criteria: ratio of sigmas between 0.25 and 4, separation of centriods gt 35 km/s, both mainlines have SN>3
	     ;If the criteria are not satisfied the fit goes back to a single component fit
	     if Nbfitlines gt 0 and firstset eq 0 and state.secondcomp_smart eq 1 then begin
		vel = [(*fit.pars)[fit.nparstart], (*fit.pars)[fit.bparstart]]
		sig = [(*fit.pars)[fit.nparstart+1], (*fit.pars)[fit.bparstart+1]]
		sn = [(*fit.pars)[fit.nparstart+2]/(*fit.sigpars)[fit.nparstart+2], (*fit.pars)[fit.bparstart+2]/(*fit.sigpars)[fit.bparstart+2]]
		if max(sig) gt abs(vel[0]-vel[1])/2. or abs(vel[0]-vel[1]) lt 75 or min(sn) lt 3 or min(finite(sn)) eq 0 then begin
		   
		   ;Remove broad lines from fit structure
		   
		   Nbfitlines = 0
		   fit.bfitpars[*] = 0
		   fit.bparstart = -1
                   fit.bparend = -1
	  	   
		   ; choose the starting values:
          	      kubeviz_linefit_getstartvals, col=col, row=row, silent=silent

          	   ; choose the parameter limits:
          	      kubeviz_linefit_setparinfo
  
          	   ; fit these parameters:
          	      if state.continuumfit_mode eq 0 then kubeviz_linefit_fit, (*fit.wave), (*fit.spec)-contsub, (*fit.weights), (*fit.gpars), pars, /savepars $
          	      else				   kubeviz_linefit_fit, (*fit.wave), (*fit.spec),	  (*fit.weights), (*fit.gpars), pars, /savepars
		endif
	     endif
	     
	     if state.domontecarlo gt 0 then begin
        	 pars_bootstrap = replicate(-999.D,state.Nmontecarlo,N_Elements(pars)) 
                 ; repeat for each bootstrap iteration:
        	 for boot=0,state.Nmontecarlo-1 do begin
        	     if state.continuumfit_mode eq 0 then spec_boot = reform((*fit.spec_bootstrap)[boot,*]) - reform(contsub_bootstrap[boot,*]) $
		                                     else spec_boot = reform((*fit.spec_bootstrap)[boot,*])
		     kubeviz_linefit_fit, (*fit.wave), spec_boot, (*fit.weights), (*fit.gpars), pars_boot
		     pars_bootstrap[boot,*] = pars_boot
        	 endfor
        	 fit.pars_bootstrap = ptr_new(pars_bootstrap)
             endif
             
             ; fix the line position and width for following linesets (firstset becomes 1 once the first valid lineset was fit):
	     if firstset eq 0 then begin 
	        firstset=1
		fit.firstset =1
	     endif else fit.firstset = 0	
	     
	     if Nbfitlines eq 0 and Nnfitlines gt 0 then linetypes = ['N']
             if Nnfitlines eq 0 and Nbfitlines gt 0 then linetypes = ['B']
             if Nnfitlines gt 0 and Nbfitlines gt 0 then linetypes = ['N','B']
             kubeviz_linefit_setfixline, linetypes
           
          ; propagate parameters to results cube if keep =1
             if keep eq 1 then kubeviz_linefit_keepfit
           
         endif else begin ; if not gaussian then fit moments:

          ; reshape mfitpars to match the results container
             fit.mfitpars = fit.nfitpars
             fit.mfitpars[0:1] = 0
             fit.mfitpars = shift(fit.mfitpars,-2)

          ; first check if any "narrow" line is selected, here narrow read as a line in general   
             if Nnfitlines gt 0 then begin
               kubeviz_linefit_moments, (*fit.wave), (*fit.spec)-contsub, (*fit.weights), pars, /savepars, col=col, row=row
               if state.domontecarlo gt 0 then begin
        	 pars_bootstrap = replicate(-999.D,state.Nmontecarlo,N_Elements(pars)) 
           ; repeat for each bootstrap iteration:
        	 for boot=0,state.Nmontecarlo-1 do begin
        	     spec_boot = reform((*fit.spec_bootstrap)[boot,*]) - reform(contsub_bootstrap[boot,*])
        	     kubeviz_linefit_moments, (*fit.wave), spec_boot, (*fit.weights), pars_boot, col=col, row=row, boot=boot
        	     pars_bootstrap[boot,*] = pars_boot
        	 endfor
        	 fit.pars_bootstrap = ptr_new(pars_bootstrap)
               endif
               if keep eq 1 then kubeviz_linefit_keepmom
             endif
           endelse
       endif
    endif
endif

end

;-----------------------------------------------------------------------------
function kubeviz_linefit_testspec, col=col, row=row
common kubeviz_state

if N_elements(col) eq 0 then col = state.col
if N_elements(row) eq 0 then row = state.row 

; get spectrum and wavelength vector
case state.specmode of
    0: begin
        spec = (*state.datacube)[col,row,*]
    end
    1: begin
        kubeviz_medianspec, /sum
        spec = *state.medspec
    end
    2: begin        
        kubeviz_medianspec
        spec = *state.medspec
    end
    3: begin        
        kubeviz_medianspec, /wavg
        spec = *state.medspec
    end
    4: begin
        ; median subtracted: not relevant for line-fitting so fit whole spectrum
        spec = (*state.datacube)[col,row,*]
    end
    5: begin        
        kubeviz_medianspec, /optimal
        spec = *state.medspec
    end
endcase

; select out part of spectrum which is not missing / bad data):
okspec = where(finite(spec) eq 1., Nokspec)

if Nokspec gt 0 then return, 1 else return, 0

end


;-----------------------------------------------------------------------------
pro kubeviz_linefit_dofit, userpar=userpar, col=col, row=row, nosave=nosave
common kubeviz_state

if N_elements(col) eq 0 then col = state.col
if N_elements(row) eq 0 then row = state.row 
if N_elements(nosave) eq 0 then keep = 1 else keep = 0 

if state.domontecarlo gt 0 then spec_bootstrap = replicate(0.,state.Nmontecarlo,state.Nwpix) else spec_bootstrap = replicate(-1.,1,state.Nwpix)

; get spectrum and wavelength vector
case state.specmode of
    0: begin
        spec = (*state.datacube)[col,row,*]
        noise = (*state.noise)[col,row,*]
        if state.domontecarlo gt 0 and total(spec, /nan) ne 0 then for i=0,state.Nmontecarlo-1 do spec_bootstrap[i,*] = (*state.montecarlocubes[i])[col,row,*]
    end
    1: begin
        kubeviz_medianspec, /sum
        spec = *state.medspec
        noise = *state.nmedspec
        if state.domontecarlo gt 0 and total(spec, /nan) ne 0 then for i=0,state.Nmontecarlo-1 do spec_bootstrap[i,*] = *state.medspec_montecarlo[i]
    end
    2: begin        
        kubeviz_medianspec
        spec = *state.medspec
        noise = *state.nmedspec
        if state.domontecarlo gt 0 and total(spec, /nan) ne 0 then for i=0,state.Nmontecarlo-1 do spec_bootstrap[i,*] = *state.medspec_montecarlo[i]
    end
    3: begin        
        kubeviz_medianspec, /wavg
        spec = *state.medspec
        noise = *state.nmedspec
        if state.domontecarlo gt 0 and total(spec, /nan) ne 0 then for i=0,state.Nmontecarlo-1 do spec_bootstrap[i,*] = *state.medspec_montecarlo[i]
    end
    4: begin
        ; median subtracted: not relevant for line-fitting so fit whole spectrum
        printf, state.log_lun, 'Median subtracted spectrum cannot be fit: Fitting whole spectrum'
        spec = (*state.datacube)[col,row,*]
        noise = (*state.noise)[col,row,*] 
        if state.domontecarlo gt 0 and total(spec, /nan) ne 0 then for i=0,state.Nmontecarlo-1 do spec_bootstrap[i,*] = (*state.montecarlocubes[i])[col,row,*] 
    end
    5: begin        
        kubeviz_medianspec, /optimal
        spec = *state.medspec
        noise = *state.nmedspec
        if state.domontecarlo gt 0 and total(spec, /nan) ne 0 then for i=0,state.Nmontecarlo-1 do spec_bootstrap[i,*] = *state.medspec_montecarlo[i]
    end
endcase

; select out part of spectrum not NaN
okspec = where(finite(spec) eq 1., Nokspec)

; FIXME! Needs to reset the results container but NOT the widgets
; unless we are definitely running in interactive mode!!
;kubeviz_linefit_reset, update=0

if Nokspec gt 0 then begin
    
  firstset=0

  ; save status of the fix button and start parameters for wavelength
  ; parameters (to be reset after fitting all lines):
  pnfix = state.pnfix
  gnfit = state.gnfit
  pbfix = state.pbfix
  gbfit = state.gbfit


; We need to fit all selected lines, but we will do so in linesets.
; First fit a lineset including line widths and offsets (e.g. the default is Halpha + [NII]). 
; IF a line width and offset has already been computed for a previous
; set, they are kept fixed.
; Note: use of okspec for the bootstrap spectrum *might* break if the
;       invalid data is different in the bootstrap cube.
  
  for set=0, state.lineset_max-1 do kubeviz_linefit_fitset, (*state.wave)[okspec], spec[okspec], noise[okspec], (*state.lineset_ind[set]),  spec_bootstrap[*,okspec], keep, firstset, col=col, row=row
  
  ; reset the fix button and start parameters:
  state.pnfix = pnfix
  state.gnfit = gnfit
  state.pbfix = pbfix
  state.gbfit = gbfit

endif

end

;----------------------------------------------------
pro kubeviz_linefit_pointerswitch, load=load, save=save, linefitmode=linefitmode, errmethod=errmethod

common kubeviz_state

; If not set use state variables
if N_elements(linefitmode) eq 0 then linefitmode = state.linefit_mode
if N_elements(errmethod)   eq 0 then errmethod   = state.domontecarlo

if N_elements(load) gt 0 then begin
   case linefitmode of
        0: begin ; spaxels
        case errmethod of
           0: begin
              state.sp_nrescube = state.noi_sp_nrescube
              state.sp_brescube = state.noi_sp_brescube
              state.sp_crescube = state.noi_sp_crescube
              state.sp_mrescube = state.noi_sp_mrescube
              state.sp_nerrrescube = state.noi_sp_nerrrescube
              state.sp_berrrescube = state.noi_sp_berrrescube
              state.sp_cerrrescube = state.noi_sp_cerrrescube
              state.sp_merrrescube = state.noi_sp_merrrescube
           end
           1: begin
              state.sp_nrescube = state.boot_sp_nrescube
              state.sp_brescube = state.boot_sp_brescube
              state.sp_crescube = state.boot_sp_crescube
              state.sp_mrescube = state.boot_sp_mrescube
              state.sp_nerrrescube = state.boot_sp_nerrrescube
              state.sp_berrrescube = state.boot_sp_berrrescube
              state.sp_cerrrescube = state.boot_sp_cerrrescube
              state.sp_merrrescube = state.boot_sp_merrrescube
           end
           2: begin
              state.sp_nrescube = state.mc1_sp_nrescube
              state.sp_brescube = state.mc1_sp_brescube
              state.sp_crescube = state.mc1_sp_crescube
              state.sp_mrescube = state.mc1_sp_mrescube
              state.sp_nerrrescube = state.mc1_sp_nerrrescube
              state.sp_berrrescube = state.mc1_sp_berrrescube
              state.sp_cerrrescube = state.mc1_sp_cerrrescube
              state.sp_merrrescube = state.mc1_sp_merrrescube
           end
           3: begin
              state.sp_nrescube = state.mc2_sp_nrescube
              state.sp_brescube = state.mc2_sp_brescube
              state.sp_crescube = state.mc2_sp_crescube
              state.sp_mrescube = state.mc2_sp_mrescube
              state.sp_nerrrescube = state.mc2_sp_nerrrescube
              state.sp_berrrescube = state.mc2_sp_berrrescube
              state.sp_cerrrescube = state.mc2_sp_cerrrescube
              state.sp_merrrescube = state.mc2_sp_merrrescube
           end 
           4: begin
              state.sp_nrescube = state.mc3_sp_nrescube
              state.sp_brescube = state.mc3_sp_brescube
              state.sp_crescube = state.mc3_sp_crescube
              state.sp_mrescube = state.mc3_sp_mrescube
              state.sp_nerrrescube = state.mc3_sp_nerrrescube
              state.sp_berrrescube = state.mc3_sp_berrrescube
              state.sp_cerrrescube = state.mc3_sp_cerrrescube
              state.sp_merrrescube = state.mc3_sp_merrrescube
           end
        endcase
        end   
        1: begin ; masks:
        case errmethod of
           0: begin
              state.nrescube = state.noi_nrescube
              state.brescube = state.noi_brescube
              state.crescube = state.noi_crescube
              state.mrescube = state.noi_mrescube
              state.nerrrescube = state.noi_nerrrescube
              state.berrrescube = state.noi_berrrescube
              state.cerrrescube = state.noi_cerrrescube
              state.merrrescube = state.noi_merrrescube
           end
           1: begin
              state.nrescube = state.boot_nrescube
              state.brescube = state.boot_brescube
              state.crescube = state.boot_crescube
              state.mrescube = state.boot_mrescube
              state.nerrrescube = state.boot_nerrrescube
              state.berrrescube = state.boot_berrrescube
              state.cerrrescube = state.boot_cerrrescube
              state.merrrescube = state.boot_merrrescube
           end
           2: begin
              state.nrescube = state.mc1_nrescube
              state.brescube = state.mc1_brescube
              state.crescube = state.mc1_crescube
              state.mrescube = state.mc1_mrescube
              state.nerrrescube = state.mc1_nerrrescube
              state.berrrescube = state.mc1_berrrescube
              state.cerrrescube = state.mc1_cerrrescube
              state.merrrescube = state.mc1_merrrescube
           end
           3: begin
              state.nrescube = state.mc2_nrescube
              state.brescube = state.mc2_brescube
              state.crescube = state.mc2_crescube
              state.mrescube = state.mc2_mrescube
              state.nerrrescube = state.mc2_nerrrescube
              state.berrrescube = state.mc2_berrrescube
              state.cerrrescube = state.mc2_cerrrescube
              state.merrrescube = state.mc2_merrrescube
           end
           4: begin
              state.nrescube = state.mc3_nrescube
              state.brescube = state.mc3_brescube
              state.crescube = state.mc3_crescube
              state.mrescube = state.mc3_mrescube
              state.nerrrescube = state.mc3_nerrrescube
              state.berrrescube = state.mc3_berrrescube
              state.cerrrescube = state.mc3_cerrrescube
              state.merrrescube = state.mc3_merrrescube
           end
          endcase
          end
   endcase
endif

if N_elements(save) gt 0 then begin
   case linefitmode of
        0: begin ; spaxels
        case errmethod of
           0: begin
             state.noi_sp_nrescube     = state.sp_nrescube         
             state.noi_sp_brescube     = state.sp_brescube         
             state.noi_sp_crescube     = state.sp_crescube
             state.noi_sp_mrescube     = state.sp_mrescube         
             state.noi_sp_nerrrescube  = state.sp_nerrrescube      
             state.noi_sp_berrrescube  = state.sp_berrrescube      
             state.noi_sp_cerrrescube  = state.sp_cerrrescube
             state.noi_sp_merrrescube  = state.sp_merrrescube      
           end
           1: begin
             state.boot_sp_nrescube     = state.sp_nrescube        
             state.boot_sp_brescube     = state.sp_brescube        
             state.boot_sp_crescube     = state.sp_crescube
             state.boot_sp_mrescube     = state.sp_mrescube        
             state.boot_sp_nerrrescube  = state.sp_nerrrescube     
             state.boot_sp_berrrescube  = state.sp_berrrescube     
             state.boot_sp_cerrrescube  = state.sp_cerrrescube  
             state.boot_sp_merrrescube  = state.sp_merrrescube     
           end
           2: begin
             state.mc1_sp_nrescube     = state.sp_nrescube         
             state.mc1_sp_brescube     = state.sp_brescube         
             state.mc1_sp_crescube     = state.sp_crescube
             state.mc1_sp_mrescube     = state.sp_mrescube         
             state.mc1_sp_nerrrescube  = state.sp_nerrrescube      
             state.mc1_sp_berrrescube  = state.sp_berrrescube      
             state.mc1_sp_cerrrescube  = state.sp_cerrrescube
             state.mc1_sp_merrrescube  = state.sp_merrrescube      
           end
           3: begin
             state.mc2_sp_nrescube     = state.sp_nrescube         
             state.mc2_sp_brescube     = state.sp_brescube         
             state.mc2_sp_crescube     = state.sp_crescube
             state.mc2_sp_mrescube     = state.sp_mrescube         
             state.mc2_sp_nerrrescube  = state.sp_nerrrescube      
             state.mc2_sp_berrrescube  = state.sp_berrrescube      
             state.mc2_sp_cerrrescube  = state.sp_cerrrescube
             state.mc2_sp_merrrescube  = state.sp_merrrescube
            end 
           4: begin
             state.mc3_sp_nrescube     = state.sp_nrescube         
             state.mc3_sp_brescube     = state.sp_brescube         
             state.mc3_sp_crescube     = state.sp_crescube
             state.mc3_sp_mrescube     = state.sp_mrescube         
             state.mc3_sp_nerrrescube  = state.sp_nerrrescube      
             state.mc3_sp_berrrescube  = state.sp_berrrescube      
             state.mc3_sp_cerrrescube  = state.sp_cerrrescube
             state.mc3_sp_merrrescube  = state.sp_merrrescube      
           end
        endcase
        end   
        1: begin ; masks:
        case errmethod of
           0: begin
             state.noi_nrescube     = state.nrescube       
             state.noi_brescube     = state.brescube       
             state.noi_crescube     = state.crescube
             state.noi_mrescube     = state.mrescube
             state.noi_nerrrescube  = state.nerrrescube    
             state.noi_berrrescube  = state.berrrescube    
             state.noi_cerrrescube  = state.cerrrescube
             state.noi_merrrescube  = state.merrrescube    
           end
           1: begin
             state.boot_nrescube     = state.nrescube      
             state.boot_brescube     = state.brescube      
             state.boot_crescube     = state.crescube
             state.boot_mrescube     = state.mrescube      
             state.boot_nerrrescube  = state.nerrrescube           
             state.boot_berrrescube  = state.berrrescube           
             state.boot_cerrrescube  = state.cerrrescube
             state.boot_merrrescube  = state.merrrescube           
           end
           2: begin
             state.mc1_nrescube     = state.nrescube       
             state.mc1_brescube     = state.brescube       
             state.mc1_crescube     = state.crescube
             state.mc1_mrescube     = state.mrescube
             state.mc1_nerrrescube  = state.nerrrescube    
             state.mc1_berrrescube  = state.berrrescube    
             state.mc1_cerrrescube  = state.cerrrescube
             state.mc1_merrrescube  = state.merrrescube    
           end
           3: begin
             state.mc2_nrescube     = state.nrescube       
             state.mc2_brescube     = state.brescube       
             state.mc2_crescube     = state.crescube
             state.mc2_mrescube     = state.mrescube       
             state.mc2_nerrrescube  = state.nerrrescube    
             state.mc2_berrrescube  = state.berrrescube    
             state.mc2_cerrrescube  = state.cerrrescube
             state.mc2_merrrescube  = state.merrrescube
            end 
           4: begin
             state.mc3_nrescube     = state.nrescube       
             state.mc3_brescube     = state.brescube       
             state.mc3_crescube     = state.crescube
             state.mc3_crescube     = state.mrescube       
             state.mc3_nerrrescube  = state.nerrrescube    
             state.mc3_berrrescube  = state.berrrescube    
             state.mc3_cerrrescube  = state.cerrrescube
             state.mc3_merrrescube  = state.merrrescube    
           end
          endcase
        end
   endcase        
endif

end

;------------------------------------------------------------
function kubeviz_linefit_flagsigclip, array, mode
common kubeviz_state

case mode of  ;mode is 1 for Line offset and 2 for Line width
  1 : okfit = where(array ne 0 and array ne state.gnlims[1,0] and array ne state.gnlims[1,1], Nok)
  2 : okfit = where(array gt 0 and array ne state.gnlims[2,0] and array ne state.gnlims[2,1], Nok)
  else: print, '[BUG]'
endcase

if Nok gt 1 then begin
   arr = array[okfit]
   percs = kubeviz_percentile(arr,[16,50,84])
   lowerlim = percs[1] - 4 *  (percs[1]-percs[0])
   upperlim = percs[1] + 4 *  (percs[2]-percs[1])
   return, [lowerlim, upperlim]
endif else return, [-1E9, 1E9]
     
end

;------------------------------------------------------------
pro kubeviz_linefit_autoflag, doflag=doflag, flagset=newflag
; flag as bad all spaxels or masks which have no line with flux S/N >=limit
; doflag is a binary mask to apply to say which spaxels to flag

; Flag = sum of binary flags:
; 16 x manual flag
;  8 x velocity and dispersion within 4 sigma of median values for
;      sigma=|16/84%ile - 50%ile|
;  4 x valid velocity and dispersion values
;  2 x errors on velocity and dispersion below threshold error
;  1 x S/N >= S/N_thresh

common kubeviz_state

; doflag is a mask to apply the autoflagging to only certain spaxels
; newflag returns new flag variable in case it is needed locally

if n_elements(doflag) eq 0 then doflag = replicate(1,state.Ncol,state.Nrow)

sn_thresh = state.mask_sn_thresh
maxvelerr = state.mask_maxvelerr

kubeviz_linefit_pointerswitch, /load

case state.linefit_type of
  0: begin ; Gaussian fit
  case state.linefit_mode of
    0: begin ; spaxels
        
        ; mean one-sigma error is half of the
        ; DJW 12-10-2016
        ; NEW: distance from -1sig ONLY relevant to determine S/N
        nline_onesig = abs(reform((*state.sp_nerrrescube)[*,*,3+indgen(state.Nlines),1]))
        bline_onesig = abs(reform((*state.sp_berrrescube)[*,*,3+indgen(state.Nlines),1]))

        sn_nline = reform((*state.sp_nrescube)[*,*,3+indgen(state.Nlines)])/nline_onesig
        sn_bline = reform((*state.sp_brescube)[*,*,3+indgen(state.Nlines)])/bline_onesig
        
	nan = where(finite(sn_nline) eq 0, Nnan)
	if Nnan gt 0 then sn_nline[nan] = 0.
	nan = where(finite(sn_bline) eq 0, Nnan)
	if Nnan gt 0 then sn_bline[nan] = 0.
	
	;If lineset is 1 consider using NII for S/N if higher than Halpha. Useful for kinematics.
	if state.selected_lineset eq 1 then begin
	    maxsnn = max(sn_nline, dimension=3)
	    maxsnb = max(sn_bline, dimension=3) 
	endif else begin
	    main_ind = where(abs(state.lines - kubeviz_getmainline()) lt 1)
	    maxsnn = sn_nline[*,*,main_ind]
	    maxsnb = sn_bline[*,*,main_ind]
	endelse
	
	maxerrokn = (abs((*state.sp_nerrrescube)[*,*,1,0]) lt maxvelerr and abs((*state.sp_nerrrescube)[*,*,1,1]) lt maxvelerr) and (abs((*state.sp_nerrrescube)[*,*,2,0]) lt maxvelerr and abs((*state.sp_nerrrescube)[*,*,2,1]) lt maxvelerr)     
	maxerrokb = (abs((*state.sp_berrrescube)[*,*,1,0]) lt maxvelerr and abs((*state.sp_berrrescube)[*,*,1,1]) lt maxvelerr) and (abs((*state.sp_berrrescube)[*,*,2,0]) lt maxvelerr and abs((*state.sp_berrrescube)[*,*,2,1]) lt maxvelerr)     
	
	if state.fitconstr then begin ;if fit constraints is active the errors can be zero while the fit is still good
	 velvalokn =   (  ((*state.sp_nerrrescube)[*,*,1,0] ne -999) and ((*state.sp_nerrrescube)[*,*,1,1] ne -999) and ((*state.sp_nerrrescube)[*,*,2,0] ne -999) and ((*state.sp_nerrrescube)[*,*,2,1] ne -999) )
	 velvalokb =   (  ((*state.sp_berrrescube)[*,*,1,0] ne -999) and ((*state.sp_berrrescube)[*,*,1,1] ne -999) and ((*state.sp_berrrescube)[*,*,2,0] ne -999) and ((*state.sp_berrrescube)[*,*,2,1] ne -999) )  
    
        endif else begin
	 velvalokn =     ( (((*state.sp_nerrrescube)[*,*,1,0] ne 0 and (*state.sp_nerrrescube)[*,*,1,0] ne -999) and ((*state.sp_nerrrescube)[*,*,1,1] ne 0 and (*state.sp_nerrrescube)[*,*,1,1] ne -999))  $  
                     and   (((*state.sp_nerrrescube)[*,*,2,0] ne 0 and (*state.sp_nerrrescube)[*,*,2,0] ne -999) and ((*state.sp_nerrrescube)[*,*,2,1] ne 0 and (*state.sp_nerrrescube)[*,*,2,1] ne -999)) )    
         velvalokb =     ( (((*state.sp_berrrescube)[*,*,1,0] ne 0 and (*state.sp_berrrescube)[*,*,1,0] ne -999) and ((*state.sp_berrrescube)[*,*,1,1] ne 0 and (*state.sp_berrrescube)[*,*,1,1] ne -999))  $ 
                     and   (((*state.sp_berrrescube)[*,*,2,0] ne 0 and (*state.sp_berrrescube)[*,*,2,0] ne -999) and ((*state.sp_berrrescube)[*,*,2,1] ne 0 and (*state.sp_berrrescube)[*,*,2,1] ne -999)) )    

	endelse
	
   	;nvelclip = kubeviz_linefit_flagsigclip((reform((*state.sp_nrescube)[*,*,1])), 1)
	;nsigclip = kubeviz_linefit_flagsigclip((reform((*state.sp_nrescube)[*,*,2])), 2)
	;bvelclip = kubeviz_linefit_flagsigclip((reform((*state.sp_brescube)[*,*,1])), 1)
	;bsigclip = kubeviz_linefit_flagsigclip((reform((*state.sp_brescube)[*,*,2])), 2)
	
	;clipokn =  (((*state.sp_nrescube)[*,*,1] gt nvelclip[0]) and ((*state.sp_nrescube)[*,*,1] lt nvelclip[1]) and ((*state.sp_nrescube)[*,*,1] ne 0) and ((*state.sp_nrescube)[*,*,2] lt nsigclip[1]) and  ((*state.sp_nrescube)[*,*,2] gt 0))
	;clipokb =  (((*state.sp_brescube)[*,*,1] gt bvelclip[0]) and ((*state.sp_brescube)[*,*,1] lt bvelclip[1]) and ((*state.sp_brescube)[*,*,1] ne 0) and ((*state.sp_brescube)[*,*,2] lt bsigclip[1]) and  ((*state.sp_brescube)[*,*,2] gt 0))
	  
	badspaxels = (*state.badpixelimg)
        setmaskn = 16B*(badspaxels eq 1) + 1B*(maxsnn lt sn_thresh) + 2B*(maxerrokn eq 0) + 4B*(velvalokn eq 0) ;+ 8B*(clipokn eq 0)  
	setmaskb = 16B*(badspaxels eq 1) + 1B*(maxsnb lt sn_thresh) + 2B*(maxerrokb eq 0) + 4B*(velvalokb eq 0) ;+ 8B*(clipokb eq 0)  
	setmask  = min( [[[setmaskn]], [[setmaskb]]]  , dimension=3)
	
	newflag = doflag*(setmaskn) + (1-doflag)*(*state.sp_nrescube)[*,*,0]
        (*state.sp_nrescube)[*,*,0] = doflag*(setmaskn) + (1-doflag)*(*state.sp_nrescube)[*,*,0]
        (*state.sp_brescube)[*,*,0] = doflag*(setmaskb) + (1-doflag)*(*state.sp_brescube)[*,*,0]
        (*state.sp_crescube)[*,*,0] = doflag*(setmask)  + (1-doflag)*(*state.sp_crescube)[*,*,0]
        
	for i=0, state.Nmontecarlo_percs-1 do begin
            (*state.sp_nerrrescube)[*,*,0,i] = doflag*(setmaskn) + (1-doflag)*(*state.sp_nrescube)[*,*,0]
            (*state.sp_berrrescube)[*,*,0,i] = doflag*(setmaskb) + (1-doflag)*(*state.sp_brescube)[*,*,0]
            (*state.sp_cerrrescube)[*,*,0,i] = doflag*(setmask)  + (1-doflag)*(*state.sp_crescube)[*,*,0]   
        endfor

    end
    
    1: begin ; masks:
      
        nrescube = (*state.nrescube)
        brescube = (*state.brescube)
        crescube = (*state.crescube)
        nerrrescube = (*state.nerrrescube)
        berrrescube = (*state.berrrescube)
        cerrrescube = (*state.cerrrescube)
        
        ; mean one-sigma error is half of the
        ; distance from +1sig to -1sig which are stored in elements 0,1:
        nline_onesig = .5*(reform(nerrrescube[*,3+indgen(state.Nlines),0])-reform(nerrrescube[*,3+indgen(state.Nlines),1]))
        bline_onesig = .5*(reform(berrrescube[*,3+indgen(state.Nlines),0])-reform(berrrescube[*,3+indgen(state.Nlines),1]))
        sn_nline = reform(nrescube[*,3+indgen(state.Nlines)])/nline_onesig
        sn_bline = reform(brescube[*,3+indgen(state.Nlines)])/bline_onesig
        for iline=0, state.Nlines-1 do begin
             zeroerr_nline = where(nerrrescube[*,3+iline,0] eq 0. and nerrrescube[*,3+iline,1] eq 0., Nzeroerr_nline)
             if (Nzeroerr_nline gt 0.) then begin
                 sn_nline_line = sn_nline[*,iline]
                 sn_nline_line[zeroerr_nline] = 0.
                 sn_nline[*,iline] = sn_nline_line
             endif
             zeroerr_bline = where(berrrescube[*,3+iline,0] eq 0. and berrrescube[*,3+iline,1] eq 0., Nzeroerr_bline)
             if (Nzeroerr_bline gt 0.) then begin
                 sn_bline_line = sn_bline[*,iline]
                 sn_bline_line[zeroerr_bline] = 0.
                 sn_bline[*,iline] = sn_bline_line
             endif
         endfor
        maxsn = max( [ [max(sn_nline, dimension=2)], [max(sn_bline, dimension=2)] ], dimension=2 ) 
        velerrok = ( ( ((nerrrescube[*,1,0] ne 0 and nerrrescube[*,1,0] ne -999) and (nerrrescube[*,1,1] ne 0 and nerrrescube[*,1,1] ne -999))   $ 
                    or ((berrrescube[*,1,0] ne 0 and berrrescube[*,1,0] ne -999) and (berrrescube[*,1,1] ne 0 and berrrescube[*,1,1] ne -999)) ) $ 
                 and ( ((nerrrescube[*,2,0] ne 0 and nerrrescube[*,2,0] ne -999) and (nerrrescube[*,2,1] ne 0 and nerrrescube[*,2,1] ne -999))   $ 
                    or ((berrrescube[*,2,0] ne 0 and berrrescube[*,2,0] ne -999) and (berrrescube[*,2,1] ne 0 and berrrescube[*,2,1] ne -999)) )   $ 
               and ((abs(nerrrescube[*,1,0]) lt maxvelerr and abs(nerrrescube[*,1,1]) lt maxvelerr) or (abs(berrrescube[*,1,0]) lt maxvelerr and abs(berrrescube[*,1,1]) lt maxvelerr)) $ 
               and ((abs(nerrrescube[*,2,0]) lt maxvelerr and abs(nerrrescube[*,2,1]) lt maxvelerr) or (abs(berrrescube[*,2,0]) lt maxvelerr and abs(berrrescube[*,2,1]) lt maxvelerr)) )
        setok = (maxsn gt sn_thresh and velerrok eq 1)
        nrescube[*,0] = (1-setok)
        brescube[*,0] = (1-setok)
        crescube[*,0] = (1-setok)
        for i=0, state.Nmontecarlo_percs-1 do begin
            nerrrescube[*,0,i] = (1-setok)
            berrrescube[*,0,i] = (1-setok)
            cerrrescube[*,0,i] = (1-setok)
        endfor
        
        ;copy back into the state variables
        
        (*state.nrescube)        =   nrescube
        (*state.brescube)        =   brescube
        (*state.crescube)        =   crescube
        (*state.nerrrescube)  =   nerrrescube
        (*state.berrrescube)  =   berrrescube
        (*state.cerrrescube)  =   cerrrescube
                
    end
    endcase
  end
  
  1: begin ;moments
  case state.linefit_mode of
    0: begin ; spaxels
         
        mrescube = (*state.sp_mrescube)
        merrrescube = (*state.sp_merrrescube)
        
	if max((*state.sp_merrrescube)[*,*,0])  gt -998 then begin
	
        for iline=0, state.Nlines-1 do begin
        
         ; mean one-sigma error is half of the
         ; distance from +1sig to -1sig which are stored in elements 0,1:
         line_onesig = .5*(reform(merrrescube[*,*,6*iline,0])-reform(merrrescube[*,*,6*iline,1]))
         sn_line = reform(mrescube[*,*,6*iline])/line_onesig
         zeroerr_line = where(merrrescube[*,*,6*iline,0] eq 0. and merrrescube[*,*,6*iline,1] eq 0., Nzeroerr_line)
             if (Nzeroerr_line gt 0.) then sn_line[zeroerr_line] = 0.
         
         velerrok = ( ((merrrescube[*,*,6*iline+1,0] ne 0 and merrrescube[*,*,6*iline+1,0] ne -999) and (merrrescube[*,*,6*iline+1,1] ne 0 and merrrescube[*,*,6*iline+1,1] ne -999))   $ 
                 and  ((merrrescube[*,*,6*iline+2,0] ne 0 and merrrescube[*,*,6*iline+2,0] ne -999) and (merrrescube[*,*,6*iline+2,1] ne 0 and merrrescube[*,*,6*iline+2,1] ne -999))   $ 
                 and (abs(merrrescube[*,*,6*iline+1,0]) lt maxvelerr and abs(merrrescube[*,*,6*iline+1,1]) lt maxvelerr)  and  $
                     (abs(merrrescube[*,*,6*iline+2,0]) lt maxvelerr and abs(merrrescube[*,*,6*iline+2,1]) lt maxvelerr))
         
         badspaxels = (*state.badpixelimg)
         setok = (sn_line gt sn_thresh and velerrok eq 1 and badspaxels eq 0)
         mrescube[*,*,6*iline+5] = (1-setok)
         for i=0, state.Nmontecarlo_percs-1 do begin
             merrrescube[*,*,6*iline+5,i] = (1-setok)
         endfor
        endfor
        
        ;copy back into the state variables
                
        (*state.sp_mrescube)     =   mrescube
        (*state.sp_merrrescube)  =   merrrescube
	
	endif
          
    end
    1: begin ; masks:
      
        mrescube = (*state.mrescube)
        merrrescube = (*state.merrrescube)
        
	if max((*state.sp_merrrescube)[*,*,0])  gt -998 then begin
         
	 for iline=0, state.Nlines-1 do begin
        
         ; mean one-sigma error is half of the
         ; distance from +1sig to -1sig which are stored in elements 0,1:
         line_onesig = .5*(reform(merrrescube[*,6*iline,0])-reform(merrrescube[*,6*iline,1]))
         sn_line = reform(mrescube[*,6*iline])/line_onesig
         zeroerr_line = where(merrrescube[*,6*iline,0] eq 0. and merrrescube[*,6*iline,1] eq 0., Nzeroerr_line)
             if (Nzeroerr_line gt 0.) then sn_line[zeroerr_line] = 0.
       
         velerrok = ( ((merrrescube[*,6*iline+1,0] ne 0 and merrrescube[*,6*iline+1,0] ne -999) and (merrrescube[*,6*iline+1,1] ne 0 and merrrescube[*,6*iline+1,1] ne -999))   $ 
                 and  ((merrrescube[*,6*iline+2,0] ne 0 and merrrescube[*,6*iline+2,0] ne -999) and (merrrescube[*,6*iline+2,1] ne 0 and merrrescube[*,6*iline+2,1] ne -999))   $ 
                 and (abs(merrrescube[*,6*iline+1,0]) lt maxvelerr and abs(merrrescube[*,6*iline+1,1]) lt maxvelerr)  and  $
                     (abs(merrrescube[*,6*iline+2,0]) lt maxvelerr and abs(merrrescube[*,6*iline+2,1]) lt maxvelerr))
         
         setok = (sn_line gt sn_thresh and velerrok eq 1)
         mrescube[*,6*iline+5] = (1-setok)
         for i=0, state.Nmontecarlo_percs-1 do begin
             merrrescube[*,6*iline+5,i] = (1-setok)
         endfor
        endfor
        
        ;copy back into the state variables
        (*state.mrescube)        =   mrescube
        (*state.merrrescube)  =   merrrescube
        
	endif        
    end
  endcase       
end
endcase

kubeviz_linefit_pointerswitch, /save

end


;------------------------------------------------------------
pro kubeviz_linefit_adjflag, xfit, yfit, retval=retval

; flag as bad all spaxels or masks which have no line with flux S/N >=limit

; Flag = sum of binary flags:
; 16 x manual flag
;  8 x velocity and dispersion within 4 sigma of median values for
;      sigma=|16/84%ile - 50%ile|
;  4 x valid velocity and dispersion values
;  2 x errors on velocity and dispersion below threshold error
;  1 x S/N >= S/N_thresh

common kubeviz_state

sn_thresh = state.mask_sn_thresh
maxvelerr = state.mask_maxvelerr

kubeviz_linefit_pointerswitch, /load

case state.linefit_type of
  0: begin ; Gaussian fit
        
        
        ; mean one-sigma error is half of the
        ; DJW 12-10-2016
        ; NEW: distance from -1sig ONLY relevant to determine S/N

        sn_nline = reform((*state.sp_nrescube)[xfit,yfit,3+indgen(state.Nlines)])/abs(reform((*state.sp_nerrrescube)[xfit,yfit,3+indgen(state.Nlines),1]))
        sn_bline = reform((*state.sp_brescube)[xfit,yfit,3+indgen(state.Nlines)])/abs(reform((*state.sp_berrrescube)[xfit,yfit,3+indgen(state.Nlines),1]))
	
	nan = where(finite(sn_nline) eq 0, Nnan)
	if Nnan gt 0 then sn_nline[nan] = 0.
	nan = where(finite(sn_bline) eq 0, Nnan)
	if Nnan gt 0 then sn_bline[nan] = 0.

	main_ind = where(abs(state.lines - kubeviz_getmainline()) lt 1)
	maxsnn = sn_nline[main_ind]
	maxsnb = sn_bline[main_ind]
	
	maxerrokn = (abs((*state.sp_nerrrescube)[xfit,yfit,1,0]) lt maxvelerr and abs((*state.sp_nerrrescube)[xfit,yfit,1,1]) lt maxvelerr) and (abs((*state.sp_nerrrescube)[xfit,yfit,2,0]) lt maxvelerr and abs((*state.sp_nerrrescube)[xfit,yfit,2,1]) lt maxvelerr)     
	maxerrokb = (abs((*state.sp_berrrescube)[xfit,yfit,1,0]) lt maxvelerr and abs((*state.sp_berrrescube)[xfit,yfit,1,1]) lt maxvelerr) and (abs((*state.sp_berrrescube)[xfit,yfit,2,0]) lt maxvelerr and abs((*state.sp_berrrescube)[xfit,yfit,2,1]) lt maxvelerr)     
	
	if state.fitconstr then begin ;if fit constraints is active the errors can be zero while the fit is still good
	 velvalokn =   (  ((*state.sp_nerrrescube)[xfit,yfit,1,0] ne -999) and ((*state.sp_nerrrescube)[xfit,yfit,1,1] ne -999) and ((*state.sp_nerrrescube)[xfit,yfit,2,0] ne -999) and ((*state.sp_nerrrescube)[xfit,yfit,2,1] ne -999) )
	 velvalokb =   (  ((*state.sp_berrrescube)[xfit,yfit,1,0] ne -999) and ((*state.sp_berrrescube)[xfit,yfit,1,1] ne -999) and ((*state.sp_berrrescube)[xfit,yfit,2,0] ne -999) and ((*state.sp_berrrescube)[xfit,yfit,2,1] ne -999) )  
    
        endif else begin
	 velvalokn =     ( (((*state.sp_nerrrescube)[xfit,yfit,1,0] ne 0 and (*state.sp_nerrrescube)[xfit,yfit,1,0] ne -999) and ((*state.sp_nerrrescube)[xfit,yfit,1,1] ne 0 and (*state.sp_nerrrescube)[xfit,yfit,1,1] ne -999))  $  
                     and   (((*state.sp_nerrrescube)[xfit,yfit,2,0] ne 0 and (*state.sp_nerrrescube)[xfit,yfit,2,0] ne -999) and ((*state.sp_nerrrescube)[xfit,yfit,2,1] ne 0 and (*state.sp_nerrrescube)[xfit,yfit,2,1] ne -999)) )    
         velvalokb =     ( (((*state.sp_berrrescube)[xfit,yfit,1,0] ne 0 and (*state.sp_berrrescube)[xfit,yfit,1,0] ne -999) and ((*state.sp_berrrescube)[xfit,yfit,1,1] ne 0 and (*state.sp_berrrescube)[xfit,yfit,1,1] ne -999))  $ 
                     and   (((*state.sp_berrrescube)[xfit,yfit,2,0] ne 0 and (*state.sp_berrrescube)[xfit,yfit,2,0] ne -999) and ((*state.sp_berrrescube)[xfit,yfit,2,1] ne 0 and (*state.sp_berrrescube)[xfit,yfit,2,1] ne -999)) )    

	endelse
	
	nvelclip = kubeviz_linefit_flagsigclip((reform((*state.sp_nrescube)[*,*,1])), 1)
	nsigclip = kubeviz_linefit_flagsigclip((reform((*state.sp_nrescube)[*,*,2])), 2)
	bvelclip = kubeviz_linefit_flagsigclip((reform((*state.sp_brescube)[*,*,1])), 1)
	bsigclip = kubeviz_linefit_flagsigclip((reform((*state.sp_brescube)[*,*,2])), 2)
	
	clipokn =  (((*state.sp_nrescube)[xfit,yfit,1] gt nvelclip[0]) and ((*state.sp_nrescube)[xfit,yfit,1] lt nvelclip[1]) and ((*state.sp_nrescube)[xfit,yfit,1] ne 0) and ((*state.sp_nrescube)[xfit,yfit,2] lt nsigclip[1]) and  ((*state.sp_nrescube)[xfit,yfit,2] gt 0))
	clipokb =  (((*state.sp_brescube)[xfit,yfit,1] gt nvelclip[0]) and ((*state.sp_brescube)[xfit,yfit,1] lt nvelclip[1]) and ((*state.sp_brescube)[xfit,yfit,1] ne 0) and ((*state.sp_brescube)[xfit,yfit,2] lt nsigclip[1]) and  ((*state.sp_brescube)[xfit,yfit,2] gt 0))
	  
        setmaskn = 1B*(maxsnn lt sn_thresh) + 2B*(maxerrokn eq 0) + 4B*(velvalokn eq 0) + 8B*(clipokn eq 0)
	setmaskb = 1B*(maxsnb lt sn_thresh) + 2B*(maxerrokb eq 0) + 4B*(velvalokb eq 0) + 8B*(clipokb eq 0)
	setmask  = min( [[[setmaskn]], [[setmaskb]]]  , dimension=3)
	
        (*state.sp_nrescube)[xfit,yfit,0] = setmaskn 
        (*state.sp_brescube)[xfit,yfit,0] = setmaskb 
        (*state.sp_crescube)[xfit,yfit,0] = setmask  
        
	for i=0, state.Nmontecarlo_percs-1 do begin
            (*state.sp_nerrrescube)[xfit,yfit,0,i] = setmaskn
            (*state.sp_berrrescube)[xfit,yfit,0,i] = setmaskb
            (*state.sp_cerrrescube)[xfit,yfit,0,i] = setmask 
        endfor

	
  end
    
  
  1: begin ;moments
         
        mrescube = (*state.sp_mrescube)
        merrrescube = (*state.sp_merrrescube)
        
	if max((*state.sp_merrrescube)[*,*,0])  gt -998 then begin
	
        for iline=0, state.Nlines-1 do begin
        
         ; mean one-sigma error is half of the
         ; distance from +1sig to -1sig which are stored in elements 0,1:
         line_onesig = .5*(reform(merrrescube[*,*,6*iline,0])-reform(merrrescube[*,*,6*iline,1]))
         sn_line = reform(mrescube[*,*,6*iline])/line_onesig
         zeroerr_line = where(merrrescube[*,*,6*iline,0] eq 0. and merrrescube[*,*,6*iline,1] eq 0., Nzeroerr_line)
             if (Nzeroerr_line gt 0.) then sn_line[zeroerr_line] = 0.
         
         velerrok = ( ((merrrescube[*,*,6*iline+1,0] ne 0 and merrrescube[*,*,6*iline+1,0] ne -999) and (merrrescube[*,*,6*iline+1,1] ne 0 and merrrescube[*,*,6*iline+1,1] ne -999))   $ 
                 and  ((merrrescube[*,*,6*iline+2,0] ne 0 and merrrescube[*,*,6*iline+2,0] ne -999) and (merrrescube[*,*,6*iline+2,1] ne 0 and merrrescube[*,*,6*iline+2,1] ne -999))   $ 
                 and (abs(merrrescube[*,*,6*iline+1,0]) lt maxvelerr and abs(merrrescube[*,*,6*iline+1,1]) lt maxvelerr)  and  $
                     (abs(merrrescube[*,*,6*iline+2,0]) lt maxvelerr and abs(merrrescube[*,*,6*iline+2,1]) lt maxvelerr))
         
         badspaxels = (*state.badpixelimg)
         setok = (sn_line gt sn_thresh and velerrok eq 1 and badspaxels eq 0)
         mrescube[*,*,6*iline+5] = (1-setok)
         for i=0, state.Nmontecarlo_percs-1 do begin
             merrrescube[*,*,6*iline+5,i] = (1-setok)
         endfor
        endfor
        
        ;copy back into the state variables
                
        (*state.sp_mrescube)     =   mrescube
        (*state.sp_merrrescube)  =   merrrescube
	
	endif
          
    end
endcase

kubeviz_linefit_pointerswitch, /save
retval = setmask

end



;------------------------------------------------------------------
pro kubeviz_linefit_fitadjall  
; fit spaxels with flag=BAD where they are adjacent to other spaxels
; with flag=OK, using the solution from the adjacent spaxels for an
; initial guess. Iterate until no more successes.
common kubeviz_state

kubeviz_linefit_pointerswitch, /load

colstart = state.col
rowstart = state.row

nrescube    = (*state.sp_nrescube)
brescube    = (*state.sp_brescube)
crescube    = (*state.sp_crescube)
nerrrescube = (*state.sp_nerrrescube)
berrrescube = (*state.sp_berrrescube)
cerrrescube = (*state.sp_cerrrescube)
badspaxels  = (*state.badpixelimg)

flag = min( [[[nrescube[*,*,0]]], [[brescube[*,*,0]]]]  , dimension=3)

;First reset fit of manually defined bad spaxels (flag=16)
manualbad_n = where(nrescube[*,*,0] eq 16, Nmanualbad)
if Nmanualbad gt 0 then begin
  indices = array_indices(flag, manualbad_n)
  for ii=0, Nmanualbad-1 do begin
    state.col = indices[0,ii]
    state.row = indices[1,ii]
    kubeviz_linefit_reset, linetype= 'N'
  endfor
endif
manualbad_b = where(brescube[*,*,0] eq 16, Nmanualbad)
if Nmanualbad gt 0 then begin
  indices = array_indices(flag, manualbad_b)
  for ii=0, Nmanualbad-1 do begin
    state.col = indices[0,ii]
    state.row = indices[1,ii]
    kubeviz_linefit_reset, linetype= 'B'
  endfor
endif

badfit = (badspaxels eq 0 and flag gt 0)
Nbad_init = total(badfit, /integer)

Hope     = replicate(1B,state.Ncol,state.Nrow)
Nokneigh = replicate(0, state.Ncol,state.Nrow)

x = indgen(state.Ncol)#replicate(1,state.Nrow)
y = replicate(1,state.Ncol)#indgen(state.Nrow)

if state.Nloopadj gt 0 then Niter = state.Nloopadj else Niter = 100000L

outrange = where((x lt state.fitallrange[0]) or (x gt state.fitallrange[1]) or (y lt state.fitallrange[2]) or (y gt state.fitallrange[3]), Nout)
if Nout gt 0 then begin
  Hope[outrange] = 0
  badfit[outrange] = 0
endif

;Display cancel button
kubeviz_linefit_cancel, /init

;Start loop
while total(badfit*Hope) gt 1 and Niter ge 0 do begin
   
   Niter = Niter-1  ; Reduce by 1 the number of remaining iterations
   
   isbadfit = where((badfit eq 1) and (x ge state.fitallrange[0]) and (x le state.fitallrange[1]) and (y ge state.fitallrange[2]) and (y le state.fitallrange[3]), Nbadfit)
   for ii=0L,Nbadfit-1 do Nokneigh[x[isbadfit[ii]],y[isbadfit[ii]]] = kubeviz_linefit_Nokadjacentspax(x[isbadfit[ii]], y[isbadfit[ii]], flag) 
   
   ;No hope for those with no good neighbors
   no_okneigh = where(Nokneigh eq 0, Nno_okneigh)
   if Nno_okneigh gt 0 then Hope[no_okneigh] = 0
   	 
   for ii=8,1,-1 do begin
      isbadfitnow = where((Hope eq 1) and (Nokneigh eq ii) and (badfit eq 1) and (x ge state.fitallrange[0]) and (x le state.fitallrange[1]) and (y ge state.fitallrange[2]) and (y le state.fitallrange[3]), Nfitnow)
      if Nfitnow gt 0 then break
   endfor
   
   ;HERE TO BE FIXED
   kubeviz_statuslinecore,  /clear
   kubeviz_statuslinecore, '[PROGRES] Running step with '+kubeviz_str(ii)+' neighbours, for '+kubeviz_str(Nfitnow)+' pixels.',0 
   
   for i=0L,Nfitnow-1 do begin
      
      x_fit = x[isbadfitnow[i]]
      y_fit = y[isbadfitnow[i]]
      
      state.col = x_fit
      state.row = y_fit
      
      ; look at results from adjacent spaxels with ok fit (if there are any)
      if kubeviz_linefit_guessfromadjacentspax(x_fit, y_fit, flag, nrescube, brescube, nerrrescube, berrrescube) then begin
	
	; using as initial guess, try to fit:			      
	kubeviz_linefit_dofit
       
	;If the cube is very big speed up by reducing GUI live updates
	if state.Ncol lt 200 and state.Nrow lt 200 then begin 
	  kubeviz_plotspax
	  kubeviz_linefit_update, /update_userpars
	endif else begin
	  kubeviz_plotspax, /fast_update 
	endelse
	
	kubeviz_plotspeczoom

	; assess individual spaxel using autoflag criteria, 
	; and put in local flag variable:
	kubeviz_linefit_adjflag, x_fit, y_fit, retval=newflag
	flag = min( [[[(*state.sp_nrescube)[*,*,0]]], [[(*state.sp_brescube)[*,*,0]]]]  , dimension=3)
	
	if newflag eq 0 then begin
     	   Hope[max([x_fit-1,0]):min([x_fit+1,state.Ncol-1]),max([y_fit-1,0]):min([y_fit+1,state.Nrow-1])] = 1
	endif else begin
     	   Hope[x_fit, y_fit] = 0
     	endelse
    	 
      endif
      
      kubeviz_linefit_cancel, /cancel, event = event
      if event then GOTO, JUMPADJ
   endfor
   
   kubeviz_plotspax
   ; re-evaluate bad fits:
   
   badfit = (badspaxels eq 0 and flag gt 0)
  
endwhile

JUMPADJ:
kubeviz_linefit_cancel, /destroy
kubeviz_linefit_resetuser, /ALL, /startonly

state.col  = colstart 
state.row  = rowstart 

Nbad_final = total((badspaxels eq 0 and flag gt 0), /integer)
kubeviz_statuslinecore, /close
printf, state.log_lun, ''
printf, state.log_lun, '[KUBEVIZ] The fit was improved for '+strcompress(Nbad_init-Nbad_final, /re)+' spaxels.'
  
end


;------------------------------------------------------------------
pro kubeviz_linefit_fitadj  
; fit single spaxel with using the solution from the adjacent spaxels for an initial guess.
common kubeviz_state

kubeviz_linefit_pointerswitch, /load

nrescube = (*state.sp_nrescube)
brescube = (*state.sp_brescube)
crescube = (*state.sp_crescube)
nerrrescube = (*state.sp_nerrrescube)
berrrescube = (*state.sp_berrrescube)
cerrrescube = (*state.sp_cerrrescube)

;save user defined initial guesses (if any)
nuser = state.gnfit
buser = state.gbfit

; look at results from adjacent spaxels with ok fit (if there are any)
if kubeviz_linefit_guessfromadjacentspax(state.col, state.row, crescube[*,*,0], nrescube, brescube, nerrrescube, berrrescube) then begin

 ; using as initial guess, try to fit:			      
 kubeviz_linefit_dofit

 state.gnfit = nuser
 state.gbfit = buser
           
 ; assess individual spaxel using autoflag criteria, 
 kubeviz_linefit_adjflag, state.col, state.row
 
endif

end
;----------------------------------------------------
pro kubeviz_linefit_saveres, fname=fname
; save the results to a FITS data cube
; add header keywords describing the planes in the cube
common kubeviz_state

printf, state.log_lun, '[KUBEVIZ] Saving FITS results file...'

; choose filename for flagged results cube output
if n_elements(fname) eq 0 then begin
    fname = dialog_pickfile(filter='*.fits', path=state.outdir, /write, get_path=path, /OVERWRITE_PROMPT)
    len = strlen(fname) - strlen(path)
    if len eq 0 then return     ; don't do anything if no filename is given
endif

; unflagged version:
rawfname = strmid(fname,0,strlen(fname)-5)+'_noflag.fits'

; If the current method is gaussian evaluate the object redshift
; TEMPORARY SOLUTION (dummy variable):
redshift_out = -1
;if state.linefit_type eq 0 and state.redshift gt 0.1 then redshift_out = kubeviz_linefit_getredshift() else redshift_out = -1

; combine results cubes:
kubeviz_linefit_pointerswitch, /load

case state.linefit_type of 
  0: begin ; gaussian fit
  case state.linefit_mode of
    0: begin ; spaxels
        rescube = replicate(0.D, state.Ncol, state.Nrow , 6*(3+state.Nlines) + 3*(1+state.Nlines) + 1 )
        
	rescube[*,*,indgen(3+state.Nlines)]		       = (*state.sp_nrescube)
        rescube[*,*,(3+state.Nlines)+indgen(3+state.Nlines)]   = (*state.sp_brescube)
        rescube[*,*,2*(3+state.Nlines)+indgen(1+state.Nlines)] = (*state.sp_crescube)
        rescube[*,*,2*(3+state.Nlines)+(1+state.Nlines)+indgen(3+state.Nlines)]   = (*state.sp_nerrrescube)[*,*,*,0] ; plus 1-sig
        rescube[*,*,3*(3+state.Nlines)+(1+state.Nlines)+indgen(3+state.Nlines)]   = (*state.sp_nerrrescube)[*,*,*,1] ; minus 1-sig
        rescube[*,*,4*(3+state.Nlines)+(1+state.Nlines)+indgen(3+state.Nlines)]   = (*state.sp_berrrescube)[*,*,*,0] ; plus 1-sig
        rescube[*,*,5*(3+state.Nlines)+(1+state.Nlines)+indgen(3+state.Nlines)]   = (*state.sp_berrrescube)[*,*,*,1] ; minus 1-sig
        rescube[*,*,6*(3+state.Nlines)+(1+state.Nlines)+indgen(1+state.Nlines)]   = (*state.sp_cerrrescube)[*,*,*,0] ; plus 1-sig
        rescube[*,*,6*(3+state.Nlines)+2*(1+state.Nlines)+indgen(1+state.Nlines)] = (*state.sp_cerrrescube)[*,*,*,1] ; minus 1-sig
        rescube[*,*,6*(3+state.Nlines)+3*(1+state.Nlines)]                        = (*state.sp_chisq) ; minus 1-sig
	
	rawrescube = rescube
	
	;Extract masks for narrow, broad and continuum
        nflag = (*state.sp_nrescube)[*,*,0] + 1
        badflag = where(nflag gt 1, Nbad)
	if Nbad gt 0 then nflag[badflag] *= !values.f_nan
	bflag = (*state.sp_brescube)[*,*,0] + 1
	badflag = where(bflag gt 1, Nbad)
	if Nbad gt 0 then bflag[badflag] *= !values.f_nan
        cflag = (*state.sp_crescube)[*,*,0] + 1
	badflag = where(cflag gt 1, Nbad)
	if Nbad gt 0 then cflag[badflag] *= !values.f_nan

	for i = 1+min(indgen(3+state.Nlines)), max(indgen(3+state.Nlines)) do  rescube[*,*,i] *= nflag
        for i = 1+min((3+state.Nlines)+indgen(3+state.Nlines)), max((3+state.Nlines)+indgen(3+state.Nlines))     do rescube[*,*,i] *= bflag
        for i = 1+min(2*(3+state.Nlines)+indgen(1+state.Nlines)),max(2*(3+state.Nlines)+indgen(1+state.Nlines))  do rescube[*,*,i] *= cflag
        for i = 1+min(2*(3+state.Nlines)+(1+state.Nlines)+indgen(2*(3+state.Nlines))), max(2*(3+state.Nlines)+(1+state.Nlines)+indgen(2*(3+state.Nlines))) do rescube[*,*,i] *= nflag 
        for i = 1+min(4*(3+state.Nlines)+(1+state.Nlines)+indgen(2*(3+state.Nlines))), max(4*(3+state.Nlines)+(1+state.Nlines)+indgen(2*(3+state.Nlines))) do rescube[*,*,i] *= bflag 
        for i = 1+min(6*(3+state.Nlines)+(1+state.Nlines)+indgen(2*(1+state.Nlines))), max(6*(3+state.Nlines)+(1+state.Nlines)+indgen(2*(1+state.Nlines))) do rescube[*,*,i] *= cflag 

    end
    1: begin ; masks
        rescube = replicate(0.D, state.maxNmask , 6*(3+state.Nlines) + 3*(1+state.Nlines) )
        flag = (*state.nrescube)[*,0] + 1
	badflag = where(flag gt 1, Nbad)
	if Nbad gt 0 then flag[badflag] *= !values.f_nan
        
	rescube[*,indgen(3+state.Nlines)]                    = (*state.nrescube)
        rescube[*,(3+state.Nlines)+indgen(3+state.Nlines)]   = (*state.brescube)
        rescube[*,2*(3+state.Nlines)+indgen(1+state.Nlines)] = (*state.crescube)
        rescube[*,2*(3+state.Nlines)+(1+state.Nlines)+indgen(3+state.Nlines)]   = (*state.nerrrescube)[*,*,0] ; plus 1-sig
        rescube[*,3*(3+state.Nlines)+(1+state.Nlines)+indgen(3+state.Nlines)]   = (*state.nerrrescube)[*,*,1] ; minus 1-sig
        rescube[*,4*(3+state.Nlines)+(1+state.Nlines)+indgen(3+state.Nlines)]   = (*state.berrrescube)[*,*,0] ; plus 1-sig
        rescube[*,5*(3+state.Nlines)+(1+state.Nlines)+indgen(3+state.Nlines)]   = (*state.berrrescube)[*,*,1] ; minus 1-sig
        rescube[*,6*(3+state.Nlines)+(1+state.Nlines)+indgen(1+state.Nlines)]   = (*state.cerrrescube)[*,*,0] ; plus 1-sig
        rescube[*,6*(3+state.Nlines)+2*(1+state.Nlines)+indgen(1+state.Nlines)] = (*state.cerrrescube)[*,*,1] ; minus 1-sig
        
	rawrescube = rescube
	
	for i=1, (6*(3+state.Nlines) + 3*(1+state.Nlines) - 1) do rescube[*,i] *= flag
    end
  endcase
  end
  
  1: begin ;moments
  
  case state.linefit_mode of
    0: begin ; spaxels
        ;All the moments are saved. The flag is applied individually for each line
        rescube    = replicate(0.D, state.Ncol, state.Nrow , (7*state.Nlines) + 6*(2*state.Nlines) )
        rawrescube = replicate(0.D, state.Ncol, state.Nrow , (7*state.Nlines) + 6*(2*state.Nlines) )
        for iline=0,state.Nlines-1 do begin
           rescube[*,*,iline*6]              = (*state.sp_mrescube)[*,*,iline*6+5]   ;Flag goes to the first plane
           rescube[*,*,iline*6+1+indgen(5)]  = (*state.sp_mrescube)[*,*,iline*6+indgen(5)]
           rescube[*,*,6*state.Nlines+iline] = (*state.sp_crescube)[*,*,iline+1]
           rescube[*,*,7*state.Nlines+10*iline+indgen(5)]    = (*state.sp_merrrescube)[*,*,iline*6+indgen(5),0]   ; plus 1-sig
           rescube[*,*,7*state.Nlines+10*iline+5+indgen(5)]  = (*state.sp_merrrescube)[*,*,iline*6+indgen(5),1]   ; minus 1-sig
           rescube[*,*,7*state.Nlines+10*state.Nlines+iline] = (*state.sp_cerrrescube)[*,*,iline+1,0]             ; plus 1-sig
           rescube[*,*,7*state.Nlines+11*state.Nlines+iline] = (*state.sp_cerrrescube)[*,*,iline+1,1]             ; minus 1-sig
        endfor
        rawrescube = rescube
        for iline=0,state.Nlines-1 do begin
           flag = (*state.sp_mrescube)[*,*,iline*6+5]
           flag_bad_col = [-1]
           flag_bad_row = [-1]
           flag_index = 0
           for icol=0,state.Ncol-1 do begin
               for irow=0,state.Nrow-1 do begin
                   if flag[icol,irow] gt 0 then begin
                       if flag_bad_col[0] eq -1 then begin
                           flag_bad_col[0] = icol
                           flag_bad_row[0] = irow
                       endif else begin
                           flag_index++
                           flag_bad_col = [flag_bad_col,icol]
                           flag_bad_row = [flag_bad_row,irow]
                       endelse
                   endif
               endfor
           endfor
           if flag_bad_col[0] ne -1 then begin
            for i=0,flag_index do begin  ;NaNs are not applied to the flag plane
             rescube[flag_bad_col[i],flag_bad_row[i],iline*6+1+indgen(5)] = !VALUES.F_NAN
             rescube[flag_bad_col[i],flag_bad_row[i],6*state.Nlines-1+iline] = !VALUES.F_NAN
             rescube[flag_bad_col[i],flag_bad_row[i],7*state.Nlines-1+10*iline+indgen(5)] = !VALUES.F_NAN
             rescube[flag_bad_col[i],flag_bad_row[i],7*state.Nlines-1+10*iline+6+indgen(5)] = !VALUES.F_NAN
             rescube[flag_bad_col[i],flag_bad_row[i],7*state.Nlines+10*state.Nlines-1+iline] = !VALUES.F_NAN
             rescube[flag_bad_col[i],flag_bad_row[i],7*state.Nlines+11*state.Nlines-1+iline] = !VALUES.F_NAN
            endfor
           endif
        endfor 
    end
    
    1: begin ; masks
        mrescube = (*state.mrescube)
        merrrescube = (*state.merrrescube)
        crescube = (*state.crescube)
        cerrrescube = (*state.cerrrescube)
        rescube    = replicate(0.D, state.maxNmask , (7*state.Nlines) + 6*(2*state.Nlines) )
        rawrescube = replicate(0.D, state.maxNmask , (7*state.Nlines) + 6*(2*state.Nlines) )
        for iline=0,state.Nlines-1 do begin
           rescube[*,iline*6] = mrescube[*,iline*6+5]
           rescube[*,iline*6+1+indgen(5)] = mrescube[*,iline*6+indgen(5)]
           rescube[*,6*state.Nlines+iline] = crescube[*,iline+1]
           rescube[*,7*state.Nlines+10*iline+indgen(5)] = merrrescube[*,iline*6+indgen(5),0]   ; plus 1-sig
           rescube[*,7*state.Nlines+10*iline+5+indgen(5)] = merrrescube[*,iline*6+indgen(5),1] ; minus 1-sig
           rescube[*,7*state.Nlines+10*state.Nlines+iline] = cerrrescube[*,iline+1,0]          ; plus 1-sig
           rescube[*,7*state.Nlines+11*state.Nlines+iline] = cerrrescube[*,iline+1,1]          ; minus 1-sig
        endfor
        rawrescube = rescube
        for iline=0,state.Nlines-1 do begin
           flag = mrescube[*,iline*6+5]
           flag_bad_mask = [-1]
           flag_index = 0
           for imask=0L,state.maxNmask-1 do begin
            if flag[imask] eq 0 then begin
                if flag_bad_mask[0] eq -1 then begin
                    flag_bad_mask[0] = imask
                endif else begin
                    flag_index++
                    flag_bad_mask = [flag_bad_mask,imask]
                endelse
            endif
           endfor
           if flag_bad_mask[0] ne -1 then begin
           for i=0,flag_index do begin ;NaNs are not applied to the flag line
            rescube[flag_bad_mask[i],iline*6+1+indgen(5)] = !VALUES.F_NAN
            rescube[flag_bad_mask[i],6*state.Nlines-1+iline] = !VALUES.F_NAN
            rescube[flag_bad_mask[i],7*state.Nlines-1+10*iline+indgen(5)] = !VALUES.F_NAN
            rescube[flag_bad_mask[i],7*state.Nlines-1+10*iline+6+indgen(5)] = !VALUES.F_NAN
            rescube[flag_bad_mask[i],7*state.Nlines+10*state.Nlines-1+iline] = !VALUES.F_NAN
            rescube[flag_bad_mask[i],7*state.Nlines+11*state.Nlines-1+iline] = !VALUES.F_NAN
           endfor
           endif
        endfor 
    end
  endcase ;Mode
  end
endcase ;Type
   
; setup header, update wcs refpixels to deal with trimmed images
xyad, (*state.indatahead), state.Startcol, state.Startrow, RA, DEC

mkhdr, hdr, rescube
sxdelpar, hdr, 'COMMENT'
sxaddpar, hdr, 'DATACUBE', state.filename           	    , 'Name of the original datacube'
sxaddpar, hdr, 'INSTRUME', strupcase(state.instr)           , 'Instrument name'
sxaddpar, hdr, 'BAND'    , strupcase(state.band)       	    , 'Band/Filter name'
sxaddpar, hdr, 'Z_IN'    , state.redshift           	    , 'Redshift used to fit the lines', format='F8.5'
sxaddpar, hdr, 'Z_OUT'   , redshift_out           	    , 'Redshift obtained from the fit', format='F8.5'
sxaddpar, hdr, 'FITTYPE' , fix(state.linefit_type)  	    , 'Type of fit: 0: gauss 1: moments'
sxaddpar, hdr, 'RESTYPE' , fix(state.linefit_mode)  	    , 'Type of result: 0: spaxels 1: masks'
sxaddpar, hdr, 'ERMETHOD', state.domontecarlo       	    , 'Err method: 0:Noise, 1:Bootstrap 2,3,4:MC1,2,3'
sxaddpar, hdr, 'SCNOISE' , fix(state.scalenoiseerrors)      , 'Scale noise Err. 0:Off 1:On (Only if ERMETH=0)'
sxaddpar, hdr, 'MCNOISE' , fix(state.useMonteCarlonoise)    , 'MonteCarlo noise 0:Off 1:On (Only if ERMETH>0)'
sxaddpar, hdr, 'SPATSMTH', state.smooth             	    , 'Spatial smoothing applied'
sxaddpar, hdr, 'SPECSMTH', state.specsmooth         	    , 'Spectral smoothing applied'
sxaddpar, hdr, 'LINESET' , state.selected_lineset   	    , 'Lineset used'

;The following makes sense only if all spaxels have been fitted (not for masks)
if state.linefit_mode eq 0 then begin 
  sxaddpar, hdr, 'CTYPE1'  , sxpar(*state.indatahead,'CTYPE1'), 'TAN projection used'
  sxaddpar, hdr, 'CTYPE2'  , sxpar(*state.indatahead,'CTYPE2'), 'TAN projection used'
  sxaddpar, hdr, 'CRPIX1'  , 1  			      , '[pix] Reference pixel in x'	     
  sxaddpar, hdr, 'CRPIX2'  , 1  			      , '[pix] Reference pixel in y'	     
  sxaddpar, hdr, 'CRVAL1'  , RA 			      , '[deg] RA at ref. pixel'
  sxaddpar, hdr, 'CRVAL2'  , DEC			      , '[deg] DEC at ref. pixel'
  sxaddpar, hdr, 'CD1_1'   , sxpar(*state.indatahead,'CD1_1') , '[] x-component of East'
  sxaddpar, hdr, 'CD2_1'   , sxpar(*state.indatahead,'CD2_1') , '[] x-component of North'
  sxaddpar, hdr, 'CD1_2'   , sxpar(*state.indatahead,'CD1_2') , '[] y-component of East'
  sxaddpar, hdr, 'CD2_2'   , sxpar(*state.indatahead,'CD2_2') , '[] y-component of North'
  sxaddpar, hdr, 'CDELT1'  , sxpar(*state.indatahead,'CDELT1'), '[deg] Pixel resolution in x'
  sxaddpar, hdr, 'CDELT2'  , sxpar(*state.indatahead,'CDELT2'), '[deg] Pixel resolution in y'
  sxaddpar, hdr, 'CUNIT1'  , sxpar(*state.indatahead,'CUNIT1'), 'Unit of x-axis'
  sxaddpar, hdr, 'CUNIT2'  , sxpar(*state.indatahead,'CUNIT2'), 'Unit of y-axis'
  sxaddpar, hdr, 'STARTCOL', state.Startcol+1		      , 'Index of the first col in the original frame'
  sxaddpar, hdr, 'STARTROW', state.Startrow+1		      , 'Index of the first row in the original frame'
endif

case state.linefit_type of
  0: begin ;Gaussian fit
      par = 0
      for restype=0, 8 do begin
          if restype ne 2 and restype lt 7 then Nrestype = (3+state.Nlines) else Nrestype = (1+state.Nlines)
          case restype of
              0: begin
                  linetype= 'narrow line'
                  datatype= ''
              end
              1:begin
                  linetype= 'broad line'
                  datatype= ''
              end
              2: begin
                  linetype= 'continuum'
                  datatype= ''
              end
              3: begin
                  linetype= 'narrow line'
                  datatype= 'plus 1-sig'
              end
              4: begin
                  linetype= 'narrow line'
                  datatype= 'minus 1-sig'
              end
              5: begin
                  linetype= 'broad line'
                  datatype= 'plus 1-sig'
              end
              6: begin
                  linetype= 'broad line'
                  datatype= 'minus 1-sig'
              end
              7: begin
                  linetype= 'continuum'
                  datatype= 'plus 1-sig'
              end
              8: begin
                    linetype= 'continuum'
                    datatype= 'minus 1-sig'          
              end
          endcase
          for i=0, Nrestype-1 do begin
              par++
              case i of
                  0: strpar='flag'
                  1: if (restype eq 2 or restype ge 7) then strpar=(*state.linenames)[i-1]+' flux' else strpar='dv'
                  2: if (restype eq 2 or restype ge 7) then strpar=(*state.linenames)[i-1]+' flux' else strpar='sig(v)'
                  else: if (restype eq 2 or restype ge 7) then strpar=(*state.linenames)[i-1]+' flux' else strpar=(*state.linenames)[i-3]+' flux'
              endcase
              sxaddpar, hdr, 'P'+kubeviz_str(par,format='(I3)'), linetype+' '+datatype+' '+strpar  
          endfor
      endfor
      sxaddpar, hdr, 'P'+kubeviz_str(par,format='(I3)'), 'Reduced chi square, (chisq/dof)' 
  end
  1: begin ;Moments
  par = 0
      for restype=0, 5 do begin
          case restype of
              0: begin
                  linetype= 'moment'
                  datatype= ''
                  startpar=0 & endpar=5
              end
              1: begin
                  linetype= 'continuum'
                  datatype= ''
                  startpar=0 & endpar=0
              end
              2: begin
                  linetype= 'moment'
                  datatype= 'plus 1-sig'
                  startpar=1 & endpar=5
              end
              3: begin
                  linetype= 'moment'
                  datatype= 'minus 1-sig'
                  startpar=1 & endpar=5
              end
              4: begin
                  linetype= 'continuum'
                  datatype= 'plus 1-sig'
                  startpar=0 & endpar=0
              end
              5: begin
                  linetype= 'continuum'
                  datatype= 'minus 1-sig'
                  startpar=0 & endpar=0
              end
          endcase
          for iline=0, state.Nlines-1 do begin
             for i=startpar, endpar  do begin
               par++
               case i of
                  0: if restype ne 1 and restype lt 4 then strpar=(*state.linenames)[iline]+' flag' else strpar = (*state.linenames)[iline]+' flux'
                  1: strpar=(*state.linenames)[iline]+' 0th (flux)'
                  2: strpar=(*state.linenames)[iline]+' 1st (dv)'
                  3: strpar=(*state.linenames)[iline]+' 2nd (sigv)'
                  4: strpar=(*state.linenames)[iline]+' 3rd (norm Skewness)'
                  5: strpar=(*state.linenames)[iline]+' 4th (norm Kurtosis)'
               endcase
               sxaddpar, hdr, 'P'+kubeviz_str(par,format='(I3)'), linetype+' '+datatype+' '+strpar 
             endfor         
          endfor
      endfor
  end
endcase    

writefits, fname, rescube, hdr
writefits, rawfname, rawrescube, hdr

end


;-------------------------------------------------------------
pro kubeviz_linefit_loadres, fname, dir=dir
; load the results from a FITS data cube
common kubeviz_state

if n_elements(dir) eq 0 then dir=''

;Read data and head
if file_test(dir+fname) eq 0 then begin
  print, '[ERROR] Results file does not exist!'
  return
endif else rescube = readfits(dir+fname, reshead, /silent)

;Perform consistency tests before loading the data into result cubes
if fxpar(reshead, 'LINESET') ne state.selected_lineset then begin
  print, '[ERROR] The current lineset is not compatible with the saved results!'
  return
endif

if fxpar(reshead, 'ERMETHOD') lt 0 or fxpar(reshead, 'ERMETHOD') gt 4 then begin
  print, '[ERROR ] Unrecognized error method!'
  return
endif

if fxpar(reshead, 'SPATSMTH') ne state.smooth or fxpar(reshead, 'SPECSMTH') ne state.specsmooth then begin
  printf, state.log_lun, '[WARNING] The current spatial or spectral smoothing do not match those for the saved results '
  printf, state.log_lun, '[WARNING] Saved results are loaded anyway. '
endif

if fxpar(reshead, 'Z_IN') - state.redshift gt 0.1 then begin
  printf, state.log_lun, '[WARNING] The current redshift significantly differs from the one stored in the results file. '
  printf, state.log_lun, '[WARNING] Saved results are loaded anyway. Use with caution! '
endif

reserror_method = fxpar(reshead, 'ERMETHOD')
reslinefit_type = fxpar(reshead, 'FITTYPE')
reslinefit_mode = fxpar(reshead, 'RESTYPE')
kubeviz_linefit_pointerswitch, /load, linefitmode=reslinefit_mode, errmethod=reserror_method


; unpack results cubes:
case reslinefit_mode of
    0: begin ; spaxels
        ;Additional check specific for the spaxel results
        if fxpar(reshead, 'NAXIS') ne 3 or fxpar(reshead, 'NAXIS1') ne state.Ncol or fxpar(reshead, 'NAXIS2') ne state.Nrow then begin
          print, '[ERROR ] Saved results cube does not match current cube dimensions!'
          return
        endif

        case reslinefit_type of 
           0: begin
             
               if fxpar(reshead, 'NAXIS3') ne 6*(3+state.Nlines) + 3*(1+state.Nlines) +1 then begin
                 print, '[ERROR ] The current number of lines does not match those in the saved results cube!'
                 return
               endif   
        
              (*state.sp_nrescube) = rescube[*,*,indgen(3+state.Nlines)] 
              (*state.sp_brescube) = rescube[*,*,(3+state.Nlines)+indgen(3+state.Nlines)] 
              (*state.sp_crescube) = rescube[*,*,2*(3+state.Nlines)+indgen(1+state.Nlines)]
              (*state.sp_nerrrescube)[*,*,*,0] = rescube[*,*,2*(3+state.Nlines)+(1+state.Nlines)+indgen(3+state.Nlines)] ; plus 1-sig
              (*state.sp_nerrrescube)[*,*,*,1] = rescube[*,*,3*(3+state.Nlines)+(1+state.Nlines)+indgen(3+state.Nlines)] ; minus 1-sig
              (*state.sp_berrrescube)[*,*,*,0] = rescube[*,*,4*(3+state.Nlines)+(1+state.Nlines)+indgen(3+state.Nlines)] ; plus 1-sig
              (*state.sp_berrrescube)[*,*,*,1] = rescube[*,*,5*(3+state.Nlines)+(1+state.Nlines)+indgen(3+state.Nlines)] ; minus 1-sig
              (*state.sp_cerrrescube)[*,*,*,0] = rescube[*,*,6*(3+state.Nlines)+(1+state.Nlines)+indgen(1+state.Nlines)] ; plus 1-sig
              (*state.sp_cerrrescube)[*,*,*,1] = rescube[*,*,6*(3+state.Nlines)+2*(1+state.Nlines)+indgen(1+state.Nlines)] ; minus 1-sig
        
              (*state.sp_nrescube) = kubeviz_remove_badvalues((*state.sp_nrescube))
              (*state.sp_brescube) = kubeviz_remove_badvalues((*state.sp_brescube))
              (*state.sp_crescube) = kubeviz_remove_badvalues((*state.sp_crescube))

              (*state.sp_nerrrescube) = kubeviz_remove_badvalues((*state.sp_nerrrescube),repval=-999)
              (*state.sp_berrrescube) = kubeviz_remove_badvalues((*state.sp_berrrescube),repval=-999)
              (*state.sp_cerrrescube) = kubeviz_remove_badvalues((*state.sp_cerrrescube),repval=-999)
          end
          1:begin
              
              if fxpar(reshead, 'NAXIS3') ne (7*state.Nlines) + 6*(2*state.Nlines) then begin
                 print, '[ERROR ] The current number of lines does not match those in the saved results cube!'
                 return
              endif  
              for iline=0,state.Nlines-1 do begin
                (*state.sp_mrescube)[*,*,iline*6+5]              = rescube[*,*,iline*6]                                   
                (*state.sp_mrescube)[*,*,iline*6+indgen(5)]      = rescube[*,*,iline*6+1+indgen(5)]                       
                (*state.sp_crescube)[*,*,iline+1]                = rescube[*,*,6*state.Nlines+iline]                      
                (*state.sp_merrrescube)[*,*,iline*6+indgen(5),0] = rescube[*,*,7*state.Nlines+10*iline+indgen(5)]                ; plus 1-sig
                (*state.sp_merrrescube)[*,*,iline*6+indgen(5),1] = rescube[*,*,7*state.Nlines+10*iline+5+indgen(5)]              ; minus 1-sig
                (*state.sp_cerrrescube)[*,*,iline+1,0]           = rescube[*,*,7*state.Nlines+10*state.Nlines+iline]             ; plus 1-sig
                (*state.sp_cerrrescube)[*,*,iline+1,1]           = rescube[*,*,7*state.Nlines+11*state.Nlines+iline]             ; minus 1-sig
              endfor
              
              (*state.sp_mrescube) = kubeviz_remove_badvalues((*state.sp_mrescube))
              (*state.sp_crescube) = kubeviz_remove_badvalues((*state.sp_crescube))
              
              (*state.sp_merrrescube) = kubeviz_remove_badvalues((*state.sp_merrrescube),repval=-999)
              (*state.sp_cerrrescube) = kubeviz_remove_badvalues((*state.sp_cerrrescube),repval=-999)
          end
        endcase      
        
    end
    1: begin ; masks
    
        ;Additional check specific for the mask results
        if fxpar(reshead, 'NAXIS') ne 2  then begin
          print, '[ERROR ] Saved results cube is corrupted!'
          return
        endif

        if fxpar(reshead, 'NAXIS2') ne 6*(3+state.Nlines) + 3*(1+state.Nlines) then begin
          print, '[ERROR ] The current number of lines does not match those in the saved results cube!'
          return
        endif   
        
        case reslinefit_type of 
           0: begin
              (*state.nrescube) = rescube[*,*,indgen(3+state.Nlines)] 
              (*state.brescube) = rescube[*,*,(3+state.Nlines)+indgen(3+state.Nlines)] 
              (*state.crescube) = rescube[*,*,2*(3+state.Nlines)+indgen(1+state.Nlines)]
              (*state.nerrrescube)[*,*,0] = rescube[*,2*(3+state.Nlines)+(1+state.Nlines)+indgen(3+state.Nlines)] ; plus 1-sig
              (*state.nerrrescube)[*,*,1] = rescube[*,3*(3+state.Nlines)+(1+state.Nlines)+indgen(3+state.Nlines)] ; minus 1-sig
              (*state.berrrescube)[*,*,0] = rescube[*,4*(3+state.Nlines)+(1+state.Nlines)+indgen(3+state.Nlines)] ; plus 1-sig
              (*state.berrrescube)[*,*,1] = rescube[*,5*(3+state.Nlines)+(1+state.Nlines)+indgen(3+state.Nlines)] ; minus 1-sig
              (*state.cerrrescube)[*,*,0] = rescube[*,6*(3+state.Nlines)+(1+state.Nlines)+indgen(1+state.Nlines)] ; plus 1-sig
              (*state.cerrrescube)[*,*,1] = rescube[*,6*(3+state.Nlines)+2*(1+state.Nlines)+indgen(1+state.Nlines)] ; minus 1-sig
        
              (*state.nrescube) = kubeviz_remove_badvalues((*state.nrescube))
              (*state.brescube) = kubeviz_remove_badvalues((*state.brescube))
              (*state.crescube) = kubeviz_remove_badvalues((*state.crescube))

              (*state.nerrrescube) = kubeviz_remove_badvalues((*state.nerrrescube),repval=-999)
              (*state.berrrescube) = kubeviz_remove_badvalues((*state.berrrescube),repval=-999)
              (*state.cerrrescube) = kubeviz_remove_badvalues((*state.cerrrescube),repval=-999)
           end
           1: begin
              
              for iline=0,state.Nlines-1 do begin
                (*state.mrescube)[*,iline*6+5]              = rescube[*,iline*6]                                  
                (*state.mrescube)[*,iline*6+indgen(5)]      = rescube[*,iline*6+1+indgen(5)]                      
                (*state.crescube)[*,iline+1]                = rescube[*,6*state.Nlines+iline]                     
                (*state.merrrescube)[*,iline*6+indgen(5),0] = rescube[*,7*state.Nlines+10*iline+indgen(5)]               ; plus 1-sig
                (*state.merrrescube)[*,iline*6+indgen(5),1] = rescube[*,7*state.Nlines+10*iline+5+indgen(5)]             ; minus 1-sig
                (*state.cerrrescube)[*,iline+1,0]           = rescube[*,7*state.Nlines+10*state.Nlines+iline]            ; plus 1-sig
                (*state.cerrrescube)[*,iline+1,1]           = rescube[*,7*state.Nlines+11*state.Nlines+iline]            ; minus 1-sig
              endfor
              
              (*state.mrescube) = kubeviz_remove_badvalues((*state.mrescube))
              (*state.crescube) = kubeviz_remove_badvalues((*state.crescube))
              
              (*state.merrrescube) = kubeviz_remove_badvalues((*state.merrrescube),repval=-999)
              (*state.cerrrescube) = kubeviz_remove_badvalues((*state.cerrrescube),repval=-999)

           end
        endcase      
    end
endcase

kubeviz_linefit_pointerswitch, /save, linefitmode=reslinefit_mode, errmethod=reserror_method
printf, state.log_lun, '[KUBEVIZ] Results file successfully loaded! ' 

if reserror_method ne state.domontecarlo then begin
   printf, state.log_lun, '[WARNING] Loaded results have been produced with an error method that is not currently set.'
   case reserror_method of
    0: printf, state.log_lun, '[WARNING] Switch to "Noise cube" errors to see the results.'
    1: printf, state.log_lun, '[WARNING] Switch to "Bootstrap" errors to see the results.'
    2: printf, state.log_lun, '[WARNING] Switch to "Monte Carlo 1" errors to see the results.'
    3: printf, state.log_lun, '[WARNING] Switch to "Monte Carlo 2" errors to see the results.'
    4: printf, state.log_lun, '[WARNING] Switch to "Monte Carlo 3" errors to see the results.'
   endcase 
endif

; This load call is necessary to reload the current error method into the main results container
kubeviz_linefit_pointerswitch, /load
state.linefit_type = reslinefit_type
end

;--------------------------------------------------------------
pro kubeviz_linefit_resetall, linetype=linetype, pars=pars, update=update

common kubeviz_state

if n_elements(update) eq 0 then update=0

if n_elements(linetype) eq 0 then begin
    linetype = ['N','B','C']
    pars=indgen(3+state.Nlines)
endif 

if (size(linetype))[0] eq 0 then linetype=[linetype]
Nlinetype = n_elements(linetype)
Npars = N_elements(pars)

; current values:
kubeviz_linefit_pointerswitch, /load

case state.linefit_mode of
    0: begin ; spaxels:
        nrescube = (*state.sp_nrescube)[state.fitallrange[0]:state.fitallrange[1],state.fitallrange[2]:state.fitallrange[3],*]
        nerrrescube = (*state.sp_nerrrescube)[state.fitallrange[0]:state.fitallrange[1],state.fitallrange[2]:state.fitallrange[3],*,*]
        brescube = (*state.sp_brescube)[state.fitallrange[0]:state.fitallrange[1],state.fitallrange[2]:state.fitallrange[3],*]
        berrrescube = (*state.sp_berrrescube)[state.fitallrange[0]:state.fitallrange[1],state.fitallrange[2]:state.fitallrange[3],*,*]
        crescube = (*state.sp_crescube)[state.fitallrange[0]:state.fitallrange[1],state.fitallrange[2]:state.fitallrange[3],*]
        cerrrescube = (*state.sp_cerrrescube)[state.fitallrange[0]:state.fitallrange[1],state.fitallrange[2]:state.fitallrange[3],*,*]
    end
    1: begin ; masks:
        nrescube = reform((*state.nrescube)[state.imask-1,*])
        nerrrescube = reform((*state.nerrrescube)[state.imask-1,*,*])
        brescube = reform((*state.brescube)[state.imask-1,*])
        berrrescube = reform((*state.berrrescube)[state.imask-1,*,*])
        crescube = reform((*state.crescube)[state.imask-1,*])
        cerrrescube = reform((*state.cerrrescube)[state.imask-1,*,*])
    end
endcase

for ilinetype = 0, Nlinetype-1 do begin

    for ipar=0, Npars-1 do begin
        par = pars[ipar]
        if par gt 0 then nullval = -999 else nullval = 16
        case linetype[ilinetype] of
            'N': begin
        	nrescube[*,*,par] = nullval
        	nerrrescube[*,*,par,*] = nullval
            end
            'B': begin
        	brescube[*,*,par] = nullval
        	berrrescube[*,*,par,*] = nullval
            end
            'C': begin
        	cpar = par-2
        	if cpar gt 0 then begin
        	    crescube[*,*,cpar] = nullval
        	    cerrrescube[*,*,cpar,*] = nullval
        	endif
            end
        endcase
    endfor
endfor

; record in state variable:
case state.linefit_mode of
    0: begin ; spaxels:
        (*state.sp_nrescube)[state.fitallrange[0]:state.fitallrange[1],state.fitallrange[2]:state.fitallrange[3],*] = nrescube
        (*state.sp_nerrrescube)[state.fitallrange[0]:state.fitallrange[1],state.fitallrange[2]:state.fitallrange[3],*,*] = nerrrescube
        (*state.sp_brescube)[state.fitallrange[0]:state.fitallrange[1],state.fitallrange[2]:state.fitallrange[3],*] = brescube
        (*state.sp_berrrescube)[state.fitallrange[0]:state.fitallrange[1],state.fitallrange[2]:state.fitallrange[3],*,*] = berrrescube
        (*state.sp_crescube)[state.fitallrange[0]:state.fitallrange[1],state.fitallrange[2]:state.fitallrange[3],*] = crescube
        (*state.sp_cerrrescube)[state.fitallrange[0]:state.fitallrange[1],state.fitallrange[2]:state.fitallrange[3],*,*] = cerrrescube
    end
    1: begin ; masks:
        (*state.nrescube)[state.imask-1,*] = nrescube
        (*state.nerrrescube)[state.imask-1,*,*] = nerrrescube
        (*state.brescube)[state.imask-1,*] = brescube
        (*state.berrrescube)[state.imask-1,*,*] = berrrescube
        (*state.crescube)[state.imask-1,*] = crescube
        (*state.cerrrescube)[state.imask-1,*,*] = cerrrescube
    end
endcase

kubeviz_linefit_pointerswitch, /save

if update then begin
 if state.cubesel gt 3 and (1+state.fitallrange[1]-state.fitallrange[0]) eq state.Ncol and (1+state.fitallrange[3]-state.fitallrange[2]) eq state.Nrow  then begin
    state.cubesel = 0
 endif

 kubeviz_linefit_update
 kubeviz_plotspeczoom
 kubeviz_plotspax 
 kubeviz_plotinfo
endif


end

;--------------------------------------------------------------
pro kubeviz_linefit_reset, linetype=linetype, pars=pars, update=update
; reset FIT results default is ALL lines 

common kubeviz_state
common kubeviz_widgetids

if n_elements(update) eq 0 then update=0

if n_elements(linetype) eq 0 then begin
    linetype = ['N','B','C']
    pars=indgen(3+state.Nlines)
endif 

if (size(linetype))[0] eq 0 then linetype=[linetype]
Nlinetype = n_elements(linetype)
Npars = N_elements(pars)

if state.debug eq 1 then begin
   msg = 'RESET pars: '
   for i=0,Npars-1 do msg+=kubeviz_str(pars[i],format='(I2)')+' '
   msg += '; linetype='
   for i=0,Nlinetype-1 do msg+=linetype[i]
   print, msg
endif

; current values:
kubeviz_linefit_pointerswitch, /load

case state.linefit_mode of
    0: begin ; spaxels:
        nrescube = reform((*state.sp_nrescube)[state.col,state.row,*])
        nerrrescube = reform((*state.sp_nerrrescube)[state.col,state.row,*,*])
        brescube = reform((*state.sp_brescube)[state.col,state.row,*])
        berrrescube = reform((*state.sp_berrrescube)[state.col,state.row,*,*])
        crescube = reform((*state.sp_crescube)[state.col,state.row,*])
        cerrrescube = reform((*state.sp_cerrrescube)[state.col,state.row,*,*])
    end
    1: begin ; masks:
        nrescube = reform((*state.nrescube)[state.imask-1,*])
        nerrrescube = reform((*state.nerrrescube)[state.imask-1,*,*])
        brescube = reform((*state.brescube)[state.imask-1,*])
        berrrescube = reform((*state.berrrescube)[state.imask-1,*,*])
        crescube = reform((*state.crescube)[state.imask-1,*])
        cerrrescube = reform((*state.cerrrescube)[state.imask-1,*,*])
    end
endcase
for ilinetype = 0, Nlinetype-1 do begin

    for ipar=0, Npars-1 do begin
        print, ipar
	par = pars[ipar]
        if par gt 0 then nullval = -999 else nullval = 16
        case linetype[ilinetype] of
            'N': begin
        	nrescube[par] = nullval
        	nerrrescube[par,*] = nullval
        	if par gt 0 then widget_control, widgetids.nresetbutton[par], set_button=0
            end
            'B': begin
        	brescube[par] = nullval
        	berrrescube[par,*] = nullval
        	if par gt 0 then widget_control, widgetids.bresetbutton[par], set_button=0
            end
            'C': begin
        	cpar = par-2
        	if cpar gt 0 then begin
        	    crescube[cpar] = nullval
        	    cerrrescube[cpar,*] = nullval
        	endif
            end
        endcase
    endfor

endfor

; record in state variable:
case state.linefit_mode of
    0: begin ; spaxels:
        (*state.sp_nrescube)[state.col,state.row,*] = nrescube
        (*state.sp_nerrrescube)[state.col,state.row,*,*] = nerrrescube
        (*state.sp_brescube)[state.col,state.row,*] = brescube
        (*state.sp_berrrescube)[state.col,state.row,*,*] = berrrescube
        (*state.sp_crescube)[state.col,state.row,*] = crescube
        (*state.sp_cerrrescube)[state.col,state.row,*,*] = cerrrescube
    end
    1: begin ; masks:
        (*state.nrescube)[state.imask-1,*] = nrescube
        (*state.nerrrescube)[state.imask-1,*,*] = nerrrescube
        (*state.brescube)[state.imask-1,*] = brescube
        (*state.berrrescube)[state.imask-1,*,*] = berrrescube
        (*state.crescube)[state.imask-1,*] = crescube
        (*state.cerrrescube)[state.imask-1,*,*] = cerrrescube
    end
endcase

kubeviz_linefit_pointerswitch, /save
if update then begin
 kubeviz_linefit_update
 kubeviz_plotspeczoom
endif

end

;----------------------------------------------------------
pro kubeviz_linefit_resetuser, linetype=linetype, pars=pars, ALL=ALL, startonly=startonly
; reset USER parameters
common kubeviz_state

if n_elements(ALL) eq 0 then ALL=0
if n_elements(startonly) eq 0 then startonly=0

if ALL eq 1 then begin
    linetype = ['N','B']
    pars=indgen(3+state.Nlines)
endif 
if (size(linetype))[0] eq 0 then linetype=[linetype]
Nlinetype = n_elements(linetype)
Npars = N_elements(pars)

if state.debug eq 1 then begin
  msg = 'RESET USER pars: '
  for i=0,Npars-1 do msg+=kubeviz_str(pars[i],format='(I2)')+' '
  msg += '; linetype='
  for i=0,Nlinetype-1 do msg+=linetype[i]
  print, msg
endif

for ilinetype = 0, Nlinetype-1 do begin

    for ipar=0, Npars-1 do begin
        par = pars[ipar]
        if par gt 0 then begin
            nullval = -999.
            case linetype[ilinetype] of
                'N': begin
                    state.gnfit[par] = nullval
                    if not startonly then begin
		      state.gnlims[par,0] = nullval
                      state.gnlims[par,1] = nullval
		    endif  
                end
                'B': begin
                    state.gbfit[par] = nullval
                    if not startonly then begin
		     state.gblims[par,0] = nullval
                     state.gblims[par,1] = nullval
		    endif 
                end
            endcase
        endif
    endfor

endfor

kubeviz_linefit_update, /update_userpars
kubeviz_plotspeczoom

end

;--------------------------------------------------------
pro kubeviz_instrres_init
; initialize the instrumental resolution values
common kubeviz_state, state

findpro, 'kubeviz', dirlist=kubeviz_dir, /noprint 

state.instrres_mode = 0

case 1 of
   strlowcase(state.instr) eq 'kmos3d': begin 
	;Resolution keyword must be present in the primary header of K3D products from the official release
	prihead = headfits(state.indir+state.filename, exten=0)
        resorder = kubeviz_fxpar_sp(prihead, 'K3D RES ORDER', /HIER)
	if resorder gt 0 then begin
	   ; Print the following message only once
	   if total(state.instrres_extpoly) eq 0 then printf, state.log_lun, '[KUBEVIZ] Instrumental resolution coefficients extracted from the primary header.'
           ; convert coefficients to work with Angstroms rather than microns:
	   for i=0, resorder do state.instrres_extpoly[i]= kubeviz_fxpar_sp(prihead, 'K3D RES COEFF'+kubeviz_str(i), /HIER)*(1.e-4)^i
	   state.instrres_mode = 1
        endif else begin
	  print, '[WARNING] Polynomial resolution coefficients must be present in the primary header of the official KMOS3D cubes'
	  return
	endelse  
   end
   strlowcase(state.instr) eq 'kmos': begin 
        ; First try and find the resolution coefficients in the primary header of the science file
	; TODO update to V2 keywords
	prihead = headfits(state.indir+state.filename, exten=0)
        resorder = kubeviz_fxpar_sp(prihead, 'RESORDER', /K3D)
	if resorder gt 0 then begin
	   ; Print the following message only once
	   if total(state.instrres_extpoly) eq 0 then printf, state.log_lun, '[KUBEVIZ] Instrumental resolution coefficients extracted from the primary header.'
           ; convert coefficients to work with Angstroms rather than microns:
	   for i=0, resorder do state.instrres_extpoly[i]= kubeviz_fxpar_sp(prihead, 'RESP'+kubeviz_str(i), /K3D)*(1.e-4)^i
	   state.instrres_mode = 1
        endif else begin
	    ; if not found search for the closest (in time) polynomial fit to the instrumental resolution:
	    spawn, 'ls '+kubeviz_dir[0]+'/templates/kmos/KMOSpoly_arc_'+strlowcase(state.band)+'*.txt' , instrres_file_list          
	    Narcs = N_elements(instrres_file_list)
	    arctime = strmid(instrres_file_list, 27, 24, /reverse_offset)
	    scitime = kubeviz_fxpar_sp(prihead, 'DATE-OBS')
	    timedif = fltarr(Narcs)
	    for i =0, Narcs-1 do timedif[i] = abs(kubeviz_timediff(scitime, arctime[i]))
	    mindif = min(timedif, minind)
	    ; and read the resolution coefficients
            readcol, instrres_file_list[minind], ifu, p0, p1, p2, p3, p4, format='I,D,D,D,D,D', comment='#', /silent 
            use = (where(ifu eq state.ifu, Nuse))[0]
            if Nuse eq 1 then begin
            	; Print the following message only once
		if total(state.instrres_extpoly) eq 0 then printf, state.log_lun, '[KUBEVIZ] Instrumental resolution coefficients extracted from the kubeviz archive'
                state.instrres_extpoly[0]=p0[use]
            	state.instrres_extpoly[1]=p1[use]*1.e-4
            	state.instrres_extpoly[2]=p2[use]*1.e-8
            	state.instrres_extpoly[3]=p3[use]*1.e-12
            	state.instrres_extpoly[4]=p4[use]*1.e-16
	    	state.instrres_mode = 1
	    endif
	endelse 
    end
    strlowcase(state.instr) eq 'sinfoni': begin ; use OH_GaussProf templates if available
        case strupcase(state.band) of
            'J' : begin
	          if state.pixscale eq 250 then instr_res_file = kubeviz_dir+'templates/sinfoni/OH_GaussProf_j250.fits' $
		  else printf, state.log_lun, '[WARNING] Pixel scale '+state.pixscale+' not yet supported for band '+state.band
            end
	    'H' : begin
	          if state.pixscale eq 250 then instr_res_file = kubeviz_dir+'templates/sinfoni/OH_GaussProf_h250.fits' $
		  else printf, state.log_lun, '[WARNING] Pixel scale '+state.pixscale+' not yet supported for band '+state.band
	    end
	    'K' : begin
	          case state.pixscale of 
		  100: instr_res_file = kubeviz_dir+'templates/sinfoni/OH_GaussProf_k100.fits'
		  250: instr_res_file = kubeviz_dir+'templates/sinfoni/OH_GaussProf_k250.fits'
		  else: printf, state.log_lun, '[WARNING] Pixel scale '+state.pixscale+' not yet supported for band '+state.band
		  endcase
            end
	    else: printf, state.log_lun, '[WARNING] Band '+state.band+' not yet supported'
	endcase
	
        if file_test(instr_res_file) eq 1 then begin
	
	  tpl_raw_data = readfits(instr_res_file, tpl_hdr, /silent)
          crvalt = fxpar(tpl_hdr, 'CRVAL1')
          cdeltt = fxpar(tpl_hdr, 'CDELT1')
          crpixt = fxpar(tpl_hdr, 'CRPIX1')
	  naxist = fxpar(tpl_hdr, 'NAXIS1')
          if ((crvalt eq 0) and (cdeltt eq 0) and (crpixt eq 0)) then begin
          crvalt = fxpar(tpl_hdr, 'CRVAL2')
          cdeltt = fxpar(tpl_hdr, 'CDELT2')
          crpixt = fxpar(tpl_hdr, 'CRPIX2')
	  naxist = fxpar(tpl_hdr, 'NAXIS2')
          endif
        
          tlambda = ((dindgen(naxist) + 1 - crpixt) * cdeltt + crvalt)
          tpl_fit = mpfitpeak(tlambda, tpl_raw_data, params, nterms=3, /positive)
          state.instrres_tplsig = params[2] * 1.e4 ;Convert to Angstrom
	  state.instrres_mode = 2
	  
	endif else begin
	  if state.debug eq 1 then printf, state.log_lun, '[WARNING] Instrumental resolution templates missing' 
	  printf, state.log_lun, '[KUBEVIZ] Fitting skylines in the noisecube.'
	  kubeviz_linefit_skylines
	endelse
    end
    strlowcase(state.instr) eq 'muse': begin 
      if total(state.instrres_extpoly) eq 0 then printf, state.log_lun, '[KUBEVIZ] Instrumental resolution obtained from the MUSE manual (Figure 11)'
      state.instrres_extpoly[0]=499.446
      state.instrres_extpoly[1]=-0.201608
      state.instrres_extpoly[2]=0.000130290
      state.instrres_extpoly[3]=-7.87479e-09
      state.instrres_extpoly[4]=0
      state.instrres_mode = 1
    end
    strlowcase(state.instr) eq 'vimos': begin
       printf, state.log_lun, '[ DEBUG ] Instrumental resolution is not properly implemented'
       printf, state.log_lun, '[KUBEVIZ] Please input the coefficients manually.'    
       ;QUICK FIX for VIMOS data  -- TODO handle different grisms
;       state.instrres_varpoly[0] = 2650
       state.instrres_extpoly[0]=2650
       state.instrres_extpoly[1]=0.
       state.instrres_extpoly[2]=0.
       state.instrres_extpoly[3]=0.
       state.instrres_extpoly[4]=0
      
       state.instrres_mode = 1
    end
    strlowcase(state.instr) eq 'sami': begin
       printf, state.log_lun, '[ DEBUG ] Instrumental resolution is not properly implemented'
       printf, state.log_lun, '[KUBEVIZ] Please input the coefficients manually.'    
       ;QUICK FIX for red grism
;       state.instrres_varpoly[0] = 4580
       state.instrres_extpoly[0]=4580.
       state.instrres_extpoly[1]=0.
       state.instrres_extpoly[2]=0.
       state.instrres_extpoly[3]=0.
       state.instrres_extpoly[4]=0.
       state.instrres_mode = 1
    end
    strlowcase(state.instr) eq 'wifes': begin
          case strupcase(state.band) of
            'B': state.instrres_extpoly[0]=3000.
	    'R': state.instrres_extpoly[0]=7000.
	  endcase 
          state.instrres_extpoly[1]=0.
          state.instrres_extpoly[2]=0.
          state.instrres_extpoly[3]=0.
          state.instrres_extpoly[4]=0.
          state.instrres_mode = 1
    end
    else: kubeviz_linefit_skylines
 endcase

end

;------------------------------------------
function kubeviz_chooselines, linesets, lines, lineset=lineset
common kubeviz_state, state

; JAG - use only a subset of the emission lines. Use global selection if none passed as argument:
if n_elements(lineset) eq 0 then selected_lineset = state.selected_lineset else selected_lineset=lineset
if state.debug eq 1 then print, '[DEBUG] selected_lineset: ', selected_lineset

; shorten to include only those with their central wavelength +/-10A of ends of the wavelength range
oklinesinspec = where(lines-(*state.wave)[0] gt 10. and (*state.wave)[state.Nwpix-1]-lines gt 10., Nlines)
oklinesinset  = where(linesets eq selected_lineset, Nlines)
; selected_lineset>0: selected lineset only:
if selected_lineset gt 0 then oklines = kubeviz_SetIntersection(oklinesinset,oklinesinspec) else oklines = oklinesinspec

return, oklines
end

;----------------------------------------------
function kubeviz_create_rescubes, type
;
; One for narrow lines, one for broad lines, one for underlying continuum
; and one each for the errors on these parameters.
; FLAG is global - i.e. identical for all results of a given mask or spaxel.

; Line results cubes:
; plane0  :  flag  flag1 (byte 1; 0: good  1: bad) 
; plane1  :  velocity offset (narrow)
; plane2  :  velocity width (narrow)
; plane3 - plane(N+2) : amplitude1 - amplitudeN

; Continuum cubes:
; plane0  :  flag  flag1 (byte 1; 0: good  1: bad) 
; plane1 - plane(N) : amplitude1 - amplitudeN

; Moment cubes:
; plane i: moment i -> Sum, mean, standard dev, skewness, kurtosis, flag

common kubeviz_state

case type of
  'm_line': begin
    line_rescube_null   = replicate(0.D, state.maxNmask , 3+state.Nlines)
    line_rescube_null[*,0] = 0 ; default flag is good
    return, line_rescube_null
  end
  'm_cont': begin
    cont_rescube_null   = replicate(0.D, state.maxNmask , 1+state.Nlines)
    cont_rescube_null[*,0] = 0
    return, cont_rescube_null
  end  
  'm_moment': begin
    moment_rescube_null = replicate(0.D, state.maxNmask , max([6*state.Nlines,1]))  ; The max function avoids a 0 dimension array if Nlines=0
    for i = 0 , state.Nlines-1 do moment_rescube_null[*,i*6+5] = 0 ; the flag for the moments is independent for each line
    return, moment_rescube_null
  end 
  'm_line_err': return, replicate(-999.D, state.maxNmask , 3+state.Nlines, state.Nmontecarlo_percs ) 
  'm_cont_err': return, replicate(-999.D, state.maxNmask , 1+state.Nlines, state.Nmontecarlo_percs )
  'm_moment_err': return, replicate(-999.D, state.maxNmask , max([6*state.Nlines,1]), state.Nmontecarlo_percs )
  
  'sp_line': begin
    sp_line_rescube_null   = replicate(0.D, state.Ncol, state.Nrow , 3+state.Nlines)
    sp_line_rescube_null[*,*,0] = 0 ; default flag is good
    return, sp_line_rescube_null
  end
  'sp_cont': begin
    sp_cont_rescube_null   = replicate(0.D, state.Ncol, state.Nrow , 1+state.Nlines)
    sp_cont_rescube_null[*,*,0] = 0
    return, sp_cont_rescube_null
  end  
  'sp_moment': begin
    sp_moment_rescube_null = replicate(0.D, state.Ncol, state.Nrow , max([6*state.Nlines,1]))
    for i = 0 , state.Nlines-1 do sp_moment_rescube_null[*,*,i*6+5] = 0 ; the flag for the moments is independent for each line
    return, sp_moment_rescube_null
  end 
  'sp_line_err': return, replicate(-999.D, state.Ncol, state.Nrow, 3+state.Nlines, state.Nmontecarlo_percs)
  'sp_cont_err': return, replicate(-999.D, state.Ncol, state.Nrow, 1+state.Nlines, state.Nmontecarlo_percs)
  'sp_moment_err': return, replicate(-999.D, state.Ncol, state.Nrow, max([6*state.Nlines,1]), state.Nmontecarlo_percs)
endcase

end

;----------------------------------------------
pro kubeviz_linefit_init

common kubeviz_state, state
common kubeviz_linesdb
common flags, flag_ok, flag_bad

;Initialise fitall range
state.fitallrange=[0,state.Ncol-1,0,state.Nrow-1]

;Initialize continuum and fitting range defaults
;Not needed anymore as the first init is done when the
;state structure is created
;kubeviz_set_autoflag_defaults
;kubeviz_set_linefitrange_defaults
;kubeviz_set_continuumfit_defaults

; observed frame wavelengths:
if state.vacuum eq 1 then begin
  ; Print the following message only once
  if total(state.lines_all) eq 0 then printf, state.log_lun, '[KUBEVIZ] The instrument works in vacuum. Rest wavelengths converted to vacuum'
  airtovac, linesdb.lines_rest, lines_rest 
endif else lines_rest = linesdb.lines_rest

lines=lines_rest*(1.+state.redshift)
ok_lines_all = where(lines-(*state.wave)[0] ge 0. and (*state.wave)[state.Nwpix-1]-lines ge 0., Nok_lines_all)

if Nok_lines_all gt 0 then begin
  
  state.Nlines_all = Nok_lines_all
  state.lines_all = lines(ok_lines_all)
  
  ; select all lines which will be fit
  oklines = kubeviz_chooselines(linesdb.linesets,lines)
  Nlines = N_Elements(oklines)

  if oklines[0] ne -1 then begin
      state.lines = lines[oklines] 
      state.linesets = ptr_new(linesdb.linesets[oklines])
      state.linenames = ptr_new(linesdb.linenames[oklines])
      state.linefancynames = ptr_new(linesdb.linefancynames[oklines])
      ; sets of lines to fit simultaneously:
      state.lineset_max = max(linesdb.linesets)
      for i=0, state.lineset_max-1 do state.lineset_ind[i] = ptr_new(where((*state.linesets) eq i+1))
  endif else begin
     printf, state.log_lun, '[WARNING] No valid selected lines in wavelength range!'
     printf, state.log_lun, '[WARNING] Check the selected lineset.'
     Nlines = 0
  endelse
  state.Nlines = Nlines
  
  ; Create bisectors of lines: this is used to define the maximum valid
  ; wavelength range within which to consider a particular line:
  ; NOTE: Has to include bisectors of ALL lines (not just those to be
  ; fit) otherwise it will not work properly.
  
  for i=0, Nlines-1 do begin
      linesabove = where(state.lines_all gt state.lines[i], Nlinesabove)
      if Nlinesabove eq 0 then nextline=-1. else nextline=min(state.lines_all[linesabove])
      linesbelow = where((state.lines_all gt 0 and state.lines_all lt state.lines[i]), Nlinesbelow)
      if Nlinesbelow eq 0 then previousline=-1. else previousline = max(state.lines_all[linesbelow])
      if previousline eq - 1 then range_min = (*state.wave)[0] else range_min = .5*(previousline+state.lines[i])
      if nextline eq -1 then range_max = (*state.wave)[state.Nwpix-1] else range_max = .5*(nextline+state.lines[i])
      state.line_bisectors[i,0] = range_min
      state.line_bisectors[i,1] = range_max
  end

endif else begin
    printf, state.log_lun, '[WARNING] No valid lines in wavelength range!'
    printf, state.log_lun, '[WARNING] Check the redshift value.'
endelse

; Create the results data cubes 
; MASK RESULTS

state.chisq = ptr_new(replicate(0.D, state.maxNmask))
state.nrescube = ptr_new(kubeviz_create_rescubes('m_line'))
state.brescube = ptr_new(kubeviz_create_rescubes('m_line'))
state.crescube = ptr_new(kubeviz_create_rescubes('m_cont'))
state.mrescube = ptr_new(kubeviz_create_rescubes('m_moment'))

state.nerrrescube = ptr_new(kubeviz_create_rescubes('m_line_err'))
state.berrrescube = ptr_new(kubeviz_create_rescubes('m_line_err'))
state.cerrrescube = ptr_new(kubeviz_create_rescubes('m_cont_err'))
state.merrrescube = ptr_new(kubeviz_create_rescubes('m_moment_err'))

state.noi_nrescube = ptr_new(kubeviz_create_rescubes('m_line'))
state.noi_brescube = ptr_new(kubeviz_create_rescubes('m_line'))
state.noi_crescube = ptr_new(kubeviz_create_rescubes('m_cont'))
state.noi_mrescube = ptr_new(kubeviz_create_rescubes('m_moment'))

state.noi_nerrrescube = ptr_new(kubeviz_create_rescubes('m_line_err'))
state.noi_berrrescube = ptr_new(kubeviz_create_rescubes('m_line_err'))
state.noi_cerrrescube = ptr_new(kubeviz_create_rescubes('m_cont_err'))
state.noi_merrrescube = ptr_new(kubeviz_create_rescubes('m_moment_err'))

; SPAXEL RESULTS:

state.sp_chisq = ptr_new(replicate(0.D, state.Ncol, state.Nrow))

state.sp_nrescube = ptr_new(kubeviz_create_rescubes('sp_line'))
state.sp_brescube = ptr_new(kubeviz_create_rescubes('sp_line'))
state.sp_crescube = ptr_new(kubeviz_create_rescubes('sp_cont'))
state.sp_mrescube = ptr_new(kubeviz_create_rescubes('sp_moment'))

state.sp_nerrrescube = ptr_new(kubeviz_create_rescubes('sp_line_err'))
state.sp_berrrescube = ptr_new(kubeviz_create_rescubes('sp_line_err'))
state.sp_cerrrescube = ptr_new(kubeviz_create_rescubes('sp_cont_err'))
state.sp_merrrescube = ptr_new(kubeviz_create_rescubes('sp_moment_err'))

state.noi_sp_nrescube = ptr_new(kubeviz_create_rescubes('sp_line'))
state.noi_sp_brescube = ptr_new(kubeviz_create_rescubes('sp_line'))
state.noi_sp_crescube = ptr_new(kubeviz_create_rescubes('sp_cont'))
state.noi_sp_mrescube = ptr_new(kubeviz_create_rescubes('sp_moment'))

state.noi_sp_nerrrescube = ptr_new(kubeviz_create_rescubes('sp_line_err'))
state.noi_sp_berrrescube = ptr_new(kubeviz_create_rescubes('sp_line_err'))
state.noi_sp_cerrrescube = ptr_new(kubeviz_create_rescubes('sp_cont_err'))
state.noi_sp_merrrescube = ptr_new(kubeviz_create_rescubes('sp_moment_err'))

; Initialize other linefit variables (only if this is the first time)

if ptr_valid(state.spaxselect) eq 0 then  state.spaxselect = ptr_new(fltarr(state.Ncol , state.Nrow, state.Nmask) )
if ptr_valid(state.medspec)    eq 0 then  state.medspec    = ptr_new(fltarr(state.Nwpix))
if ptr_valid(state.nmedspec)   eq 0 then  state.nmedspec   = ptr_new(fltarr(state.Nwpix))

if ptr_valid(state.lineresimg)    eq 0 then  begin
     state.lineresimg        = ptr_new(replicate(0.D, state.Ncol, state.Nrow))
     state.unm_lineresimg    = ptr_new(replicate(0.D, state.Ncol, state.Nrow))
     state.lineerrresimg     = ptr_new(replicate(0.D, state.Ncol, state.Nrow))
     state.unm_lineerrresimg = ptr_new(replicate(0.D, state.Ncol, state.Nrow))
endif else begin  
     (*state.lineresimg)        = replicate(0.D, state.Ncol, state.Nrow)
     (*state.unm_lineresimg)    = replicate(0.D, state.Ncol, state.Nrow)
     (*state.lineerrresimg)     = replicate(0.D, state.Ncol, state.Nrow)
     (*state.unm_lineerrresimg) = replicate(0.D, state.Ncol, state.Nrow)
  endelse

; Initialize the instrumental resolution values only if required
case state.instrres_mode of   0: if total(state.instrres_varpoly) eq 0 then kubeviz_instrres_init 
   1: if total(state.instrres_extpoly) eq 0 then kubeviz_instrres_init 
   2: if total(state.instrres_tplsig)  eq 0 then kubeviz_instrres_init 
endcase
end

;-----------------------------------------------------------
pro kubeviz_linefit_typeswitch

common kubeviz_state
common kubeviz_widgetids

case state.linefit_type of
  0:  begin
    str_linefit_type  = "GAUSS"
    case state.secondcomp_mode of
       0: begin
        str_linefit_1stcomp = '1st Comp.'
	str_linefit_2ndcomp = '2nd Comp.'
       end
       1: begin
        str_linefit_1stcomp = 'Blue'
	str_linefit_2ndcomp = 'Red'
       end
       2: begin
        str_linefit_1stcomp = 'Narrow'
	str_linefit_2ndcomp = 'Broad'
       end
    endcase
    for i = 1,2 do begin
       widget_control, widgetids.ncomp_id[i], set_value= str_linefit_1stcomp
       widget_control, widgetids.bcomp_id[i], set_value= str_linefit_2ndcomp
    endfor
    for i = 3,state.Nlines+2 do begin
      widget_control, widgetids.ncomp_id[i], set_value= str_linefit_1stcomp
      widget_control, widgetids.bcomp_id[i], set_value= str_linefit_2ndcomp
    endfor
    for i = 1,state.Nlines+2 do begin
       widget_control, widgetids.gnfit_id[i],      map=1
       widget_control, widgetids.gbfit_id[i],      map=1
       widget_control, widgetids.gnlims_id[i],     map=1
       widget_control, widgetids.gblims_id[i],     map=1
       widget_control, widgetids.nfixbutton[i],    map=1
       widget_control, widgetids.bfixbutton[i],    map=1
       widget_control, widgetids.bimagebutton[i],  map=1
    endfor
    widget_control, widgetids.columns_id[3] ,   set_value='start'
    widget_control, widgetids.columns_id[5] ,   set_value='min'
    widget_control, widgetids.columns_id[6] ,   set_value='max'
    widget_control, widgetids.columns_id[7] ,   set_value='fix'
    widget_control, widgetids.columns_id[8] ,   set_value='reset'
    widget_control, widgetids.fitconstrlabel ,  set_value='Fit with constraints?'
    widget_control, widgetids.fitconstrbutton,  map=1
    widget_control, widgetids.fixratioslabel ,  set_value='Fix line ratios?'
    widget_control, widgetids.fixratiosbutton,  map=1
    widget_control, widgetids.secondcomplabel , set_value='2nd component: '
    widget_control, widgetids.secondcompbuttons[0], map=1
    widget_control, widgetids.secondcompbuttons[1], map=1
    widget_control, widgetids.secondcompbuttons[2], map=1
    widget_control, widgetids.mom_thresh_id  ,  map=0
    ;widget_control, widgetids.fitoption_id[2],  sensitive=1
  end
  1:  begin
    str_linefit_type = "MOMENTS"
    for i = 1,2 do begin
      widget_control, widgetids.ncomp_id[i], set_value= 'Main Line'
      widget_control, widgetids.bcomp_id[i], set_value= ''
    endfor
    for i = 3,state.Nlines+2 do begin
      widget_control, widgetids.ncomp_id[i], set_value= 'Flux'
      widget_control, widgetids.bcomp_id[i], set_value= '1st/2nd'
    endfor
    for i = 1,state.Nlines+2 do begin
       widget_control, widgetids.gnfit_id[i],      map=0
       widget_control, widgetids.gbfit_id[i],      map=0
       widget_control, widgetids.gnlims_id[i],     map=0
       widget_control, widgetids.gblims_id[i],     map=0
       widget_control, widgetids.nfixbutton[i],    map=0
       widget_control, widgetids.bfixbutton[i],    map=0
       widget_control, widgetids.bimagebutton[i],  map=0
    endfor
    widget_control, widgetids.columns_id[3] ,   set_value=''
    widget_control, widgetids.columns_id[5] ,   set_value=''
    widget_control, widgetids.columns_id[6] ,   set_value=''
    widget_control, widgetids.columns_id[7] ,   set_value=''
    widget_control, widgetids.columns_id[8] ,   set_value=''
    widget_control, widgetids.fitconstrlabel ,  set_value='Moments thresh:'
    widget_control, widgetids.fitconstrbutton,  map=0
    widget_control, widgetids.fixratioslabel ,  set_value=''
    widget_control, widgetids.fixratiosbutton,  map=0
    widget_control, widgetids.secondcomplabel , set_value=''
    widget_control, widgetids.secondcompbuttons[0], map=0
    widget_control, widgetids.secondcompbuttons[1], map=0
    widget_control, widgetids.secondcompbuttons[2], map=0
    widget_control, widgetids.mom_thresh_id   ,  map=1
    ;widget_control, widgetids.fitoption_id[2],  sensitive=0
  end
endcase

widget_control, widgetids.linefit_type_id, set_value=str_linefit_type

end

;-----------------------------------------------------------
pro kubeviz_linefit_update, update_userpars=update_userpars, fast_update=fast_update
; update the linefit window
common kubeviz_state
common kubeviz_widgetids
common flags, flag_ok, flag_bad

if n_elements(fast_update) eq 0 then fast_update=0 

; previous results from the selected mask or spaxel:


kubeviz_linefit_pointerswitch, /load

case state.linefit_mode of
    0: begin ; spaxels:
          nrescube = reform((*state.sp_nrescube)[state.col,state.row,*])
          brescube = reform((*state.sp_brescube)[state.col,state.row,*])
          crescube = reform((*state.sp_crescube)[state.col,state.row,*])
          mrescube = reform((*state.sp_mrescube)[state.col,state.row,*])
          nerrrescube = reform((*state.sp_nerrrescube)[state.col,state.row,*,*])
          berrrescube = reform((*state.sp_berrrescube)[state.col,state.row,*,*])
          cerrrescube = reform((*state.sp_cerrrescube)[state.col,state.row,*,*])
          merrrescube = reform((*state.sp_merrrescube)[state.col,state.row,*,*])
    end
    1: begin ; masks:
          nrescube = reform((*state.nrescube)[state.imask-1,*])
          brescube = reform((*state.brescube)[state.imask-1,*])
          crescube = reform((*state.crescube)[state.imask-1,*])
          mrescube = reform((*state.mrescube)[state.imask-1,*])
          nerrrescube = reform((*state.nerrrescube)[state.imask-1,*,*])
          berrrescube = reform((*state.berrrescube)[state.imask-1,*,*])
          cerrrescube = reform((*state.cerrrescube)[state.imask-1,*,*])
          merrrescube = reform((*state.merrrescube)[state.imask-1,*,*])
    end
endcase


; mask selection:
case state.linefit_mode of
    0: begin
	str_linefit_mode  = "SPAXEL"
	str_linefit_sel   = "["+kubeviz_str(state.col)+","+kubeviz_str(state.row)+"]"
	str_linefit_chisq = kubeviz_str((*state.sp_chisq)[state.col,state.row],format='(f10.4)')
	widget_control, widgetids.stepmask_id, map=0
    end
    1: begin
	str_linefit_mode = "MASK"
	str_linefit_sel = kubeviz_str(state.imask)+'/'+kubeviz_str(state.Nmask)
	str_linefit_chisq = kubeviz_str((*state.chisq)[state.imask-1],format='(f10.4)')
	widget_control, widgetids.stepmask_id, map=1
    end
endcase
widget_control, widgetids.linefit_mode_id, set_value=str_linefit_mode
widget_control, widgetids.imask_lf_id, set_value=str_linefit_sel

; Results:
if state.linefit_type eq 0 then begin ; GAUSSIAN fit
   
   case strmid(state.par_imagebutton,0,1) of
      'N': flag =  nrescube[0]
      'B': flag =  brescube[0]
      'C': flag =  crescube[0]
      else: flag = nrescube[0]
   endcase    
   
   for par = 1, 2 do begin
       widget_control, widgetids.pnfit_id[par], set_value=kubeviz_resultsstring(nrescube[par],[nerrrescube[par,0],nerrrescube[par,1]])
       widget_control, widgetids.pbfit_id[par], set_value=kubeviz_resultsstring(brescube[par],[berrrescube[par,0],berrrescube[par,1]])
   endfor
   for iline = 0, state.Nlines-1 do begin
       par = iline+3
       cpar = iline+1
       widget_control, widgetids.pnfit_id[par], set_value=kubeviz_resultsstring(nrescube[par],[nerrrescube[par,0],nerrrescube[par,1]])
       widget_control, widgetids.pbfit_id[par], set_value=kubeviz_resultsstring(brescube[par],[berrrescube[par,0],berrrescube[par,1]])
       widget_control, widgetids.pcfit_id[cpar], set_value=kubeviz_resultsstring(crescube[cpar],[cerrrescube[cpar,0],cerrrescube[cpar,1]], format='(F10.4)')
   endfor
   widget_control, widgetids.chisq_id, set_value=str_linefit_chisq


endif else begin ;Moments
   
   if state.Nlines gt 0 then begin
     main_ind = where(abs(state.lines - kubeviz_getmainline()) lt 1) ; Very dangerous piece of code, it works unless there are lines closer than 1A
     flag = mrescube[6*main_ind+5]
     for par = 1, 2 do begin
       widget_control, widgetids.pnfit_id[par], set_value=kubeviz_resultsstring(mrescube[6*main_ind+par],[merrrescube[6*main_ind+par,0],merrrescube[6*main_ind+par,1]])
       widget_control, widgetids.pbfit_id[par], set_value=''
     endfor
     for iline = 0, state.Nlines-1 do begin
       respar = 6*iline
       par = iline+3
       cpar = iline+1
       widget_control, widgetids.pnfit_id[par], set_value=kubeviz_resultsstring(mrescube[respar],[merrrescube[respar,0],merrrescube[respar,1]])
       widget_control, widgetids.blamb_id[par], set_value=kubeviz_resultsstring(mrescube[respar+1],[merrrescube[respar+1,0],merrrescube[respar+1,1]])
       widget_control, widgetids.pbfit_id[par], set_value=kubeviz_resultsstring(mrescube[respar+2],[merrrescube[respar+2,0],merrrescube[respar+2,1]])
       widget_control, widgetids.pcfit_id[cpar], set_value=kubeviz_resultsstring(crescube[cpar],[cerrrescube[cpar,0],cerrrescube[cpar,1]], format='(F10.4)')
    endfor
   endif 
   widget_control, widgetids.chisq_id, set_value='--' 

endelse

; update line position:
for iline=0, state.Nlines-1 do begin
    if state.linefit_type eq 0 then begin
     ; narrow line:
     if nerrrescube[1,0] eq -999. then ndv = 0. else ndv = nrescube[1]
     lambda_best_narrow = state.lines[iline]*(1.+ndv/state.ckms)
     widget_control, widgetids.nlamb_id[iline+3], set_value=kubeviz_str(lambda_best_narrow,format='(F8.2)')
     ; broad line:
     if berrrescube[1,0] eq -999. then bdv = 0. else bdv = brescube[1]
     lambda_best_broad = state.lines[iline]*(1.+bdv/state.ckms)
     widget_control, widgetids.blamb_id[iline+3], set_value=kubeviz_str(lambda_best_broad,format='(F8.2)')
    endif else begin
     lambda_best_narrow = state.lines[iline]*(1.+mrescube[6*iline+1]/state.ckms)
     widget_control, widgetids.nlamb_id[iline+3], set_value=kubeviz_str(lambda_best_narrow,format='(F8.2)')
    endelse 
endfor

if fast_update eq 0 then begin
   
   ; Flag widget
   if flag gt 0 then widget_control, widgetids.flag_button, set_value=flag_bad else widget_control, widgetids.flag_button, set_value=flag_ok

   ; redshift:
   widget_control, widgetids.redshift_id, set_value=kubeviz_str(state.redshift, format='(F9.6)')

   ; instrumental resolution:
   widget_control, widgetids.instrres_id, set_value=kubeviz_str(kubeviz_getinstrres(), format='(f12.2)')
   case state.instrres_mode of
       0: str_instrres_mode = '(Polynomial fit to cube variance)'
       1: str_instrres_mode = '(Polynomial fit to external arcs)'
       2: str_instrres_mode = '(Spectral Templates)'
   endcase
   widget_control, widgetids.instrres_mode_id, set_value=str_instrres_mode

   ; Error method control
   case state.domontecarlo of
       0: str_errormethod = "Noise Cube"
       1: str_errormethod = "Bootstrap"
       2: str_errormethod = "MonteCarlo 1"
       3: str_errormethod = "MonteCarlo 2"
       4: str_errormethod = "MonteCarlo 3"
   endcase
   widget_control,  widgetids.errormethod_id, set_value=str_errormethod

   ;Various Options plots control
   if state.plotMonteCarlodistrib eq 0 then str_mcplot = "OFF" else str_mcplot = "ON"
   if state.saveMonteCarlodistrib eq 0 then str_mcsave = "OFF" else str_mcsave = "ON"
   if state.useMonteCarlonoise    eq 0 then str_mcnois = "OFF" else str_mcnois = "ON"
   if state.scaleNoiseerrors      eq 0 then str_scnois = "OFF" else str_scnois = "ON"
   case state.linefit_type of 
      0: if ptr_valid(state.gauss_initmap) then str_inputmap_flag = "ON" else str_inputmap_flag = "OFF"
      1: if ptr_valid(state.mom_windowmap) then str_inputmap_flag = "ON" else str_inputmap_flag = "OFF"
   endcase
  
   widget_control, widgetids.montecarloplot_id,  set_value=str_mcplot
   widget_control, widgetids.montecarlosave_id,  set_value=str_mcsave
   widget_control, widgetids.montecarlonoise_id, set_value=str_mcnois
   widget_control, widgetids.scalenoise_id,      set_value=str_scnois
   widget_control, widgetids.inputmap_flag,      set_value=str_inputmap_flag
   
   ;Fit range and continuum parameters control
   widget_control, widgetids.maxwoffb_id, set_value=string(state.maxwoffb,format='(F8.2)')
   widget_control, widgetids.maxwoffr_id, set_value=string(state.maxwoffr,format='(F8.2)')
   widget_control, widgetids.mom_thresh_id, set_value=string(state.mom_thresh,format='(F8.2)')
   widget_control, widgetids.masksnthresh_id,  set_value=string(state.mask_sn_thresh,format='(F8.2)')
   widget_control, widgetids.maskmaxvelerr_id, set_value=string(state.mask_maxvelerr,format='(F8.2)')

   widget_control, widgetids.continuummode_button, set_button=state.continuumfit_mode
   widget_control, widgetids.continuumfit_minoff_id, set_value=string(state.continuumfit_minoff,format='(F8.2)')
   widget_control, widgetids.continuumfit_maxoff_id, set_value=string(state.continuumfit_maxoff,format='(F8.2)')
   widget_control, widgetids.continuumfit_minperc_id, set_value=string(state.continuumfit_minperc,format='(F8.2)')
   widget_control, widgetids.continuumfit_maxperc_id, set_value=string(state.continuumfit_maxperc,format='(F8.2)')
   widget_control, widgetids.continuumfit_order_id, set_value=string(state.continuumfit_order,format='(I1)')
      
   ; Also update user-selection parameters to the stored ones?
   ; default is off (i.e. stick to those parameters currently set)
   if N_elements(update_userpars) eq 0 then update_userpars = 0
   if update_userpars eq 1 then begin

       for par=1, state.Nlines+2 do begin 
   	   widget_control, widgetids.gnfit_id[par],    set_value=kubeviz_userparstring(state.gnfit[par])  
   	   widget_control, widgetids.gbfit_id[par],    set_value=kubeviz_userparstring(state.gbfit[par])  
   	   widget_control, widgetids.gnlims_id[par,0], set_value=kubeviz_userparstring(state.gnlims[par,0])
   	   widget_control, widgetids.gblims_id[par,0], set_value=kubeviz_userparstring(state.gblims[par,0])
   	   widget_control, widgetids.gnlims_id[par,1], set_value=kubeviz_userparstring(state.gnlims[par,1])
   	   widget_control, widgetids.gblims_id[par,1], set_value=kubeviz_userparstring(state.gblims[par,1])
   	   ; buttons:
   	   widget_control, widgetids.nfixbutton[par], set_button=state.pnfix[par]
   	   widget_control, widgetids.bfixbutton[par], set_button=state.pbfix[par]
   	   if par gt 2 then begin
   	       line = par-3
   	       widget_control, widgetids.pndofitbutton[par], set_button=state.pndofit[line]
   	       widget_control, widgetids.pbdofitbutton[par], set_button=state.pbdofit[line]
   	       widget_control, widgetids.pcdofitbutton[par], set_button=state.pcdofit[line]
   	       widget_control, widgetids.nshowbutton[par], set_button=state.nshow[line]
   	       widget_control, widgetids.bshowbutton[par], set_button=state.bshow[line]     
   	       widget_control, widgetids.cshowbutton[par], set_button=state.cshow[line] 
   	   endif
       endfor

       widget_control, widgetids.fitconstrbutton, set_button=state.fitconstr
       widget_control, widgetids.fixratiosbutton, set_button=state.fitfixratios
       widget_control, widgetids.secondcompbuttons[state.secondcomp_mode], set_button=1
       widget_control, widgetids.secondcompsmartbutton, set_button=state.secondcomp_smart
       
       case state.instrres_mode of
         0: for i=0, state.max_polycoeff_instrres do widget_control, widgetids.polycoeffres_id[i], set_value=kubeviz_str(state.instrres_varpoly[i]) 
         1: for i=0, state.max_polycoeff_instrres do widget_control, widgetids.polycoeffres_id[i], set_value=kubeviz_str(state.instrres_extpoly[i]) 
         else:
       endcase 

    endif
 endif
end

;-----------------------------------------------------------------------
pro kubeviz_linefit_image_update, newpar

common kubeviz_state
common kubeviz_widgetids
common flags, flag_ok, flag_bad

; Save the old par value
oldpar = state.par_imagebutton

;Assign the new value to the state value
state.par_imagebutton = newpar

if newpar ne '' then begin 
  
  ;Read current values
  kubeviz_linefit_pointerswitch, /load

  if state.linefit_type eq 0 then begin
   nrescube    = (*state.sp_nrescube)
   brescube    = (*state.sp_brescube)
   crescube    = (*state.sp_crescube)
   nerrrescube = (*state.sp_nerrrescube)
   berrrescube = (*state.sp_berrrescube)
   cerrrescube = (*state.sp_cerrrescube)
  endif else begin
   crescube    = (*state.sp_crescube)
   cerrrescube = (*state.sp_cerrrescube)
   brescube    = kubeviz_moment_reshape((*state.sp_mrescube))
   berrrescube = kubeviz_moment_reshape((*state.sp_merrrescube), /error)
   nrescube    = kubeviz_moment_reshape((*state.sp_mrescube))
   nerrrescube = kubeviz_moment_reshape((*state.sp_merrrescube), /error)
  endelse
    
  case newpar of
    'CHISQ': begin
         
	 okflag = where(crescube[*,*,0] eq 0, Nokflag)
         (*state.lineresimg) = replicate(!VALUES.F_NAN,state.Ncol,state.Nrow)
         (*state.lineerrresimg) = replicate(!VALUES.F_NAN,state.Ncol,state.Nrow)
	 (*state.lineresimg)[okflag] = (*state.sp_chisq)[okflag]
	 (*state.lineerrresimg)[okflag] = (*state.sp_chisq)[okflag] 
  	 (*state.unm_lineresimg) = (*state.sp_chisq)
  	 (*state.unm_lineerrresimg) =  (*state.sp_chisq) 
    end
    'FLAG': begin

  	 (*state.lineresimg) = nrescube[*,*,0]
  	 (*state.unm_lineresimg) = nrescube[*,*,0]
  	 (*state.lineerrresimg) = nerrrescube[*,*,0]
  	 (*state.unm_lineerrresimg) = nerrrescube[*,*,0]

    end
    else: begin
      linetype = strmid(newpar,0,1)
      par = fix(strmid(newpar,1))  
          
      case linetype of
         'N': begin
            okflag = where(nrescube[*,*,0] eq 0, Nokflag)
            parplane = nrescube[*,*,par]
            ; error plans shows average of +/- errors:
            parerrplane = 0.5*(reform(nerrrescube[*,*,par,0])-reform(nerrrescube[*,*,par,1]))
         end
         'B': begin
            okflag = where(brescube[*,*,0] eq 0, Nokflag)
            parplane = brescube[*,*,par]
            parerrplane = 0.5*(reform(berrrescube[*,*,par,0])-reform(berrrescube[*,*,par,1]))
         end
         'C': begin
            okflag = where(crescube[*,*,0] eq 0, Nokflag)
            cpar = par-2
            parplane = crescube[*,*,cpar]
            parerrplane = 0.5*(reform(cerrrescube[*,*,cpar,0])-reform(cerrrescube[*,*,cpar,1]))
         end
      endcase	
    
      if Nokflag eq 0 then begin
         ; NO valid results: revert to previous selection if FLAG is ON
         if state.flagmode eq 1 then begin
            kubeviz_imagebutton_update, newpar, 0
            state.par_imagebutton = oldpar 
         endif
         (*state.unm_lineresimg) = parplane
         (*state.unm_lineerrresimg) = parerrplane
         print, '[WARNING] No spaxels with OK flag.'
      endif else begin
          
           if ((max(abs(parerrplane)) eq 0. or max(abs(parerrplane)) eq 999.) and min(abs(nerrrescube[*,*,par,0]))) gt 998 then begin 
               printf, state.log_lun, '[WARNING] No valid results in this plane as yet.'
               if state.flagmode eq 1 then begin
                  kubeviz_imagebutton_update, newpar, 0
                  state.par_imagebutton = oldpar 
               endif
           endif else begin
               (*state.lineresimg) = replicate(!VALUES.F_NAN,state.Ncol,state.Nrow)
               (*state.lineerrresimg) = replicate(!VALUES.F_NAN,state.Ncol,state.Nrow)
               (*state.lineresimg)[okflag] = parplane[okflag]
               (*state.lineerrresimg)[okflag] = parerrplane[okflag]
               (*state.unm_lineresimg) = parplane
               (*state.unm_lineerrresimg) = parerrplane
               badspaxelvector = where((*state.badpixelimg) eq 1, Nbad)
               if Nbad gt 0 then begin
        	 (*state.unm_lineresimg)[badspaxelvector] = !values.f_nan
        	 (*state.unm_lineerrresimg)[badspaxelvector] = !values.f_nan
               endif 
           endelse
      endelse
     end
  endcase  
    
  ;Keep the button pressed, reverse default behaviour
  ;Unset only if the old button was different to avoid blinking effects
  if state.par_imagebutton ne oldpar then kubeviz_imagebutton_update, oldpar, 0
  kubeviz_imagebutton_update, state.par_imagebutton, 1

endif

end

;---------------------------------------------------------------------
pro kubeviz_imagebutton_update, par, setval

common kubeviz_state
common kubeviz_widgetids

case par of
   'FLAG' : widget_control, widgetids.flag_imagebutton, set_button=setval
   'CHISQ' : widget_control, widgetids.chisq_imagebutton, set_button=setval 
   else: begin
      linetype = strmid(par,0,1)
      linepar  = fix(strmid(par,1))  
      case linetype of
        'N': widget_control, widgetids.nimagebutton[linepar], set_button=setval
        'B': widget_control, widgetids.bimagebutton[linepar], set_button=setval
        'C': widget_control, widgetids.cimagebutton[linepar], set_button=setval
        else:
      endcase
  end
endcase    

end

;-----------------------------------------------------------------------
pro kubeviz_change_redshift, redshift
common kubeviz_state
common kubeviz_widgetids

; reset all parameters with new redshift + Force reset of all buttons
kubeviz_linefit_reset 
state.par_imagebutton = ''
for par=1, state.Nlines+2 do begin ;0 is for the mask
   state.pnfix[par] = 0
   state.pbfix[par] = 0
endfor   
for iline=0, state.Nlines-1 do begin
   state.pndofit[iline] = 0
   state.pbdofit[iline] = 0
   state.pcdofit[iline] = 0
   state.nshow[iline] = 0
   state.bshow[iline] = 0
   state.cshow[iline] = 0
endfor                 
state.fitconstr = 0
state.fitfixratios = 0

; kill current fitwindow, if linefit image is shown revert to datacube
widget_control, widgetids.base5, /destroy 
if state.cubesel gt 3 then begin
 state.cubesel = 0
 kubeviz_plotspax 
 kubeviz_plotinfo
endif

state.redshift = redshift
kubeviz_linefit_init    ;reinitialise lines
kubeviz_setup_linefit   ;realize fitwindow again
xmanager, 'kubeviz_linefit',  widgetids.base5, /no_block
kubeviz_linefit_update, /update_userpars ; and update linefit window

end

;-----------------------------------------------------------------------
pro kubeviz_set_autoflag_defaults
common kubeviz_state
; define the default values for the autoflag 

state.mask_sn_thresh = 3
state.mask_maxvelerr = 50.

end

;-----------------------------------------------------------------------
pro kubeviz_set_linefitrange_defaults
common kubeviz_state
; define the default values for the line fitting region (wrt the main
; line in each set):

state.maxwoffb = 80.
state.maxwoffr = 80.
state.mom_thresh = 0.5

end

;-----------------------------------------------------------------------
pro kubeviz_set_continuumfit_defaults
common kubeviz_state
; define the default values for the continuum fitting regions:
; as for SDSS: default to use the mean value of the spectrum within the
; 40-60%iles in the continuum region defined as the interval within
; +/-[200,500] Angstroms of the lineset.

state.continuumfit_minoff = 200.
state.continuumfit_maxoff = 500.
state.continuumfit_minperc = 40.
state.continuumfit_maxperc = 60.
state.continuumfit_order = 0.

end

;-----------------------------------------------------------------------
pro kubeviz_set_user_fitpars, fitpars
common kubeviz_state

case n_elements(fitpars) of
  6: begin ;Only contpars are given
    state.continuumfit_minoff = fitpars[0]
    state.continuumfit_maxoff = fitpars[1]
    state.continuumfit_minperc = fitpars[2]
    state.continuumfit_maxperc = fitpars[3]
    state.continuumfit_order = fitpars[4]
    state.continuumfit_mode = fitpars[5]
  end
  8: begin ;Also the fitrange is specified
    state.maxwoffb = fitpars[0]
    state.maxwoffr = fitpars[1]
    state.continuumfit_minoff = fitpars[2]
    state.continuumfit_maxoff = fitpars[3]
    state.continuumfit_minperc = fitpars[4]
    state.continuumfit_maxperc = fitpars[5]
    state.continuumfit_order = fitpars[6]
    state.continuumfit_mode = fitpars[7]
  end	     
endcase

end


;-----------------------------------------------------------------------
pro kubeviz_setup_viewers, menubar

common kubeviz_state
common kubeviz_widgetids
common kubeviz_winsize
common flags, flag_ok, flag_bad

; Define the widget bases
base1 = widget_base(title='spaxel viewer: '+state.filename, /col, xoffset=10, mbar=menubar, /TLB_KILL_REQUEST_EVENTS)
base2 = widget_base(title='spectrum viewer: '+state.filename, group_leader=base1, /col, xoffset=winsize.xwin1+30, /TLB_SIZE_EVENTS,  /TLB_KILL_REQUEST_EVENTS)
base3 = widget_base(title='spectral zoom', group_leader=base1, /col, xoffset=winsize.xwin1+30, yoffset=350, map=0, /TLB_KILL_REQUEST_EVENTS, /TLB_SIZE_EVENTS)

widgetids.base1 = base1
widgetids.base2 = base2
widgetids.base3 = base3

;-----------------------
; base 1 - spaxel window

draw1 = widget_draw(base1, xsize=state.Ncol * state.zoomfac, ysize=state.Nrow * state.zoomfac, $
                    /button_events, /keyboard_events, /align_center, /motion_events)

rowcol = widget_base(base1, /row, frame=1)
widgetids.colmin_id = widget_label(rowcol, xsize=65, /align_right, value=' ')
colbar = widget_draw(rowcol, xsize=winsize.xwin1-140, ysize=20)
widgetids.colmax_id = widget_label(rowcol, xsize=65, /align_left, value=' ')

rowcol2 = widget_base(base1, /row, frame=1)
dummy = widget_label(rowcol2, xsize=40, value='Wave:', frame=0, /align_center)
widgetids.slid_id  = widget_slider(rowcol2, /drag, xsize=winsize.xwin1-170, $
                      maximum=state.Nwpix-1, value=state.wpix, /suppress_value)
dummy = widget_label(rowcol2, xsize=10, value='', frame=0, /align_center)
widgetids.spatzoom_id = widget_label(rowcol2, xsize=60, value='Zoom: '+string(state.zoomspax, format='(i2)'), frame=0, /align_left)
dummy = widget_button(rowcol2, value='-', uvalue='LESSSPATZOOM')
dummy = widget_button(rowcol2, value='+', uvalue='MORESPATZOOM')


;-----------------------------------------------------------------------
filemenu = widget_button(menubar, value='File', /menu)  
dummy =  widget_button(filemenu, value='Open', uval='Open')
dummy =  widget_button(filemenu, value='Save Image as EPS', uval='Print')
dummy =  widget_button(filemenu, value='Save Image as FITS', uval='SaveImage')
dummy =  widget_button(filemenu, value='Save Cube as FITS',  uval='SaveCube')
dummy =  widget_button(filemenu, value='Display FITS Header', uval='DisplayHeader')
dummy =  widget_button(filemenu, value='Save Session', uval='SaveSession')
dummy =  widget_button(filemenu, value='Load Session', uval='LoadSession')
dummy =  widget_button(filemenu, value='Quit', uval='Quit')

cubemenu = widget_button(menubar, value='Cube', /menu) 
dummy =  widget_button(cubemenu, value='Data', uval='Data')
dummy =  widget_button(cubemenu, value='Noise', uval='Noise')
dummy =  widget_button(cubemenu, value='BadPixels', uval='BadPixels')
dummy =  widget_button(cubemenu, value='S/N', uval='SN')
dummy =  widget_button(cubemenu, value='-------------')
dummy =  widget_button(cubemenu, value='Linefit', uval='Linefit')
dummy =  widget_button(cubemenu, value='Linefit errors', uval='LineErrors')
dummy =  widget_button(cubemenu, value='Linefit S/N', uval='LineSN')

zcutmenu = widget_button(menubar, value='Zcut', /menu)  
dummy =  widget_button(zcutmenu, value='HistEq', uval='HistEq')
dummy =  widget_button(zcutmenu, value='Zscale', uval='Zscale')
dummy =  widget_button(zcutmenu, value='MinMax', uval='MinMax')
dummy =  widget_button(zcutmenu, value='99.5%', uval='99.5')
dummy =  widget_button(zcutmenu, value='99%',   uval='99.0')
dummy =  widget_button(zcutmenu, value='97%',   uval='97.0')
dummy =  widget_button(zcutmenu, value='95%',   uval='95.0')
dummy =  widget_button(zcutmenu, value='-------------')
dummy =  widget_button(zcutmenu, value='User linear', uval='UserLin')
dummy =  widget_button(zcutmenu, value='User sqrt',  uval='UserSqrt')
dummy =  widget_button(zcutmenu, value='User log10', uval='UserLog')
dummy =  widget_button(zcutmenu, value='User Parameters', uval='UserPars')

colormenu = widget_button(menubar, value='Colour', /menu)  
dummy =  widget_button(colormenu, value='Grey', uval='Grey')
dummy =  widget_button(colormenu, value='BlueWhite', uval='BlueWhite')
dummy =  widget_button(colormenu, value='Heat', uval='Heat')
dummy =  widget_button(colormenu, value='STD GAMMA-II', uval='STD GAMMA-II')
dummy =  widget_button(colormenu, value='Rainbow', uval='Rainbow')
dummy =  widget_button(colormenu, value='Invert', uval='Invert')

cursormenu = widget_button(menubar, value='Cursor', /menu)  
dummy =  widget_button(cursormenu, value='Crosshair', uval='Crosshair')
dummy =  widget_button(cursormenu, value='None', uval='None')

spaxelmenu = widget_button(menubar, value='Spax.Masks', /menu)  
dummy =  widget_button(spaxelmenu, value='Select', uval='Select')
dummy =  widget_button(spaxelmenu, value='Deselect', uval='Deselect')
dummy =  widget_button(spaxelmenu, value='Clear', uval='Clear')
dummy =  widget_button(spaxelmenu, value='Optimal Mask', uval='OptimalMask')
dummy =  widget_button(spaxelmenu, value='---------------')
dummy =  widget_button(spaxelmenu, value='New Mask', uval='NewMask')
dummy =  widget_button(spaxelmenu, value='Save', uval='SAVE')
dummy =  widget_button(spaxelmenu, value='Load', uval='Load')
dummy =  widget_button(spaxelmenu, value='Delete Mask', uval='DeleteMask')
dummy =  widget_button(spaxelmenu, value='---------------')
dummy =  widget_button(spaxelmenu, value='Go to Mask',    uval='GoToMask')
dummy =  widget_button(spaxelmenu, value='Previous Mask', uval='PrevMask')
dummy =  widget_button(spaxelmenu, value='Next Mask',     uval='NextMask')
dummy =  widget_button(spaxelmenu, value='---------------')
dummy =  widget_button(spaxelmenu, value='Mask Parameters', uval='MaskPars')



modemenu = widget_button(menubar, value='Mode', /menu)  
dummy =  widget_button(modemenu, value='Slice', uval='Slice')
dummy =  widget_button(modemenu, value='-------------')
dummy =  widget_button(modemenu, value='Sum1', uval='Sum1')
dummy =  widget_button(modemenu, value='Median1', uval='Median1')
dummy =  widget_button(modemenu, value='Weighted Avg1', uval='WeightedAvg1') 
dummy =  widget_button(modemenu, value='Weighted Med1', uval='WeightedMed1')
dummy =  widget_button(modemenu, value='MedSub1', uval='MedSub1')
dummy =  widget_button(modemenu, value='-------------')
dummy =  widget_button(modemenu, value='Sum2', uval='Sum2')
dummy =  widget_button(modemenu, value='Median2', uval='Median2')
dummy =  widget_button(modemenu, value='Weighted Avg2', uval='WeightedAvg2') 
dummy =  widget_button(modemenu, value='Weighted Med2', uval='WeightedMed2')
dummy =  widget_button(modemenu, value='MedSub2', uval='MedSub2')
dummy =  widget_button(modemenu, value='-------------')
dummy =  widget_button(modemenu, value='Med2-Med1', uval='Med2-Med1')
dummy =  widget_button(modemenu, value='Med1-Med2', uval='Med1-Med2')

errorsmenu = widget_button(menubar, value='Errors', /menu)  
dummy =  widget_button(errorsmenu, value='Use Noise-cube', uval='NoiseErrors')
dummy =  widget_button(errorsmenu, value='Use Bootstraps', uval='BootstrapErrors')
dummy =  widget_button(errorsmenu, value='Use Monte Carlo 1', uval='Mc1Errors')
dummy =  widget_button(errorsmenu, value='Use Monte Carlo 2', uval='Mc2Errors')
dummy =  widget_button(errorsmenu, value='Use Monte Carlo 3', uval='Mc3Errors')

optionsmenu = widget_button(menubar, value='Options', /menu)  
dummy =  widget_button(optionsmenu, value='Smooth parameters', uval='SmoothPars')
dummy =  widget_button(optionsmenu, value='FITALL range', uval='FitallRange')
dummy =  widget_button(optionsmenu, value='Use  Montecarlo Noise (on/off)', uval='MonteCarloNoise')
dummy =  widget_button(optionsmenu, value='Save Montecarlo Plots (on/off)', uval='MonteCarloPlot')
dummy =  widget_button(optionsmenu, value='Save Montecarlo PDFs  (on/off)', uval='MonteCarloSave')
dummy =  widget_button(optionsmenu, value='Scale Noisecube error (on/off)', uval='NoiseCubeErrScale')
dummy =  widget_button(optionsmenu, value='Load Results File', uval='LoadResultFile')

helpmenu = widget_button(menubar, value='Help', /menu)  
dummy    = widget_button(helpmenu, value='What''s new', uval='HelpWhatIsNew')  
dummy    = widget_button(helpmenu, value='Instructions', uval='HelpInstructions')  
dummy    = widget_button(helpmenu, value='Keyboard shortcuts', uval='HelpShortcuts')  

;-----------------------------------------------------------------------

infobase = widget_base(base1, /col, /align_center, frame=1, xsize=winsize.xwin1)
rowbase1 = widget_base(infobase, /row, /align_center, frame=0)
dummy   = widget_label(rowbase1, /align_right, xsize=40, value='Image:')
widgetids.pixima_id  = widget_label(rowbase1, /align_left, xsize=55, value=' ')
dummy   = widget_label(rowbase1, /align_right, xsize=40, value='Phys:')
widgetids.pixphys_id = widget_label(rowbase1, /align_left, xsize=55, value=' ')
widgetids.pixval_id = widget_label(rowbase1, xsize=92, /align_center, value=' ')
dummy   = widget_label(rowbase1, xsize=43, /align_left, value='Mask:')
widgetids.imask_id = widget_label(rowbase1, /align_center, xsize=45, value=' ')

widgetids.cube_id = widget_label(rowbase1, xsize=90, value=' ', frame=1)

rowbase2 = widget_base(infobase, /row, /align_center, frame=0)
dummy = widget_label(rowbase2, /align_right,  xsize=40, value='WCS:')
widgetids.wcs_id = widget_label(rowbase2,/align_left, xsize=235, value=' ')

dummy   = widget_label(rowbase2, xsize=50, /align_right, value='Smooth:')
widgetids.smooth_id = widget_label(rowbase2, /align_center, xsize=55, value=' ')
widgetids.imgmode_id = widget_label(rowbase2, xsize=90, value=' ', frame=1)

;-----------------------------------------------------------------------

; base 2 - spectrum window
draw2 = widget_draw(base2, xsize=winsize.xwin2, ysize=winsize.ywin2, /button_events, /keyboard_events, /motion_events)
rowbase = widget_base(base2, /row)

dummy = widget_label(rowbase, value='Min')
widgetids.zminsp_id = widget_text(rowbase, /editable, xsize=5, ysize=1, uvalue='ZMINSP')
dummy = string(state.zmin_spec, format='(f7.1)')
widget_control, widgetids.zminsp_id, set_value=strtrim(dummy, 2)
dummy = widget_label(rowbase, value='Max')
widgetids.zmaxsp_id = widget_text(rowbase, /editable, xsize=5, ysize=1, uvalue='ZMAXSP')
dummy = string(state.zmax_spec, format='(f7.1)')
widget_control, widgetids.zmaxsp_id, set_value=strtrim(dummy, 2)

button = widget_button(rowbase, value='Fix Scale', uvalue='SCALE')
button = widget_button(rowbase, value='Toggle Zoom', uvalue='ZOOM')
button = widget_button(rowbase, value='Save', uval='SAVESPEC')

kubeviz_xpdmenu, [ '/Sel Range/ {', '/Range1/', '/Range2/', '/Reset/', '}' ], rowbase
kubeviz_xpdmenu, [ '/Mode/ {', '/Spaxel/', '/Sum/', '/Median/', '/W. Avg/', '/MedSub/', '/Optimal/', '}' ], rowbase  

;------------------------------

; base 3 - spectrum zoom window
draw3 = widget_draw(base3, xsize=winsize.xwin3, ysize=winsize.ywin3, /button_events, /keyboard_events)
rowbase0 = widget_base(base3, /row)
button = widget_button(rowbase0, value='  +  ', uvalue='ZOOMIN')
button = widget_button(rowbase0, value='  -  ', uvalue='ZOOMOUT')
widgetids.zoom_id  = widget_label(rowbase0, value=string(2*state.zoomrange, format='(i4)')+' pixel')

widget_control, base1, /realize
widget_control, base2, /realize
widget_control, base3, /realize

; Get their IDs
widget_control, draw1, get_value=wid1
widget_control, draw2, get_value=wid2
widget_control, draw3, get_value=wid3
widget_control, colbar, get_value=wid4

; assign IDs
widgetids.wid1 = wid1
widgetids.wid2 = wid2
widgetids.wid3 = wid3
widgetids.wid4 = wid4
widgetids.draw2 = draw2
widgetids.draw3 = draw3

end

;------------------------------

pro kubeviz_setup_linefit

common kubeviz_state
common kubeviz_widgetids
common kubeviz_winsize
common flags, flag_ok, flag_bad

dimensions = GET_SCREEN_SIZE(RESOLUTION=resolution)
if state.debug eq 1 then print, '[DEBUG] screen size: ', dimensions
def_xsize = 890 
def_ysize = 560 + state.Nlines * 120

state.scroll=0 ;Reset scroll value

if .8*dimensions[0] lt def_xsize then begin
    xsize = .8*dimensions[0]
    state.scroll = 1
endif else xsize = def_xsize
if .9*dimensions[1] lt def_ysize then begin
    xsize = def_xsize + 20  ;20 is the xsize of the scrollbar
    ysize = .9*dimensions[1]
    state.scroll = 1
endif else ysize = def_ysize

if state.scroll eq 1 then $
  base5 = widget_base(title='linefit: '+state.filename, group_leader=base1, /col, xoffset=winsize.xwin1+700, $
                      /TLB_ICONIFY_EVENTS, /TLB_KILL_REQUEST_EVENTS, /TLB_SIZE_EVENTS, /scroll, scr_xsize=xsize, scr_ysize=ysize) else $
  base5 = widget_base(title='linefit: '+state.filename, group_leader=base1, /col, xoffset=winsize.xwin1+700, $
                      /TLB_ICONIFY_EVENTS, /TLB_KILL_REQUEST_EVENTS, /TLB_SIZE_EVENTS)

widgetids.base5 = base5

case state.linefit_mode of
    0: begin ; spaxels:
        nrescube = reform((*state.sp_nrescube)[state.col,state.row,*])
        nerrrescube = reform((*state.sp_nerrrescube)[state.col,state.row,*,*])
        brescube = reform((*state.sp_brescube)[state.col,state.row,*])
        berrrescube = reform((*state.sp_berrrescube)[state.col,state.row,*,*])
        crescube = reform((*state.sp_crescube)[state.col,state.row,*])
        cerrrescube = reform((*state.sp_cerrrescube)[state.col,state.row,*,*])
    end
    1: begin ; masks:
        nrescube = reform((*state.nrescube)[state.imask-1,*])
        brescube = reform((*state.brescube)[state.imask-1,*])
        crescube = reform((*state.crescube)[state.imask-1,*])
        nerrrescube = reform((*state.nerrrescube)[state.imask-1,*,*])
        berrrescube = reform((*state.berrrescube)[state.imask-1,*,*])
        cerrrescube = reform((*state.cerrrescube)[state.imask-1,*,*])
    end
endcase

case state.linefit_mode of
    0: begin
        str_linefit_mode = "SPAXEL"
        str_linefit_sel = "["+kubeviz_str(state.col)+","+kubeviz_str(state.row)+"]"
    end
    1: begin
        str_linefit_mode = "MASK"
        str_linefit_sel = kubeviz_str(state.imask)+'/'+kubeviz_str(state.Nmask)
    end
endcase

case state.linefit_type of 
  0: str_linefit_type = "GAUSS"
  1: str_linefit_type = "MOMENTS"
endcase

case state.domontecarlo of
    0: str_errormethod = "Noise Cube"
    1: str_errormethod = "Bootstrap"
    2: str_errormethod = "MonteCarlo 1"
    3: str_errormethod = "MonteCarlo 2"
    4: str_errormethod = "MonteCarlo 3"
endcase

header_row = widget_base(base5, /row)
dummy = widget_label(header_row, xsize=60, /align_left, value='Fit Type: ')
widgetids.linefit_type_id = widget_label(header_row, xsize=50, /align_left, value=str_linefit_type)
dummy = widget_label(header_row, /align_right, xsize=90, value='Error method: ')
widgetids.errormethod_id = widget_label(header_row, xsize = 90, /align_left, value=str_errormethod)
dummy = widget_label(header_row, xsize=30, /align_right, value='z =')
widgetids.redshift_id = widget_text(header_row, /editable, xsize=8, ysize=1, uvalue='REDSHIFT', value=kubeviz_str(state.redshift, format='(F9.6)'))
dummy = widget_label(header_row, xsize=20, /align_left,   value='')
widgetids.linefit_mode_id = widget_label(header_row, xsize=40, /align_left, value=str_linefit_mode)
widgetids.imask_lf_id = widget_label(header_row, xsize=55, /align_left, value=str_linefit_sel)
widgetids.stepmask_id = widget_base(header_row, /row, map=0)
button = widget_button(widgetids.stepmask_id, value='-', uvalue='PREVMASK')
button = widget_button(widgetids.stepmask_id, value='+', uvalue='NEXTMASK')

dummy = widget_label(header_row, xsize=80, /align_center, value='LIMITS')
dummy = widget_label(header_row, xsize=190, /align_center, value='CONTROL')

columns_row = widget_base(base5, /row)
widgetids.columns_id[0]  = widget_label(columns_row, xsize=90,  /align_left,   value='')
widgetids.columns_id[1]  = widget_label(columns_row, xsize=70,  /align_left,   value='component')
widgetids.columns_id[2]  = widget_label(columns_row, xsize=130, /align_left,   value='lambda_cen')
widgetids.columns_id[3]  = widget_label(columns_row, xsize=70,  /align_center, value='start')
widgetids.columns_id[4]  = widget_label(columns_row, xsize=180, /align_left,   value='best')
widgetids.columns_id[5]  = widget_label(columns_row, xsize=54,  /align_left,   value='min')
widgetids.columns_id[6]  = widget_label(columns_row, xsize=54,  /align_center, value='max')
widgetids.columns_id[7]  = widget_label(columns_row, xsize=35,  /align_right,  value='fix')
widgetids.columns_id[8]  = widget_label(columns_row, xsize=45,  /align_center, value='reset')
widgetids.columns_id[9]  = widget_label(columns_row, xsize=30,  /align_center, value='image')
widgetids.columns_id[10] = widget_label(columns_row, xsize=35,  /align_center, value='fit')
widgetids.columns_id[11] = widget_label(columns_row, xsize=30,  /align_center, value='show')

; initialize par to 1 (par=0 is the flag) central wavelength:
par = 1

; narrow:
strpar = 'N'+kubeviz_str(par,format='(I2)')
row_base = widget_base(base5, /row, ysize=30)
dummy = widget_label(row_base, xsize=90, /align_left, value='Line offset')
widgetids.ncomp_id[par] = widget_label(row_base, xsize=70, /align_left, value='Narrow')
widgetids.nlamb_id[par] = widget_label(row_base, xsize=130, /align_left, value='')
widgetids.gnfit_id[par] = widget_text(widget_base(row_base, /row) , /editable, xsize=8, ysize=1, uvalue='SP'+strpar, value=kubeviz_userparstring(state.gnfit[par]))
widgetids.pnfit_id[par] = widget_label(row_base, xsize=160, /align_left, value=kubeviz_resultsstring(nrescube[par],[nerrrescube[par,0],nerrrescube[par,1]]))
lims_base = widget_base(row_base, /row)
widgetids.gnlims_id[par,0] = widget_text(lims_base, /editable, xsize=8, ysize=1, uvalue='MINP'+strpar, value=kubeviz_userparstring(state.gnlims[par,0]))
widgetids.gnlims_id[par,1] = widget_text(lims_base, /editable, xsize=8, ysize=1, uvalue='MAXP'+strpar, value=kubeviz_userparstring(state.gnlims[par,1]))
button_base1 = widget_base(row_base, /row, /nonexclusive)
button_base2 = widget_base(row_base, /row, /nonexclusive)
widgetids.nfixbutton[par]   = widget_button(button_base1, xsize=35, value=' ', uval='FIX'+strpar)
widgetids.nresetbutton[par] = widget_button(button_base1, xsize=30, value=' ', uval='RESET'+strpar)
widgetids.nimagebutton[par] = widget_button(button_base2, xsize=35, value=' ', uval='IMAGE'+strpar)

; broad:
strpar = 'B'+kubeviz_str(par,format='(I2)')
row_base = widget_base(base5, /row, ysize=30)
dummy = widget_label(row_base, xsize=90, /align_center, value='(km/s)')
widgetids.bcomp_id[par]  = widget_label(row_base, xsize=70, /align_left, value='Broad')
widgetids.blamb_id[par]  = widget_label(row_base, xsize=130, /align_left, value='')
widgetids.gbfit_id[par]  = widget_text(widget_base(row_base, /row), /editable, xsize=8, ysize=1, uvalue='SP'+strpar, value=kubeviz_userparstring(state.gbfit[par]))
widgetids.pbfit_id[par]  = widget_label(row_base, xsize=160, /align_left, value=kubeviz_resultsstring(brescube[par],[berrrescube[par,0],berrrescube[par,1]]))
lims_base = widget_base(row_base, /row)
widgetids.gblims_id[par,0] = widget_text(lims_base , /editable, xsize=8, ysize=1, uvalue='MINP'+strpar, value=kubeviz_userparstring(state.gblims[par,0]))
widgetids.gblims_id[par,1] = widget_text(lims_base, /editable, xsize=8, ysize=1, uvalue='MAXP'+strpar, value=kubeviz_userparstring(state.gblims[par,1]))
button_base1 = widget_base(row_base, /row, /nonexclusive)
button_base2 = widget_base(row_base, /row, /nonexclusive)
widgetids.bfixbutton[par]   = widget_button(button_base1, xsize=35, value=' ', uval='FIX'+strpar)
widgetids.bresetbutton[par] = widget_button(button_base1, xsize=30, value=' ', uval='RESET'+strpar)
widgetids.bimagebutton[par] = widget_button(button_base2, xsize=35, value=' ', uval='IMAGE'+strpar)

; width:
par++

strpar = 'N'+kubeviz_str(par,format='(I2)')
row_base = widget_base(base5, /row, ysize=30)
dummy = widget_label(row_base, xsize=90, /align_left, value='Line width')
widgetids.ncomp_id[par]  = widget_label(row_base, xsize=70, /align_left, value='Narrow')
widgetids.nlamb_id[par]  = widget_label(row_base, xsize=130, /align_left, value='')
widgetids.gnfit_id[par]  = widget_text(widget_base(row_base, /row), /editable, xsize=8, ysize=1, uvalue='SP'+strpar, value=kubeviz_userparstring(state.gnfit[par]))
widgetids.pnfit_id[par]  = widget_label(row_base, xsize=160, /align_left, value=kubeviz_resultsstring(nrescube[par],[nerrrescube[par,0],nerrrescube[par,1]]))
lims_base = widget_base(row_base, /row)
widgetids.gnlims_id[par,0] = widget_text(lims_base, /editable, xsize=8, ysize=1, uvalue='MINP'+strpar, value=kubeviz_userparstring(state.gnlims[par,0]))
widgetids.gnlims_id[par,1] = widget_text(lims_base, /editable, xsize=8, ysize=1, uvalue='MAXP'+strpar, value=kubeviz_userparstring(state.gnlims[par,1]))
button_base1 = widget_base(row_base, /row, /nonexclusive)
button_base2 = widget_base(row_base, /row, /nonexclusive)
widgetids.nfixbutton[par]   = widget_button(button_base1, xsize=35, value=' ', uval='FIX'+strpar)
widgetids.nresetbutton[par] = widget_button(button_base1, xsize=30, value=' ', uval='RESET'+strpar)
widgetids.nimagebutton[par] = widget_button(button_base2, xsize=35, value=' ', uval='IMAGE'+strpar)

; broad:
strpar = 'B'+kubeviz_str(par,format='(I2)')
row_base = widget_base(base5, /row, ysize=30)
dummy = widget_label(row_base, xsize=90, /align_center, value='(km/s)')
widgetids.bcomp_id[par]  = widget_label(row_base, xsize=70, /align_left, value='Broad')
widgetids.blamb_id[par]  = widget_label(row_base, xsize=130, /align_left, value='')
widgetids.gbfit_id[par]  = widget_text(widget_base(row_base, /row) , /editable, xsize=8, ysize=1, uvalue='SP'+strpar, value=kubeviz_userparstring(state.gbfit[par]))
widgetids.pbfit_id[par]  = widget_label(row_base, xsize=160, /align_left, value=kubeviz_resultsstring(brescube[par],[berrrescube[par,0],berrrescube[par,1]]))
lims_base = widget_base(row_base, /row)
widgetids.gblims_id[par,0] = widget_text(lims_base, /editable, xsize=8, ysize=1, uvalue='MINP'+strpar, value=kubeviz_userparstring(state.gblims[par,0]))
widgetids.gblims_id[par,1] = widget_text(lims_base, /editable, xsize=8, ysize=1, uvalue='MAXP'+strpar, value=kubeviz_userparstring(state.gblims[par,1]))
button_base1 = widget_base(row_base, /row, /nonexclusive)
button_base2 = widget_base(row_base, /row, /nonexclusive)
widgetids.bfixbutton[par]   = widget_button(button_base1, xsize=35, value=' ', uval='FIX'+strpar)
widgetids.bresetbutton[par] = widget_button(button_base1, xsize=30, value=' ', uval='RESET'+strpar)
widgetids.bimagebutton[par] = widget_button(button_base2, xsize=35, value=' ', uval='IMAGE'+strpar)

for iline=0, state.Nlines-1 do begin
  par++
  ; narrow:
  strpar = 'N'+kubeviz_str(par,format='(I2)')
  nline_row = widget_base(base5, /row, ysize=30)
  dummy = widget_label(nline_row, xsize=90, /align_left, value=(*state.linefancynames)[iline])
  widgetids.ncomp_id[par]      = widget_label(nline_row, xsize=70, /align_left, value='Narrow')
  lambda_best_narrow = state.lines[iline]*(1.+(*state.nrescube)[state.imask-1,1]/state.ckms)
  widgetids.nlamb_id[par]      = widget_label(nline_row, xsize=130, /align_left, value=kubeviz_str(lambda_best_narrow,format='(F8.2)'))
  widgetids.gnfit_id[par]      = widget_text(widget_base(nline_row, /row), /editable, xsize=8, ysize=1, uvalue='SP'+strpar, value=kubeviz_userparstring(state.gnfit[par]))
  widgetids.pnfit_id[par]      = widget_label(nline_row, xsize=160, /align_left, value=kubeviz_resultsstring(nrescube[par],[nerrrescube[par,0],nerrrescube[par,1]]))
  lims_base		       = widget_base(nline_row, /row)
  widgetids.gnlims_id[par,0]   = widget_text(lims_base, /editable, xsize=8, ysize=1, uvalue='MINP'+strpar, value=kubeviz_userparstring(state.gnlims[par,0]))
  widgetids.gnlims_id[par,1]   = widget_text(lims_base, /editable, xsize=8, ysize=1, uvalue='MAXP'+strpar, value=kubeviz_userparstring(state.gnlims[par,1]))
  nline_buttonbase1	       = widget_base(nline_row, /row, /nonexclusive)
  nline_buttonbase2	       = widget_base(nline_row, /row, /nonexclusive)
  widgetids.nfixbutton[par]    = widget_button(nline_buttonbase1, xsize=35, value=' ', uval='FIX'+strpar)
  widgetids.nresetbutton[par]  = widget_button(nline_buttonbase1, xsize=30, value=' ', uval='RESET'+strpar)
  widgetids.nimagebutton[par]  = widget_button(nline_buttonbase2, xsize=35, value=' ', uval='IMAGE'+strpar)
  widgetids.pndofitbutton[par] = widget_button(nline_buttonbase2, xsize=35, value=' ', uval='FIT'+strpar)
  widgetids.nshowbutton[par]   = widget_button(nline_buttonbase2, xsize=25, value=' ', uval='SHOW'+strpar)
  

   ; broad:
  strpar = 'B'+kubeviz_str(par,format='(I2)')
  bline_row = widget_base(base5, /row, ysize=30)
  dummy 		   = widget_label(bline_row, xsize=90, ysize=30, /align_left, value='')
  widgetids.bcomp_id[par]      = widget_label(bline_row, xsize=70, /align_left, value='Broad')
  lambda_best_broad = state.lines[iline]*(1.+(*state.brescube)[state.imask-1,1]/state.ckms)
  widgetids.blamb_id[par]      = widget_label(bline_row, xsize=130, /align_left, value=kubeviz_str(lambda_best_broad,format='(F8.2)'))
  widgetids.gbfit_id[par]      = widget_text(widget_base(bline_row,  /row), /editable, xsize=8, ysize=1, uvalue='SP'+strpar, value=kubeviz_userparstring(state.gbfit[par]))
  widgetids.pbfit_id[par]      = widget_label(bline_row, xsize=160, /align_left, value=kubeviz_resultsstring(brescube[par],[berrrescube[par,0],berrrescube[par,1]]))
  lims_base		   = widget_base(bline_row, /row)
  widgetids.gblims_id[par,0]   = widget_text(lims_base, /editable, xsize=8, ysize=1, uvalue='MINP'+strpar, value=kubeviz_userparstring(state.gblims[par,0]))
  widgetids.gblims_id[par,1]   = widget_text(lims_base, /editable, xsize=8, ysize=1, uvalue='MAXP'+strpar, value=kubeviz_userparstring(state.gblims[par,1]))
  bline_buttonbase1	   = widget_base(bline_row, /row, /nonexclusive)
  bline_buttonbase2	   = widget_base(bline_row, /row, /nonexclusive)
  widgetids.bfixbutton[par]    = widget_button(bline_buttonbase1, xsize=35, value=' ', uval='FIX'+strpar)
  widgetids.bresetbutton[par]  = widget_button(bline_buttonbase1, xsize=30, value=' ', uval='RESET'+strpar)
  widgetids.bimagebutton[par]  = widget_button(bline_buttonbase2, xsize=35, value=' ', uval='IMAGE'+strpar)
  widgetids.pbdofitbutton[par] = widget_button(bline_buttonbase2, xsize=35, value=' ', uval='FIT'+strpar)
  widgetids.bshowbutton[par]   = widget_button(bline_buttonbase2, xsize=25, value=' ', uval='SHOW'+strpar)

  ; continuum:
  cpar = iline+1
  strpar = 'C'+kubeviz_str(par,format='(I2)')
  cline_row = widget_base(base5, /row, ysize=30)
  dummy = widget_label(cline_row, xsize=90, /align_left, value='')
  widgetids.ccomp_id[par] = widget_label(cline_row, xsize=70, /align_left, value='Continuum')
  dummy = widget_label(cline_row, xsize=130, /align_left, value='')
  dummy_base = widget_base(cline_row,  /row, map=0) 
  dummy = widget_text(dummy_base, xsize=8, /align_left, /editable, value='')
  widgetids.pcfit_id[cpar] = widget_label(cline_row, xsize=160, /align_left, value=kubeviz_resultsstring(crescube[cpar],[cerrrescube[cpar,0],cerrrescube[cpar,1]],format='(F10.4)'))
  dummy_base = widget_base(cline_row,  /row, map=0) 
  dummy = widget_text(dummy_base, xsize=8, /align_left, /editable, value='')
  dummy = widget_text(dummy_base, xsize=8, /align_left, /editable, value='')
  dummy = widget_label(cline_row, xsize=73, /align_left, value='')
  cline_buttonbase = widget_base(cline_row, /row, /nonexclusive)
  widgetids.cimagebutton[par]  = widget_button(cline_buttonbase, xsize=35, value=' ', uval='IMAGE'+strpar)
  widgetids.pcdofitbutton[par] = widget_button(cline_buttonbase, xsize=35, value=' ', uval='FIT'+strpar)
  widgetids.cshowbutton[par]   = widget_button(cline_buttonbase, xsize=25, value=' ', uval='SHOW'+strpar)

endfor

; CHISQ and FLAGS 
fitdetails_row = widget_base(base5, /row, frame=1)
dummy = widget_label(fitdetails_row, value='Red-chi-sq:')
widgetids.chisq_id = widget_label(fitdetails_row, /align_right, xsize=70, value=kubeviz_str(0.D,format='(f10.4)'))
dummy = widget_label(fitdetails_row, /align_left, value=' image:')
flag_buttonbase = widget_base(fitdetails_row, /row, /nonexclusive)
widgetids.chisq_imagebutton = widget_button(flag_buttonbase,  xsize=30, value=' ', uval='IMAGECHISQ')

kubeviz_linefit_define_flags
dummy = widget_label(fitdetails_row, /align_right, xsize=40,  value='Flag: ')
widgetids.flag_button = widget_button(fitdetails_row, value=flag_ok, /BITMAP, uvalue='FLAG', tooltip='Toggle fit status good or bad')
dummy = widget_label(fitdetails_row, /align_left, value=' image:')
flag_buttonbase = widget_base(fitdetails_row, /row, /nonexclusive)
widgetids.flag_imagebutton = widget_button(flag_buttonbase,  xsize=30, value=' ', uval='IMAGEFLAG')

dummy = widget_label(fitdetails_row, xsize=135, /align_left, value='Autoflag S/N thresh:')
widgetids.masksnthresh_id  = widget_text(fitdetails_row, /editable, xsize=8, ysize=1, uvalue='MASKSNTHRESH',  value=string(state.mask_sn_thresh,format='(F8.2)'))
dummy = widget_label(fitdetails_row, xsize=160, /align_center, value=' Max Vel/Sig error (km/s):')
widgetids.maskmaxvelerr_id = widget_text(fitdetails_row, /editable, xsize=8, ysize=1, uvalue='MASKMAXVELERR', value=string(state.mask_maxvelerr,format='(F8.2)'))

; DATA, CONTINUUM Fitting ranges
fitparsbase = widget_base(base5, /col, frame=1)
continuumfit_headerrow = widget_base(fitparsbase, /row)
dummy = widget_label(continuumfit_headerrow, xsize=200, /align_left, value='Fitting Range (l/r side, A)')
dummy = widget_label(continuumfit_headerrow, xsize=120, /align_left, value='Continuum Fit:')
dummy = widget_label(continuumfit_headerrow, xsize=190, /align_center, value='Range of Offset (min/max, A)')
dummy = widget_label(continuumfit_headerrow, xsize=190, /align_center, value='Range of Percentiles included')
dummy = widget_label(continuumfit_headerrow, xsize=100, /align_center, value='Poly. Order')

continuumfit_row = widget_base(fitparsbase, /row)
dummy = widget_label(continuumfit_row, xsize=10, /align_left, value='')
widgetids.maxwoffb_id = widget_text(continuumfit_row, /editable, xsize=8, uvalue='LINEFITMAXOFFB', value=string(state.maxwoffb,format='(F8.2)'))
dummy = widget_label(continuumfit_row, xsize=10, /align_left, value='')
widgetids.maxwoffr_id = widget_text(continuumfit_row, /editable, xsize=8, uvalue='LINEFITMAXOFFR', value=string(state.maxwoffr,format='(F8.2)'))
dummy = widget_label(continuumfit_row, xsize=40, /align_left, value='')
contmode_buttonbase  = widget_base(continuumfit_row, /row, /nonexclusive)
widgetids.continuummode_button = widget_button(contmode_buttonbase, value='MPFIT CONT', uval='CONTMODE')
dummy = widget_label(continuumfit_row, xsize=40, /align_left, value='')
widgetids.continuumfit_minoff_id = widget_text(continuumfit_row, /editable, xsize=8, ysize=1, uvalue='CONTMINOFF', value=string(state.continuumfit_minoff,format='(F8.2)'))
dummy = widget_label(continuumfit_row, xsize=10, /align_left, value='')
widgetids.continuumfit_maxoff_id = widget_text(continuumfit_row, /editable, xsize=8, ysize=1, uvalue='CONTMAXOFF', value=string(state.continuumfit_maxoff,format='(F8.2)'))
dummy = widget_label(continuumfit_row, xsize=40, /align_left, value='')
widgetids.continuumfit_minperc_id = widget_text(continuumfit_row, /editable, xsize=8, ysize=1, uvalue='CONTMINPERC', value=string(state.continuumfit_minperc,format='(F8.2)'))
dummy = widget_label(continuumfit_row, xsize=10, /align_left, value='')
widgetids.continuumfit_maxperc_id = widget_text(continuumfit_row, /editable, xsize=8, ysize=1, uvalue='CONTMAXPERC', value=string(state.continuumfit_maxperc,format='(F8.2)'))
dummy = widget_label(continuumfit_row, xsize=50, /align_left, value='')
widgetids.continuumfit_order_id = widget_text(continuumfit_row, /editable, xsize=1, ysize=1, uvalue='CONTORDER', value=string(state.continuumfit_order,format='(I1)'))

; Line fitting region:
fitoptionsbase = widget_base(base5, /col, frame=1)
row1 = widget_base(fitoptionsbase, /row)
widgetids.fitconstrlabel        = widget_label(row1, xsize=135, /align_left, value='Fit with constraints?')
fitconstraints_buttonbase       = widget_base(row1, /row, /nonexclusive)
widgetids.fitconstrbutton       = widget_button(fitconstraints_buttonbase, value=' ', uval='FITCONSTRAINTS')
momthresh_base                  = widget_base(row1, /row, map=0)
widgetids.mom_thresh_id         = widget_text(momthresh_base, /editable, xsize=8, ysize=1, uvalue='MOMTHRESH', value=string(state.mom_thresh,format='(F8.2)'))
widgetids.fixratioslabel        = widget_label(row1, xsize=100, /align_right, value='Fix line ratios?')
fixratios_buttonbase            = widget_base(row1, /row, /nonexclusive)
widgetids.fixratiosbutton       = widget_button(fixratios_buttonbase, value=' ', uval='FIXRATIOS')
widgetids.secondcomplabel       = widget_label(row1, xsize=90, /align_left, value=' 2nd component:')
secondcomp_buttonbase           = widget_base(row1, /row, /exclusive)
widgetids.secondcompbuttons[0]  = widget_button(secondcomp_buttonbase, value='FAINTER' , uvalue='FAINT_2ND')
widgetids.secondcompbuttons[1]  = widget_button(secondcomp_buttonbase, value='LARGER OFFSET', uvalue='HIVEL_2ND')
widgetids.secondcompbuttons[2]  = widget_button(secondcomp_buttonbase, value='LARGER WIDTH', uvalue='BROAD_2ND')
widgetids.secondcompsmartlabel  = widget_label(row1, xsize=45, /align_left, value='Smart?')
secondcompsmart_buttonbase      = widget_base(row1, /row, /nonexclusive)
widgetids.secondcompsmartbutton = widget_button(secondcompsmart_buttonbase, value=' ', uval='SMART_2ND')

row2 = widget_base(fitoptionsbase, /row)

if state.plotMonteCarlodistrib eq 0 then str_mcplot  = "OFF" else str_mcplot = "ON"
if state.saveMonteCarlodistrib eq 0 then str_mcsave  = "OFF" else str_mcsave = "ON"
if state.useMonteCarlonoise    eq 0 then str_mcnoise = "OFF" else str_mcnoise = "ON"
if state.scaleNoiseerrors      eq 0 then str_scnoise = "OFF" else str_scnoise = "ON"
case state.linefit_type of 
 0: if ptr_valid(state.gauss_initmap) then str_inputmap_flag = "ON" else str_inputmap_flag = "OFF"
 1: if ptr_valid(state.mom_windowmap) then str_inputmap_flag = "ON" else str_inputmap_flag = "OFF"
endcase

dummy = widget_label(row2, xsize=170, /align_left, value='Montecarlo PDFs plot/save: ')
widgetids.montecarloplot_id = widget_label(row2, value=str_mcplot)
dummy = widget_label(row2, xsize=10, /align_center, value='/')
widgetids.montecarlosave_id = widget_label(row2, value=str_mcsave)
dummy = widget_label(row2, xsize=130, /align_right, value=' Montecarlo noise: ')
widgetids.montecarlonoise_id = widget_label(row2, value=str_mcnoise)
dummy = widget_label(row2, xsize=150, /align_right,value=' Scale Noisecube error: ')
widgetids.scalenoise_id = widget_label(row2, value=str_scnoise)
dummy = widget_label(row2, xsize=150, /align_right,value=' Input Startvals File: ')
widgetids.inputmap_flag = widget_label(row2, value=str_inputmap_flag)

;INSTRUMENTAL RESOLUTION ROWS (2)
instrresbase = widget_base(base5, /col, frame=1)
instrres_row1 = widget_base(instrresbase, /row)
dummy = widget_label(instrres_row1, value='Compute Instrumental Resolution: ')
button = widget_button(instrres_row1, value='FIT TO VARIANCE', uvalue='FITSKY', tooltip='Fit to lines in variance spectrum')
findpro, 'kubeviz', dirlist=kubeviz_dir, /noprint 
case 1 of
 (state.instr eq 'kmos') or (state.instr eq 'kmos3d') : button = widget_button(instrres_row1, value='USE POLYNOMIAL', uvalue='POLYSKY', $
                              tooltip='Use value from polynomial')
 state.instr eq 'sinfoni' : button = widget_button(instrres_row1, value='USE TEMPLATES', uvalue='TPLSKY', $
                              tooltip='Fit to instrument templates', $ 
			      sensitive=file_test(kubeviz_dir+'templates/sinfoni/OH_GaussProf*.fits'))
 else:       button = widget_button(instrres_row1, value='USE TEMPLATES', uvalue='TPLSKY', $
                              tooltip='Fit to instrument templates', $
			      sensitive=0)			      
endcase
dummy = widget_label(instrres_row1, /align_right, xsize=125, value='R (at main line):')
widgetids.instrres_id = widget_label(instrres_row1, xsize=60, /align_left, value=kubeviz_str(kubeviz_getinstrres(), format='(f10.2)'))

case state.instrres_mode of
    0: str_instrres_mode = '(Polynomial fit to cube variance)'
    1: str_instrres_mode = '(Polynomial fit to external arcs)'
    2: str_instrres_mode = '(Spectral Templates)'
endcase
widgetids.instrres_mode_id = widget_label(instrres_row1, /align_left, xsize=250, value=str_instrres_mode)
			      
instrres_row2 = widget_base(instrresbase, /row)
dummy = widget_label(instrres_row2, value='Polynomial coefficients: ')
case state.instrres_mode of
 0: for i=0, state.max_polycoeff_instrres do widgetids.polycoeffres_id[i] = widget_text(instrres_row2, /editable, xsize=14, ysize=1, $
                                             uvalue='POLYCOEFFPAR'+kubeviz_str(i,format='(I1)'), value=kubeviz_str(state.instrres_varpoly[i]))
 1: for i=0, state.max_polycoeff_instrres do widgetids.polycoeffres_id[i] = widget_text(instrres_row2, /editable, xsize=14, ysize=1, $
                                             uvalue='POLYCOEFFPAR'+kubeviz_str(i,format='(I1)'), value=kubeviz_str(state.instrres_extpoly[i]))
 2: for i=0, state.max_polycoeff_instrres do widgetids.polycoeffres_id[i] = widget_text(instrres_row2, /editable, xsize=14, ysize=1, $
                                             uvalue='POLYCOEFFPAR'+kubeviz_str(i,format='(I1)'), value=kubeviz_str(state.instrres_varpoly[i]))					     
endcase

optionbase = widget_base(base5, /col, frame=1, ypad=1)
buttonsrow1 = widget_base(optionbase, /row, ypad=0)
widgetids.fitoption_id[0] = widget_button(buttonsrow1, value='     FIT     ', uvalue='FIT', $
                        tooltip='Fit lines')
widgetids.fitoption_id[1] = widget_button(buttonsrow1, value='   FIT ADJ   ', uvalue='FITADJ', $
                        tooltip='Fit lines')
widgetids.fitoption_id[2] = widget_button(buttonsrow1, value='   RESET FIT   ', uvalue='RESETFIT', $
                        tooltip='Reset fit results')

dummy                     = widget_label(buttonsrow1, xsize=40, /align_left, value='')
widgetids.fitoption_id[3] = widget_button(buttonsrow1, value='  AUTO FLAG  ', uvalue='FLAGALL', $
                        tooltip='Auto-Flag all masks / spaxels')
dummy                     = widget_label(buttonsrow1, xsize=40, /align_left, value='')
widgetids.fitoption_id[4] = widget_button(buttonsrow1, value='  RESET ALL PARS  ', uvalue='RESETALL', $
                        tooltip='Reset all parameters (fit, user, settings)')
widgetids.fitoption_id[5] = widget_button(buttonsrow1, value='  RESET USER PARS  ', uvalue='RESETUSER', $
                        tooltip='Reset user defined constraints')

buttonsrow2 = widget_base(optionbase, /row, ypad=0)
widgetids.fitoption_id[10] = widget_button(buttonsrow2, value='   FIT ALL   ', uvalue='FITALL', $
                        tooltip='Fit all masks / spaxels')
widgetids.fitoption_id[11] = widget_button(buttonsrow2, value=' FIT ADJ ALL ', uvalue='FITADJALL', $
                        tooltip='Fit bad spaxels using initial guess from adjacent spaxels')
widgetids.fitoption_id[12] = widget_button(buttonsrow2, value=' RESET FIT ALL ', uvalue='RESETFITALL', $
                        tooltip='Reset fit results')
			
dummy                     = widget_label(buttonsrow2, xsize=40, /align_left, value='')
widgetids.fitoption_id[13] = widget_button(buttonsrow2, value=' FLAG ON/OFF ', uvalue='RESIMAMASK', $
                        tooltip='Switch masking of results on/off')
dummy                     = widget_label(buttonsrow2, xsize=40, /align_left, value='')
widgetids.fitoption_id[14] = widget_button(buttonsrow2, value='   SAVE   ', uvalue='SAVE', $
                        tooltip='Save fit results in FITS format')
widgetids.fitoption_id[15] = widget_button(buttonsrow2, value='  MASK/SPAXEL  ', uvalue='MODE', $
                        tooltip='Switch MASK/SPAXEL fitting modes')
widgetids.fitoption_id[16] = widget_button(buttonsrow2, value=' GAUSS/MOMENTS ', uvalue='TYPE', $
                        tooltip='Switch GAUSSIAN/MOMENTS fitting types')

; Realize the window
widget_control, base5, /realize

end

;-----------------------------------------------------------------------
pro kubeviz_spax_event, ev

; handle events generated by the spaxel viewer window
common kubeviz_state
common kubeviz_widgetids
common kubeviz_winsize

update = 0   ; reset update flag 
nomove = 0

widget_control, ev.id, get_uvalue=value
if (n_elements(value) eq 0) then value = ''

name = strmid(tag_names(ev, /structure_name), 7, 4)
case (name) of

    "KILL": kubeviz_destroy

    "DRAW": begin
        	
	; ev.type=0  is button press, ev.type=1 is button release
        ; ev.press=1 is left button, ev.press=4  is right button.
        if (ev.type eq 0 and ev.press   eq 1 and state.drag_spax eq 0) then state.drag_spax=1
        if (ev.type eq 0 and ev.press   eq 4 and state.drag_spax eq 0) then state.drag_spax=2
	if (ev.type eq 1 and ev.release eq 4 and state.drag_spax ne 0) then state.drag_spax=0
	if (ev.type eq 1 and ev.release eq 1 and state.drag_spax ne 0) then begin 
	  state.drag_spax=0
          if state.cursormode eq 2 or state.cursormode eq 3 then begin
	    kubeviz_medianspec
	    kubeviz_linefit_update
    	    kubeviz_plotspec
    	    kubeviz_plotspeczoom
	  endif  
 	endif  

        if state.drag_spax eq 2 then begin  ; change contrast/brighness levels     
            update = 1

            x0 = -1
            if ev.X ge 0 and ev.X lt winsize.xwin1 then $
              x0 = float(ev.X) / float(winsize.xwin1)

            y0 = -1
            if ev.Y ge 0 and ev.Y lt winsize.ywin1 then $
              y0 = float(ev.Y) / float(winsize.ywin1)

            if (x0 ne -1 and y0 ne -1) then begin
                ; from ATV.pro - needs improvement
                ; try something like DF's windowimage
                ;print, x0, y0

                brightness = 0.5
                contrast = 0.5
                x = brightness * 255 * x0
                y = contrast * 255 * y0 > 2

                high = x+y & low = x-y
                diff = (high-low) > 1
                slope = float(255) / diff
                intercept = -slope * low
                ;print, x, y, high, low, diff, slope, intercept
                p = long(findgen(255)*slope+intercept)
                kubeviz_modify_colour, ci=p
            endif
        endif

        if state.drag_spax eq 1 then begin  ; update cursor position
	    
	    drawXsize = float(state.zoomfac*state.Ncol)
	    drawYsize = float(state.zoomfac*state.Nrow)
	    
	    state.col = ev.X / (drawXsize/(state.spaxrange[1]-state.spaxrange[0]+1) ) + state.spaxrange[0] 
            state.row = ev.Y / (drawYsize/(state.spaxrange[3]-state.spaxrange[2]+1) ) + state.spaxrange[2] 
	    	    
	    case 1 of
	      (state.spaxrange[0] eq 0): $
	          if ev.X lt drawXsize*0.85 and ev.Y gt drawYsize*0.15 and ev.Y lt drawYsize*0.85 then nomove=1
	      (state.spaxrange[1] eq state.Ncol-1): $
	          if ev.X gt drawXsize*0.15 and ev.Y gt drawYsize*0.15 and ev.Y lt drawYsize*0.85 then nomove=1
	      (state.spaxrange[2] eq 0): $ 
	          if ev.X gt drawXsize*0.15 and ev.X lt drawXsize*0.85 and ev.Y lt drawYsize*0.85 then nomove=1
	      (state.spaxrange[3] eq state.Nrow-1): $
	          if ev.X gt drawXsize*0.15 and ev.X lt drawXsize*0.85 and ev.Y gt drawYsize*0.15 then nomove=1
	      else: $
	          if ev.X gt drawXsize*0.15 and ev.X lt drawXsize*0.85 and ev.Y gt drawYsize*0.15 and ev.Y lt drawYsize*0.85 then nomove=1
	    endcase
	   	    	    	    
            if state.col le 0 then state.col = 0
            if state.col ge state.Ncol then state.col = state.Ncol-1
            if state.row le 0 then state.row = 0
            if state.row ge state.Nrow then state.row = state.Nrow-1
	    
            if state.cursormode eq 2 then begin
	       case state.maskmode of
	         0: (*state.spaxselect)[state.col, state.row, state.imask-1] = 1
		 1: begin
		    tempmask = (*state.spaxselect)[state.col-state.maskradius-1:state.col+state.maskradius+1,$
		                                   state.row-state.maskradius-1:state.row+state.maskradius+1, state.imask-1]
		    dist_ellipse, distima, [(size(tempmask))[1:2]], ((size(tempmask))[1]-1)/2, ((size(tempmask))[1]-1)/2, 1, 0
		    okspax = where(distima le state.maskradius, Nspax)
		    if Nspax gt 0 then begin
		      tempmask[okspax] = 1
		      (*state.spaxselect)[state.col-state.maskradius-1:state.col+state.maskradius+1,$
		                          state.row-state.maskradius-1:state.row+state.maskradius+1, state.imask-1] = tempmask
		    endif		    
		 end
		 2: begin
		    halfsize = floor(state.maskradius/2.)
		    if state.col-halfsize ge 0 then x0 = state.col-halfsize else x0 = 0
		    if state.col+halfsize lt state.Ncol then x1 = state.col+halfsize else x1 = state.Ncol-1
		    if state.row-halfsize ge 0 then y0 = state.row-halfsize else y0 = 0
		    if state.row+halfsize lt state.Nrow then y1 = state.row+halfsize else y1 = state.Nrow-1
		    if halfsize mod 2 eq 0 then x1 -= 1 & y1 -= 1
		    (*state.spaxselect)[x0:x1,y0:y1, state.imask-1] = 1
		 end
	       endcase
            endif
	    
	    if state.cursormode eq 3 then begin
	       case state.maskmode of
	         0: (*state.spaxselect)[state.col, state.row, state.imask-1] = 0
		 1: begin
		    tempmask = (*state.spaxselect)[state.col-state.maskradius-1:state.col+state.maskradius+1,$
		                                   state.row-state.maskradius-1:state.row+state.maskradius+1, state.imask-1]
		    dist_ellipse, distima, [(size(tempmask))[1:2]], ((size(tempmask))[1]-1)/2, ((size(tempmask))[1]-1)/2, 1, 0
		    okspax = where(distima le state.maskradius, Nspax)
		    if Nspax gt 0 then begin
		      tempmask[okspax] = 0
		      (*state.spaxselect)[state.col-state.maskradius-1:state.col+state.maskradius+1,$
		                          state.row-state.maskradius-1:state.row+state.maskradius+1, state.imask-1] = tempmask
		    endif		    
		 end
		 2: begin
		    halfsize = floor(state.maskradius/2.)
		    if halfsize mod 2 eq 0 $
		     then (*state.spaxselect)[state.col-halfsize:state.col+halfsize-1, state.row-halfsize:state.row+halfsize-1, state.imask-1] = 0 $
		     else (*state.spaxselect)[state.col-halfsize:state.col+halfsize, state.row-halfsize:state.row+halfsize, state.imask-1] = 0 
		 end
	       endcase
	    endif
	    
            if state.cursormode lt 2 then update = 2 else begin
	     update = 0
	     kubeviz_plotspax, /fast_update, nomove=nomove
             kubeviz_plotinfo
            endelse
        endif

        ; here we deal with keyboard press. update only on key press not on key release
        if ev.press eq 1 then kubeviz_keyboard_handler, ev.key, ev.ch, update, nomove

    end  ; end draw

    "SLID": begin
        update = 1
        state.wpix = ev.value
    end
    
    "BUTT": begin
        ;print, '[DEBUG] button value: ', value
        update=1
        
	if (value eq "MORESPATZOOM") then begin
            update=1
	    if state.zoomspax ge 1 then state.zoomspax += 1
	    widget_control, widgetids.spatzoom_id, set_value='Zoom: '+string(state.zoomspax, format='(i2)')

        endif
	
	if (value eq "LESSSPATZOOM") then begin
            update=1
	    if state.zoomspax gt 1 then state.zoomspax -= 1
	    widget_control, widgetids.spatzoom_id, set_value='Zoom: '+string(state.zoomspax, format='(i2)')
        endif
	
        if (value eq "Quit") then begin
            update=0
            kubeviz_destroy
            return
        endif

        if (value eq "Open") then begin
            
            smooth_orig = state.smooth
            specsmooth_orig = state.specsmooth
            lineset_orig = state.selected_lineset
	    ctab_orig = state.ctab
	    zcuts_orig = state.zcuts
            
            ; get new cube 
            kubeviz_selectcube, datafile, dir, noisefile, noisedir
            if n_elements(datafile) eq 0 then break      ; no file selected 

            ;save active session
            kubeviz_savesession, fname=state.cwdir+'lastsession.sav'
           
            kubeviz_destroy                           ;kill previous session
            kubeviz_common                            ;reset public variables
            
            kubeviz_setcolour, ctab_orig              
            state.smooth = smooth_orig
            state.specsmooth = specsmooth_orig
            state.selected_lineset = lineset_orig
            state.filename = datafile
            state.indir = dir
            state.zcuts = zcuts_orig       
            state.zoomrange = 32            ; start half-range of spectral zoom window
                 
            ;Create the file where to store fit err msg
            errmsg_file_info = file_info(state.cwdir+'fit_errmsg_'+strmid(state.filename,0,strlen(state.filename)-5)+'.txt')
            if errmsg_file_info.exists eq 0 then spawn, 'touch '+state.cwdir+'fit_errmsg_'+strmid(state.filename,0,strlen(state.filename)-5)+'.txt'
                    
            kubeviz_getdata, dir, datafile, noisedir=noisedir, noisefile=noisefile       ; load selected file
            kubeviz_linefit_init
            kubeviz_create, fname
           
        endif

        if (value eq "Print") then begin
            kubeviz_plotspax, /ps
            printf, state.log_lun, '[KUBEVIZ] Created postscript file: '+strmid(state.filename,0,strlen(state.filename)-5)+'.ps'
        endif
	
        if (value eq "SaveCube") then begin
            update=0
            kubeviz_savecube
        endif
	
	if (value eq "SaveImage") then begin
            update=0
            kubeviz_saveimage
        endif
	
        if (value eq "DisplayHeader") then begin
            update=0
            kubeviz_selecthead
        endif

        if (value eq "SaveSession") then begin
            update=0
            kubeviz_savesession
        endif

        if (value eq "LoadSession") then begin
            kubeviz_loadsession
        endif       

        if (value eq "Data") then begin
            state.cubesel=0
        endif

        if (value eq "Noise") then begin
            state.cubesel=1
        endif      

        if (value eq "BadPixels") then begin
            state.cubesel=2
        endif 
        
        if (value eq "SN") then begin
            state.cubesel=3
        endif      

        if (value eq "Linefit") then begin
	    state.cubesel=4
            if state.par_imagebutton eq '' then begin
                state.par_imagebutton='FLAG'
                state.lineresimg = ptr_new((*state.sp_nrescube)[*,*,0])
            endif
        endif
        
        if (value eq "LineErrors") then begin
            state.cubesel=5
            if state.par_imagebutton eq '' then begin
                state.par_imagebutton='FLAG'
                state.lineerrresimg = ptr_new((*state.sp_nerrrescube)[*,*,0])
            endif
        endif
        
        if (value eq "LineSN") then begin
            state.cubesel=6
            if state.par_imagebutton eq '' then begin
                state.par_imagebutton='FLAG'
                state.lineresimg    = ptr_new((*state.sp_nrescube)[*,*,0])
                state.lineerrresimg = ptr_new(replicate(1.D,state.Ncol, state.Nrow))
            endif
        endif  

        if value eq "HelpWhatIsNew"    then kubeviz_help_whatisnew
        if value eq "HelpInstructions" then kubeviz_help_instructions
        if value eq "HelpShortcuts"    then kubeviz_help_shortcuts


        if value eq "Grey"         then kubeviz_setcolour, 0
        if value eq "BlueWhite"    then kubeviz_setcolour, 1
        if value eq "Heat"         then kubeviz_setcolour, 3
        if value eq "STD GAMMA-II" then kubeviz_setcolour, 5
        if value eq "Rainbow"      then kubeviz_setcolour, 13

        if value eq "Invert" then begin 
            state.invert = -1 * state.invert 
            kubeviz_modify_colour
        endif

        if value eq "MinMax"   then state.zcuts  = 2
        if value eq "Zscale"   then state.zcuts  = 3
        if value eq "HistEq"   then state.zcuts  = 4
        if value eq "99.5"     then state.zcuts  = 11
	if value eq "99.0"     then state.zcuts  = 12
        if value eq "97.0"     then state.zcuts  = 13
        if value eq "95.0"     then state.zcuts  = 14
        if value eq "UserLin"  then state.zcuts  = 1
        if value eq "UserSqrt" then state.zcuts  = 5
	if value eq "UserLog"  then state.zcuts  = 6
	if value eq "UserPars" then begin 
	   kubeviz_zcut_parameters
	   update=0
	endif   

        if value eq "None"      then begin
	   update=2
	   state.cursormode = 0
	endif   
        if value eq "Crosshair" then begin
           update = 2
	   state.cursormode = 1
           state.linefit_mode = 0
           state.specmode = 0
        endif   
        if value eq "Select"    then begin 
           update = 2
	   state.cursormode = 2
           state.linefit_mode = 1
           widget_control, widgetids.stepmask_id, map=1
        endif   
	if value eq "Deselect"    then begin
	   update=2
	   state.cursormode = 3
	endif   
        if value eq "Clear"       then begin
	   update=2
	   kubeviz_clear_select
	endif   
        if value eq "OptimalMask" then begin
	   update=2
	   kubeviz_optimal_mask
	endif   
        if value eq "SAVE"        then begin
	   update=0
	   kubeviz_save_mask
	endif   
        if value eq "Load"        then begin
	   update=2
	   kubeviz_load_mask
	endif  
        if value eq "NewMask"     then begin
	   update=2
	   kubeviz_newmask
	endif   
        if value eq "DeleteMask"  then begin
	   update=2
	   kubeviz_deletemask
        endif
	if value eq "GoToMask"  then begin
	   update=2
	   kubeviz_gotomask
        endif
	if value eq "PrevMask"  then begin
	   state.imask = state.imask - 1
           if state.imask le 0 then state.imask = 1
	   update=2
        endif
	if value eq "NextMask"  then begin
	   state.imask = state.imask + 1
	   if state.imask gt state.Nmask then state.imask = state.Nmask
	   update=2
        endif
	if value eq "MaskPars"  then begin
	   update=2
	   kubeviz_maskpars
        endif
	
        if value eq "Slice" then begin
            state.imgmode = 0
            kubeviz_spec_reset_range
        endif

        if value eq "Sum1" then begin 
            state.imgmode = 1
            update=3
        endif

        if value eq "Median1" then begin 
            state.imgmode = 2
            update=3
        endif

        if value eq "WeightedAvg1" then begin 
            state.imgmode = 3
            update=3
        endif

        if value eq "WeightedMed1" then begin 
            state.imgmode = 4
            update=3
        endif
        
        if value eq "MedSub1" then begin 
            state.imgmode = 5
            update=3
        endif

        if value eq "Sum2" then begin 
            state.imgmode = 6
            update=3
        endif

        if value eq "Median2" then begin 
            state.imgmode = 7
            update=3
        endif

        if value eq "WeightedAvg2" then begin 
            state.imgmode = 8
            update=3
        endif

        if value eq "WeightedMed2" then begin 
            state.imgmode = 9
            update=3
        endif
        
        if value eq "MedSub2" then begin  
            state.imgmode = 10
            update=3
        endif

        if value eq "Med2-Med1" then begin
            state.imgmode = 11
            kubeviz_medsum_image_update, imgmode = 11
	    kubeviz_medsum_image_update, imgmode = 12
	    update=1
        endif
	
	if value eq "Med1-Med2" then begin
            state.imgmode = 12
            kubeviz_medsum_image_update, imgmode = 11
	    kubeviz_medsum_image_update, imgmode = 12
	    update=1
        endif
        	
	if value eq "SmoothPars" then begin
	     kubeviz_smooth_parameters
	     update=0
	endif
	
	if value eq "FitallRange" then begin
	     kubeviz_fitall_range
	     update=0
	endif     
	
        if value eq "BootstrapErrors" then begin
            if state.domontecarlo ne 1 then begin
                if state.Nbootstrap eq 0 then begin
                    fname = dialog_pickfile(filter='*.fits', /read, /must_exist, get_path=path)
                    len = strlen(fname) - strlen(path)
                    if len gt 0 then begin
                        state.domontecarlo=1
                        kubeviz_readbootstrapcubes, fname
                        kubeviz_setupmontecarlocubes
                        if state.cubesel gt 3 then kubeviz_linefit_image_reset
                    endif else printf, state.log_lun, '[WARNING] No bootstrap file selected'
                endif else begin
                        state.domontecarlo=1
                        state.montecarlocubes = state.bootstrapcubes
                        state.Nmontecarlo = state.Nbootstrap
                        if state.useMonteCarlonoise eq 1 then state.noise = state.bootstrapnoise
                        if (state.cubesel gt 3 and max(abs((*state.boot_sp_nrescube)[*,*,1])) eq 0.) then kubeviz_linefit_image_reset
                endelse
                update = 3
            endif
        endif
        
        if value eq "Mc1Errors" then begin
            if state.domontecarlo ne 2 then begin
                if state.Nmc1 eq 0 then begin
                    state.domontecarlo=2
                    if state.Nmc2 gt 0 then state.Nmc1=state.Nmc2 else state.Nmc1 = 100
                    kubeviz_createmc1cubes
                    kubeviz_setupmontecarlocubes
                    if state.cubesel gt 3 then kubeviz_linefit_image_reset
                endif else begin
                    state.domontecarlo=2
                    state.montecarlocubes = state.mc1cubes
                    state.Nmontecarlo = state.Nmc1
                    if state.useMonteCarlonoise eq 1 then state.noise = state.mc1noise
                    if (state.cubesel gt 3 and max(abs((*state.mc1_sp_nrescube)[*,*,1])) eq 0.) then kubeviz_linefit_image_reset
                endelse 
                update = 3
             endif
        endif
        
        if value eq "Mc2Errors" then begin
            if state.domontecarlo ne 3 then begin
                if state.Nmc2 eq 0 then begin
                    state.domontecarlo=3
                    if state.Nmc1 gt 0 then state.Nmc2=state.Nmc1 else state.Nmc2 = 100
                    kubeviz_createmc2cubes
                    kubeviz_setupmontecarlocubes
                    if state.cubesel gt 3 then kubeviz_linefit_image_reset
                endif else begin
                    state.domontecarlo=3
                    state.montecarlocubes = state.mc2cubes
                    state.Nmontecarlo = state.Nmc2
                    if state.useMonteCarlonoise eq 1 then state.noise = state.mc2noise
                    if (state.cubesel gt 3 and max(abs((*state.mc2_sp_nrescube)[*,*,1])) eq 0.) then kubeviz_linefit_image_reset
                endelse 
                update = 3
            endif
        endif
        
        if value eq "Mc3Errors" then begin
            if state.domontecarlo ne 4 then begin
                if state.Nmc3 eq 0 then begin
                    if state.redshift gt 0 then begin
                      state.domontecarlo=4
                      if state.Nmc1 gt 0 then state.Nmc3=state.Nmc1 else begin
                      if state.Nmc2 gt 0 then state.Nmc3=state.Nmc2 else state.Nmc3 = 100 
                      endelse
                      kubeviz_createmc3cubes
                      kubeviz_setupmontecarlocubes
                    if state.cubesel gt 3 then kubeviz_linefit_image_reset
                    endif else printf, state.log_lun, '[WARNING] Set the correct redshift before calling this method.'
                endif else begin
                    state.domontecarlo=4
                    state.montecarlocubes = state.mc3cubes
                    state.Nmontecarlo = state.Nmc3
                    if state.useMonteCarlonoise eq 1 then state.noise = state.mc3noise
                    if (state.cubesel gt 3 and max(abs((*state.mc3_sp_nrescube)[*,*,1])) eq 0.) then kubeviz_linefit_image_reset
                endelse 
                update = 3
            endif
        endif
        
        if value eq "NoiseErrors" then begin
            if state.domontecarlo gt 0 then begin
             state.domontecarlo=0
             state.noise = state.noisecube
             if (state.cubesel gt 3 and max(abs((*state.noi_sp_nrescube)[*,*,1])) eq 0.) then kubeviz_linefit_image_reset
             update = 3 
            endif 
        endif        
        
	if value eq "MonteCarloNoise" then begin
            ; switch bootstrap noise cubes on/off
            state.useMonteCarlonoise = 1-state.useMonteCarlonoise
            case state.useMonteCarlonoise of
                0: begin
		  printf, state.log_lun, "[KUBEVIZ] Montecarlo Noise switched OFF"
		  state.noise = state.noisecube
                end
		1: begin 
		  printf, state.log_lun, "[KUBEVIZ] Montecarlo Noise switched ON"
		  case state.domontecarlo of
		   0: state.noise = state.noisecube
		   1: state.noise = state.bootstrapnoise
		   2: state.noise = state.mc1noise
		   3: state.noise = state.mc2noise
		   4: state.noise = state.mc3noise
		  endcase
                end
	    endcase
            update = 3
        endif
	
        if value eq "MonteCarloPlot" then begin
            ; switch bootstrap plotting on or off
            state.plotMonteCarlodistrib = 1-state.plotMonteCarlodistrib
            case state.plotMonteCarlodistrib of
                0: printf, state.log_lun, "[KUBEVIZ] Montecarlo plotting switched OFF"
                1: printf, state.log_lun, "[KUBEVIZ] Montecarlo plotting switched ON"
            endcase
            update = 1
        endif
	
	if value eq "MonteCarloSave" then begin
            ; switch bootstrap plotting on or off
            state.saveMonteCarlodistrib = 1-state.saveMonteCarlodistrib
            case state.saveMonteCarlodistrib of
                0: printf, state.log_lun, "[KUBEVIZ] Montecarlo PDF saving switched OFF"
                1: printf, state.log_lun, "[KUBEVIZ] Montecarlo PDF saving switched ON"
            endcase
            update = 1
        endif
	
	if value eq "NoiseCubeErrScale" then begin
            ; switch the scaling of noisecube errors based on the chi-sq on or off
           state.scaleNoiseerrors = 1-state.scaleNoiseerrors
	   case state.scaleNoiseerrors of
	       0: printf, state.log_lun, "[KUBEVIZ] Noisecube errors scaling switched OFF"
	       1: printf, state.log_lun, "[KUBEVIZ] Noisecube errors scaling switched ON"
	   endcase
	   update = 1
        endif
        
        if value eq "LoadResultFile" then begin
            ; load an external file where the fit results are saved
            ffile = dialog_pickfile(filter='*.fits', path=state.cwdir, /read, /must_exist, get_path=path)
            kubeviz_splitpath, ffile, fname, path
            if strlen(fname) gt 0 then kubeviz_linefit_loadres, fname, dir=path else printf, state.log_lun, '[WARNING] No result file selected'
	    kubeviz_linefit_typeswitch
        endif

    end

endcase
if update eq 3 then begin
    update = 1
    if state.imgmode  gt 0 then kubeviz_medsum_image_update
    if state.specmode gt 0 then kubeviz_medianspec
endif    
if update eq 2 then begin
    kubeviz_linefit_update
    kubeviz_plotspax, /fast_update , nomove=nomove
    kubeviz_plotinfo
    kubeviz_plotspec
    kubeviz_plotspeczoom
endif
if update eq 1 then begin
    kubeviz_linefit_update
    kubeviz_plotspax, nomove=nomove
    kubeviz_plotinfo
    kubeviz_plotspec
    kubeviz_plotspeczoom
endif

end

;-----------------------------------------------------------------------
pro kubeviz_linefit_image_reset

common kubeviz_state

state.par_imagebutton = ''
state.cubesel = 0
printf, state.log_lun, '[WARNING] No valid results in this plane for this error method as yet.'
printf, state.log_lun, '[KUBEVIZ] Reverting to data cube'

end

;-----------------------------------------------------------------------
pro kubeviz_spec_event, ev
; handle events generated by the spectrum window
common kubeviz_state
common kubeviz_widgetids
common kubeviz_winsize

update = 0 ; reset update flag

widget_control, ev.id, get_uvalue=value
if (n_elements(value) eq 0) then value = ''

name = strmid(tag_names(ev, /structure_name), 7, 4)
case (name) of
    
    "BASE" : begin ; resize the spec window 
        ; get the new window size
        widget_control, ev.top, tlb_get_size=resize
        ; update winsize-state vector 
        if resize[0] lt 500 then resize[0] = 500
	if resize[1] lt 250 then resize[1] = 250
	winsize.xwin2 = resize[0]
        winsize.ywin2 = resize[1]
        ; update the draw widgetsize
        widget_control, widgetids.draw2, draw_xsize=winsize.xwin2, draw_ysize=winsize.ywin2
        ; update the plot
        kubeviz_plotspec
    end

    "BUTT": begin

        if value eq "SCALE"   then state.scale  = -1 * state.scale
        if value eq "Spaxel"  then begin
            state.specmode = 0
            state.linefit_mode = 0
        endif
        if value eq "Sum"  then begin
            state.specmode = 1
            state.linefit_mode = 1
	    kubeviz_medianspec
        endif
        if value eq "Median"  then begin
            state.specmode = 2
            state.linefit_mode = 1
	    kubeviz_medianspec
        endif
        if value eq "W. Avg"  then begin
            state.specmode = 3
            state.linefit_mode = 1
	    kubeviz_medianspec
        endif
	if value eq "MedSub" then begin
            state.specmode = 4
            state.linefit_mode = 0
	    kubeviz_medianspec
        endif
        if value eq "Optimal"  then begin
            state.specmode = 5
            state.linefit_mode = 1
	    kubeviz_medianspec
        endif

        if value eq "Range1"  then state.wavsel = 1
        if value eq "Range2"  then state.wavsel = 2
        if value eq "Reset"   then kubeviz_spec_reset_range

        ;if value eq "MARKER"  then state.marker = -1 * state.marker

        if value eq "ZOOM" then begin
            if state.zoommap eq 0 then state.zoommap=1 else state.zoommap=0
            widget_control, widgetids.base3, Map=state.zoommap
        endif

        if value eq "SAVESPEC" then kubeviz_savespec

        update = 1

    end

    "TEXT": begin
        update = 1
        if (value eq "ZMINSP") then begin
            widget_control, ev.id, get_value=value
            state.zmin_spec = value
            ; signal return-press by moving cursor to the left
            widget_control, ev.id, set_value=strtrim(value, 2)
        endif
        if (value eq "ZMAXSP") then begin 
            widget_control, ev.id, get_value=value
            state.zmax_spec = value
            widget_control, ev.id, set_value=strtrim(value, 2)
        endif
    end

    
    "DRAW": begin
        ; ev.press=1 is left button, ev.press=4  is right button
        if (ev.type eq 0 and ev.press eq 1 and state.drag_spec eq 0) then state.drag_spec=1
        if (ev.type eq 0 and ev.press eq 4 and state.drag_spec eq 0) then state.drag_spec=2
        if (ev.type eq 1 and state.drag_spec ne 0) then state.drag_spec=0

        if state.drag_spec eq 1 then begin     
            ; positional numbers '40' and 'winsize.xwin2-30 are hard coded in kubeviz_plotspec 
	    ; if you change here then change there as well
            frac = (ev.X - 40.) / (winsize.xwin2-30 - 40.)
            pix = fix(frac*state.Nwpix)
            if pix ge 0 and  pix le state.Nwpix-1 then begin
                state.wpix = pix
                widget_control, widgetids.slid_id, set_value=pix
                if state.imgmode eq 0 or state.imgmode eq 5 or state.imgmode eq 10 then update = 1 else update=2
            endif
        endif

        if ev.press eq 1 then begin
            kubeviz_keyboard_handler, ev.key, ev.ch, update
	    
            case string(ev.ch) of
              "s": begin     ; select wavelength range
                    if state.wavsel eq 0 then state.wavsel=1 ; default to range1
                    if state.npress eq 0 then begin
                        if state.wavsel eq 1 then begin
                            state.wavrange1[0] = state.wpix
                            state.wavrange1[1] = state.wpix
                        endif
                        if state.wavsel eq 2 then begin
                            state.wavrange2[0] = state.wpix
                            state.wavrange2[1] = state.wpix
                        endif
                        state.npress = 1
                        update = 1
                    endif else begin
                        if state.wavsel  eq 1 then state.wavrange1[1] = state.wpix
                        if state.wavsel  eq 2 then state.wavrange2[1] = state.wpix
                        if state.imgmode gt 0 then kubeviz_medsum_image_update
                        state.npress = 0
                        update = 1
                    endelse
                end
               "z": begin
                    new_redshift = ( (*state.wave)[state.wpix] / kubeviz_getmainline(redshift=0) ) - 1
		    kubeviz_change_redshift, new_redshift  
                end
                else:
            endcase
        endif
    end
    else:
endcase

if update eq 2 then begin
    kubeviz_plotspax, /fast_update
    kubeviz_plotinfo
    kubeviz_plotspec
    kubeviz_plotspeczoom
    kubeviz_linefit_update
endif

if update eq 1 then begin
    kubeviz_plotspax
    kubeviz_plotinfo
    kubeviz_plotspec
    kubeviz_plotspeczoom
    kubeviz_linefit_update
endif


end


;-----------------------------------------------------------------------
pro kubeviz_speczoom_event, ev
; handle events generated by the zoom spectrum window
common kubeviz_state
common kubeviz_widgetids
common kubeviz_winsize

update = 0 ; reset update flag

widget_control, ev.id, get_uvalue=value
if (n_elements(value) eq 0) then value = ''

name = strmid(tag_names(ev, /structure_name), 7, 4)
case (name) of

    "BASE" : begin ; resize the speczoom window 
        ; get the new window size
        widget_control, ev.top, tlb_get_size=resize
        ; update winsize-state vector 
        if resize[0] lt 250 then resize[0] = 250
	if resize[1] lt 125 then resize[1] = 125
	winsize.xwin3 = resize[0]
        winsize.ywin3 = resize[1]
        ; update the draw widgetsize
        widget_control, widgetids.draw3, draw_xsize=winsize.xwin3, draw_ysize=winsize.ywin3
        ; update the plot
        kubeviz_plotspeczoom
    end

    "KILL": begin  ; close window
        state.zoommap=0
        widget_Control, widgetids.base3, Map=state.zoommap
        return
    end

    "BUTT": begin
        if value eq "ZOOMOUT" or value eq "ZOOMIN" then begin ; zoom

            if value eq "ZOOMOUT" then state.zoomrange = state.zoomrange * 2
            if value eq "ZOOMIN" then  state.zoomrange = state.zoomrange / 2
            
            npix = state.Nwpix / 2
            if state.zoomrange gt npix then state.zoomrange = npix
            if state.zoomrange lt 1 then state.zoomrange =  1
            
            widget_control, widgetids.zoom_id, $
                            set_value=string(2*state.zoomrange,format='(i4)')+' pixel'
                                ;set_value=string(2*state.zoomrange)+' pixel'
            kubeviz_plotspeczoom
        endif
        
    end

    "DRAW": begin
        if ev.press eq 1 then begin

            if ev.type eq 0 then begin
                ; 40 and 10 in xleft and xright respectively are hardcoded in kubeviz_plotzoomspec 
                xleft  = 40.
                xright = winsize.xwin3 - 10.
                frac = (ev.X - xleft) / (xright - xleft)
                z0 = (*state.wave)[state.wpix]
                x1 = state.wpix - state.zoomrange
                x2 = state.wpix + state.zoomrange
                if x1 lt 0 then x1 = 0
                if x2 gt state.Nwpix-1 then x2 = state.Nwpix-1
                pix = fix(x1+frac*(x2-x1))
                if pix ge 0 and  pix le state.Nwpix-1 then begin
                    widget_control, widgetids.slid_id, set_value=pix
                    state.wpix = pix
                    if state.imgmode eq 0 or state.imgmode eq 5 or state.imgmode eq 10 then update = 1 else update=2 
		endif      
            endif

            kubeviz_keyboard_handler, ev.key, ev.ch, update

            case string(ev.ch) of

                "s": begin     ; select wavelength range
                    if state.wavsel eq 0 then state.wavsel=1 ; default to range1
                    if state.npress eq 0 then begin
                        pix = state.wpix
                        if state.wavsel eq 1 then begin
                            state.wavrange1[0] = pix
                            state.wavrange1[1] = pix
                        endif
                        if state.wavsel eq 2 then begin
                            state.wavrange2[0] = pix
                            state.wavrange2[1] = pix
                        endif
                        update = 1
                        state.npress = 1
                    endif else begin
                        pix = state.wpix
                        if state.wavsel eq 1 then state.wavrange1[1] = pix
                        if state.wavsel eq 2 then state.wavrange2[1] = pix
                        if state.imgmode gt 0 then kubeviz_medsum_image_update
                        update = 1
                        state.npress = 0
                    endelse
                end

                "g": begin
                    case state.cubesel of
		      0: begin
		        case state.specmode of
		           0 : spectrum =  (*state.datacube)[state.col,state.row,*]
                           4 : spectrum =  (*state.datacube)[state.col,state.row,*] - (*state.medspec)
                           else : spectrum =  (*state.medspec) ; Includes med, sum, wavg and optimal
			 endcase
		       end
		      1: begin
		        case state.specmode of
		           0 : spectrum =  (*state.noisecube)[state.col,state.row,*]
                           4 : spectrum =  (*state.noisecube)[state.col,state.row,*] - (*state.nmedspec)
                           else : spectrum =  (*state.nmedspec)
			 endcase
		       end
		      else: begin
		        case state.specmode of
		           0 : spectrum =  (*state.datacube)[state.col,state.row,*]
                           4 : spectrum =  (*state.datacube)[state.col,state.row,*] - (*state.medspec)
                           else : spectrum =  (*state.medspec)
			 endcase
		       end 	 	   
                    endcase
		    
                    bin = 40 ;half-range in pixels over which to fit - hard coded 
                    z0 = state.wpix
                    xg = (*state.wave)[z0-bin:z0+bin]
                    yg = fltarr(1+2*bin)
                    yg[*]=spectrum[z0-bin:z0+bin]
                    gfit = mpfitpeak(xg, yg, a, nterms=4, perror=aerr, chisq=chisq)
                    printf, state.log_lun, '[KUBEVIZ]  Best-fit Gauss parameters: amplitude, centroid, width, y-offset:'
                    printf, state.log_lun, '[KUBEVIZ] ', a
		    printf, state.log_lun, '[KUBEVIZ] ', aerr
                    kubeviz_plotspeczoom
                    if state.cubesel gt 1 then printf, state.log_lun, '[WARNING] Fit applied to data, but alternative cube currently selected!'
		    oplot, xg, gfit, color=253, thick=2 
                end
                
                "z": begin
                    new_redshift = ( (*state.wave)[state.wpix] / kubeviz_getmainline(redshift=0) ) - 1
                    kubeviz_change_redshift, new_redshift  
                end

                else:
            endcase
        endif
    end
    else:     ; no matching case found
endcase


if update eq 2 then begin
    kubeviz_plotspax, /fast_update
    kubeviz_plotinfo
    kubeviz_plotspec
    kubeviz_plotspeczoom
endif

if update eq 1 then begin
    kubeviz_plotspax
    kubeviz_plotinfo
    kubeviz_plotspec
    kubeviz_plotspeczoom
endif

end

;-----------------------------------------------
pro kubeviz_linefit_event, ev

common kubeviz_state
common kubeviz_widgetids
common flags, flag_ok, flag_bad

widget_control, ev.id, get_uvalue=value
if (n_elements(value) eq 0) then value = ''

name = strmid(tag_names(ev, /structure_name), 7, 4)

case (name) of

    "TLB_": begin
      if  strmid(tag_names(ev, /structure_name), 11, 4) eq "ICON" then begin
       state.linefitmap = 1-ev.iconified
      endif
    end

    "TEXT": begin 
; record user entered parameter values
        widget_control, ev.id, get_value=parval
        widget_control, ev.id, set_value=strtrim(parval, 2)

        if value eq 'REDSHIFT' then kubeviz_change_redshift, parval
        
        if strmid(value,0,2) eq 'SP' then begin
            linetype = strmid(value,2,1)
            case linetype of
                'N': begin
                    par = fix(strmid(value,3))
                    state.gnfit[par] = parval
                end
                'B':begin
                    par = fix(strmid(value,3))
                    state.gbfit[par] = parval
                end
            endcase
        endif
        if strmid(value,0,4) eq 'MINP' then begin
            linetype = strmid(value,4,1)
            case linetype of
                'N': begin
                    par = fix(strmid(value,5))
                    state.gnlims[par,0] = parval
                end
                'B': begin
                    par = fix(strmid(value,5))
                    state.gblims[par,0] = parval
                end
            endcase
        endif
        if strmid(value,0,4) eq 'MAXP' then begin
            linetype = strmid(value,4,1)
            case linetype of
                'N': begin
                    par = fix(strmid(value,5))
                    state.gnlims[par,1] = parval
                end
                'B': begin
                    par = fix(strmid(value,5))
                    state.gblims[par,1] = parval
                end
            endcase
        endif
        
        if value eq 'LINEFITMAXOFFB' then begin
          if parval ge 0. then state.maxwoffb = parval else printf, state.log_lun, '[WARNING] Range should be > 0!'
        endif  
        if value eq 'LINEFITMAXOFFR' then begin
          if parval ge 0. then state.maxwoffr = parval else printf, state.log_lun, '[WARNING] Range should be > 0!'
        endif
        if value eq 'MASKSNTHRESH' then begin
          if parval ge 0. then state.mask_sn_thresh = parval else printf, state.log_lun, '[WARNING] S/N threshold should be > 0!'
        endif  
        if value eq 'MASKMAXVELERR' then begin
          if parval ge 0. then state.mask_maxvelerr = parval else printf, state.log_lun, '[WARNING] Max velocity error should be > 0!'
        endif
	if value eq 'MOMTHRESH' then begin
          state.mom_thresh = parval 
        endif
        
        if strmid(value,0,4) eq 'CONT' then begin
            contfitpar = strmid(value,4)
            case contfitpar of
                'MINOFF': if parval ge 0. then state.continuumfit_minoff = parval else printf, state.log_lun, '[WARNING] continuum min. offset should be > 0.!'                    
                'MAXOFF': state.continuumfit_maxoff = parval
                'MINPERC': if parval ge 0. then state.continuumfit_minperc = parval else printf, state.log_lun, '[WARNING] continuum min percentile should be >=0!'
                'MAXPERC': if parval le 100. then state.continuumfit_maxperc = parval else printf, state.log_lun, '[WARNING] continuum max percentile should be <=100!'
                'ORDER': state.continuumfit_order = parval
                else: print, '[BUG] No such paramater ',value
            endcase
        endif
        if strmid(value,0,12) eq 'POLYCOEFFPAR' then begin
            coeff_index = fix(strmid(value,12,1))
            case state.instrres_mode of
	      0: state.instrres_varpoly[coeff_index] = parval
	      1: state.instrres_extpoly[coeff_index] = parval
	      else:
	    endcase
            kubeviz_linefit_update
        endif
        
    kubeviz_linefit_update
    end

    "BUTT": begin ; process button events
        case (value) of

            'FLAG': begin
               case state.linefit_type of 
	    	   0: begin
	    	      case state.linefit_mode of
            	   	   0: flag = fix( (*state.sp_nrescube)[state.col, state.row, 0] )
            	   	   1: flag = fix( (*state.nrescube)[state.imask-1, 0] )
            	      endcase
	    	   end
	    	   1: begin
	    	      case state.linefit_mode of
            	   	   0: flag = fix( (*state.sp_mrescube)[state.col, state.row, 5] )
            	   	   1: flag = fix( (*state.mrescube)[state.imask-1, 5] )
            	      endcase
	    	   end  
	        endcase 
                f1 = flag ;and 1B 
                if f1 eq 0 then f1=32 else f1=0
                if state.debug then print, 'DEBUG: flag1: ', flag, f1
                case state.linefit_type of
	   	   0: begin
	   	    case state.linefit_mode of
           	   	0: begin
           	   	 (*state.sp_nrescube)[state.col, state.row, 0] = f1
           	   	 (*state.sp_brescube)[state.col, state.row, 0] = f1
           	   	 (*state.sp_crescube)[state.col, state.row, 0] = f1
           	   	 (*state.sp_nerrrescube)[state.col, state.row, 0, *] = f1
           	   	 (*state.sp_berrrescube)[state.col, state.row, 0, *] = f1
           	   	 (*state.sp_cerrrescube)[state.col, state.row, 0, *] = f1
           	   	end 
           	   	1: begin
           	   	 (*state.nrescube)[state.imask-1, 0] = f1
           	   	 (*state.brescube)[state.imask-1, 0] = f1
           	   	 (*state.crescube)[state.imask-1, 0] = f1
           	   	 (*state.nerrrescube)[state.imask-1, 0, *] = f1
           	   	 (*state.berrrescube)[state.imask-1, 0, *] = f1
           	   	 (*state.cerrrescube)[state.imask-1, 0, *] = f1
           	   	end
           	    endcase
	   	   end
	   	   1: begin
	   	    case state.linefit_mode of
           	   	0: begin
           	   	 (*state.sp_mrescube)[state.col, state.row, 5] = f1
           	   	 (*state.sp_merrrescube)[state.col, state.row, 5, *] = f1
           	   	end 
           	   	1: begin
           	   	 (*state.mrescube)[state.imask-1, 5] = f1
           	   	 (*state.merrrescube)[state.imask-1, 5, *] = f1
           	   	end
           	    endcase
	   	   end
	        endcase    
                if f1 gt 0 then widget_control, widgetids.flag_button, set_value=flag_bad
                if f1 eq 0 then widget_control, widgetids.flag_button, set_value=flag_ok
                if state.cubesel gt 3 then begin
                  kubeviz_plotspax ;linefit result are updated if they are currently displayed
                  kubeviz_plotinfo
                endif  
            end

            'FIT': begin
                case state.fitconstr of
                    0: kubeviz_linefit_dofit ; determine starting param values automatically
                    1: kubeviz_linefit_dofit, /userpar  ; use user-supplied starting param values
                endcase
                doflag = replicate(0, state.Ncol, state.Nrow)
                doflag[state.col, state.row] = 1
                kubeviz_linefit_autoflag, doflag=doflag
		kubeviz_linefit_update
                kubeviz_plotspeczoom
                if state.cubesel gt 3 then begin
                     kubeviz_plotspax ;linefit result are updated if they are currently displayed
                     kubeviz_plotinfo
                endif     
            end
            'FITALL': begin ; fit ALL spaxels or masks             
                case state.linefit_mode of
                    0: kubeviz_linefit_fitall, /spaxel
                    1: kubeviz_linefit_fitall, /mask
                endcase
                ; auto-flag after fitting all:
                kubeviz_linefit_autoflag 
                kubeviz_plotspec
                kubeviz_plotspax
                kubeviz_plotinfo
                kubeviz_plotspeczoom
                kubeviz_linefit_update
            end
            'FLAGALL': begin ; for all spaxels or masks which have been fit, flag as bad those with S/N(lineflux)<3. in all lines.
                kubeviz_linefit_autoflag
                kubeviz_linefit_update
                kubeviz_plotspax
                kubeviz_plotinfo
                
             end
	    'FITADJ': begin
	       case state.linefit_mode of
                    0: begin      ; spaxels
		      kubeviz_linefit_fitadj
                      doflag = replicate(0, state.Ncol, state.Nrow)
                      doflag[state.col, state.row] = 1
                      kubeviz_linefit_autoflag, doflag=doflag
	    	      kubeviz_linefit_update
		      kubeviz_plotspeczoom 
		      if state.cubesel gt 3 then begin
		         kubeviz_plotspax
			 kubeviz_plotinfo
		      endif	 
		    end
		    1: printf, state.log_lun, '[WARNING] Not applicable to masks'
		endcase    
	    end 
            'FITADJALL': begin ; try and fit spaxels with flag='BAD' using solutions from adjacent 'OK' spaxels as first guess solutions.
               state.zoomspax = 1
	       widget_control, /hourglass
               case state.linefit_mode of
                  0: begin      ; spaxels
                     cur_col = state.col
                     cur_row = state.row
                     kubeviz_linefit_fitadjall 
                     state.col = cur_col
                     state.row = cur_row
                     kubeviz_plotspax
		     kubeviz_plotspeczoom
                     kubeviz_plotinfo
                     kubeviz_linefit_update
                  end
                  1: printf, state.log_lun, '[WARNING] Not applicable to masks'
               endcase
            end
            'RESIMAMASK': begin ; switch on/off the masking of results in the display image
                if state.flagmode eq 0 then state.flagmode = 1 else state.flagmode = 0
                ;;; DJW UPDATE 10.2016
                if state.flagmode eq 1 then begin
                   kubeviz_imagebutton_update, state.par_imagebutton, 0
                   state.par_imagebutton = 'FLAG' ; resets to flag image display
                   kubeviz_imagebutton_update, state.par_imagebutton, 1
                endif
                ;;;
                kubeviz_plotspax
                kubeviz_plotinfo
            end
            'SAVE': kubeviz_linefit_saveres

            'RESETALL': begin
	         kubeviz_linefit_reset 
		 kubeviz_linefit_resetuser, /ALL  
		 kubeviz_set_autoflag_defaults
 		 kubeviz_set_linefitrange_defaults
 		 kubeviz_set_continuumfit_defaults
             end

            'RESETUSER': kubeviz_linefit_resetuser, /ALL  
	    
	    'RESETFIT': kubeviz_linefit_reset, /update
	    
	    'RESETFITALL': kubeviz_linefit_resetall, /update
            
            'MODE': begin
                if state.linefit_mode eq 0 then begin
                    state.linefit_mode = 1
                    state.specmode   = 1 ; sum
                    state.cursormode = 2
                endif else begin
                    state.linefit_mode = 0
                    state.specmode   = 0 ; slice
                    state.cursormode = 1
                endelse
                kubeviz_linefit_update
                kubeviz_plotspax
                kubeviz_plotspec
                kubeviz_plotspeczoom
            end

            'PREVMASK':  begin
                state.imask = state.imask - 1
                if state.imask le 0 then state.imask = 1
                kubeviz_medianspec
		kubeviz_linefit_update
                kubeviz_plotspax 
                kubeviz_plotinfo
                kubeviz_plotspec
                kubeviz_plotspeczoom
            end
            'NEXTMASK': begin
            	state.imask = state.imask + 1
            	if state.imask gt state.Nmask then state.imask = state.Nmask
            	kubeviz_medianspec
		kubeviz_linefit_update
            	kubeviz_plotspax 
            	kubeviz_plotinfo
            	kubeviz_plotspec
            	kubeviz_plotspeczoom
            end
            'FITSKY': begin
                ; choose the instrumental resolution fit mode:
                state.instrres_mode = 0
                if total(state.instrres_varpoly) eq 0 then kubeviz_linefit_skylines
                kubeviz_linefit_update, /update_userpars
            end
            'POLYSKY' : begin
                state.instrres_mode = 1
                kubeviz_linefit_update, /update_userpars
            end
	    'TPLSKY' : begin
		state.instrres_mode = 2
                kubeviz_linefit_update
            end
            'TYPE': begin
                if state.linefit_type eq 0 then begin
                   state.linefit_type = 1
                endif else begin
                    state.linefit_type = 0
                endelse
                kubeviz_linefit_typeswitch
                kubeviz_linefit_update
                kubeviz_plotspax
                kubeviz_plotspec
                kubeviz_plotspeczoom
            end

            'FITCONSTRAINTS': state.fitconstr = 1-state.fitconstr
	    
	    'FIXRATIOS': state.fitfixratios = 1-state.fitfixratios
	    
	    'CONTMODE': state.continuumfit_mode = 1-state.continuumfit_mode
            
	    'FREE_2ND': begin
	         state.secondcomp_mode = 0
		 kubeviz_linefit_typeswitch
	    end
	    
	    'HIVEL_2ND': begin
	         state.secondcomp_mode = 1
		 kubeviz_linefit_typeswitch
	    end
	    
	    'BROAD_2ND': begin
	         state.secondcomp_mode = 2
		 kubeviz_linefit_typeswitch
	    end
	    
	    'SMART_2ND': state.secondcomp_smart = 1- state.secondcomp_smart
	    
            else: begin

                if strmid(value,0,3) eq 'FIX' then begin

                    linetype = strmid(value,3,1)
                    case linetype of
                        'N': begin
                            par = fix(strmid(value,4))
                            state.pnfix[par] = 1-state.pnfix[par]
                        end
                        'B': begin
                            par = fix(strmid(value,4))
                            state.pbfix[par] = 1-state.pbfix[par]
                        end
                    endcase

                endif
                if strmid(value,0,3) eq 'FIT' then begin

                    linetype = strmid(value,3,1)
                    par = fix(strmid(value,4))
                    line = par - 3
                    case linetype of
                        'N': state.pndofit[line] = 1-state.pndofit[line]
                        'B': state.pbdofit[line] = 1-state.pbdofit[line]
                        'C': begin
                           ; fit continuum for all lines in set (i.e. fit continuum region).
                           lineset = (*state.linesets)[line] ;kubeviz_linesetforline((*state.linenames)[line])
                           lines_cfit = kubeviz_chooselines((*state.linesets), state.lines, lineset=lineset)
                           if lines_cfit[0] eq -1 then Nlines_cfit=0 else Nlines_cfit = N_Elements(lines_cfit)
                           if Nlines_cfit gt 0 then begin
                              printf, state.log_lun, '[KUBEVIZ] Setting continuum fitting for all lines in set.'
                              for iline=0, Nlines_cfit-1 do begin
                                state.pcdofit[lines_cfit[iline]] = 1-state.pcdofit[lines_cfit[iline]]
                                ipar = lines_cfit[iline] + 3
                                widget_control, widgetids.pcdofitbutton[ipar], set_button=state.pcdofit[lines_cfit[iline]]
                              endfor
                           endif   
                        end
                     endcase

                endif        
                if strmid(value,0,5) eq 'RESET' then begin

                    linetype = strmid(value,5,1)
                    par = fix(strmid(value,6))  
                    kubeviz_linefit_reset, linetype=linetype, pars=[par]

                endif
                if strmid(value,0,4) eq 'SHOW' then begin

                    linetype = strmid(value,4,1)
                    par = fix(strmid(value,5))  
                    line = par-3
                    case linetype of
                        'N': state.nshow[line]=1-state.nshow[line]
                        'B': state.bshow[line]=1-state.bshow[line]
                        'C': begin 
                           ; show continuum for all lines in set (i.e. show continuum region).
                           lineset = (*state.linesets)[line];kubeviz_linesetforline((*state.linenames)[line])
                           lines_cfit = kubeviz_chooselines((*state.linesets), state.lines, lineset=lineset)
                           if lines_cfit[0] eq -1 then Nlines_cfit=0 else Nlines_cfit = N_Elements(lines_cfit)
                           if Nlines_cfit gt 0 then begin
                              for iline=0, Nlines_cfit-1 do begin
                                state.cshow[lines_cfit[iline]] = 1-state.cshow[lines_cfit[iline]]
                                ipar = lines_cfit[iline] + 3
                                widget_control, widgetids.cshowbutton[ipar], set_button=state.cshow[lines_cfit[iline]]
                              endfor
                           endif
                        end
                    endcase
                    kubeviz_plotspeczoom

                endif
                if strmid(value,0,5) eq 'IMAGE' then begin
		    kubeviz_plotspax, fitima_update=strmid(value,5)
                    kubeviz_plotinfo
                    kubeviz_linefit_update
                endif
                
            endelse

        endcase
    
    
    end ; end button events

    else: ; no matching events found

end ; end case loop over all events

end

;-----------------------------------------------------------------------
pro kubeviz_linefit_fitall, spaxel=spaxel, mask=mask

common kubeviz_state

widget_control, /hourglass
state.zoomspax = 1
kubeviz_plotspax

if keyword_set(spaxel) then begin
   ;Initialise time counter
   t_start = systime(/seconds)
   step_count = 0L
   abs_step = (1L*(state.fitallrange[1]-state.fitallrange[0]+1) * (state.fitallrange[3]-state.fitallrange[2]+1) * state.percent_step)/100
   ;Save current selection
   cur_col = state.col
   cur_row = state.row
   ;Display cancel button
   kubeviz_linefit_cancel, /init
   ;Start loop
   for col=state.fitallrange[0],state.fitallrange[1] do begin
       state.col = col
       for row=state.fitallrange[2], state.fitallrange[3] do begin
   	   state.row = row
   	   kubeviz_statusline, step_count, abs_step, t_start
   	   ;Is the spectrum good at all?
   	   if kubeviz_linefit_testspec() then begin
   	     kubeviz_plotspax, /fast_update
   	     case state.fitconstr of
   		 0: kubeviz_linefit_dofit ; determine starting param values automatically
   		 1: kubeviz_linefit_dofit, /userpar ; use user-supplied starting param values
   	     endcase
   	     ;if state.linefitmap eq 1 then kubeviz_linefit_update, /fast_update
   	     if state.zoommap eq 1 then kubeviz_plotspeczoom
   	   endif
	 kubeviz_linefit_cancel, /cancel, event = event
   	 if event then GOTO, JUMPS
   	endfor
   endfor
   JUMPS:
   kubeviz_statusline, /close
   kubeviz_linefit_cancel, /destroy
   state.col = cur_col
   state.row = cur_row
endif

if keyword_set(mask) then begin
     cur_mask = state.imask
     ;Display cancel button
     kubeviz_linefit_cancel, /init
     ;Start loop
     for mask = 1L, state.Nmask do begin
	 state.imask = mask
     	 kubeviz_statuslinecore, '[PROGRES] Fitting mask '+kubeviz_str(mask), 0
     	 if kubeviz_linefit_testspec() then begin
	   kubeviz_plotspax, /fast_update
     	   case state.fitconstr of
     	       0: kubeviz_linefit_dofit ; determine starting param values automatically
     	       1: kubeviz_linefit_dofit, /userpar ; use user-supplied starting param values
     	   endcase
     	   ;if state.linefitmap eq 1 then kubeviz_linefit_update, /fast_update
     	   if state.zoommap eq 1 then kubeviz_plotspeczoom
         endif
	 kubeviz_linefit_cancel, /cancel, event = event
   	 if event then GOTO, JUMPM
     endfor
     JUMPM:
     kubeviz_statusline, /close
     kubeviz_linefit_cancel, /destroy
     state.imask = cur_mask
endif

end

;-----------------------------------------------------------------------
pro kubeviz_spec_reset_range

common kubeviz_state

state.wavrange1[0] = 0 & state.wavrange1[1] = 0
state.wavrange2[0] = 0 & state.wavrange2[1] = 0
(*state.img1)  = fltarr(state.Ncol, state.Nrow)  ; zero image
(*state.img2)  = fltarr(state.Ncol, state.Nrow)  ; zero image
(*state.nimg1) = fltarr(state.Ncol, state.Nrow) ; zero image
(*state.nimg2) = fltarr(state.Ncol, state.Nrow) ; zero image
state.wavsel = 1

end

;-------------------------------------------------------------------------------
function kubeviz_montecarlonoise, mode, operation=operation, range=range

common kubeviz_state

case state.domontecarlo of
   1: begin
     Nmc = state.Nbootstrap
     mccubes = state.bootstrapcubes
   end
   2: begin
     Nmc = state.Nmc1
     mccubes = state.mc1cubes
   end
   3: begin
     Nmc = state.Nmc2
     mccubes = state.mc2cubes
   end
   4: begin
     Nmc = state.Nmc3
     mccubes = state.mc3cubes
   end
endcase

case mode of
   'cube': begin
      mcnoise = fltarr(state.Ncol, state.Nrow, state.Nwpix)
      for col=0,state.Ncol-1 do begin
         for row=0,state.Nrow-1 do begin
           for ik=0, state.Nwpix-1 do begin
             temp = fltarr(Nmc)
             for mccube = 0, Nmc-1 do temp[mccube] = (*mccubes[mccube])[col,row,ik]
             mcnoise[col,row,ik] =.5*(kubeviz_percentile(temp,(*state.montecarlo_percs)[0])-kubeviz_percentile(temp,(*state.montecarlo_percs)[1]))
          endfor
         endfor
      endfor   
   end
   'image': begin
      mcnoise = fltarr(state.Ncol, state.Nrow)
      for col=0,state.Ncol-1 do begin
        for row=0,state.Nrow-1 do begin
            temp = fltarr(Nmc)
            case operation of
              'sum':   for mccube = 0, Nmc-1 do temp[mccube] = total((*mccubes[mccube])[col,row,range[0]:range[1]])
              'med':   for mccube = 0, Nmc-1 do temp[mccube] = median((*mccubes[mccube])[col,row,range[0]:range[1]])
              'w_med': for mccube = 0, Nmc-1 do temp[mccube] = kubeviz_weighted_median((*mccubes[mccube])[col,row,range[0]:range[1]],(1./((*state.noise)[col,row, range[0]:range[1]])^2))
              'w_avg': for mccube = 0, Nmc-1 do temp[mccube] = total((*mccubes[mccube])[col,row,range[0]:range[1]]/((*state.noise)[col,row, range[0]:range[1]])^2) / total(1./((*state.noise)[col,row, range[0]:range[1]])^2)
            endcase
            mcnoise[col,row] = .5*(kubeviz_percentile(temp,(*state.montecarlo_percs)[0])-kubeviz_percentile(temp,(*state.montecarlo_percs)[1]))
        endfor
     endfor
   end
   'spec' : begin
     mcnoise = fltarr(state.Nwpix)
     for ik=0, state.Nwpix-1 do begin
       temp = fltarr(Nmc)
       for mccube = 0, Nmc-1 do temp[mccube] = (*state.medspec_montecarlo[mccube])[ik]
       mcnoise[ik] = .5*(kubeviz_percentile(temp,(*state.montecarlo_percs)[0])-kubeviz_percentile(temp,(*state.montecarlo_percs)[1]))
     endfor
   end

endcase   

return, mcnoise      
      
end



;-----------------------------------------------------------------------
pro kubeviz_medsum_image_update, imgmode=imgmode, lim=lim
common kubeviz_state

if n_elements(imgmode) eq 0 then imgmode=state.imgmode

;setup the correct mode, range and weight option based on imgmode
case 1 of 
 (imgmode eq 0) : mode='none'
 (imgmode eq 1) or (imgmode eq 3)  or (imgmode eq 6)  or (imgmode eq 8)  : mode='sum'
 (imgmode eq 2) or (imgmode eq 4)  or (imgmode eq 5)  or (imgmode eq 7) or $
 (imgmode eq 9) or (imgmode eq 10) or (imgmode eq 11) or (imgmode eq 12) : mode='med'
endcase

case 1 of
 (imgmode eq 1) or (imgmode eq 2) or (imgmode eq 3) or (imgmode eq 4) or (imgmode eq 5) or (imgmode eq 11) : begin 
   range = 1
   if n_elements(lim) eq 0 then lo = state.wavrange1[0] else lo = lim[0]
   if n_elements(lim) eq 0 then hi = state.wavrange1[1] else hi = lim[1]
   if hi lt lo then kubeviz_util_swap, lo, hi
   (*state.img1)  = fltarr(state.Ncol, state.Nrow) ; zero image
   (*state.nimg1) = fltarr(state.Ncol, state.Nrow) ; zero noise image
   (*state.bimg1) = fltarr(state.Ncol, state.Nrow) ; zero noise image
   if hi lt lo then kubeviz_util_swap, lo, hi
 end   
 (imgmode eq 6) or (imgmode eq 7) or (imgmode eq 8) or (imgmode eq 9) or (imgmode eq 10) or (imgmode eq 12) : begin
   range = 2
   if n_elements(lim) eq 0 then lo = state.wavrange2[0] else lo = lim[0]
   if n_elements(lim) eq 0 then hi = state.wavrange2[1] else hi = lim[1]
   if hi lt lo then kubeviz_util_swap, lo, hi
   (*state.img2) = fltarr(state.Ncol, state.Nrow) ; zero image
   (*state.nimg2) = fltarr(state.Ncol, state.Nrow) ; zero noise image
   (*state.bimg2) = fltarr(state.Ncol, state.Nrow) ; zero noise image
  end
  else:   
endcase

if (imgmode eq 3) or (imgmode eq 4) or (imgmode eq 8) or (imgmode eq 9) then weighted = 1 else weighted = 0

case mode of
  'none' : ;Dont do anything if imgmode is slice
  'sum'  : begin
    if range eq  1 then begin
       if lo ne hi then begin
          if weighted eq 0 then begin
             (*state.img1) = total( (*state.datacube)[*,*, lo:hi], 3 ) 
             ; noise image: assumes uncorrelated pixels in wavelength axis for noisecube method otherwise use montecarlo noise
             if (state.domontecarlo eq 0) or (state.useMonteCarlonoise eq 0) then (*state.nimg1) = sqrt( (total( ((*state.noise)[*,*, lo:hi])^2 , 3))) $
             else  (*state.nimg1) = kubeviz_montecarlonoise('image',operation='sum', range=[lo,hi])
          endif else begin               ;Although inconsistent the weighted average is here
             (*state.img1) = total( (*state.datacube)[*,*, lo:hi]/((*state.noise)[*,*, lo:hi])^2, 3 ) / total( 1./((*state.noise)[*,*, lo:hi])^2, 3 )
             ; noise image: assumes uncorrelated pixels in wavelength axis
             if (state.domontecarlo eq 0) or (state.useMonteCarlonoise eq 0) then (*state.nimg1) = sqrt(1./(total( ((*state.noise)[*,*, lo:hi])^(-2) , 3))) $
             else  (*state.nimg1) = kubeviz_montecarlonoise('image',operation='w_avg', range=[lo,hi])  
          endelse
	  ;Setup badpixel image (MF: ASSUMES THAT ALL THE PIXELS IN RANGE ARE BAD TO BE CONSIDERED BAD)
	  (*state.bimg1)  = fix(total( (*state.badpixelmask)[*,*, lo:hi], 3 ,/integer) / (hi-lo+1))
       endif
       (*state.img1)  = kubeviz_remove_badvalues((*state.img1)) 
       (*state.nimg1) = kubeviz_remove_badvalues((*state.nimg1))
    endif
    if range eq 2 then begin  
       if lo ne hi then begin
          if weighted eq 0 then begin
             (*state.img2) = total( (*state.datacube)[*,*, lo:hi], 3 ) 
             ; noise image: assumes uncorrelated pixels in wavelength axis
             if (state.domontecarlo eq 0) or (state.useMonteCarlonoise eq 0) then (*state.nimg2) = sqrt( (total( ((*state.noise)[*,*, lo:hi])^2 , 3)))  $
             else  (*state.nimg2) = kubeviz_montecarlonoise('image',operation='sum', range=[lo,hi])
          endif else begin
             (*state.img2) = total( (*state.datacube)[*,*, lo:hi]/((*state.noise)[*,*, lo:hi])^2, 3 ) / total( 1./((*state.noise)[*,*, lo:hi])^2, 3 )
             ; noise image: assumes uncorrelated pixels in wavelength axis
             if (state.domontecarlo eq 0) or (state.useMonteCarlonoise eq 0) then (*state.nimg2) = sqrt(1./(total( ((*state.noise)[*,*, lo:hi])^(-2) , 3))) $
             else  (*state.nimg2) = kubeviz_montecarlonoise('image',operation='w_avg', range=[lo,hi])  
          endelse
	  ;Setup badpixel image
	  (*state.bimg2)  = fix(total( (*state.badpixelmask)[*,*, lo:hi], 3 ,/integer) / (hi-lo+1))
       endif
       (*state.img2) = kubeviz_remove_badvalues((*state.img2)) 
       (*state.nimg2) = kubeviz_remove_badvalues((*state.nimg2))
    endif
   end
  
  'med': begin
  if range eq 1 then begin
      if lo ne hi then begin
          if weighted eq 0 then begin
             (*state.img1) = median( (*state.datacube)[*,*, lo:hi], dimension=3)
              ; noise for median image is the same as for the mean image 
             if (state.domontecarlo eq 0) or (state.useMonteCarlonoise eq 0) then (*state.nimg1) = sqrt( (total( ((*state.noise)[*,*, lo:hi])^2 , 3)))  / (1+hi-lo) $
             else  (*state.nimg1) = kubeviz_montecarlonoise('image',operation='med', range=[lo,hi])
          endif else begin
             for col = 0, state.Ncol-1 do begin
               for row = 0, state.Nrow-1 do begin
                   (*state.img1)[col,row] = kubeviz_weighted_median( (*state.datacube)[col,row, lo:hi], (1./((*state.noise)[col,row, lo:hi])^2))
               endfor
             endfor 
             ; noise for weighted median image is the same as for the weigthed mean image 
             if (state.domontecarlo eq 0) or (state.useMonteCarlonoise eq 0) then (*state.nimg1) = sqrt(1./(total( ((*state.noise)[*,*, lo:hi])^(-2) , 3)))  $
             else  (*state.nimg1) = kubeviz_montecarlonoise('image',operation='w_med', range=[lo,hi])
          endelse
	   ;Setup badpixel image
	  (*state.bimg1)  = fix(total( (*state.badpixelmask)[*,*, lo:hi], 3 ,/integer) / (hi-lo+1))
      endif
      (*state.img1) = kubeviz_remove_badvalues((*state.img1)) 
      (*state.nimg1) = kubeviz_remove_badvalues((*state.nimg1))
    endif
    if range eq 2 then begin 
      if lo ne hi then begin
          if weighted eq 0 then begin
             (*state.img2) = median( (*state.datacube)[*,*, lo:hi], dimension=3)
             ; noise for median image is the same as for the mean image 
             if (state.domontecarlo eq 0) or (state.useMonteCarlonoise eq 0) then (*state.nimg2) = sqrt( (total( ((*state.noise)[*,*, lo:hi])^2 , 3)))  / (1+hi-lo) $
             else  (*state.nimg2) = kubeviz_montecarlonoise('image',operation='med', range=[lo,hi])
          endif else begin
            for col = 0, state.Ncol-1 do begin
               for row = 0, state.Nrow-1 do begin
                   (*state.img2)[col,row] = kubeviz_weighted_median( (*state.datacube)[col,row, lo:hi], (1./((*state.noise)[col,row, lo:hi])^2))
               endfor
             endfor 
             ; noise for weighted median image is the same as for the weighted mean image 
             if (state.domontecarlo eq 0) or (state.useMonteCarlonoise eq 0) then (*state.nimg2) = sqrt(1./(total( ((*state.noise)[*,*, lo:hi])^(-2) , 3)))  $
             else  (*state.nimg2) = kubeviz_montecarlonoise('image',operation='w_med', range=[lo,hi])
          endelse
	   ;Setup badpixel image
	  (*state.bimg2)  = fix(total( (*state.badpixelmask)[*,*, lo:hi], 3 ,/integer) / (hi-lo+1))
      endif
      (*state.img2) = kubeviz_remove_badvalues((*state.img2)) 
      (*state.nimg2) = kubeviz_remove_badvalues((*state.nimg2))
    endif
  end 
endcase

end 

 
;-----------------------------------------------------------------------
pro kubeviz_plotspax, ps=ps, fast_update=fast_update, fitima_update=fitima_update, nomove=nomove
; plot the spaxel plane
common kubeviz_state
common kubeviz_widgetids
common kubeviz_winsize

wset, widgetids.wid1 

if n_elements(ps) gt 0            then set_plot,'ps'
if n_elements(fast_update) eq 0   then fast_update = 0
if n_elements(nomove) eq 0        then nomove = 1
if n_elements(fitima_update) gt 0 then begin
  state.cubesel = 4
  par_image     = fitima_update
endif else par_image  = state.par_imagebutton


if fast_update eq 0 then begin

   ; summed and median images etc only make sense for data, noise and S/N
   ; cubes. Fit results cubes are treated independently of imgmode.
   
   if state.cubesel ge 4 and state.cubesel le 6 then begin
   	 kubeviz_linefit_image_update, par_image
	 badimage = (*state.badpixelimg)
	 case state.cubesel of
      	    4: if state.flagmode eq 1 then image = *state.lineresimg else image = *state.unm_lineresimg
      	    5: if state.flagmode eq 1 then image = *state.lineerrresimg else image = *state.unm_lineerrresimg
      	    6: if state.flagmode eq 1 then  image = kubeviz_getsnimage(*state.lineresimg, *state.lineerrresimg, /abs) $
   	   	  else image = kubeviz_getsnimage(*state.unm_lineresimg, *state.unm_lineerrresimg, /abs)
	 endcase   
   endif else begin
   	   case state.imgmode of
   	       0: begin
	           badimage = (*state.badpixelmask)[*,*, state.wpix]
   		   case state.cubesel of 
   		       0: image = (*state.datacube)[*,*, state.wpix] 
   		       1: image = (*state.noise)[*,*, state.wpix] 
   		       2: image = badimage
		       3: image = kubeviz_getsnimage((*state.datacube)[*,*, state.wpix],(*state.noise)[*,*, state.wpix]) 
   		   endcase
   	       end
   	       1: begin
	           badimage = *state.bimg1
   		   case state.cubesel of 
   		       0: image = *state.img1 ; summed image1
   		       1: image = *state.nimg1 ; noise for combined image1
		       2: image = badimage
   		       3: image = kubeviz_getsnimage(*state.img1, *state.nimg1)
   		   endcase
   	       end
   	       2: begin
	           badimage = *state.bimg1
   		   case state.cubesel of 
   		       0: image = *state.img1 ; median image1
   		       1: image = *state.nimg1 ; median noise image1
		       2: image = badimage
   		       3: image = kubeviz_getsnimage(*state.img1, *state.nimg1)
   		   endcase
   	       end
   	       3: begin
	           badimage = *state.bimg1
   		   case state.cubesel of 
   		       0: image = *state.img1 ; weighted summed image1
   		       1: image = *state.nimg1 ; noise for combined image1
		       2: image = badimage
   		       3: image = kubeviz_getsnimage(*state.img1, *state.nimg1)
   		   endcase
   	       end
   	       4: begin
	           badimage = *state.bimg1
   		   case state.cubesel of 
   		       0: image = *state.img1 ; weighted median image1
   		       1: image = *state.nimg1 ; weigthed median noise image1
		       2: image = badimage
   		       3: image = kubeviz_getsnimage(*state.img1, *state.nimg1)
   		   endcase
   	       end
   	       5: begin
	           badimage = fix((*state.bimg1 + (*state.badpixelmask)[*,*, state.wpix]) / 2)
   		   case state.cubesel of 
   		       0: image = (*state.datacube)[*,*, state.wpix] - (*state.img1) ; difference image
   		       1: image = sqrt( ((*state.noise)[*,*, state.wpix])^2 + (*state.nimg1)^2 ) ; noise for difference image
		       2: image = badimage
   		       3: image = kubeviz_getsnimage(( (*state.datacube)[*,*, state.wpix] - (*state.img1) ) , ( sqrt( ((*state.noise)[*,*, state.wpix])^2 + (*state.nimg1)^2 ) ))
   		   endcase
   	       end
   	       6: begin
	           badimage = *state.bimg2
   		   case state.cubesel of 
   		       0: image = *state.img2 ; summed image2
   		       1: image = *state.nimg2 ; noise for combined image2
		       2: image = badimage
   		       3: image = kubeviz_getsnimage(*state.img2, *state.nimg2)
   		   endcase
   	       end
   	       7: begin
	           badimage = *state.bimg2
   		   case state.cubesel of 
   		       0: image = *state.img2 ; median image2
   		       1: image = *state.nimg2 ; median noise image2
		       2: image = badimage
   		       3: image = kubeviz_getsnimage(*state.img2, *state.nimg2)
   		   endcase
   	       end
   	       8: begin
	           badimage = *state.bimg2
   		   case state.cubesel of 
   		       0: image = *state.img2 ; weighted summed image2
   		       1: image = *state.nimg2 ; noise for combined image2
   		       2: image = badimage
		       3: image = kubeviz_getsnimage(*state.img2, *state.nimg2)
   		   endcase
   	       end
   	       9: begin
	           badimage = *state.bimg2
   		   case state.cubesel of 
   		       0: image = *state.img2 ; weighted median image2
   		       1: image = *state.nimg2 ; median noise image2
		       2: image = badimage
   		       3: image = kubeviz_getsnimage(*state.img2, *state.nimg2)
   		   endcase
   	       end
   	       10: begin
	       	   badimage = fix((*state.bimg2 + (*state.badpixelmask)[*,*, state.wpix]) / 2)
   		   case state.cubesel of 
   		       0: image = (*state.datacube)[*,*, state.wpix] - (*state.img2) ; difference image
   		       1: image = sqrt( ((*state.noise)[*,*, state.wpix])^2 + (*state.nimg2)^2 ) ; noise for difference image
		       2: image = badimage
   		       3: image = kubeviz_getsnimage( ( (*state.datacube)[*,*, state.wpix] - (*state.img2) ) , ( sqrt( ((*state.noise)[*,*, state.wpix])^2 + (*state.nimg2)^2 ) ) )
   		   endcase
   	       end
   	       11: begin
	           badimage = fix((*state.bimg2 + *state.bimg1) / 2)
   		   case state.cubesel of 
   		       0: image = (*state.img2) - (*state.img1)
   		       1: image = sqrt( (*state.nimg1)^2 + (*state.nimg2)^2 ) 
		       2: image = badimage
   		       3: image = kubeviz_getsnimage( ( (*state.img2) - (*state.img1) ) , ( sqrt( (*state.nimg1)^2 + (*state.nimg2)^2 )  ) )
   		   endcase
   	       end
	       12: begin
	           badimage = fix((*state.bimg2 + *state.bimg1) / 2)
   		   case state.cubesel of 
   		       0: image = (*state.img1) - (*state.img2)
   		       1: image = sqrt( (*state.nimg1)^2 + (*state.nimg2)^2 ) 
		       2: image = badimage
   		       3: image = kubeviz_getsnimage( ( (*state.img1) - (*state.img2) ) , ( sqrt( (*state.nimg1)^2 + (*state.nimg2)^2 )  ) )
   		   endcase
   	       end
   	   endcase
   endelse
  
   ; Replace badpixels with NaNs
   badpixels = float(badimage+1)
   replace = where(badpixels eq 2, Nreplace)
   if Nreplace gt 0 then badpixels[replace] = !values.f_nan
   image *= badpixels

   case state.zcuts of
       1: begin
   	   range = [state.zmin_ima, state.zmax_ima]
	   frame = bytscl(image, range[0], range[1], TOP=state.maxcol)
       end
       2: begin
   	   range = [min(image,/nan), max(image,/nan)]
   	   frame = bytscl(image, range[0], range[1], TOP=state.maxcol)
       end
       3: begin
   	   range = kubeviz_zscale_range(image)
   	   frame = bytscl(image, range[0], range[1], TOP=state.maxcol)
       end
       4:  begin 
  	     if (min(image,/nan) ne max(image,/nan)) and (finite(min(image,/nan)) eq 1 ) then begin
   	      frame = hist_equal(image, OMIN= valmin, OMAX=valmax, TOP=state.maxcol) 
   	      range = [valmin, valmax]
   	     endif else begin
   	       frame = bytscl(image, 0, 0, TOP=state.maxcol)
	       range = [0,0]
   	    endelse
       end
       5: begin
   	   frame = bytarr(state.Ncol, state.Nrow)
   	   range = [state.zmin_ima, state.zmax_ima]
    	   ok = where(image gt 0) 
	   if range[0] le 0 then range[0] = min(image[ok], /nan) 
   	   frame[ok] = bytscl( sqrt(image[ok]), sqrt(range[0]), sqrt(range[1]), TOP=state.maxcol)
       end
       6: begin
   	   frame = bytarr(state.Ncol, state.Nrow)
	   range = [state.zmin_ima, state.zmax_ima]
    	   ok = where(image gt 0)
	   if range[0] le 0 then range[0] = min(image[ok], /nan) 
   	   frame[ok] = bytscl( alog10(image[ok]), alog10(range[0]), alog10(range[1]), TOP=state.maxcol)
       end
       
       11: begin ; 99.5% cut
   	   range = kubeviz_dataclip(image, percentage=99.5 )
   	   frame = bytscl(image, range[0], range[1], TOP=state.maxcol)
       end
       12: begin ; 99% cut
   	   range = kubeviz_dataclip(image, percentage=99.0 )
   	   frame = bytscl(image, range[0], range[1], TOP=state.maxcol)
       end
       13: begin ; 97% cut
   	   range = kubeviz_dataclip(image, percentage=97.0 )
   	   frame = bytscl(image, range[0], range[1], TOP=state.maxcol)
       end
       14: begin ; 95% cut
   	   range = kubeviz_dataclip(image, percentage=95.0 )
   	   frame = bytscl(image, range[0], range[1], TOP=state.maxcol)
       end

    endcase

   frame = frame +255*badimage ;Plot NaNs in white
   (*state.curr_ima) = image
   (*state.curr_byteima) = frame

endif 

if(!d.name eq 'PS') then begin ;FIXME implement zoom on the postscript as well
    device, filename=strmid(state.filename,0,strlen(state.filename)-5)+'.ps', /color, $ 
		 encapsulated=0, xsize=6, ysize=6, /inch, bits_per_pixel=8
    displimage = congrid(frame, 500, 500) ;Rebin to a higher resolution
    kubeviz_plotimage, bytscl(displimage), imgxrange=[0, state.ncol-1], imgyrange=[0,state.nrow-1], xstyle=1, ystyle=1, /iso, $
	       title=state.filename, charsize=1.1, charthick=3.0
    device, /close  ; close postscript file
endif else begin
    oldrange = state.spaxrange
    if state.zoomspax gt 1 then begin
       Nzoomedcol = state.Ncol/(state.zoomspax)
       Nzoomedrow = state.Nrow/(state.zoomspax)
       if not nomove then begin
          state.spaxrange = [state.col-(Nzoomedcol/2),state.col+(Nzoomedcol/2)-1,state.row-(Nzoomedrow/2),state.row+(Nzoomedrow/2)-1]
          if state.spaxrange[0] lt 0 then begin
            state.spaxrange[0] = 0
	    state.spaxrange[1] = Nzoomedcol-1
          endif
          if state.spaxrange[1] gt state.Ncol-1 then begin
            state.spaxrange[0] = state.Ncol-Nzoomedcol 
	    state.spaxrange[1] = state.Ncol-1
          endif
          if state.spaxrange[2] lt 0 then begin
            state.spaxrange[2] = 0 
	    state.spaxrange[3] = Nzoomedrow-1
          endif
          if state.spaxrange[3] gt state.Nrow-1 then begin
            state.spaxrange[2] = state.Nrow-Nzoomedrow 
	    state.spaxrange[3] = state.Nrow-1
          endif 
       endif
     endif else state.spaxrange = [0,state.Ncol-1,0,state.Nrow-1]
     
     ;If range changes then redo the display gridding
     if not array_equal(oldrange, state.spaxrange) or fast_update eq 0 then (*state.curr_byteima_grid) = congrid((*state.curr_byteima)[state.spaxrange[0]:state.spaxrange[1],state.spaxrange[2]:state.spaxrange[3]], state.Ncol * state.zoomfac, state.Nrow * state.zoomfac, /center)

endelse 

tv, (*state.curr_byteima_grid)

if state.cursormode eq 1 then kubeviz_plotcrosshair
if state.cursormode gt 1 then begin
    kubeviz_markspaxel
    kubeviz_plotcrosshair
endif

if fast_update eq 0 then begin
   set_plot,'x'
   kubeviz_plotcolourbar, range
   ; update the zcut histogram if open
   if (xregistered('kubeviz_zcutpars', /noshow)) then kubeviz_plotzcuthist
endif else wait,0.000001 ;No idea why but this is necessary to have a full display of the crosshair in fast update mode

end

;-----------------------------------------------------------------------
pro kubeviz_plotcolourbar, range
common kubeviz_state
common kubeviz_widgetids
common kubeviz_winsize

; plot the colour bar
xsize = winsize.xwin1-100
colbar = fltarr(xsize, 20)
for i = 0, xsize-1 do colbar[i,*] = (1.*state.maxcol/xsize) * i
wset, widgetids.wid4
tv, bytscl(colbar, TOP=state.maxcol)
if abs(range[0]) lt 1e5 then formatmin = '(f8.1)' else formatmin = '(E8.2)'
if abs(range[1]) lt 1e5 then formatmax = '(f8.1)' else formatmax = '(E8.2)'
widget_control, widgetids.colmin_id, set_value=kubeviz_str(range[0], format=formatmin)
widget_control, widgetids.colmax_id, set_value=kubeviz_str(range[1], format=formatmax)

end

;-----------------------------------------------------------------------
pro kubeviz_getspec, spec, nspec, zspec
common kubeviz_state
; get the appropriate (selected) 1d-spectrum

case state.specmode of
    0: begin
        spec = reform((*state.datacube)[state.col,state.row,*])
        nspec = reform((*state.noise)[state.col,state.row,*])
        zspec = replicate(0.,state.nwpix)
    end
    1: begin
        spec  = reform(*state.medspec)
        nspec = reform(*state.nmedspec)
        zspec = replicate(0.,state.nwpix)
    end
    2: begin
        spec  = reform(*state.medspec)
        nspec = reform(*state.nmedspec)
        zspec = replicate(0.,state.nwpix)
    end
    3: begin        
        spec = reform(*state.medspec)
        nspec = reform(*state.nmedspec)
        zspec = replicate(0.,state.nwpix)
    end
     4: begin
        spec = reform((*state.datacube)[state.col,state.row,*]) - reform((*state.medspec))
        nspec = sqrt( reform((*state.noise)[state.col,state.row,*])^2 + reform((*state.nmedspec))^2 ) 
        zspec = replicate(0.,state.nwpix)
    end
   5: begin        
        spec = reform(*state.medspec)
        nspec = reform(*state.nmedspec)
        zspec = replicate(0.,state.nwpix)
    end
endcase

end

;-----------------------------------------------------------------------
pro kubeviz_plotspec
; plot a spectrum (in the spectrum window)
common kubeviz_state
common kubeviz_widgetids
common kubeviz_winsize

wset, widgetids.wid2

z0 = (*state.wave)[state.wpix]

kubeviz_getspec, spec, nspec, zspec

case state.cubesel of
    0: showspec = spec
    1: showspec = nspec
    2: showspec = zspec
    3: begin
        showspec = spec
        showspec2 = nspec
    end
    else: showspec = spec
endcase

;MFumagalli 2015 set buffer according to cube width
if(n_elements(showspec) gt 400) then tmpbfr=100 else tmpbfr=1

if state.scale eq 1 then range = [state.zmin_spec,state.zmax_spec] $
else begin
 ymin = min(showspec[tmpbfr:state.Nwpix-tmpbfr])
 ymax = max(showspec[tmpbfr:state.Nwpix-tmpbfr])
 if finite(ymin) and finite(ymax) then range=[ymin, ymax] else range=[0,1]
endelse

case state.specmode of
    0: title="Spaxel: " + string(state.col, format='(i4)') + "," + string(state.row, format='(i4)')
    1: title='Sum of selected spaxels' 
    2: title='Median of selected spaxels' 
    3: title='Weigh.avg. of selected spaxels' 
    4: title='Median subtracted'
    5: title='Robertson extracted'
endcase 

;Add wavelength information

if (*state.wave)[0] gt 0 then $
       lamrangestr = "   Slice: " + string(state.wpix, format='(i4)')+ $
                     "   !7k!3: " + string(z0,format='(f9.2)') + $
	             "   Range: " + string(state.wavsel, format='(i1)') $
else   lamrangestr = "   Slice: " + string(state.wpix, format='(i4)')+ $
	             "   Range: " + string(state.wavsel, format='(i1)') 		     

title += lamrangestr  

; plot panel. If you change the position values, modify them as well in spec_event DRAW
position = [40, 20, winsize.xwin2-30, winsize.ywin2-22]

erase, 255
plot, (*state.wave), showspec, /nodata, /xstyle, /ystyle, yrange=range, title=title, $
  color=254, background=255, /device, position=position, charsize=1.05 

; highlight selected wavelength range
if state.wavrange1[0] ne state.wavrange1[1] then kubeviz_markwrange, showspec, 1
if state.wavrange2[0] ne state.wavrange2[1] then kubeviz_markwrange, showspec, 2
; overplot the data
oplot, (*state.wave), showspec, color=254

; overplot second set of data?
if state.cubesel eq 3 then oplot, (*state.wave), showspec2, color=253

if state.marker eq 1 then $
    oplot, [z0, z0], [-1e10, 1e10], thick=1.0, color=253

; overplot marker to highlight starting wavelength when selecting a range
if state.npress eq 1 and state.wavsel eq 1 then begin
    z0 = (*state.wave)[state.wavrange1[0]]
    oplot, [z0, z0], [-1E10,1E10], color=254, thick=1.0
endif
if state.npress eq 1 and state.wavsel eq 2 then begin
    z0 = (*state.wave)[state.wavrange2[0]]
    oplot, [z0, z0], [-1E10,1E10], color=254, thick=1.0
endif

plot, (*state.wave), showspec, /nodata, /xstyle, /ystyle, yrange=range, $
  color=254, background=255, /device, position=position, charsize=1.05, /noerase


end

;-----------------------------------------------------------------------
pro kubeviz_plotspeczoom
; plot the zoom-in spectrum 
common kubeviz_state
common kubeviz_widgetids
common kubeviz_winsize

sigtofwhm = 2.35482

wset, widgetids.wid3

z0 = (*state.wave)[state.wpix]
kubeviz_getspec, spec, nspec, zspec

x1 = state.wpix - state.zoomrange
x2 = state.wpix + state.zoomrange
if x1 lt 0 then x1 = 0
if x2 gt state.Nwpix-1 then x2 = state.Nwpix-1
;print, x1, x2

;x1=0
;x2=state.Nwpix-1

case state.cubesel of
    0: showspec = spec
    1: showspec = nspec
    2: showspec = zspec
    3: begin
        showspec = spec
        showspec2 = nspec
    end
    else: showspec = spec
endcase

position = [40, 20, winsize.xwin3-10, winsize.ywin3-20]

if state.scale eq 1 then range = [state.zmin_spec,state.zmax_spec] else begin
    ymax = 1.1*max(showspec[x1:x2])
    ymin = min([0.,1.1*min(showspec[x1:x2])])
    if finite(ymax) and finite(ymin) then range = [ymin,ymax] else range=[0,1]
endelse

; define plot panel
erase, 255
plot, (*state.wave)[x1:x2], showspec[x1:x2], xstyle=5, /ystyle, /nodata, $
  color=254, background=255, /device, position=position, yr=range
; highlight selected wavelength ranges
if state.wavrange1[0] ne state.wavrange1[1] then kubeviz_markwrange, showspec[x1:x2], 1
if state.wavrange2[0] ne state.wavrange2[1] then kubeviz_markwrange, showspec[x1:x2], 2
; plot pixel x-axis (top)
axis, xaxis=1, /xstyle, color=254, xrange=[x1,x2]
; plot wcs x-axis (bottom)
axis, xaxis=0, /xstyle, color=254, xrange=[ (*state.wave)[x1], (*state.wave)[x2] ]
; overplot the data, set psym=10 to plot as a histogram
xx = (*state.wave)[x1:x2]
oplot, xx, showspec[x1:x2], color=254, psym=10
; overplot second set of data?
if state.cubesel eq 3 then oplot, (*state.wave)[x1:x2], showspec2[x1:x2], color=253

; overplot any fitted gaussians from the linefit window:
kubeviz_linefit_pointerswitch, /load

case state.linefit_mode of
    0: begin ; spaxels:
        nrescube = reform((*state.sp_nrescube)[state.col,state.row,*])
        brescube = reform((*state.sp_brescube)[state.col,state.row,*])
        crescube = reform((*state.sp_crescube)[state.col,state.row,*])
        mrescube = reform((*state.sp_mrescube)[state.col,state.row,*])
    end
    1: begin ; masks:
        nrescube = reform((*state.nrescube)[state.imask-1,*])
        brescube = reform((*state.brescube)[state.imask-1,*])
        crescube = reform((*state.crescube)[state.imask-1,*])
        mrescube = reform((*state.mrescube)[state.imask-1,*])
    end
endcase

if state.Nlines gt 0 and state.cubesel ne 1 and state.cubesel ne 2 then begin
  nfitlines = nrescube[3:state.Nlines+2]
  bfitlines = brescube[3:state.Nlines+2]
  cfitlines = crescube[1:state.Nlines]
  mfitlines = mrescube[0:6*state.Nlines-1:6]

  nshowline = where(nfitlines gt 0. and state.nshow eq 1, Nnshowline)
  bshowline = where(bfitlines gt 0. and state.bshow eq 1, Nbshowline)
  mshowline = where(                    state.nshow eq 1, Nmshowline)       ;mfitlines gt 0. and

  ; sum of lines:
  totfit = replicate(0., x2-x1+1)
  
  if state.linefit_type eq 0 then begin ;GAUSSIAN
     if Nnshowline gt 0 or Nbshowline gt 0 then begin
      	  cont = replicate(!values.f_nan, x2-x1+1,state.lineset_max+1)
      	  for lineset=1,state.lineset_max do begin
      	     lines_cfit = kubeviz_chooselines((*state.linesets), state.lines, lineset=lineset)
      	     if lines_cfit[0] ne -1 then begin 
      	       centline   = median([state.lines[lines_cfit]],/even)
      	       xxset= where(xx ge centline-state.maxwoffb and xx le centline+state.maxwoffr, Nxxset)
      	       if Nxxset gt 0 then begin
	          contvals = cfitlines[lines_cfit]
      	          okcontvals = where(contvals ne -999 and contvals ne 0, Nok)
		  if Nok gt 0 then cont[xxset,lineset] = mean(contvals[okcontvals]) else cont[xxset,lineset]=0.00
      	       endif
      	       xxbluecont= where(xx ge centline-state.continuumfit_maxoff and xx le centline-state.continuumfit_minoff, Nxxbluecont)
      	       if Nxxbluecont gt 0 then begin
      	          oplot, replicate(xx[xxbluecont[0]],2), [-1E10,1E10], color=252, linestyle=3, thick=1.0 
      	          oplot, replicate(xx[xxbluecont[Nxxbluecont-1]],2), [-1E10,1E10], color=252, linestyle=3, thick=1.0 
      	       endif
      	       xxredcont= where(xx ge centline+state.continuumfit_minoff and xx le centline+state.continuumfit_maxoff, Nxxredcont)
      	       if Nxxredcont gt 0 then begin
      	          oplot, replicate(xx[xxredcont[0]],2), [-1E10,1E10], color=252, linestyle=3, thick=1.0 
      	          oplot, replicate(xx[xxredcont[Nxxredcont-1]],2), [-1E10,1E10], color=252, linestyle=3, thick=1.0 
      	       endif
      	     endif 
      	  endfor 
	  totfit += median(cont, dimension=2, /even) 
     endif
     
     if Nnshowline gt 0 then begin 	    
     	  for i=0,Nnshowline-1 do begin
              line = nshowline[i]
	      set  = (*state.linesets)[line] ;kubeviz_linesetforline((*state.linenames)[line])
              dv = nrescube[1]
              if dv eq -999. then dv = 0.
              pos = state.lines[line]*(1.+dv/state.ckms)
              sigv = nrescube[2]
              if sigv eq -999. then sigv = 0.
              width = state.lines[line]*(sigv/state.ckms)
              instrres_A = pos/(sigtofwhm*kubeviz_getinstrres(lambda=pos)) ; convert from R to Angstroms
              totwidthsq = width^2 + instrres_A^2
              norm = nfitlines[line]
              gn = (norm/sqrt(2.*!DPI*totwidthsq))*exp(-((xx-pos)^2/totwidthsq)/2.)
	      totfit +=  gn
	      if total(gn) gt 0 then begin ; Is the line visible at all in the current view?
	       if state.cshow[nshowline[i]] eq 1 then oplot, xx, (cont[*,set] + gn), color=253, linestyle=0, thick=2. $
	       else oplot, xx, gn,  color=253, linestyle=0, thick=2.    
	      endif
          endfor
      endif
      
      if Nbshowline gt 0 then begin
          for i=0,Nbshowline-1 do begin
              line = bshowline[i]
              set  = (*state.linesets)[line] ;kubeviz_linesetforline((*state.linenames)[line])
	      dv = brescube[1]
              if dv eq -999. then dv = 0.
              pos = state.lines[line]*(1.+dv/state.ckms)
              sigv = brescube[2]
              if sigv eq -999. then sigv = 0.
              width = state.lines[line]*(sigv/state.ckms)
              instrres_A = pos/(sigtofwhm*kubeviz_getinstrres(lambda=pos)) ; convert from R to Angstroms
              totwidthsq = width^2 + instrres_A^2
              norm = bfitlines[line]
              gb = (norm/sqrt(2.*!DPI*totwidthsq))*exp(-((xx-pos)^2/totwidthsq)/2.)
              totfit += gb
              if total(gb) gt 0 then begin ; Is the line visible at all in the current view?
	       if state.cshow[bshowline[i]] eq 1 then oplot, xx, (cont[*,set] + gb), color=252, linestyle=0, thick=2. $
	       else oplot, xx, gb,  color=252, linestyle=0, thick=2.
	      endif 
          endfor
      endif
      if Nnshowline+Nbshowline gt 0 then begin
        okregion = where(abs(totfit) gt 1E-30, Nok)
        if Nok gt 0 then oplot, xx[okregion], totfit[okregion], color=249, linestyle=2, thick=2.  
      endif	
  endif else begin ;MOMENTS
      if Nmshowline gt 0 then begin
      
       cont = replicate(!values.f_nan, x2-x1+1,state.lineset_max+1)
       for lineset=1,state.lineset_max do begin
          lines_cfit = kubeviz_chooselines((*state.linesets), state.lines, lineset=lineset)
          if lines_cfit[0] ne -1 then begin 
            centline   = median([state.lines[lines_cfit]],/even)
            xxset= where(xx ge centline-state.maxwoffb and xx le centline+state.maxwoffr, Nxxset)
            if Nxxset gt 0 then begin
               cont[xxset,lineset] = mean(cfitlines[lines_cfit])
            endif
            xxbluecont= where(xx ge centline-state.continuumfit_maxoff and xx le centline-state.continuumfit_minoff, Nxxbluecont)
            if Nxxbluecont gt 0 then begin
                oplot, replicate(xx[xxbluecont[0]],2), [-1E10,1E10], color=252, linestyle=3, thick=1.0 
      	        oplot, replicate(xx[xxbluecont[Nxxbluecont-1]],2), [-1E10,1E10], color=252, linestyle=3, thick=1.0 
      	    endif
            xxredcont= where(xx ge centline+state.continuumfit_minoff and xx le centline+state.continuumfit_maxoff, Nxxredcont)
            if Nxxredcont gt 0 then begin
                oplot, replicate(xx[xxredcont[0]],2), [-1E10,1E10], color=252, linestyle=3, thick=1.0 
      	        oplot, replicate(xx[xxredcont[Nxxredcont-1]],2), [-1E10,1E10], color=252, linestyle=3, thick=1.0 
            endif
          endif 
       endfor 
      
       totfit += median(cont, dimension=2, /even) 
      
          for i=0,Nmshowline-1 do begin
              line = mshowline[i]
              set  = (*state.linesets)[line]
	      dv = mrescube[6*line+1]
              if dv eq -999. then dv = 0.
              pos = state.lines[line]*(1.+dv/state.ckms)
              sigv = mrescube[6*line+2]
              if sigv eq -999. then sigv = 0.
              width = state.lines[line]*(sigv/state.ckms)
              instrres_A = pos/(sigtofwhm*kubeviz_getinstrres(lambda=pos)) ; convert from R to Angstroms
              totwidthsq = width^2 + instrres_A^2
              norm = mfitlines[line]
              gn = (norm/sqrt(2.*!DPI*totwidthsq))*exp(-((xx-pos)^2/totwidthsq)/2.)
	      totfit += gn
              if total(gn) gt 0 then begin ; Is the line visible at all in the current view?
		if state.cshow[mshowline[i]] eq 1 then oplot, xx, (cont[*,set] + gn), color=253, linestyle=0, thick=2.  $
	        else oplot, xx, gn,  color=253, linestyle=0, thick=2.
              endif
	  endfor
       endif
       if Nmshowline gt 0 then begin
          okregion = where(abs(totfit) gt 1E-30, Nok)
          if Nok gt 0 then oplot, xx[okregion], totfit[okregion], color=249, linestyle=2, thick=2.  
       endif
  endelse

endif  

; overplot marker to highlight selected wavelength
oplot, [z0, z0], [-1E10,1E10], color=253, linestyle=2, thick=1.0

; overplot marker to highlight starting wavelength when selecting a range
if state.npress eq 1 and state.wavsel eq 1 then begin
    z0 = (*state.wave)[state.wavrange1[0]]
    oplot, [z0, z0], [-1E10,1E10], color=254, thick=2.0
endif
if state.npress eq 1 and state.wavsel eq 2 then begin
    z0 = (*state.wave)[state.wavrange2[0]]
    oplot, [z0, z0], [-1E10,1E10], color=254, thick=2.0
endif



end

;-----------------------------------------------------------------------
pro kubeviz_plotinfo
common kubeviz_state
common kubeviz_widgetids

if state.col lt 1000 then formatcol ='(i3)' else formatcol='(i4)'
if state.row lt 1000 then formatrow ='(i3)' else formatrow='(i4)'

imacoords  = '(' + string(state.col, format=formatcol) + ',' $
                + string(state.row, format=formatrow) + ')' 

physcoords = '(' + string(state.col+state.Startcol, format=formatcol) + ',' $
                 + string(state.row+state.Startrow, format=formatrow) + ')' 

widget_control, widgetids.pixima_id,  set_value=imacoords
widget_control, widgetids.pixphys_id, set_value=physcoords

case state.cubesel of
    0: begin ; data
        case state.imgmode of
          0: val = (*state.badpixelmask)[state.col, state.row, state.wpix] eq 1 ? !values.f_nan : (*state.datacube)[state.col, state.row, state.wpix]
          1: val = (*state.bimg1)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.img1)[state.col, state.row]		       
          2: val = (*state.bimg1)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.img1)[state.col, state.row]		       
          3: val = (*state.bimg1)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.img1)[state.col, state.row]		       
          4: val = (*state.bimg1)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.img1)[state.col, state.row]		       
          5: val = (*state.bimg1)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.datacube)[state.col, state.row, state.wpix] - (*state.img1)[state.col, state.row]
          6: val = (*state.bimg2)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.img2)[state.col, state.row]
          7: val = (*state.bimg2)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.img2)[state.col, state.row]
          8: val = (*state.bimg2)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.img2)[state.col, state.row]
          9: val = (*state.bimg2)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.img2)[state.col, state.row]
         10: val = (*state.bimg2)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.datacube)[state.col, state.row, state.wpix] - (*state.img2)[state.col, state.row]
         11: val = (*state.bimg2)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.img2)[state.col, state.row] - (*state.img1)[state.col, state.row]
	 12: val = (*state.bimg2)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.img1)[state.col, state.row] - (*state.img2)[state.col, state.row]
        end                                                                       
        widget_control, widgetids.cube_id, set_value='DATA'
    end
    1: begin ; noise
        case state.imgmode of
          0: val = (*state.badpixelmask)[state.col, state.row, state.wpix] eq 1 ? !values.f_nan : (*state.noise)[state.col, state.row, state.wpix]
          1: val = (*state.bimg1)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.nimg1)[state.col, state.row]
          2: val = (*state.bimg1)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.nimg1)[state.col, state.row]
          3: val = (*state.bimg1)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.nimg1)[state.col, state.row]
          4: val = (*state.bimg1)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.nimg1)[state.col, state.row]
          5: val = (*state.bimg1)[state.col, state.row] 		   eq 1 ? !values.f_nan : sqrt(((*state.noise)[state.col, state.row, state.wpix])^2 + ((*state.nimg1)[state.col, state.row])^2)
          6: val = (*state.bimg2)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.nimg2)[state.col, state.row]
          7: val = (*state.bimg2)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.nimg2)[state.col, state.row]
          8: val = (*state.bimg2)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.nimg2)[state.col, state.row]
          9: val = (*state.bimg2)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.nimg2)[state.col, state.row]
         10: val = (*state.bimg2)[state.col, state.row] 		   eq 1 ? !values.f_nan : sqrt(((*state.noise)[state.col, state.row, state.wpix])^2 + ((*state.nimg2)[state.col, state.row])^2)
         11: val = (*state.bimg2)[state.col, state.row] 		   eq 1 ? !values.f_nan : sqrt(((*state.nimg1)[state.col, state.row])^2 + ((*state.nimg2)[state.col, state.row])^2)
         12: val = (*state.bimg2)[state.col, state.row] 		   eq 1 ? !values.f_nan : sqrt(((*state.nimg1)[state.col, state.row])^2 + ((*state.nimg2)[state.col, state.row])^2)
	end
        widget_control, widgetids.cube_id, set_value='NOISE'
    end
    2: begin ; bad pixels
        val = (*state.badpixelmask)[state.col, state.row, state.wpix]
        widget_control, widgetids.cube_id, set_value='BAD'
    end
    3: begin ;S/N
        case state.imgmode of
          0: val = (*state.badpixelmask)[state.col, state.row, state.wpix] eq 1 ? !values.f_nan : (*state.datacube)[state.col, state.row, state.wpix]/(*state.noise)[state.col, state.row, state.wpix]
          1: val = (*state.bimg1)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.img1)[state.col, state.row]/(*state.nimg1)[state.col, state.row]
          2: val = (*state.bimg1)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.img1)[state.col, state.row]/(*state.nimg1)[state.col, state.row]
          3: val = (*state.bimg1)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.img1)[state.col, state.row]/(*state.nimg1)[state.col, state.row]
          4: val = (*state.bimg1)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.img1)[state.col, state.row]/(*state.nimg1)[state.col, state.row]
          5: val = (*state.bimg1)[state.col, state.row] 		   eq 1 ? !values.f_nan : ((*state.datacube)[state.col, state.row, state.wpix] - (*state.img1)[state.col, state.row])/(sqrt(((*state.noise)[state.col, state.row, state.wpix])^2 + ((*state.nimg1)[state.col, state.row])^2))
          6: val = (*state.bimg2)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.img2)[state.col, state.row]/(*state.nimg2)[state.col, state.row]
          7: val = (*state.bimg2)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.img2)[state.col, state.row]/(*state.nimg2)[state.col, state.row]
          8: val = (*state.bimg2)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.img2)[state.col, state.row]/(*state.nimg2)[state.col, state.row]
          9: val = (*state.bimg2)[state.col, state.row] 		   eq 1 ? !values.f_nan : (*state.img2)[state.col, state.row]/(*state.nimg2)[state.col, state.row]
         10: val = (*state.bimg2)[state.col, state.row] 		   eq 1 ? !values.f_nan : ((*state.datacube)[state.col, state.row, state.wpix] - (*state.img2)[state.col, state.row])/(sqrt(((*state.noise)[state.col, state.row, state.wpix])^2 + ((*state.nimg2)[state.col, state.row])^2))
         11: val = (*state.bimg2)[state.col, state.row] 		   eq 1 ? !values.f_nan : ((*state.img2)[state.col, state.row] - (*state.img1)[state.col, state.row])/sqrt(((*state.nimg1)[state.col, state.row])^2 + ((*state.nimg2)[state.col, state.row])^2)
         12: val = (*state.bimg2)[state.col, state.row] 		   eq 1 ? !values.f_nan : ((*state.img1)[state.col, state.row] - (*state.img2)[state.col, state.row])/sqrt(((*state.nimg1)[state.col, state.row])^2 + ((*state.nimg2)[state.col, state.row])^2)
	end
        widget_control, widgetids.cube_id, set_value='S/N'
    end
    4: begin ; linefit
        if state.flagmode eq 1 then val = (*state.lineresimg)[state.col, state.row] else $
        val = (*state.unm_lineresimg)[state.col, state.row]
        widget_control, widgetids.cube_id, set_value='FIT'
    end
    5: begin ; linefit errors
        if state.flagmode eq 1 then val = (*state.lineerrresimg)[state.col, state.row] else $
        val = (*state.unm_lineerrresimg)[state.col, state.row]
        widget_control, widgetids.cube_id, set_value='FIT ERR'
    end
    6: begin ; linefit S/N
        if state.flagmode eq 1 then $ 
        val = ((*state.lineerrresimg)[state.col, state.row] gt 0.) ? abs((*state.lineresimg)[state.col,state.row])/(*state.lineerrresimg)[state.col, state.row] : 0.D else $
        val = ((*state.unm_lineerrresimg)[state.col, state.row] gt 0.) ? abs((*state.unm_lineresimg)[state.col,state.row])/(*state.unm_lineerrresimg)[state.col, state.row] : 0.D
        widget_control, widgetids.cube_id, set_value='FIT S/N'
    end
endcase

widget_control, widgetids.spatzoom_id, set_value='Zoom: '+string(state.zoomspax, format='(i2)')

widget_control, widgetids.pixval_id, set_value=string(val, format='(f12.4)')

smooth_string = string(state.smooth, format='(i2)') + 'x' + strtrim(string(state.smooth, format='(i2)'),2) + 'x' + strtrim(string(state.specsmooth, format='(i3)'),2)
widget_control, widgetids.smooth_id, set_value=smooth_string

mask_string = string(state.imask, format='(i3)') + '/' + strtrim(string(state.Nmask, format='(i3)'), 2)
widget_control, widgetids.imask_id, set_value=mask_string


; Use an external module to compute wcs solution
extast, (*state.indatahead),  astr, noparams   
z0 = string((*state.wave)[state.wpix], format='(f9.2)')
case 1 of 
   noparams eq 2 and (*state.wave)[0] gt 0: begin
     xyad, (*state.indatahead), state.col+state.Startcol, state.row+state.Startrow ,x0,y0
     radec = strsplit(adstring(x0, y0, 1), /Extract)
     wcscoords = '(' + radec[0] + ':' + radec[1] + ':' + radec[2] +', '+ radec[3] + ':' + radec[4] + ':' + radec[5] + ', ' +z0 + ')'

     ;wcscoords = '(' + string(x0, format='(f10.5)') + ',' $
     ;                + string(y0, format='(f10.5)') + ',' $
     ;                + string(z0, format='(f9.2)') + ')'
   end		     
   noparams ne 2 and (*state.wave)[0] gt 0: wcscoords = '(SPATIAL WCS NOT FOUND,' + string(z0, format='(f9.2)') + ')' 
   else: wcscoords = '(WCS NOT FOUND)'   
endcase

widget_control, widgetids.wcs_id, set_value=wcscoords

case state.imgmode of
    0: widget_control, widgetids.imgmode_id, set_value="Slice "+kubeviz_str(state.wpix)
    1: widget_control, widgetids.imgmode_id, set_value="Sum1"
    2: widget_control, widgetids.imgmode_id, set_value="Median1"
    3: widget_control, widgetids.imgmode_id, set_value="Weighted Avg1"
    4: widget_control, widgetids.imgmode_id, set_value="Weighted Med1"
    5: widget_control, widgetids.imgmode_id, set_value="MedSub1"
    6: widget_control, widgetids.imgmode_id, set_value="Sum2"
    7: widget_control, widgetids.imgmode_id, set_value="Median2"
    8: widget_control, widgetids.imgmode_id, set_value="Weighted Avg2"
    9: widget_control, widgetids.imgmode_id, set_value="Weighted Med2"
   10: widget_control, widgetids.imgmode_id, set_value="MedSub2"
   11: widget_control, widgetids.imgmode_id, set_value="Med2-Med1"
   12: widget_control, widgetids.imgmode_id, set_value="Med1-Med2"
   else: print, '[  BUG  ] error in imgmode widget update'
endcase

end

;-----------------------------------------------------------------------
pro kubeviz_keyboard_handler, key, ch, update, nomove
; process keyboard selections
common kubeviz_state
common kubeviz_widgetids

case key of
    6: begin
        update = 2
        state.col = state.col + 1
        if state.col ge state.Ncol then state.col = state.Ncol-1
    end
    5: begin
        update = 2
        state.col = state.col - 1
        if state.col le 0 then state.col = 0
    end
    
    7: begin
        update = 2
        state.row = state.row + 1
        if state.row ge state.Nrow then state.row = state.Nrow-1
    end
    8: begin
        update = 2
        state.row = state.row - 1
        if state.row le 0 then state.row = 0
    end
    else:
endcase

case string(ch) of
    "q": kubeviz_destroy
    ",": begin
        update = 1
        state.wpix = state.wpix - 1
        widget_control, widgetids.slid_id, set_value=state.wpix
        if state.wpix lt 0 then state.wpix = 0
    end
    ".": begin
        update = 1
        state.wpix = state.wpix + 1
        widget_control, widgetids.slid_id, set_value=state.wpix
        if state.wpix ge state.Nwpix then state.wpix = state.Nwpix-1
    end
    "m": begin  ; select spaxel
        (*state.spaxselect)[state.col, state.row, state.imask-1] = 1
        update = 1
    end
    "n": begin  ; de-select spaxel
        (*state.spaxselect)[state.col, state.row, state.imask-1] = 0
        update = 1
    end
    "r": begin  ; clear spaxel mask
        kubeviz_clear_select
        update = 1
    end
    "k": begin  ; decrease imask index
        state.imask -= 1
        if state.imask le 0 then state.imask = 1
        update = 1
    end
    "l": begin  ; increase imask index
        state.imask += 1
        if state.imask gt state.Nmask then state.imask = state.Nmask
        update = 1
    end
    "c":  begin ;Fit a gaussian to the current image
        kubeviz_findcentroid 
	update=0
    end	
    "f": begin ; fit using linefit:
        case state.fitconstr of
            0: kubeviz_linefit_dofit ; determine starting param values automatically
            1: kubeviz_linefit_dofit, /userpar ; use user-supplied starting param values
        endcase
	doflag = replicate(0, state.Ncol, state.Nrow)
        doflag[state.col, state.row] = 1
        kubeviz_linefit_autoflag, doflag=doflag

        update=1
	nomove=1
    end
    "a": begin ;fit using adjacent spaxels initial guesses
        case state.linefit_mode of
           0: begin	 ; spaxels
	      kubeviz_linefit_fitadj
	      doflag = replicate(0, state.Ncol, state.Nrow)
              doflag[state.col, state.row] = 1
              kubeviz_linefit_autoflag, doflag=doflag
	      kubeviz_plotspeczoom 
              kubeviz_linefit_update
	      if state.cubesel gt 3 then begin
	         kubeviz_plotspax
		 kubeviz_plotinfo
	      endif	 

	   end
	   1: printf, state.log_lun, '[WARNING] Not applicable to masks'
       endcase    
    end
    "b": begin ; toggle ok/bad fitting flag 
	
	;Extract flag here
	case state.linefit_type of 
	    0: begin
	       case state.linefit_mode of
                    0: begin
		         case strmid(state.par_imagebutton,0,1) of
                            'N': flag  = fix( (*state.sp_nrescube)[state.col, state.row, 0] )
	    	    	    'B': flag  = fix( (*state.sp_brescube)[state.col, state.row, 0] )
	    	    	    'C': flag  = fix( (*state.sp_crescube)[state.col, state.row, 0] )
	    	    	    else: flag = fix( (*state.sp_nrescube)[state.col, state.row, 0] )
	                 endcase
                    end
		    1: begin
		         case strmid(state.par_imagebutton,0,1) of
			    'N': flag = fix( (*state.nrescube)[state.imask-1, 0] )
			    'B': flag = fix( (*state.brescube)[state.imask-1, 0] )
			    'C': flag = fix( (*state.crescube)[state.imask-1, 0] )
			    else: flag = fix( (*state.nrescube)[state.imask-1, 0] )
			 endcase 
		    end	   
               endcase
	    end
	    1: begin
	       case state.linefit_mode of
                    0: flag = fix( (*state.sp_mrescube)[state.col, state.row, 5] )
                    1: flag = fix( (*state.mrescube)[state.imask-1, 5] )
               endcase
	    end  
	 endcase     
        
	;Modify flag here
	f1 = flag
        if f1 eq 0 then f1=16 else f1=0
        if state.debug then print, '[DEBUG] flag: ', flag, f1
        
	;Save new flag here
	case state.linefit_type of
	   0: begin
	    case state.linefit_mode of
            	0: begin
                 case strmid(state.par_imagebutton,0,1) of
		   'N': begin
		      (*state.sp_nrescube)[state.col, state.row, 0] = f1
                      (*state.sp_nerrrescube)[state.col, state.row, 0, *] = f1
		   end
		   'B': begin
		      (*state.sp_brescube)[state.col, state.row, 0] = f1
                      (*state.sp_berrrescube)[state.col, state.row, 0, *] = f1
		   end
		   'C': begin
		      (*state.sp_crescube)[state.col, state.row, 0] = f1
                      (*state.sp_cerrrescube)[state.col, state.row, 0, *] = f1
                   end
		   else: begin
		      (*state.sp_nrescube)[state.col, state.row, 0] = f1
                      (*state.sp_nerrrescube)[state.col, state.row, 0, *] = f1
		   end
		  endcase    
               	end 
            	1: begin
                 case strmid(state.par_imagebutton,0,1) of
		   'N': begin
		      (*state.nrescube)[state.imask-1, 0] = f1
                      (*state.nerrrescube)[state.imask-1, 0, *] = f1
		    end
		   'B': begin
                      (*state.brescube)[state.imask-1, 0] = f1
                      (*state.berrrescube)[state.imask-1, 0, *] = f1
		   end
		   'C': begin
                      (*state.crescube)[state.imask-1, 0] = f1
                      (*state.cerrrescube)[state.imask-1, 0, *] = f1
		   end
		   else: begin
		      (*state.nrescube)[state.imask-1, 0] = f1
                      (*state.nerrrescube)[state.imask-1, 0, *] = f1
		   end
		 endcase   
            	end
            endcase
	   end
	   1: begin
	    case state.linefit_mode of
            	0: begin
                 (*state.sp_mrescube)[state.col, state.row, 5] = f1
                 (*state.sp_merrrescube)[state.col, state.row, 5, *] = f1
            	end 
            	1: begin
                 (*state.mrescube)[state.imask-1, 5] = f1
                 (*state.merrrescube)[state.imask-1, 5, *] = f1
            	end
            endcase
	   end
	endcase    
        update=1
	nomove=1
    end
    else:
endcase

end

;-----------------------------------------------------------------------
pro kubeviz_destroy, abort=abort

common kubeviz_state
common kubeviz_winsize
common kubeviz_widgetids
common kubeviz_fit
common kubeviz_linesdb

if n_elements(abort) eq 0 then printf, state.log_lun, '[KUBEVIZ] Quitting...' else print, '[ ERROR ] Aborting...'

; Test if err_msg.txt is empty and if so removes it
openr, lun, state.cwdir+'fit_errmsg_'+strmid(state.filename,0,strlen(state.filename)-5)+'.txt', /get_lun
stat = fstat(lun)
if stat.size eq 0 then spawn, 'rm fit_errmsg_'+strmid(state.filename,0,strlen(state.filename)-5)+'.txt'
free_lun, lun

;destroy existent graphic windows
if xregistered('kubeviz_spax',/noshow) gt 0 then begin 
   widget_control, /destroy, widgetids.base1         
   widget_control, /destroy, widgetids.base5                

   ; save the session into default save-file:
   kubeviz_savesession, fname=state.cwdir+'lastsession.sav'
endif

;recursively free all heap variables
heap_free, state
heap_free, widgetids
heap_free, winsize
heap_free, fit
heap_free, linesdb

;Additional garbage collector for heap variables not
;found by heap_free. Only to be on the safe side.
;Also free all the open logic units
heap_gc 

if n_elements(abort) eq 0 then printf, state.log_lun, '[KUBEVIZ] Quit successful.' else retall

end


;-----------------------------------------------------------------------
pro kubeviz_selectcube, fname, path, fnoisename, noisepath, dataonly=dataonly, noiseonly=noiseonly
; select a data and/or noise FITS cube

if n_elements(dataonly)  gt 0 then readnoise = 0 else readnoise = 1
if n_elements(noiseonly) gt 0 then readdata = 0  else readdata  = 1

if readdata eq 1 then begin
 filename = dialog_pickfile(filter='*.fits', /read, /must_exist, get_path=path, title='Select a cube for reading')
 if (strlen(filename) - strlen(path)) eq 0 then return ; don't do anything if no filename is given
 fname = file_basename(filename)
endif

fits_open, path+fname, fcb
if fcb.nextend lt 2 and readnoise eq 1 then begin
  noisefilename = dialog_pickfile(filter='*.fits', /read, /must_exist, get_path=noisepath, title='Select the noise cube for reading')
  if (strlen(noisefilename) - strlen(noisepath)) eq 0 then return ; don't do anything if no filename is given
  fnoisename = file_basename(noisefilename)
endif else begin
  fnoisename = fname
  noisepath = path
endelse
fits_close, fcb

end

;-----------------------------------------------------------------------
pro kubeviz_splitpath, fullfname, fname, path
; splits file name into the name of the file and its path

fname = file_basename(fullfname)
if fname eq fullfname then path ='' else path = file_dirname(fullfname,/mark_directory) 

end

;--------------------------------------------------------
pro kubeviz_savesession, fname=fname
; save full session (state variable)
common kubeviz_state

if n_elements(fname) eq 0 then begin
    fname = dialog_pickfile(filter='*.sav', path=state.cwdir, /write, get_path=path, /OVERWRITE_PROMPT, title='Save session as:')
    len = strlen(fname) - strlen(path)
    if len eq 0 then return     ; don't do anything if no filename is given
endif

save, state, filename=fname

printf, state.log_lun, '[KUBEVIZ] Session saved.'

end

;------------------------------------------
pro kubeviz_loadsession, fname, path, spmask=spmask

; load full session (state variable)
common kubeviz_state
common kubeviz_widgetids

; Save name and size of the current cube
current_file = state.filename
current_z    = state.redshift
current_ncol = state.Ncol
current_nrow = state.Nrow
current_cwdir = state.cwdir
current_outdir = state.outdir
dummy = where(state.lines gt 0, current_nlines)

if N_elements(fname) eq 0 then begin
   fname = dialog_pickfile(filter='*.sav', path=state.cwdir, /read, /must_exist, get_path=path, title='Load session file:')
   len = strlen(fname) - strlen(path)
   if len eq 0 then return ; don't do anything if no filename is given
endif else fname=path+fname

printf, state.log_lun, '[KUBEVIZ] Loading session: ',fname

; take a copy of the existing state structure:
current_state = state

; restore state variable
heap_free, state
restore, fname, /RELAXED_STRUCTURE_ASSIGNMENT
state.cwdir = current_cwdir
; DJW 13-10-16. Check existence of outdir.
if file_test(state.outdir,/directory) eq 0 then state.outdir = current_outdir

; copy contents into old_state and rename to state (to ensure new fields are maintained!)
struct_assign, state, current_state
state = temporary(current_state)

; If the new session is for a different object, the widgets have to be re-created
; The same applies if Npix or Nlines changes
dummy = where(state.lines gt 0, nlines)

if ((current_file ne state.filename) or (current_ncol ne state.Ncol or $ 
     current_Nrow ne state.Nrow)     or (current_z eq 0) or $
     current_nlines ne nlines) then begin
      if xregistered('kubeviz_spax',/noshow) gt 0 then begin 
        widget_control, /destroy, widgetids.base1     
        widget_control, /destroy, widgetids.base5
      endif
      kubeviz_create, state.filename, scroll=state.scroll
endif

;Test the status of the errmsg logfile
errmsg_file_info = file_info(state.cwdir+'fit_errmsg_'+strmid(state.filename,0,strlen(state.filename)-5)+'.txt')
if errmsg_file_info.exists eq 0 then spawn, 'touch '+state.cwdir+'fit_errmsg_'+strmid(state.filename,0,strlen(state.filename)-5)+'.txt'

;Reset log unit to terminal
state.log_lun = -1

;read spaxels mask if provided
if n_elements(spmask) gt 0 then kubeviz_load_mask, fname=spmask

; Display zoomed map?
widget_Control, widgetids.base3, Map=state.zoommap

; Update the display for the image and spectra
kubeviz_setcolour, state.ctab
kubeviz_medianspec
kubeviz_plotspax
kubeviz_plotinfo
kubeviz_plotspec
kubeviz_plotspeczoom
kubeviz_linefit_typeswitch
kubeviz_linefit_update, /update_userpars

end

;-----------------------------------------------------------------------
pro kubeviz_findcentroid
common kubeviz_state
     
case 1 of
  state.imgmode eq 0:  begin 
      image = (*state.datacube)[*,*, state.wpix]
      weights = (1/(*state.datacube)[*,*, state.wpix])^2
  end   
  state.imgmode eq 1 or state.imgmode eq 2 or state.imgmode eq 3 or state.imgmode eq 4: begin
      image = (*state.img1)
      weights = (1/(*state.nimg1))^2
  end    
  state.imgmode eq 5: begin 
      image =  (*state.datacube)[*,*, state.wpix] - (*state.img1) 
      weights = 1/(((*state.noise)[*,*, state.wpix])^2 + (*state.nimg1)^2 )
  end    
  state.imgmode eq 6 or state.imgmode eq 7 or state.imgmode eq 8 or state.imgmode eq 9: begin 
      image = (*state.img2)
      weights = (1/(*state.nimg2))^2
  end    
  state.imgmode eq 10: begin 
      image =  (*state.datacube)[*,*, state.wpix] - (*state.img2) 
      weights = 1/(((*state.noise)[*,*, state.wpix])^2 + (*state.nimg2)^2 )
  end    
  state.imgmode eq 11: begin 
      image =  (*state.img2) - (*state.img1)
      weights = 1/((*state.nimg1)^2 + (*state.nimg2)^2 )
  end
  state.imgmode eq 12: begin 
      image =  (*state.img1) - (*state.img2)
      weights = 1/((*state.nimg1)^2 + (*state.nimg2)^2 )
  end    
endcase 

bad = where(finite(weights) eq 0, Nbad)
if Nbad gt 0 then weights[bad] = 0 

;initial guesses
baseline =  median(image[where(image ne 0)])
peak = image[state.col, state.row]
halfwidth = 2
centroidX = state.col
centroidY = state.row
PA = 0

yfit = mpfit2dpeak(image, A ,estimates=[baseline, peak, halfwidth, halfwidth, centroidX, centroidY, PA], weights=weights, /TILT)
printf, state.log_lun, '[KUBEVIZ]  Best-fit Gauss parameters: baseline, amplitude, sigmaX, sigmaY, centX, centY:'
printf, state.log_lun, '[KUBEVIZ]', A[0:5]

kubeviz_plotellipse, A[4], A[5], A[2], A[3], A[6]

return 
end

;-----------------------------------------------------------------------
pro kubeviz_optimal_mask

common kubeviz_state
sn_thresh = 5  ;hardcoded

kubeviz_clear_select
     
case state.imgmode of
  0: sn = kubeviz_getsnimage((*state.datacube)[*,*, state.wpix],(*state.noise)[*,*, state.wpix]) 
  1: sn = kubeviz_getsnimage(*state.img1, *state.nimg1)
  2: sn = kubeviz_getsnimage(*state.img1, *state.nimg1)
  3: sn = kubeviz_getsnimage(*state.img1, *state.nimg1)
  4: sn = kubeviz_getsnimage(*state.img1, *state.nimg1)
  5: sn = kubeviz_getsnimage(( (*state.datacube)[*,*, state.wpix] - (*state.img1) ) , ( sqrt( ((*state.noise)[*,*, state.wpix])^2 + (*state.nimg1)^2 ) ))
  6: sn = kubeviz_getsnimage(*state.img2, *state.nimg2)
  7: sn = kubeviz_getsnimage(*state.img2, *state.nimg2)
  8: sn = kubeviz_getsnimage(*state.img2, *state.nimg2)
  9: sn = kubeviz_getsnimage(*state.img2, *state.nimg2)
 10: sn = kubeviz_getsnimage( ( (*state.datacube)[*,*, state.wpix] - (*state.img2) ) , ( sqrt( ((*state.noise)[*,*, state.wpix])^2 + (*state.nimg2)^2 ) ) )
 11: sn = kubeviz_getsnimage( ( (*state.img2) - (*state.img1) ) , ( sqrt( (*state.nimg1)^2 + (*state.nimg2)^2 )  ) )
 11: sn = kubeviz_getsnimage( ( (*state.img1) - (*state.img2) ) , ( sqrt( (*state.nimg1)^2 + (*state.nimg2)^2 )  ) )
endcase 

for col=0,state.Ncol-1 do begin
   for row=0,state.Nrow-1 do begin
     if sn[col,row] gt sn_thresh then (*state.spaxselect)[col, row, state.imask-1] = 1
   endfor
endfor

if total(*state.spaxselect) gt 0 then begin

state.cursormode=2
state.linefit_mode = 1
kubeviz_plotspax

endif else return ;don't do anything if no spaxels reach the sn threshold

end


;-----------------------------------------------------------------------
pro kubeviz_smooth, montecarlo=montecarlo, datacube=datacube
; smooth the datacube (and corresponding noisecube) by smooth pixels.
; or smooth the bootstrapcubes

common kubeviz_state

if n_elements(montecarlo) gt 0 then montecarlo=1 else montecarlo=0
if n_elements(datacube)   gt 0 then datacube=1 else datacube=0

if state.smooth lt 1 then state.smooth=1
if state.specsmooth lt 1 then state.specsmooth=1

if state.smooth eq 1 and state.specsmooth eq 1 then begin
; no smoothing:
    if datacube eq 1 then begin 
        if ptr_valid(state.datacube) eq 1 then begin
	  (*state.datacube) = (*state.indatacube)
          (*state.noisecube) = (*state.innoisecube)
	endif else begin
	  state.datacube = state.indatacube
	  state.noisecube = state.innoisecube
	endelse  
    endif
    if montecarlo eq 1 then begin
       case state.domontecarlo of
       1: if ptr_valid(state.bootstrapcubes[0]) eq 1 then $
          for boot=0,state.Nbootstrap-1 do (*state.bootstrapcubes[boot]) = (*state.inbootstrapcubes[boot]) else $
          for boot=0,state.Nbootstrap-1 do state.bootstrapcubes[boot]    = ptr_new((*state.inbootstrapcubes[boot]))
       2: if ptr_valid(state.mc1cubes[0]) eq 1 then $
          for mc1=0,state.Nmc1-1 do (*state.mc1cubes[mc1]) = (*state.inmc1cubes[mc1]) else $
	  for mc1=0,state.Nmc1-1 do state.mc1cubes[mc1] = ptr_new((*state.inmc1cubes[mc1]))
       else:
       endcase
    endif
    
endif else begin

; find the dimension of the input datacube
    state.Ncol  = (size(*state.indatacube))[1]
    state.Nrow  = (size(*state.indatacube))[2]
    state.Nwpix = (size(*state.indatacube))[3]
    
    Ncol  = state.Ncol
    Nrow  = state.Nrow
    Nwpix = state.Nwpix
    
; Do spatial smoothing first.

; check smooth is a positive integer!
    if fix(state.smooth) ne state.smooth or state.smooth lt 1 then state.smooth = max([abs(fix(state.smooth)),1])
    printf, state.log_lun, '[KUBEVIZ] Smoothing datacubes by '+strtrim(state.smooth,2)+' spatial pixels, and '+strtrim(state.specsmooth,2)+' spectral pixels.'
    
; if smooth is an even number, reduce size of resultant cube by 1 in
; each axis:
    if state.smooth/2. eq fix(state.smooth/2.) then begin
        Ncol--
        Nrow--
        even = 1
    endif else even = 0
    maxoff = (state.smooth-1)/2.
; same for spectral axis:
    if state.specsmooth/2. eq fix(state.specsmooth/2.) then begin
        Nwpix--
        even_spec = 1
    endif else even_spec = 0        
    maxoff_spec = (state.specsmooth-1)/2.

; output cubes:
    if datacube eq 1 then begin 
      smth_datacube = replicate(0.,Ncol,Nrow,Nwpix)
      smth_noisecube = replicate(0.,Ncol,Nrow,Nwpix)
    endif
    if montecarlo eq 1 then begin 
      if state.domontecarlo eq 1 then smth_bootstrapcubes = replicate(0.,Ncol,Nrow,Nwpix,state.Nbootstrap)
      if state.domontecarlo eq 2 then smth_mc1cubes = replicate(0.,Ncol,Nrow,Nwpix,state.Nmc1)
    endif
      
; pixel by pixel (slow, but so far found no vector-based method which works for all smoothing lengths)
; It can take a long time for large cubes. better print some progress stats
    step_count = 0L
    abs_step = 1L*(Nwpix * state.percent_step)/100
    t_start = systime(/seconds)
    for k = 0, Nwpix-1 do begin
        kubeviz_statusline, step_count, abs_step, t_start
	for i = 0, Ncol-1 do begin
            for j = 0, Nrow - 1 do begin
                if i-maxoff gt -1 then i0=i-maxoff else i0=0
                if i+maxoff lt Ncol then i1=i+maxoff else i1=Ncol-1
                if j-maxoff gt -1 then j0=j-maxoff else j0=0
                if j+maxoff lt Nrow then j1=j+maxoff else j1=Nrow-1
                if even eq 1 then begin ; 
                    i0 += 0.5
                    i1 += 0.5
                    j0 += 0.5
                    j1 += 0.5
                endif
                if k-maxoff_spec gt -1 then k0=k-maxoff_spec else k0=0
                if k+maxoff_spec lt Nwpix then k1=k+maxoff_spec else k1=Nwpix-1
                if even_spec eq 1 then begin ; 
                    k0 += 0.5
                    k1 += 0.5
                endif
                if datacube eq 1 then begin
		    ;MF We need to be careful with the NaNs, median is safe, total needs the /nan keyword.
                    smth_datacube[i,j,k] = median((*state.indatacube)[i0:i1,j0:j1,k0:k1],/even)
                    smth_noisecube[i,j,k] = sqrt(total((*state.innoisecube)[i0:i1,j0:j1,k0:k1]^2, /nan))/total(finite((*state.innoisecube)[i0:i1,j0:j1,k0:k1]))
                endif
                if montecarlo eq 1 then begin
                  if state.domontecarlo eq 1 then for boot=0,state.Nbootstrap-1 do smth_bootstrapcubes[i,j,k,boot] = median((*state.inbootstrapcubes[boot])[i0:i1,j0:j1,k0:k1],/even)
                  if state.domontecarlo eq 2 then for mc1=0,state.Nmc1-1 do smth_mc1cubes[i,j,k,mc1] = median((*state.inmc1cubes[mc1])[i0:i1,j0:j1,k0:k1],/even)
                endif
            endfor
        endfor
    endfor
    kubeviz_statusline, /close
; in the even case: embed the cube into one of the original dimensions,
; with NaN in the first row/column:
    if even or even_spec then begin
        if even then begin
            start_spatial = 1
        endif else start_spatial = 0
        if even_spec then begin
            start_spec = 1
        endif else start_spec = 0
        if datacube eq 1 then begin
          smth_datacube_embed = !values.f_nan*(*state.indatacube)
          smth_datacube_embed[start_spatial:state.Ncol-1,start_spatial:state.Nrow-1,start_spec:state.Nwpix-1] = smth_datacube
          smth_datacube = smth_datacube_embed
          smth_noisecube_embed = !values.f_nan*(*state.innoisecube)
          smth_noisecube_embed[start_spatial:state.Ncol-1,start_spatial:state.Nrow-1,start_spec:state.Nwpix-1] = smth_noisecube
          smth_noisecube = smth_noisecube_embed
        endif
        if montecarlo eq 1 then begin
          if state.domontecarlo eq 1 then begin
              smth_bootstrapcubes_embed = replicate(!values.f_nan,state.Ncol,state.Nrow,state.Nwpix, state.Nbootstrap)
              for boot=0,state.Nbootstrap-1 do smth_bootstrapcubes_embed[start_spatial:state.Ncol-1,start_spatial:state.Nrow-1,start_spec:state.Nwpix-1,boot] = smth_bootstrapcubes[*,*,*,boot]
              smth_bootstrapcubes = smth_bootstrapcubes_embed
          endif
          if state.domontecarlo eq 2 then begin
              smth_mc1cubes_embed = replicate(!values.f_nan,state.Ncol,state.Nrow,state.Nwpix, state.Nmc1)
              for mc1=0,state.Nmc1-1 do smth_mc1cubes_embed[start_spatial:state.Ncol-1,start_spatial:state.Nrow-1,start_spec:state.Nwpix-1,mc1] = smth_mc1cubes[*,*,*,mc1]
              smth_mc1cubes = smth_mc1cubes_embed
          endif
        endif
    endif
    
; output, smoothed cube:
    if datacube eq 1 then begin 
        state.datacube  = ptr_new(smth_datacube, /no_copy)
        state.noisecube = ptr_new(smth_noisecube, /no_copy)
    endif
    if montecarlo eq 1 then begin
      if state.domontecarlo eq 1 then for boot=0,state.Nbootstrap-1 do state.bootstrapcubes[boot] = ptr_new(reform(smth_bootstrapcubes[*,*,*,boot]))
      if state.domontecarlo eq 2 then for mc1=0,state.Nmc1-1 do state.mc1cubes[mc1] = ptr_new(reform(smth_mc1cubes[*,*,*,mc1]))
    endif
endelse

; generate/update the bad pixel mask (1 = BAD 0 = GOOD):
if datacube eq 1 then begin 
    (*state.badpixelmask) = 1-(finite(*state.datacube)*finite(*state.noisecube))
    
    ;Additional fix for KMOS data. Where the expmask is 1 the noise is wrong. Assign to badpixel mask
    if state.instr eq 'kmos' then begin
       badkmos = where(*state.noisecube gt 1E5, Nbadkmos)
       if Nbadkmos gt 0 then (*state.badpixelmask)[badkmos] = 1
    endif
    
    (*state.badpixelimg)  = fix(total((*state.badpixelmask), 3 ,/integer)/state.Nwpix)
    ;badpixelvector = where((*state.badpixelmask) eq 1, Nbadpixelvector)
    ;if Nbadpixelvector gt 0 then begin
    ;	(*state.datacube)[badpixelvector] = 0.
    ;	(*state.noisecube)[badpixelvector] = 0.
    ;endif
    state.noise = state.noisecube ; copy into current noise container   
endif
end

;--------------------------------------------------------------------------------
pro kubeviz_readbootstrapcubes, bootstrap_fname, dir=dir, momwindowfile=momwindowfile
; read in bootstrap cubes
common kubeviz_state

if n_elements(dir) eq 0 then dir=''

bootstrap_file = file_info(dir+bootstrap_fname)
if bootstrap_file.exists eq 0 then begin
   printf, state.log_lun, '[WARNING] Unable to locate file: '+dir+bootstrap_fname
   state.domontecarlo = 0
   return
endif

primary = readfits(dir+bootstrap_fname, phdr, ext=0, /silent)
Next = strcompress(sxpar(phdr,'NEXT'),/REMOVE_ALL)
state.Nbootstrap = Next

if Next eq 0 then begin
   printf, state.log_lun, '[WARNING] The bootstrap file: '+dir+bootstrap_fname+' appears to be corrupted. Return.'
   state.domontecarlo = 0
   return
endif

; bootstrapcubes contains an array of pointers to each cube.
state.inbootstrapcubes = replicate(ptr_new(),Next)

printf, state.log_lun, '[KUBEVIZ] Loading bootstrap data cube: '+dir+bootstrap_fname
P = [state.Startcol, state.Startcol+state.Ncol-1, state.Startrow, state.Startrow+state.Nrow-1, state.Startwpix, state.Startwpix+state.Nwpix-1]
for ext=1,Next do begin
    bootstrapcube = kubeviz_fxread_cube(dir+bootstrap_fname, hdr, P[0], P[1], P[2], P[3], P[4], P[5],  ext=ext, /silent)
    bootstrapcube *= state.fluxfac
    state.inbootstrapcubes[ext-1] = ptr_new(bootstrapcube)
endfor

; If the variables momwindowfile and / or gaussinitfile exist, read in the bootstrap extensions, and record them to the state variable.
if n_elements(momwindowfile) gt 0 then begin
    state.mom_windowmap_bootstrap = replicate(ptr_new(),Next)
    for ext=1,Next do begin
    	mom_windowmap = readfits(momwindowfile,ext=ext)
      	mom_Ncol  = (size(mom_windowmap))[1]	
      	mom_Nrow  = (size(mom_windowmap))[2]  
      	mom_Nplane = (size(mom_windowmap))[3]
	if mom_Nplane eq 2 then begin
	  if mom_Ncol - state.Ncol ne 0 or mom_Nrow - state.Nrow ne 0 then begin	   
	    if mom_Ncol - state.Ncol ge 0 and mom_Nrow - state.Nrow ge 0 then mom_windowmap = mom_windowmap[indgen(state.Ncol),indgen(state.Nrow),*] else begin    
	      temp=replicate(-999,state.Ncol,state.Nrow,2)											      	          
	      if mom_Ncol - state.Ncol le 0 and mom_Nrow - state.Nrow le 0 then temp[indgen(mom_Ncol),indgen(mom_Nrow),*] = mom_windowmap	      	          
	      if mom_Ncol - state.Ncol le 0 and mom_Nrow - state.Nrow gt 0 then temp[indgen(mom_Ncol),*,*] = mom_windowmap[*,indgen(state.Nrow),*]    	          
	      if mom_Ncol - state.Ncol gt 0 and mom_Nrow - state.Nrow le 0 then temp[*,indgen(mom_Nrow),*] = mom_windowmap[indgen(state.Ncol),*,*]    	          
	      mom_windowmap = temp														      	          
	    endelse																		  
    	  endif    																		   
	  state.mom_windowmap_bootstrap[ext-1] = ptr_new(mom_windowmap) 											   
	endif
    endfor
endif
        
end

;--------------------------------------------------------------------------------
pro kubeviz_setupmontecarlocubes
; mask for bad pixels, and setup percentiles
common kubeviz_state

; Also setup the choice of percentiles to store if we are bootstrapping:
; first 2 have to be +/-1-sigma. 
med_i = 50.
sig1_i = 50.*0.6826894850
sig2_i = 50.*0.9544997241
montecarlo_percs = [med_i+sig1_i,med_i-sig1_i,med_i+sig2_i,med_i-sig2_i,med_i,10.,90.,25.,75.]
state.montecarlo_percs = ptr_new(montecarlo_percs)
state.Nmontecarlo_percs = N_elements(montecarlo_percs)

badpixelvector = where((*state.badpixelmask) eq 1, Nbadpixelvector)

if state.domontecarlo eq 1 then begin

   printf, state.log_lun, '[KUBEVIZ] Setting up bootstrap data cubes...'
   kubeviz_smooth, /montecarlo
   for boot=0,state.Nbootstrap-1 do begin
       if Nbadpixelvector gt 0 then (*state.bootstrapcubes[boot])[badpixelvector] = 0.
       ; also replace any other NaNs in the bootstrap cube 
       (*state.bootstrapcubes[boot]) = kubeviz_remove_badvalues((*state.bootstrapcubes[boot]))
   endfor
   state.bootstrapnoise = ptr_new(kubeviz_montecarlonoise('cube'))
   state.montecarlocubes = state.bootstrapcubes
   if state.useMonteCarlonoise eq 1 then state.noise = state.bootstrapnoise
   state.Nmontecarlo = state.Nbootstrap
   
   ; Create the results data cubes 
   state.boot_nrescube = ptr_new(kubeviz_create_rescubes('m_line'))
   state.boot_brescube = ptr_new(kubeviz_create_rescubes('m_line'))
   state.boot_crescube = ptr_new(kubeviz_create_rescubes('m_cont'))
   state.boot_mrescube = ptr_new(kubeviz_create_rescubes('m_moment'))

   state.boot_nerrrescube = ptr_new(kubeviz_create_rescubes('m_line_err'))
   state.boot_berrrescube = ptr_new(kubeviz_create_rescubes('m_line_err'))
   state.boot_cerrrescube = ptr_new(kubeviz_create_rescubes('m_cont_err'))
   state.boot_merrrescube = ptr_new(kubeviz_create_rescubes('m_moment_err'))
   
   state.boot_sp_nrescube = ptr_new(kubeviz_create_rescubes('sp_line'))
   state.boot_sp_brescube = ptr_new(kubeviz_create_rescubes('sp_line'))
   state.boot_sp_crescube = ptr_new(kubeviz_create_rescubes('sp_cont'))
   state.boot_sp_mrescube = ptr_new(kubeviz_create_rescubes('sp_moment'))

   state.boot_sp_nerrrescube = ptr_new(kubeviz_create_rescubes('sp_line_err'))
   state.boot_sp_berrrescube = ptr_new(kubeviz_create_rescubes('sp_line_err'))
   state.boot_sp_cerrrescube = ptr_new(kubeviz_create_rescubes('sp_cont_err'))
   state.boot_sp_merrrescube = ptr_new(kubeviz_create_rescubes('sp_moment_err'))

endif

if state.domontecarlo eq 2 then begin

   printf, state.log_lun, '[KUBEVIZ] Setting up MonteCarlo 1 data cubes...'
   kubeviz_smooth, /montecarlo
   for mc1=0,state.Nmc1-1 do begin
       mc1cube = (*state.mc1cubes[mc1])
       if Nbadpixelvector gt 0 then mc1cube[badpixelvector] = 0.
       ; also replace any other NaNs in the mc1 cube and save into state var
       (*state.mc1cubes[mc1]) = kubeviz_remove_badvalues(mc1cube)
   endfor
   state.mc1noise = ptr_new(kubeviz_montecarlonoise('cube'))
   state.montecarlocubes = state.mc1cubes
   if state.useMonteCarlonoise eq 1 then  state.noise = state.mc1noise
   state.Nmontecarlo = state.Nmc1
   
   ; Create the results data cubes 
   state.mc1_nrescube = ptr_new(kubeviz_create_rescubes('m_line'))
   state.mc1_brescube = ptr_new(kubeviz_create_rescubes('m_line'))
   state.mc1_crescube = ptr_new(kubeviz_create_rescubes('m_cont'))
   state.mc1_mrescube = ptr_new(kubeviz_create_rescubes('m_moment'))

   state.mc1_nerrrescube = ptr_new(kubeviz_create_rescubes('m_line_err'))
   state.mc1_berrrescube = ptr_new(kubeviz_create_rescubes('m_line_err'))
   state.mc1_cerrrescube = ptr_new(kubeviz_create_rescubes('m_cont_err'))
   state.mc1_merrrescube = ptr_new(kubeviz_create_rescubes('m_moment_err'))

   state.mc1_sp_nrescube = ptr_new(kubeviz_create_rescubes('sp_line'))
   state.mc1_sp_brescube = ptr_new(kubeviz_create_rescubes('sp_line'))
   state.mc1_sp_crescube = ptr_new(kubeviz_create_rescubes('sp_cont'))
   state.mc1_sp_mrescube = ptr_new(kubeviz_create_rescubes('sp_moment'))

   state.mc1_sp_nerrrescube = ptr_new(kubeviz_create_rescubes('sp_line_err'))
   state.mc1_sp_berrrescube = ptr_new(kubeviz_create_rescubes('sp_line_err'))
   state.mc1_sp_cerrrescube = ptr_new(kubeviz_create_rescubes('sp_cont_err'))
   state.mc1_sp_merrrescube = ptr_new(kubeviz_create_rescubes('sp_moment_err'))

endif

if state.domontecarlo eq 3 then begin

   printf, state.log_lun, '[KUBEVIZ] Setting up MonteCarlo 2 data cubes...'

   for mc2=0,state.Nmc2-1 do begin
       mc2cube = (*state.mc2cubes[mc2])
       if Nbadpixelvector gt 0 then mc2cube[badpixelvector] = 0.
       ; also replace any other NaNs in the mc2 cube and save into state var
       (*state.mc2cubes[mc2]) = kubeviz_remove_badvalues(mc2cube)
    endfor
    state.mc2noise = ptr_new(kubeviz_montecarlonoise('cube'))
    state.montecarlocubes = state.mc2cubes
    if state.useMonteCarlonoise eq 1 then state.noise = state.mc2noise
    state.Nmontecarlo = state.Nmc2
   
    ; Create the results data cubes 
    state.mc2_nrescube = ptr_new(kubeviz_create_rescubes('m_line'))
    state.mc2_brescube = ptr_new(kubeviz_create_rescubes('m_line'))
    state.mc2_crescube = ptr_new(kubeviz_create_rescubes('m_cont'))
    state.mc2_mrescube = ptr_new(kubeviz_create_rescubes('m_moment'))

    state.mc2_nerrrescube = ptr_new(kubeviz_create_rescubes('m_line_err'))
    state.mc2_berrrescube = ptr_new(kubeviz_create_rescubes('m_line_err'))
    state.mc2_cerrrescube = ptr_new(kubeviz_create_rescubes('m_cont_err'))
    state.mc2_merrrescube = ptr_new(kubeviz_create_rescubes('m_moment_err'))

    state.mc2_sp_nrescube = ptr_new(kubeviz_create_rescubes('sp_line'))
    state.mc2_sp_brescube = ptr_new(kubeviz_create_rescubes('sp_line'))
    state.mc2_sp_crescube = ptr_new(kubeviz_create_rescubes('sp_cont'))
    state.mc2_sp_mrescube = ptr_new(kubeviz_create_rescubes('sp_moment'))

    state.mc2_sp_nerrrescube = ptr_new(kubeviz_create_rescubes('sp_line_err'))
    state.mc2_sp_berrrescube = ptr_new(kubeviz_create_rescubes('sp_line_err'))
    state.mc2_sp_cerrrescube = ptr_new(kubeviz_create_rescubes('sp_cont_err'))
    state.mc2_sp_merrrescube = ptr_new(kubeviz_create_rescubes('sp_moment_err'))
   

endif

if state.domontecarlo eq 4 then begin

   printf, state.log_lun, '[KUBEVIZ] Setting up MonteCarlo 3 data cubes...'

   for mc3=0,state.Nmc3-1 do begin
       mc3cube = (*state.mc3cubes[mc3])
       if Nbadpixelvector gt 0 then mc3cube[badpixelvector] = 0.
       ; also replace any other NaNs in the mc3 cube and dave into state var
       (*state.mc3cubes[mc3]) = kubeviz_remove_badvalues(mc3cube)
   endfor
   state.mc3noise = ptr_new(kubeviz_montecarlonoise('cube'))
   state.montecarlocubes = state.mc3cubes
   if state.useMonteCarlonoise eq 1 then  state.noise = state.mc3noise
   state.Nmontecarlo = state.Nmc3
   
   ; Create the results data cubes 
   state.mc3_nrescube = ptr_new(kubeviz_create_rescubes('m_line'))
   state.mc3_brescube = ptr_new(kubeviz_create_rescubes('m_line'))
   state.mc3_crescube = ptr_new(kubeviz_create_rescubes('m_cont'))
   state.mc3_mrescube = ptr_new(kubeviz_create_rescubes('m_moment'))

   state.mc3_nerrrescube = ptr_new(kubeviz_create_rescubes('m_line_err'))
   state.mc3_berrrescube = ptr_new(kubeviz_create_rescubes('m_line_err'))
   state.mc3_cerrrescube = ptr_new(kubeviz_create_rescubes('m_cont_err'))
   state.mc3_merrrescube = ptr_new(kubeviz_create_rescubes('m_moment_err'))

   state.mc3_sp_nrescube = ptr_new(kubeviz_create_rescubes('sp_line'))
   state.mc3_sp_brescube = ptr_new(kubeviz_create_rescubes('sp_line'))
   state.mc3_sp_crescube = ptr_new(kubeviz_create_rescubes('sp_cont'))
   state.mc3_sp_mrescube = ptr_new(kubeviz_create_rescubes('sp_moment'))

   state.mc3_sp_nerrrescube = ptr_new(kubeviz_create_rescubes('sp_line_err'))
   state.mc3_sp_berrrescube = ptr_new(kubeviz_create_rescubes('sp_line_err'))
   state.mc3_sp_cerrrescube = ptr_new(kubeviz_create_rescubes('sp_cont_err'))
   state.mc3_sp_merrrescube = ptr_new(kubeviz_create_rescubes('sp_moment_err'))
      

endif


end


;--------------------------------------------------------------------------------
pro kubeviz_createmc1cubes
common kubeviz_state

printf, state.log_lun, '[KUBEVIZ] Creating MonteCarlo 1 cubes...'

; bootstrapcubes contains an array of pointers to each cube.
state.inmc1cubes = replicate(ptr_new(),state.Nmc1)


for mc1=0,state.Nmc1-1 do begin
   seedm = 353L + mc1*3
   
   mc1cube = (*state.indatacube) + randomn(seedm,(size(*state.indatacube))[1],(size(*state.indatacube))[2],(size(*state.indatacube))[3], /NORMAL) * (*state.innoisecube)
   state.inmc1cubes[mc1] = ptr_new(mc1cube)

endfor

end

;--------------------------------------------------------------------------------
pro kubeviz_createmc2cubes
common kubeviz_state

printf, state.log_lun, '[KUBEVIZ] Creating MonteCarlo 2 cubes...'

; bootstrapcubes contains an array of pointers to each cube.
state.mc2cubes = replicate(ptr_new(),state.Nmc2)

for mc2=0,state.Nmc2-1 do begin
   seedm = 353L + mc2*3
   
   mc2cube = (*state.datacube) +  randomn(seedm, (size(*state.datacube))[1], (size(*state.datacube))[2], (size(*state.datacube))[3], /NORMAL) * (*state.noisecube)
   state.mc2cubes[mc2] = ptr_new(mc2cube)

endfor

end

;--------------------------------------------------------------------------------
pro kubeviz_createmc3cubes
common kubeviz_state

printf, state.log_lun, '[KUBEVIZ] Creating MonteCarlo 3 cubes...'

; Save state of some fitting variables
montecarlo = state.domontecarlo
pndofit = state.pndofit
pcdofit = state.pcdofit
maxoffwb = state.maxwoffb
maxoffwr = state.maxwoffr
continuumfit_minoff = state.continuumfit_minoff 
continuumfit_maxoff = state.continuumfit_maxoff  
continuumfit_minperc = state.continuumfit_minperc   
continuumfit_maxperc = state.continuumfit_maxperc   
continuumfit_order  = state.continuumfit_order   

; Set the default parameters for this quick and dirty fit (all and only narrow lines)
state.domontecarlo = 0
Nlines = N_elements(kubeviz_chooselines((*state.linesets), state.lines, lineset=lineset))
for iline=0, Nlines-1 do begin
  state.pndofit[iline] = 1
  state.pbdofit[iline] = 0
  state.pcdofit[iline] = 1
endfor                       
kubeviz_set_continuumfit_defaults
kubeviz_set_linefitrange_defaults

; Compute residual spectra here
residual = 0.D * (*state.datacube)
for col=0, state.Ncol-1 do begin
  for row=0, state.Nrow-1 do begin
    residual[col,row,*] = kubeviz_get_residual_spec(col,row)
  endfor
endfor

; bootstrapcubes contains an array of pointers to each cube.
state.mc3cubes = replicate(ptr_new(),state.Nmc3)

for mc3=0,state.Nmc3-1 do begin
   seed1 = 424L + mc3*3
   seed2 = 385L + mc3*3
 
   lshift = fix(randomn(seed1,/normal))
   mc3cube = (*state.datacube) +  randomn(seed2, state.Ncol, state.Nrow, state.Nwpix, /NORMAL) * SHIFT(residual,0,0,lshift)
   state.mc3cubes[mc3] = ptr_new(mc3cube)

endfor

;Set original values of the parameters
state.domontecarlo = montecarlo 
state.pndofit =  pndofit 
state.pcdofit =  pcdofit 
state.maxwoffb = maxoffwb
state.maxwoffr = maxoffwr
state.continuumfit_minoff = continuumfit_minoff 
state.continuumfit_maxoff = continuumfit_maxoff 
state.continuumfit_minperc = continuumfit_minperc
state.continuumfit_maxperc = continuumfit_maxperc
state.continuumfit_order = continuumfit_order  

end

;-----------------------------------------------------------------------
pro kubeviz_getdata, dir, datafile, noisedir=noisedir, noisefile=noisefile, ext=ext, noise_ext=noise_ext,$
                     trim=trim, logarithmic=logarithmic, waveunit=waveunit, $
                     bootstrap_fname=bootstrap_fname, momwindowfile=momwindowfile, $
		     gaussinitfile=gaussinitfile

common kubeviz_state

; open a FITS file datacube (x,y,lambda) and obtain its wavelength sol.
; if fnoisename is specified the noise is obtained from a second cube.
; if the datacube has only the primary extension ext=0 otherwise (KMOS) ext=1
; if the noisecube is separate and has one extension noise_ext=0 otherwise (KMOS) noise_ext=2

if n_elements(noisefile) eq 0 then noisefile = datafile   
fits_open, dir+datafile, fcb
fits_open, noisedir+noisefile, fcbnoi

; Figure out some parameters from the primary (or first ext) header
prihead = headfits(dir+datafile, exten=0)
if strlen(state.instr) eq 0 then begin
  case 1 of 
   strtrim(fxpar(prihead, 'INSTRUME'),2) eq 'KMOS' and strtrim(fxpar(prihead, 'VERSION'),2) gt 0 : begin
      state.instr='kmos3d'
      state.vacuum = 1
      printf, state.log_lun, '[KUBEVIZ] Detected instrument: KMOS'
      printf, state.log_lun, '[KUBEVIZ] Detected data product from KMOS3D data release'+string(fxpar(prihead, 'VERSION'), format='(F4.1)')
   end
   strtrim(fxpar(prihead, 'INSTRUME'),2) eq 'KMOS': begin
      state.instr='kmos'
      state.vacuum = 1
      printf, state.log_lun, '[KUBEVIZ] Detected instrument: KMOS'
   end
   strtrim(fxpar(prihead, 'INSTRUME'),2) eq 'VIMOS': begin
      state.instr='vimos'
      printf, state.log_lun, '[KUBEVIZ] Detected instrument: VIMOS'
   end
   strtrim(fxpar(prihead, 'INSTRUME'),2) eq 'SINFONI': begin
      state.instr='sinfoni'
      state.vacuum = 1
      printf, state.log_lun, '[KUBEVIZ] Detected instrument: SINFONI'
   end
   strtrim(fxpar(prihead, 'INSTRUME'),2) eq 'MUSE': begin
      state.instr='muse'
      printf, state.log_lun, '[KUBEVIZ] Detected instrument: MUSE'
   end
   strtrim(fxpar(prihead, 'INSTRUME'),2) eq 'SAMI': begin
       state.instr='sami'
       printf, state.log_lun, '[KUBEVIZ] Detected instrument: SAMI'
   end
   strtrim(fxpar(prihead, 'INSTRMNT'),2) eq 'WiFeS': begin
       state.instr='wifes'
       printf, state.log_lun, '[KUBEVIZ] Detected instrument: WiFeS'
   end
   else: state.instr = strtrim(fxpar(prihead, 'INSTRUME'),2)
  endcase 
endif else state.instr = strlowcase(state.instr)

if state.instr eq 'wifes' then begin
  ext=0
  noise_ext=1
endif
if state.instr eq 'kmos3d' then begin
  ext=1
  noise_ext=2
endif

if n_elements(ext) eq 0 then begin 
   case fcb.nextend of 
     0: ext=0
     1: ext=1
     2: ext=1
     else: begin
        print, '[ ERROR ] FITS data format unrecognized. Please specify the data extension.'
	kubeviz_destroy, /abort
	end
   endcase  
endif     
if n_elements(noise_ext) eq 0 then begin
   case fcbnoi.nextend of 
     0: noise_ext=0
     1: noise_ext=1 
     2: noise_ext=2 
     else: begin
        print, '[ ERROR ] FITS data format unrecognized. Please specify the noise extension.'
	kubeviz_destroy, /abort
	end
   endcase  
endif     

if n_elements(trim) gt 0 then P = kubeviz_decodetrimstr(trim, fcb.axis[*,ext]) else P = replicate(-1,6)
fits_close, fcb
fits_close, fcbnoi

;Check that the input data is a cube (3 dimensions)
datahdr  = headfits(dir+datafile,  ext=ext)
noisehdr = headfits(dir+noisefile, ext=noise_ext)
Ndim_data = fxpar(datahdr, 'NAXIS') 
Ndim_noi  = fxpar(noisehdr, 'NAXIS')

; Read the FITS files
case 1 of
 Ndim_data eq 3 and Ndim_noi eq 3: begin
   ; WARNING: The use of fxread_cube is experimental and allows to save a lot of memory when reading
   ; small chunks of very large datacubes. If you find any bug or you think the loaded data is 
   ; corrupted please report it immediately to the authors and use the following (commented) readfits calls.
   ; cube      = readfits(dir+datafile, hdr, ext=ext, /silent)
   ; noisecube = readfits(noisedir+noisefile, noisehdr, ext=noise_ext, /silent)

   printf, state.log_lun, '[KUBEVIZ] Loading data cube: '+dir+datafile
   start_mem = MEMORY(/CURRENT) 
   cube      = kubeviz_fxread_cube(dir+datafile, hdr, P[0],P[1],P[2],P[3],P[4],p[5], ext=ext, errmsg=dataerr, /silent)
   noisecube = kubeviz_fxread_cube(noisedir+noisefile, noisehdr, P[0],P[1],P[2],P[3],P[4],p[5], ext=noise_ext, errmsg=noiseerr, /silent)
   if dataerr ne ''  ||  noiseerr ne '' then  kubeviz_destroy, /abort 
   end_mem   = (MEMORY(/CURRENT) - start_mem)/(1048576)  
 end
 Ndim_data eq 1 and Ndim_noi eq 1: begin ;Read a single spectrum
  cube      = readfits(dir+datafile, hdr, ext=ext, /silent)
  noisecube = readfits(noisedir+noisefile, noisehdr, ext=noise_ext, /silent)
  cube = reform(cube, 1, 1, n_elements(cube))
  noisecube = reform(noisecube, 1, 1, n_elements(noisecube))
 end
 else: begin
   print, '[ ERROR ] Either the input data or the noise is not a 3D datacube. Abort.'
   kubeviz_destroy, /abort
 end
endcase 


case state.instr of 
   'kmos3d': begin
      
      if state.band eq '' then begin
        state.band = strtrim(kubeviz_fxpar_sp(prihead, 'OBSBAND'),2)
        printf, state.log_lun, '[KUBEVIZ] Detected band: ', state.band 
      endif 

      state.ifu = -1 
      
      if state.fluxfac eq 1 then begin
           state.fluxfac=0.1 ; 1e3 (erg/cm^2)/(W/m^2) * 1e-4 Angstrom/micron
           printf, state.log_lun, '[KUBEVIZ] Input unit: 1E-17 W/m^2/um -> Output unit: 1E-17 erg/cm^2/s/A '
      endif
    end
    'kmos': begin
      
      if state.band eq '' then begin
        state.band = strtrim(kubeviz_fxpar_sp(prihead, 'INS.FILT1.NAME',/HIER),2)
        printf, state.log_lun, '[KUBEVIZ] Detected band: ', state.band 
      endif 

      if state.ifu eq 0 then begin  ;Detect IFU from the data header (not primary header) 
       for j=1,24 do begin 
         if kubeviz_fxpar_sp(hdr,'OCS.ARM'+strtrim(j,2)+'.NAME',/HIER) gt '0' then begin
           state.ifu = j
           printf, state.log_lun, '[KUBEVIZ] Detected ifu: ', strtrim(state.ifu,2)
           break
         endif
       endfor     
      endif else state.ifu = -1 
      
      if state.fluxfac eq 1 then begin
           state.fluxfac=0.1 ; 1e3 (erg/cm^2)/(W/m^2) * 1e-4 Angstrom/micron
           printf, state.log_lun, '[KUBEVIZ] Input unit: W/m^2/um -> Output unit: erg/cm^2/s/A '
      endif
    end
    'sinfoni': begin
      
      if state.band eq '' then begin
        state.band = strtrim(kubeviz_fxpar_sp(prihead, 'INS.FILT1.NAME',/HIER),2)
        printf, state.log_lun, '[KUBEVIZ] Detected band: ', state.band 
      endif 
      
      if state.pixscale eq 0 then begin ; 2 is because the slitlets are rectangular
          state.pixscale = round (2. * 3600. * 1000. * abs(fxpar(prihead, 'CDELT1')))
	  printf, state.log_lun, '[KUBEVIZ] Detected pixel scale: ', strtrim(state.pixscale,2),' mas'
      endif
      	  
      if state.fluxfac eq 1 then begin
           state.fluxfac=0.1 ; 1e3 (erg/cm^2)/(W/m^2) * 1e-4 Angstrom/micron
           printf, state.log_lun, '[KUBEVIZ] Input unit: W/m^2/um -> Output unit: erg/cm^2/s '
      endif
      
    end
    'muse': begin
      
      if Ndim_data eq 3 then printf, state.log_lun, '[KUBEVIZ] Memory required in reading the data: '+kubeviz_str(end_mem)+' Mb'
      
      state.noiseisvar = 1 ;The STAT extension is the variance
      
      if state.band eq '' then begin
        state.band = strtrim(kubeviz_fxpar_sp(prihead, 'INS.MODE',/HIER),2)
        printf, state.log_lun, '[KUBEVIZ] Detected mode: ', state.band 
      endif 
      
      if state.fluxfac eq 1 then begin
           state.fluxfac=1 ;Instrument works in Angstroms, flux in units of 1E-20 erg/cm/s
           printf, state.log_lun, '[KUBEVIZ] Output unit: 1E-20 erg/cm^2/s '
      endif
      
    end
    'sami': begin
        state.noiseisvar = 1 ;The STAT extension is the variance
    end
    'wifes': begin
         
	 state.noiseisvar = 1 ;The STAT extension is the variance
	 
	 if state.band eq '' then begin
           state.band = strtrim(kubeviz_fxpar_sp(prihead, 'WIFESARM'),2)
           printf, state.log_lun, '[KUBEVIZ] Detected arm: ', state.band 
         endif 
	 
	 if state.fluxfac eq 1 then begin
           state.fluxfac=1e18 ;Instrument works in Angstroms, flux in units of 1E-20 erg/cm/s
           printf, state.log_lun, '[KUBEVIZ] Output unit: 1E-18 erg/cm^2/s '
         endif
    end

    else: printf, state.log_lun, '[KUBEVIZ] Instrument: '+state.instr+' requires manual input of band and fluxfac parameters, if required.'
endcase

if state.noiseisvar eq 1 then noisecube = sqrt(TEMPORARY(noisecube)) ;Take the sqrt if the noise extension is the variance
state.indatacube  = ptr_new(temporary(cube) * state.fluxfac)
state.innoisecube = ptr_new(temporary(noisecube) * state.fluxfac)
state.indatahead  = ptr_new(hdr)
state.innoisehead = ptr_new(noisehdr)
if n_elements(prihead) gt 0 and ext ne 0 then state.inprihead   = ptr_new(prihead) ;Save primary header if it exists

; find its dimension 
state.Ncol  = (size(*state.indatacube))[1]
state.Nrow  = (size(*state.indatacube))[2]
state.Nwpix = (size(*state.indatacube))[3]

; setup a bad pixel mask for bad spaxel values:
state.badpixelmask = ptr_new(replicate(0B, state.Ncol, state.Nrow, state.Nwpix))
state.badpixelimg  = ptr_new(replicate(0B, state.Ncol, state.Nrow))

; Smooth the data as required:
kubeviz_smooth, /datacube

; determine units of wavelength to convert to Angstroms:
unit    = strcompress(sxpar(hdr,'CUNIT3'),/REMOVE_ALL)
; use keyword in header if available:
case unit of
    'Angstrom': waveunit='ANGSTROMS'
    'um':       waveunit='MICRONS'
    'MICRON':   waveunit='MICRONS'
    else:
endcase

; conversion factor:
wave_conv = 1.
if N_Elements(waveunit) gt 0 then begin
    case strupcase(waveunit) of
        'ANGSTROMS': wave_conv = 1.
        'MICRONS': wave_conv = 1.e4
        'NM': wave_conv = 10.
        else: print, 'No unit '+waveunit+'defined'
    endcase
    if state.debug eq 1 then print, '[DEBUG] Units: '+waveunit+'; conversion '+string(wave_conv)
endif

;Determine the wavelength solution
if Ndim_data eq 1 then waveaxis = '1' else waveaxis='3'
if state.debug eq 1 then print, '[DEBUG] Waveaxis :'+waveaxis

state.pix0    = sxpar(hdr,'CRPIX'+waveaxis)
state.lambda0 = sxpar(hdr,'CRVAL'+waveaxis)*wave_conv
if sxpar(hdr,'CDELT'+waveaxis) gt 0 $   ;Try different keywords for the wave resolution
   then state.dlambda = sxpar(hdr,'CDELT'+waveaxis)*wave_conv $
   else state.dlambda = sxpar(hdr,'CD'+waveaxis+'_'+waveaxis)*wave_conv

if state.lambda0 eq 0. then begin
    printf, state.log_lun, '[WARNING] Cannot find wavelength solution in header.'
    printf, state.log_lun, '[WARNING] Will use pixel number instead of lambda.'
    state.wave = ptr_new( findgen(state.Nwpix) )
endif else begin
    if keyword_set(logarithmic) $ 
      then state.wave = ptr_new(exp(state.lambda0 + (state.dlambda)*((findgen(state.Nwpix)+1)-state.pix0+state.Startwpix) )) $
      else state.wave = ptr_new(state.lambda0 + state.dlambda*((findgen(state.Nwpix)+1)-state.pix0+state.Startwpix) )
endelse

if state.transpose ne 0 then begin
    printf, state.log_lun, '[KUBEVIZ] Transposing the datacube...'
    for k = 0, state.nwpix-1 do begin
        (*state.datacube)[*,*,k] = transpose( (*state.datacube)[*,*,k] )
        (*state.noisecube)[*,*,k] = transpose( (*state.noisecube)[*,*,k] )        
    endfor
endif

;Copy the noisecube into the current noise pointer
state.noise = state.noisecube

; If the variables momwindowfile and / or gaussinitfile exist, read these in and record them to the state variable.
if n_elements(momwindowfile) gt 0 then begin
  mom_windowmap = readfits(momwindowfile,ext=0)
  mom_Ncol  = (size(mom_windowmap))[1]	
  mom_Nrow  = (size(mom_windowmap))[2]
  mom_Nplane = (size(mom_windowmap))[3]
  if mom_Nplane ne 2 then printf, state.log_lun, '[WARNING] File '+momwindowfile+' contains wrong number of image planes' else begin
    if mom_Ncol - state.Ncol ne 0 or mom_Nrow - state.Nrow ne 0 then begin
      printf, state.log_lun, '[WARNING] File '+momwindowfile+' image size does not correspond to spatial cube size! Assuming [0,0] is rigid'
      if mom_Ncol - state.Ncol ge 0 and mom_Nrow - state.Nrow ge 0 then mom_windowmap = mom_windowmap[indgen(state.Ncol),indgen(state.Nrow),*] else begin
	temp=replicate(-999,state.Ncol,state.Nrow,2)
	if mom_Ncol - state.Ncol le 0 and mom_Nrow - state.Nrow le 0 then temp[indgen(mom_Ncol),indgen(mom_Nrow),*] = mom_windowmap
	if mom_Ncol - state.Ncol le 0 and mom_Nrow - state.Nrow gt 0 then temp[indgen(mom_Ncol),*,*] = mom_windowmap[*,indgen(state.Nrow),*]
	if mom_Ncol - state.Ncol gt 0 and mom_Nrow - state.Nrow le 0 then temp[*,indgen(mom_Nrow),*] = mom_windowmap[indgen(state.Ncol),*,*]
	mom_windowmap = temp
      endelse
    endif
    state.mom_windowmap = ptr_new(mom_windowmap)
  endelse  
endif	

if n_elements(gaussinitfile) gt 0 then begin
  gauss_initmap = readfits(gaussinitfile,ext=0)
  gauss_Ncol  = (size(gauss_initmap))[1]	
  gauss_Nrow  = (size(gauss_initmap))[2]
  gauss_Nplane = (size(gauss_initmap))[3]
  ; [UPDATE!] WARNING: HARDCODED FOR 5 PLANES!!!
  if gauss_Nplane ne 5 then printf, state.log_lun, '[WARNING] File '+gaussinitfile+' contains wrong number of image planes' else begin
    if gauss_Ncol - state.Ncol ne 0 or gauss_Nrow - state.Nrow ne 0 then begin
      printf, state.log_lun, '[WARNING] File '+momwindowfile+' image size does not correspond to spatial cube size! Assuming [0,0] is rigid'
      if gauss_Ncol - state.Ncol ge 0 and gauss_Nrow - state.Nrow ge 0 then gauss_initmap = gauss_initmap[indgen(state.Ncol),indgen(state.Nrow),*] else begin
      	temp=replicate(-999,state.Ncol,state.Nrow,2)
	if gauss_Ncol - state.Ncol le 0 and gauss_Nrow - state.Nrow le 0 then temp[indgen(gauss_Ncol),indgen(gauss_Nrow),*] = gauss_initmap
	if gauss_Ncol - state.Ncol le 0 and gauss_Nrow - state.Nrow gt 0 then temp[indgen(gauss_Ncol),*,*] = gauss_initmap[*,indgen(state.Nrow),*]
	if gauss_Ncol - state.Ncol gt 0 and gauss_Nrow - state.Nrow le 0 then temp[*,indgen(gauss_Nrow),*] = gauss_initmap[indgen(state.Ncol),*,*]
	gauss_initmap = temp
      endelse
    endif
    state.gauss_initmap = ptr_new(gauss_initmap)
  endelse  
endif

end

;-----------------------------------------------------------------------
pro kubeviz_plotcrosshair
; plot a crosshair in the spaxel viewer
common kubeviz_state
common kubeviz_winsize

if state.zoomfac ge 10 then thick = 2 else thick = 1
case state.ctab of 
  0 : col = 252
  1 : col = 249
  3 : col = 252
  5 : col = 252
  13: col = 255
endcase

xspaxsize = (1.D*state.zoomfac*state.Ncol) / (state.spaxrange[1]-state.spaxrange[0]+1)
yspaxsize = (1.D*state.zoomfac*state.Nrow) / (state.spaxrange[3]-state.spaxrange[2]+1)
x0 = (state.col + 0.5 - state.spaxrange[0]) * xspaxsize
y0 = (state.row + 0.5 - state.spaxrange[2]) * yspaxsize
plots, [x0, x0], [0, winsize.ywin1], /device, thick=thick, color=col
plots, [0, winsize.xwin1], [y0, y0], /device, thick=thick, color=col

end

;-----------------------------------------------------------------------
pro kubeviz_plotellipse, xc, yc, rx, ry, paX
; plot a crosshair in the spaxel viewer
common kubeviz_state
common kubeviz_winsize
common kubeviz_widgetids

wset, widgetids.wid1

if state.zoomfac ge 10 then thick = 3 else thick = 2
case state.ctab of 
  0 : col = 253
  1 : col = 253
  3 : col = 252
  5 : col = 253
  13: col = 255
endcase

x0 = (xc) * state.zoomfac ;winsize.xwin1 / state.Ncol
y0 = (yc) * state.zoomfac ;winsize.ywin1 / state.Nrow
rmax = max([rx,ry]) * winsize.xwin1 / state.Ncol
rmin = min([rx,ry]) * winsize.xwin1 / state.Ncol
PAellipse = 180-(paX*!RADEG)

tvellipse, rmax,rmin, x0,y0, /device,  thick=thick, color=col

end

;-----------------------------------------------------------------------
pro kubeviz_markspaxel
; highlight a selected spaxel by drawing an open box around it.
common kubeviz_state

if state.cursormode lt 2 then return
mask = (*state.spaxselect)[*,*, state.imask-1]
ok = where(mask gt 0, Nok ) 

if Nok eq 0 then return
pos = array_indices(mask,ok)
for i = 0L, Nok-1L do kubeviz_plotbox, pos[0,i], pos[1,i]

end

;-----------------------------------------------------------------------
pro kubeviz_markwrange, spec, wavsel
; highlight the selected wavelength range
common kubeviz_state

if wavsel eq 1 then begin
    orient =  45
    xmin = (*state.wave)[state.wavrange1[0]]
    xmax = (*state.wave)[state.wavrange1[1]]
    klr = 251
endif

if wavsel eq 2 then begin
    orient = -45
    xmin = (*state.wave)[state.wavrange2[0]]
    xmax = (*state.wave)[state.wavrange2[1]]
    klr = 250
endif
x0 = [xmin, xmin, xmax, xmax]

ymin = min([0,state.zmin_spec,min(spec)])
ymax = max([state.zmax_spec,2.0 * max(spec)])
if finite(min(spec)) and finite(max(spec)) then y0 = [ymin, ymax, ymax, ymin] $
    	    	    	    	    	   else y0 = [-10,10,10,-10] 


polyfill, x0, y0, color=klr, noclip=0  ; use solid fill 

end

;-----------------------------------------------------------------------
pro kubeviz_plotbox, x0, y0
; over plot a box with the dimensions of a spaxel
common kubeviz_state
common kubeviz_winsize

if state.zoomfac ge 10 then thick = 2 else thick = 1
case state.ctab of 
  0 : col = 253
  1 : col = 253
  3 : col = 252
  5 : col = 253
  13: col = 255
endcase

x0 -= state.spaxrange[0]
y0 -= state.spaxrange[2]

xspaxsize = (1.D*state.zoomfac*state.Ncol) / (state.spaxrange[1]-state.spaxrange[0]+1)
yspaxsize = (1.D*state.zoomfac*state.Nrow) / (state.spaxrange[3]-state.spaxrange[2]+1)

xbox = fltarr(6) & ybox = fltarr(6)

xbox[0] = ( x0 )       * xspaxsize   
xbox[1] = ( x0 + 1.0)  * xspaxsize
xbox[2] = ( x0 + 1.0 ) * xspaxsize
xbox[3] = ( x0 )       * xspaxsize
xbox[4] = xbox[0]
xbox[5] = xbox[2]

ybox[0] = ( y0 )       * yspaxsize   
ybox[1] = ( y0 )       * yspaxsize
ybox[2] = ( y0 + 1.0 ) * yspaxsize
ybox[3] = ( y0 + 1.0 ) * yspaxsize
ybox[4] = ybox[0]
ybox[5] = ybox[2]
        
plots, xbox, ybox, /device, thick=thick, color=col

end

;-----------------------------------------------------------------------
pro kubeviz_clear_select
; clear the array that enumerates the selected spaxels
common kubeviz_state

(*state.spaxselect)[*,*,state.imask-1] = 0
(*state.medspec)  = fltarr(state.Nwpix)
(*state.nmedspec) = fltarr(state.Nwpix)

end

;-----------------------------------------------------------------------
pro kubeviz_save_mask
; save the array that enumerates the selected spaxels to a fits file
common kubeviz_state

fname = dialog_pickfile(filter='*.fits', path=state.cwdir, /write, get_path=path, /OVERWRITE_PROMPT)

len = strlen(fname) - strlen(path)
if len eq 0 then return ; don't do anything if no filename is given

writefits, fname, (*state.spaxselect)
printf, state.log_lun, '[KUBEVIZ] Wrote selected spaxel positions to file: ', fname

end

;-----------------------------------------------------------------------
pro kubeviz_load_mask, fname=fname
; load the array that enumerates the selected spaxels from a fits file
; OVERWRITES current set of masks
common kubeviz_state

if n_elements(fname) eq 0 then begin 
 fname = dialog_pickfile(filter='*.fits', /read, /must_exist, get_path=path)
 len = strlen(fname) - strlen(path)
 if len eq 0 then return ; don't do anything if no filename is given
endif

printf, state.log_lun, '[KUBEVIZ] Reading spaxels mask file: '+fname
(*state.spaxselect) = readfits(fname, /silent)
sz = size((*state.spaxselect))
if sz[0] eq 3 then state.Nmask = sz[3] else state.Nmask = 1

state.imask = 1
state.linefit_mode = 1
state.specmode   = 1 ; SUM
state.cursormode = 2
kubeviz_medianspec

end

;----------------------------------------------------------------------------
pro kubeviz_newmask
; create new spaxel mask: add to end of mask array
common kubeviz_state

state.Nmask++
; goto new mask:
state.imask = state.Nmask
; initialize spaxel array:
spaxselect = ptr_new( intarr(state.Ncol , state.Nrow, state.Nmask) )
for i=1L, state.Nmask-1 do begin
    (*spaxselect)[*,*, i-1] = (*state.spaxselect)[*,*,i-1]
endfor
heap_free, state.spaxselect 
state.spaxselect = spaxselect

end

;-----------------------------------------------------------------------------
pro kubeviz_deletemask
; delete current mask from array, move to previous mask if it exists
common kubeviz_state

if state.Nmask eq 1 then printf, state.log_lun, '[WARNING] Last mask: Do not delete (can clear spaxels instead).' else begin
    ;reform necessary to keep the second dimension where Nmask=1
    spaxselect = ptr_new( reform(intarr(state.Ncol , state.Nrow, state.Nmask-1),[state.Ncol , state.Nrow, state.Nmask-1]) ) 
    for i=1, state.Nmask-1 do begin
        if i lt state.imask then (*spaxselect)[*,*, i-1] = (*state.spaxselect)[*,*,i-1]
        if i gt state.imask then (*spaxselect)[*,*, i-1] = (*state.spaxselect)[*,*,i]
    endfor
    state.Nmask--
    state.spaxselect = spaxselect
    if state.imask gt 1 then state.imask-- 
endelse

end

;-----------------------------------------------------------------------
pro kubeviz_medianspec, med=med, sum=sum, wavg=wavg, optimal=optimal, all=all, redsh=redsh
; median (or sum) the spectrum over the selected spaxels
; optimal uses a modified version of the robertson method

common kubeviz_state

; reformat mask in 2d:
case 1 of 
   N_elements(all) gt 0 : mask = replicate(1,state.Ncol,state.Nrow)
   N_elements(redsh) gt 0 : begin 
     maskcen = replicate(0,state.Ncol,state.Nrow)
     maskcen[fix(state.Ncol/3):fix(state.Ncol/3*2),fix(state.Nrow/3):fix(state.Nrow/3*2)] = 1
     mask = maskcen
   end   
   else:  begin
     
     mask = (*state.spaxselect)[*,*,state.imask-1]
     
     ;Now split the mask into a binary mask and a mask of weights (if weights are saved in the mask)
     ok = where(mask gt 0, Nok)
     ones = where(mask eq 1, Nones) 
     if Nones ne Nok then begin 
       wswitch = 'wmask'
       weights = mask
       if Nok gt 0 then mask[ok] = 1
     endif else wswitch = 'noise'
   
   end
endcase   

case 1 of 
  keyword_set(sum)     : mode = 0
  keyword_set(med)     : mode = 1
  keyword_set(wavg)    : mode = 2
  keyword_set(optimal) : mode = 3
  else: begin
   case state.specmode of
     1: mode=0
     2: mode=1
     3: mode=2
     4: mode=1
     5: mode=3
     else: return
   endcase
  end
endcase  

if total(mask) gt 0 then begin
    
    if mode eq 3 then begin
    
     case state.imgmode of
       0: image  = (*state.datacube)[*,*, state.wpix]
       1: image = (*state.img1)
       2: image = (*state.img1)
       3: image = (*state.img1)
       4: image = (*state.img1)
       5: image = ( (*state.datacube)[*,*, state.wpix] - (*state.img1) ) 
       6: image = (*state.img2)
       7: image = (*state.img2)
       8: image = (*state.img2)
       9: image = (*state.img2)
      10: image = ( (*state.datacube)[*,*, state.wpix] - (*state.img2) )
      11: image = (*state.img2) - (*state.img1) 
      12: image = (*state.img1) - (*state.img2) 
     endcase 

     if min(mask*image, /nan) lt 0 then begin
      printf, state.log_lun, '[WARNING] The current mask includes spaxels with negative values in the selected image.'
      printf, state.log_lun, '[WARNING] The optimal method is not reliable. Sum method will be used instead.'
      mode = 1
      state.specmode=1
     endif
    endif
    
    
    case mode of 
     0: begin ;SUM
       
        sumspec = fltarr(state.Nwpix)
        for ik=0, state.Nwpix-1 do sumspec[ik] = total( mask * (*state.datacube)[*,*,ik])

        if state.domontecarlo gt 0 then begin
            state.medspec_montecarlo = replicate(ptr_new(),state.Nmontecarlo)
            for boot=0,state.Nmontecarlo-1 do begin
                sumspec_montecarlo = replicate(0.,state.Nwpix)
                for ik=0, state.Nwpix-1 do sumspec_montecarlo[ik] = total( mask * (*state.montecarlocubes[boot])[*,*,ik]) 
                state.medspec_montecarlo[boot] = ptr_new(sumspec_montecarlo)
            endfor
        endif
        
	; noise computation (still slow)
        if (state.domontecarlo eq 0) or (state.useMonteCarlonoise eq 0) then begin
          nsumspecsq = replicate(0.D,state.Nwpix) 
          for ik=0, state.Nwpix-1 do nsumspecsq[ik] = total( mask * (*state.noise)[*,*,ik]^2 )
          nsumspec = sqrt(nsumspecsq)
        endif else nsumspec = kubeviz_montecarlonoise('spec')
        
        *state.medspec = sumspec
        *state.nmedspec = nsumspec

    end
    1: begin ; MEDIAN
        
        medspec = fltarr(state.Nwpix)
        nmedspecsq = replicate(0.D,state.Nwpix) 
        
        for ik=0, state.Nwpix-1 do medspec[ik] = median( (reform((*state.datacube)[*,*,ik],state.Ncol*state.Nrow))[where(mask eq 1)] ) 
      
        if state.domontecarlo gt 0 then begin
            state.medspec_montecarlo = replicate(ptr_new(),state.Nmontecarlo)
            for boot =0,state.Nmontecarlo-1 do begin
                medspec_montecarlo = replicate(0.,state.Nwpix)
                for ik=0, state.Nwpix-1 do medspec_montecarlo[ik] =  median ( (reform((*state.montecarlocubes[boot])[*,*,ik],state.Ncol*state.Nrow))[where(mask eq 1)] )
                state.medspec_montecarlo[boot] = ptr_new(medspec_montecarlo)
            endfor
        endif
                
        if (state.domontecarlo eq 0) or (state.useMonteCarlonoise eq 0) then begin
          for ik=0, state.Nwpix-1 do nmedspecsq[ik] = total( mask * (*state.noise)[*,*,ik]^2 )
          nmedspec = sqrt(nmedspecsq)/total(mask)
        endif else nmedspec = kubeviz_montecarlonoise('spec')
        
        *state.medspec = medspec
        *state.nmedspec = nmedspec
    end
    
    2: begin ; WEIGHTED AVERAGE
        avgspec = fltarr(state.Nwpix)
        for ik=0, state.Nwpix-1 do begin
	   if wswitch eq 'noise' then weights = mask / (*state.noise)[*,*,ik]^2
	   avgspec[ik] = total( weights * (*state.datacube)[*,*,ik], /nan )/ total( weights, /nan)
        endfor
	
        if state.domontecarlo gt 0 then begin
            state.medspec_montecarlo = replicate(ptr_new(),state.Nmontecarlo)
            for boot=0,state.Nmontecarlo-1 do begin
                avgspec_montecarlo = replicate(0.,state.Nwpix)
                for ik=0, state.Nwpix-1 do begin 
		    if wswitch eq 'noise' then weights = mask / (*state.noise)[*,*,ik]^2
		    avgspec_montecarlo[ik] = total( weights * (*state.montecarlocubes[boot])[*,*,ik], /nan) / total( weights, /nan)
                endfor
		state.medspec_montecarlo[boot] = ptr_new(avgspec_montecarlo)
            endfor
        endif
    
	; noise computation (still slow)
        if (state.domontecarlo eq 0) or (state.useMonteCarlonoise eq 0) then begin
          navgspecsq = replicate(0.D,state.Nwpix) 
          for ik=0, state.Nwpix-1 do begin
	    if wswitch eq 'noise' then navgspecsq[ik] = 1./ total(  mask * ((*state.noise)[*,*,ik])^(-2), /nan ) 
	    if wswitch eq 'wmask' then navgspecsq[ik] = total( weights * (*state.noise)[*,*,ik]^2 ,/nan)/total(weights, /nan)
         endfor
	  navgspec = sqrt(navgspecsq)
        endif else navgspec = kubeviz_montecarlonoise('spec')
        
        *state.medspec = avgspec
        *state.nmedspec = navgspec
    end
    
    3: begin ;OPTIMAL
        optspec = fltarr(state.Nwpix)
        noptspec = fltarr(state.Nwpix)
        
        frac_image = mask *image / total( mask * image, /nan)
        
        for ik=0, state.Nwpix-1 do begin
               scaled_noi = (*state.noise)[*,*,ik]/frac_image
               scaled_img = (*state.datacube)[*,*,ik]/frac_image
               optspec[ik] = total( mask * scaled_img * (1./scaled_noi)^2, /nan) / total( mask * (1./scaled_noi)^2, /nan)
               noptspec[ik] = sqrt(1./total( (mask * scaled_noi)^(-2), /nan ))
        endfor
        optspec = kubeviz_remove_badvalues(optspec)
        noptspec = kubeviz_remove_badvalues(noptspec)
                
        if (state.domontecarlo gt 0) and (state.useMonteCarlonoise gt 0) then begin
            state.medspec_montecarlo = replicate(ptr_new(),state.Nmontecarlo)
            for boot =0,state.Nmontecarlo-1 do begin
                medspec_montecarlo = replicate(0.,state.Nwpix)
                for ik=0, state.Nwpix-1 do begin
                    scaled_noi = (*state.noise)[*,*,ik]/frac_image
                    scaled_img = (*state.montecarlocubes[boot])[*,*,ik]/frac_image
                    medspec_montecarlo[ik] = total( mask * scaled_img * (1./scaled_noi)^2, /nan) / total( mask * (1./scaled_noi)^2, /nan)
                endfor
                state.medspec_montecarlo[boot] = ptr_new(kubeviz_remove_badvalues(medspec_montecarlo))
            endfor
	    noptspec = kubeviz_montecarlonoise('spec')
        endif
        *state.medspec = optspec
        *state.nmedspec = noptspec
     end
    endcase

endif else begin
    *state.medspec = fltarr(state.Nwpix)
    *state.nmedspec = fltarr(state.Nwpix)
endelse

end

;-----------------------------------------------------------------------
pro kubeviz_savespec
; save the spectrum and noise vectors to ext 1 and 2 of an output file.
common kubeviz_state

fname = dialog_pickfile(filter='*.1dspec.fits', path=state.cwdir, /write, get_path=path, /OVERWRITE_PROMPT)

len = strlen(fname) - strlen(path)
if len eq 0 then return ; don't do anything if no filename is given

kubeviz_getspec, spec, nspec, zspec

; put spec in ext 0 and prepare the header
mkhdr, hdr_1d, spec, /extend
sxaddpar, hdr_1d, 'EXT_TY', 'DATA'                           , 'The type of the spectrum in this extension'
sxaddpar, hdr_1d, 'CRVAL1', (*state.wave)[0]                 , '[um] Wavelength at ref. pixel'
sxaddpar, hdr_1d, 'CRPIX1', 1                                , '[pix] Reference pixel in x'
sxaddpar, hdr_1d, 'CDELT1', state.dlambda                    , '[um] Spectral resolution'
sxaddpar, hdr_1d, 'CTYPE1', sxpar(*state.indatahead,'CTYPE3'), 'Coordinate system of x-axis'
sxaddpar, hdr_1d, 'CUNIT1', sxpar(*state.indatahead,'CUNIT3'), 'Unit of x-axis'
case state.specmode of
 0: begin
    sxaddpar, hdr_1d, 'SPTYPE', 'SPAXEL'                     , 'Single spaxel spectrum from the cube'
    sxaddpar, hdr_1d, 'X_ORIG', state.col+1                  , 'X position of the spectrum in the cube (1 based)'
    sxaddpar, hdr_1d, 'Y_ORIG', state.row+1                  , 'Y position of the spectrum in the cube (1 based)'
 end   
 1: sxaddpar, hdr_1d, 'SPTYPE', 'SUM'                        , 'Sum spectrum from the spaxel mask'
 2: sxaddpar, hdr_1d, 'SPTYPE', 'MEDIAN'                     , 'Median spectrum from the spaxel mask'
 3: sxaddpar, hdr_1d, 'SPTYPE', 'WAVERAGE'                   , 'W.average spectrum from the spaxel mask'
 4: begin
    sxaddpar, hdr_1d, 'SPTYPE', 'MEDSUB'                     , 'Median subtracted spaxel spectrum'
    sxaddpar, hdr_1d, 'X_ORIG', state.col+1                  , 'X position of the spectrum in the cube (1 based)'
    sxaddpar, hdr_1d, 'Y_ORIG', state.row+1                  , 'Y position of the spectrum in the cube (1 based)'
 end   
 5: sxaddpar, hdr_1d, 'SPTYPE', 'OPTIMAL'                     , 'Optimal extraction from the spaxel mask'
endcase
mwrfits, spec, fname, hdr_1d, /CREATE

sxaddpar, hdr_1d, 'EXT_TY', 'NOISE'                          , 'The type of the spectrum in this extension'
mwrfits, nspec, fname, hdr_1d, /silent

printf, state.log_lun, '[KUBEVIZ] Wrote spectrum to file: ', fname

end

;-----------------------------------------------------------------------
pro kubeviz_savecube
; save the spectrum and noise vectors to ext 1 and 2 of an output file.
common kubeviz_state

fname = dialog_pickfile(filter='*.fits', path=state.cwdir, /write, get_path=path, /OVERWRITE_PROMPT)

len = strlen(fname) - strlen(path)
if len eq 0 then return ; don't do anything if no filename is given

; Prepare a dummy primary header and the astro solution
mkhdr, prihdr, '', /extend
xyad, (*state.indatahead), state.Startcol, state.Startrow, RA, DEC

datahdr = (*state.indatahead)
noisehdr = (*state.innoisehead)
sxaddpar, datahdr, 'XTENSION', 'IMAGE'			       ,  'IMAGE extension', before='BITPIX'
sxaddpar, datahdr, 'EXTNAME' , 'DATA'			       ,  'This extension contains data values'
sxaddpar, datahdr, 'DATACUBE', state.filename			, 'Name of the original datacube'
sxaddpar, datahdr, 'INSTRUME', strupcase(state.instr)		, 'Instrument name'
sxaddpar, datahdr, 'BAND'    , strupcase(state.band)		, 'Band/Filter name'
sxaddpar, datahdr, 'SPATSMTH', state.smooth			, 'Spatial smoothing applied'
sxaddpar, datahdr, 'SPECSMTH', state.specsmooth 		, 'Spectral smoothing applied'
sxaddpar, datahdr, 'CTYPE1'  , sxpar(*state.indatahead,'CTYPE1'), 'TAN projection used'
sxaddpar, datahdr, 'CTYPE2'  , sxpar(*state.indatahead,'CTYPE2'), 'TAN projection used'
sxaddpar, datahdr, 'CTYPE3'  , sxpar(*state.indatahead,'CTYPE3'), 'Coordinate system of z-axis'
sxaddpar, datahdr, 'CRPIX1'  , 1				, '[pix] Reference pixel in x'         
sxaddpar, datahdr, 'CRPIX2'  , 1				, '[pix] Reference pixel in y'         
sxaddpar, datahdr, 'CRPIX3'  , 1				, '[pix] Reference pixel in z'
sxaddpar, datahdr, 'CRVAL1'  , RA				, '[deg] RA at ref. pixel'
sxaddpar, datahdr, 'CRVAL2'  , DEC				, '[deg] DEC at ref. pixel'
sxaddpar, datahdr, 'CRVAL3'  , (*state.wave)[0] 		, '[um] Wavelength at ref. pixel'
sxaddpar, datahdr, 'CD1_1'   , sxpar(*state.indatahead,'CD1_1') , '[] x-component of East'
sxaddpar, datahdr, 'CD2_1'   , sxpar(*state.indatahead,'CD2_1') , '[] x-component of North'
sxaddpar, datahdr, 'CD1_2'   , sxpar(*state.indatahead,'CD1_2') , '[] y-component of East'
sxaddpar, datahdr, 'CD2_2'   , sxpar(*state.indatahead,'CD2_2') , '[] y-component of North'
sxaddpar, datahdr, 'CDELT1'  , sxpar(*state.indatahead,'CDELT1'), '[deg] Pixel resolution in x'
sxaddpar, datahdr, 'CDELT2'  , sxpar(*state.indatahead,'CDELT2'), '[deg] Pixel resolution in y'
sxaddpar, datahdr, 'CDELT3'  , state.dlambda			, '[um] Spectral resolution'
sxaddpar, datahdr, 'CUNIT1'  , sxpar(*state.indatahead,'CUNIT1'), 'Unit of x-axis'
sxaddpar, datahdr, 'CUNIT2'  , sxpar(*state.indatahead,'CUNIT2'), 'Unit of y-axis'
sxaddpar, datahdr, 'CUNIT3'  , sxpar(*state.indatahead,'CUNIT3'), 'Unit of z-axis'
sxaddpar, datahdr, 'STARTCOL', state.Startcol+1 		, 'Index of the first col in the original frame'
sxaddpar, datahdr, 'STARTROW', state.Startrow+1 		, 'Index of the first row in the original frame'

sxaddpar, noisehdr, 'XTENSION', 'IMAGE'			       ,  'IMAGE extension', before='BITPIX'
sxaddpar, noisehdr, 'EXTNAME' , 'DATA'			       ,  'This extension contains data values'
sxaddpar, noisehdr, 'DATACUBE', state.filename			, 'Name of the original datacube'
sxaddpar, noisehdr, 'INSTRUME', strupcase(state.instr)		, 'Instrument name'
sxaddpar, noisehdr, 'BAND'    , strupcase(state.band)		, 'Band/Filter name'
sxaddpar, noisehdr, 'SPATSMTH', state.smooth			, 'Spatial smoothing applied'
sxaddpar, noisehdr, 'SPECSMTH', state.specsmooth 		, 'Spectral smoothing applied'
sxaddpar, noisehdr, 'CTYPE1'  , sxpar(*state.indatahead,'CTYPE1'), 'TAN projection used'
sxaddpar, noisehdr, 'CTYPE2'  , sxpar(*state.indatahead,'CTYPE2'), 'TAN projection used'
sxaddpar, noisehdr, 'CTYPE3'  , sxpar(*state.indatahead,'CTYPE3'), 'Coordinate system of z-axis'
sxaddpar, noisehdr, 'CRPIX1'  , 1				, '[pix] Reference pixel in x'         
sxaddpar, noisehdr, 'CRPIX2'  , 1				, '[pix] Reference pixel in y'         
sxaddpar, noisehdr, 'CRPIX3'  , 1				, '[pix] Reference pixel in z'
sxaddpar, noisehdr, 'CRVAL1'  , RA				, '[deg] RA at ref. pixel'
sxaddpar, noisehdr, 'CRVAL2'  , DEC				, '[deg] DEC at ref. pixel'
sxaddpar, noisehdr, 'CRVAL3'  , (*state.wave)[0] 		, '[um] Wavelength at ref. pixel'
sxaddpar, noisehdr, 'CD1_1'   , sxpar(*state.indatahead,'CD1_1') , '[] x-component of East'
sxaddpar, noisehdr, 'CD2_1'   , sxpar(*state.indatahead,'CD2_1') , '[] x-component of North'
sxaddpar, noisehdr, 'CD1_2'   , sxpar(*state.indatahead,'CD1_2') , '[] y-component of East'
sxaddpar, noisehdr, 'CD2_2'   , sxpar(*state.indatahead,'CD2_2') , '[] y-component of North'
sxaddpar, noisehdr, 'CDELT1'  , sxpar(*state.indatahead,'CDELT1'), '[deg] Pixel resolution in x'
sxaddpar, noisehdr, 'CDELT2'  , sxpar(*state.indatahead,'CDELT2'), '[deg] Pixel resolution in y'
sxaddpar, noisehdr, 'CDELT3'  , state.dlambda			, '[um] Spectral resolution'
sxaddpar, noisehdr, 'CUNIT1'  , sxpar(*state.indatahead,'CUNIT1'), 'Unit of x-axis'
sxaddpar, noisehdr, 'CUNIT2'  , sxpar(*state.indatahead,'CUNIT2'), 'Unit of y-axis'
sxaddpar, noisehdr, 'CUNIT3'  , sxpar(*state.indatahead,'CUNIT3'), 'Unit of z-axis'
sxaddpar, noisehdr, 'STARTCOL', state.Startcol+1 		, 'Index of the first col in the original frame'
sxaddpar, noisehdr, 'STARTROW', state.Startrow+1 		, 'Index of the first row in the original frame'

mwrfits, 0, fname, prihdr, /CREATE
case  fxpar(*state.indatahead, 'BITPIX') of
  -32: begin
       mwrfits, float((*state.datacube)),  fname, datahdr, /silent
       if state.noiseisvar eq 0 then begin
          mwrfits, float((*state.noisecube)),	fname, noisehdr, /silent
       endif else begin
          mwrfits, float((*state.noisecube)^2), fname, noisehdr, /silent
       endelse
  end
  -64: begin
       mwrfits, (*state.datacube),  fname, datahdr, /silent
       if state.noiseisvar eq 0 then begin
          mwrfits, (*state.noisecube),	fname, noisehdr, /silent
       endif else begin
          mwrfits, (*state.noisecube)^2, fname, noisehdr, /silent
       endelse
  end
endcase

printf, state.log_lun, '[KUBEVIZ] Wrote cube to file: ', fname 
 
end

;-----------------------------------------------------------------------
pro kubeviz_saveimage
; save the spectrum and noise vectors to ext 1 and 2 of an output file.
common kubeviz_state

fname = dialog_pickfile(filter='*.fits', path=state.cwdir, /write, get_path=path, /OVERWRITE_PROMPT)

len = strlen(fname) - strlen(path)
if len eq 0 then return ; don't do anything if no filename is given


;Prepare appropriate data and noise
if state.cubesel le 3 then begin ;Data or noise or SN
   case 1 of
     state.imgmode eq 0: begin
       data  = (*state.datacube)[*,*,state.wpix]
       noise = (*state.noise)[*,*,state.wpix]
     end  
     state.imgmode ge 1 and state.imgmode le 4 : begin
       data = (*state.img1)	     
       noise =  (*state.nimg1)  	 
     end
     state.imgmode eq 5: begin
       data = (*state.datacube)[*,*, state.wpix] - (*state.img1)
       noise = sqrt(((*state.noise)[*,*, state.wpix])^2 + ((*state.nimg1))^2)	       
     end
     state.imgmode ge 6 and state.imgmode le 9 : begin
       data = (*state.img2)	     
       noise =  (*state.nimg2)  	 
     end
     state.imgmode eq 10: begin
       data = (*state.datacube)[*,*, state.wpix] - (*state.img2)
       noise = sqrt(((*state.noise)[*,*, state.wpix])^2 + ((*state.nimg2))^2)	       
     end
     state.imgmode eq 11: begin
       data = (*state.img2) - (*state.img1)
       noise = sqrt(((*state.nimg1)^2 + ((*state.nimg2))^2))		  
     end
     state.imgmode eq 12: begin
       data = (*state.img1) - (*state.img2)
       noise = sqrt(((*state.nimg1)^2 + ((*state.nimg2))^2))		  
     end
   end  								     
endif else begin ;Linefit or lineerrors or lineSN
  if state.flagmode eq 0 then begin
     data = *state.unm_lineresimg
     noise = *state.unm_lineerrresimg
   endif else begin
     data = *state.lineresimg
     noise = *state.lineerrresimg
  endelse
endelse

; Prepare a dummy primary header and the astro solution
mkhdr, prihdr, '', /extend
xyad, (*state.indatahead), state.Startcol, state.Startrow, RA, DEC

datahdr = (*state.indatahead)
noisehdr = (*state.innoisehead)
sxaddpar, datahdr, 'XTENSION', 'IMAGE'			       ,  'IMAGE extension', before='BITPIX'
sxaddpar, datahdr, 'EXTNAME' , 'DATA'			       ,  'This extension contains data values'
sxaddpar, datahdr, 'NAXIS'   ,  2				
sxaddpar, datahdr, 'DATACUBE', state.filename			, 'Name of the original datacube'
sxaddpar, datahdr, 'INSTRUME', strupcase(state.instr)		, 'Instrument name'
sxaddpar, datahdr, 'BAND'    , strupcase(state.band)		, 'Band/Filter name'
sxaddpar, datahdr, 'SPATSMTH', state.smooth			, 'Spatial smoothing applied'
sxaddpar, datahdr, 'SPECSMTH', state.specsmooth 		, 'Spectral smoothing applied'
sxaddpar, datahdr, 'CTYPE1'  , sxpar(*state.indatahead,'CTYPE1'), 'TAN projection used'
sxaddpar, datahdr, 'CTYPE2'  , sxpar(*state.indatahead,'CTYPE2'), 'TAN projection used'
sxaddpar, datahdr, 'CRPIX1'  , 1				, '[pix] Reference pixel in x'         
sxaddpar, datahdr, 'CRPIX2'  , 1				, '[pix] Reference pixel in y'         
sxaddpar, datahdr, 'CRVAL1'  , RA				, '[deg] RA at ref. pixel'
sxaddpar, datahdr, 'CRVAL2'  , DEC				, '[deg] DEC at ref. pixel'
sxaddpar, datahdr, 'CD1_1'   , sxpar(*state.indatahead,'CD1_1') , '[] x-component of East'
sxaddpar, datahdr, 'CD2_1'   , sxpar(*state.indatahead,'CD2_1') , '[] x-component of North'
sxaddpar, datahdr, 'CD1_2'   , sxpar(*state.indatahead,'CD1_2') , '[] y-component of East'
sxaddpar, datahdr, 'CD2_2'   , sxpar(*state.indatahead,'CD2_2') , '[] y-component of North'
sxaddpar, datahdr, 'CDELT1'  , sxpar(*state.indatahead,'CDELT1'), '[deg] Pixel resolution in x'
sxaddpar, datahdr, 'CDELT2'  , sxpar(*state.indatahead,'CDELT2'), '[deg] Pixel resolution in y'
sxaddpar, datahdr, 'CUNIT1'  , sxpar(*state.indatahead,'CUNIT1'), 'Unit of x-axis'
sxaddpar, datahdr, 'CUNIT2'  , sxpar(*state.indatahead,'CUNIT2'), 'Unit of y-axis'
sxaddpar, datahdr, 'STARTCOL', state.Startcol+1 		, 'Index of the first col in the original frame'
sxaddpar, datahdr, 'STARTROW', state.Startrow+1 		, 'Index of the first row in the original frame'
sxdelpar, datahdr, 'NAXIS3'
sxdelpar, datahdr, 'CRPIX3' 
sxdelpar, datahdr, 'CRVAL3' 
sxdelpar, datahdr, 'CDELT3' 
sxdelpar, datahdr, 'CTYPE3' 
sxdelpar, datahdr, 'CUNIT3' 

sxaddpar, noisehdr, 'XTENSION', 'IMAGE'			       ,  'IMAGE extension', before='BITPIX'
sxaddpar, noisehdr, 'EXTNAME' , 'DATA'			       ,  'This extension contains data values'
sxaddpar, noisehdr, 'NAXIS'   ,  2				
sxaddpar, noisehdr, 'DATACUBE', state.filename			, 'Name of the original datacube'
sxaddpar, noisehdr, 'INSTRUME', strupcase(state.instr)		, 'Instrument name'
sxaddpar, noisehdr, 'BAND'    , strupcase(state.band)		, 'Band/Filter name'
sxaddpar, noisehdr, 'SPATSMTH', state.smooth			, 'Spatial smoothing applied'
sxaddpar, noisehdr, 'SPECSMTH', state.specsmooth 		, 'Spectral smoothing applied'
sxaddpar, noisehdr, 'CTYPE1'  , sxpar(*state.indatahead,'CTYPE1'), 'TAN projection used'
sxaddpar, noisehdr, 'CTYPE2'  , sxpar(*state.indatahead,'CTYPE2'), 'TAN projection used'
sxaddpar, noisehdr, 'CRPIX1'  , 1				, '[pix] Reference pixel in x'         
sxaddpar, noisehdr, 'CRPIX2'  , 1				, '[pix] Reference pixel in y'         
sxaddpar, noisehdr, 'CRVAL1'  , RA				, '[deg] RA at ref. pixel'
sxaddpar, noisehdr, 'CRVAL2'  , DEC				, '[deg] DEC at ref. pixel'
sxaddpar, noisehdr, 'CD1_1'   , sxpar(*state.indatahead,'CD1_1') , '[] x-component of East'
sxaddpar, noisehdr, 'CD2_1'   , sxpar(*state.indatahead,'CD2_1') , '[] x-component of North'
sxaddpar, noisehdr, 'CD1_2'   , sxpar(*state.indatahead,'CD1_2') , '[] y-component of East'
sxaddpar, noisehdr, 'CD2_2'   , sxpar(*state.indatahead,'CD2_2') , '[] y-component of North'
sxaddpar, noisehdr, 'CDELT1'  , sxpar(*state.indatahead,'CDELT1'), '[deg] Pixel resolution in x'
sxaddpar, noisehdr, 'CDELT2'  , sxpar(*state.indatahead,'CDELT2'), '[deg] Pixel resolution in y'
sxaddpar, noisehdr, 'CUNIT1'  , sxpar(*state.indatahead,'CUNIT1'), 'Unit of x-axis'
sxaddpar, noisehdr, 'CUNIT2'  , sxpar(*state.indatahead,'CUNIT2'), 'Unit of y-axis'
sxaddpar, noisehdr, 'STARTCOL', state.Startcol+1 		, 'Index of the first col in the original frame'
sxaddpar, noisehdr, 'STARTROW', state.Startrow+1 		, 'Index of the first row in the original frame'
sxdelpar, noisehdr, 'NAXIS3'
sxdelpar, noisehdr, 'CRPIX3'  
sxdelpar, noisehdr, 'CRVAL3'  
sxdelpar, noisehdr, 'CDELT3'  
sxdelpar, noisehdr, 'CTYPE3'  
sxdelpar, noisehdr, 'CUNIT3'  

mwrfits, 0, fname, prihdr, /CREATE
case  fxpar(*state.indatahead, 'BITPIX') of
  -32: begin
       mwrfits, float(data),  fname, datahdr, /silent
       if state.noiseisvar eq 0 then begin
          mwrfits, float(noise),	fname, noisehdr, /silent
       endif else begin
          mwrfits, float(noise^2), fname, noisehdr, /silent
       endelse
  end
  -64: begin
       mwrfits, data,  fname, datahdr, /silent
       if state.noiseisvar eq 0 then begin
          mwrfits, noise,	fname, noisehdr, /silent
       endif else begin
          mwrfits, noise^2, fname, noisehdr, /silent
       endelse
  end
endcase

printf, state.log_lun, '[KUBEVIZ] Wrote data and noise images to file: ', fname 
 
end
;-----------------------------------------------------------------------
pro kubeviz_util_swap, lo, hi
; swap the values of the two input variables

tmp = lo
lo  = hi
hi  = tmp

end

;-----------------------------------------------------------------------
pro kubeviz_setcolour, ctab
common kubeviz_state
; load selected colour table 

loadct, ctab, /silent

tvlct, r, g, b, /get  

state.red = r
state.green = g
state.blue = b

kubeviz_modify_colour

state.ctab = ctab   ;; store selected value

; set default background to white:
!p.background = 255

end

;-----------------------------------------------------------------------
pro kubeviz_modify_colour, ci=ci
; define last few colour indices 
; 
; ci  - colour index, if set use this instead of the state colour indices

common kubeviz_state

state.maxcol = 247

r = state.red
g = state.green
b = state.blue

if state.invert eq 1 then begin
    r = reverse(r)
    g = reverse(g)
    b = reverse(b)
endif

if n_elements(ci) ne 0 then begin
; this is clumsy, change it! 
; it's done like this as ci length is 255 instead of 256
    r = ci
    g = ci
    b = ci
    if state.invert eq 1 then begin
        r = reverse(r)
        g = reverse(g)
        b = reverse(b)
    endif
    r[253] = 255  &  g[253] = 0  &  b[253] = 0
    r[254] = 0  &  g[254] = 0  &  b[254] = 0
    r[252] = 175  &  g[252] = 175  &  b[252] = 175
    r[251] = 225  &  g[251] = 225  &  b[251] = 225
    r[250] = 0  &  g[250] = 0  &  b[250] = 255
    r[249] = 0  &  g[249] = 255  &  b[249] = 0
    r[248] = 255  &  g[248] = 204  &  b[248] = 229
    tvlct, state.red[r], state.green[g], state.blue[b]
    return
endif

; set index 255 to white - to make background in spectra panels white 
r[255] = 255  &  g[255] = 255  &  b[255] = 255
; set index 254 to black  - to make the foreground lines black
r[254] = 0    &  g[254] = 0    &  b[254] = 0
; set index 253 to red  - to make the marker lines red
r[253] = 255  &  g[253] = 0    &  b[253] = 0
; set index 252 to blue  
r[252] =   0  &  g[252] =   0  &  b[252] = 255   ; blue
; set index 251 to light grey  - to highlight 2nd selected lambda range
r[251] = 225  &  g[251] = 225  &  b[251] = 225
; set index 250 to dark grey  - to highlight 1st selected lambda range
r[250] = 150  &  g[250] = 150  &  b[250] = 150
; set index 249 to green
r[249] = 0    &  g[249] = 255  &  b[249] = 0
; set index 248 to pink
r[248] = 255  &  g[248] = 204  &  b[248] = 229

tvlct, r, g, b

end

;------------------------------------------------
pro kubeviz_linefit_cancel, init=init, destroy=destroy, cancel=cancel, event=event
common kubeviz_widgetids

if keyword_set(init) then begin
  
  dimensions = GET_SCREEN_SIZE()

  fitcancelbase  = widget_base(title='Fit', group_leader = widgetids.base1,  /col, /base_align_center, $
                 scr_xsize = 150, scr_ysize=40, xoffset=dimensions[0]/2, yoffset=dimensions[1]/2)
  widgetids.cancelfit_id = widget_button(fitcancelbase,  value=' Interrupt Fit ', uvalue='FITCANCEL')
  widget_control, fitcancelbase, /realize
  
endif

if keyword_set(destroy) then begin
  widget_control, widget_info(widgetids.cancelfit_id,/PARENT),/DESTROY 
endif

if keyword_set(cancel) then begin
    ev = widget_event(widgetids.cancelfit_id, /nowait, /save_hourglass) 
    if strmid(tag_names(ev, /structure_name), 7, 4) eq "BUTT" then event = 1 else event = 0
endif

end

;----------------------------------------------------------------------
pro kubeviz_smooth_parameters
common kubeviz_pars
common kubeviz_state
common kubeviz_widgetids

if (not (xregistered('kubeviz_smoothpars', /noshow))) then begin

pars.spatsmooth = state.smooth
pars.specsmooth = state.specsmooth

smoothparsbase = widget_base(title='Smooth Parameters',group_leader = widgetids.base1,  /col, /base_align_center, $
                 scr_xsize = 300, scr_ysize=190)
dummy = widget_label(smoothparsbase, value='WARNING: The current linefit results are reset ')
dummy = widget_label(smoothparsbase, value='if the smooth parameters are changed!')
dummy = widget_base(smoothparsbase, ysize=5,  /row)
row1 = widget_base(smoothparsbase, /row)
dummy = widget_label(row1, xsize=200, /align_left, value='Spatial smoothing kernel: ')
dummy = widget_text(row1, /editable, xsize=4, ysize=1, uvalue='SPATSMOOTH', value=kubeviz_str(pars.spatsmooth, format='(I2)'))
row2 = widget_base(smoothparsbase, /row)
dummy = widget_label(row2, xsize=200, /align_left, value='Spectral smoothing kernel: ')
dummy = widget_text(row2, /editable, xsize=4, ysize=1, uvalue='SPECSMOOTH', value=kubeviz_str(pars.specsmooth, format='(I2)'))
dummy = widget_base(smoothparsbase, ysize=5,  /row)
row3 = widget_base(smoothparsbase, /row)
dummy = widget_button(row3,  value=' Apply ', uvalue='APPLY')
dummy = widget_label(row3, xsize=75, /align_left, value='')
dummy = widget_button(row3,  value=' Cancel ', uvalue='CANCEL')

widget_control, smoothparsbase, /realize
xmanager, 'kubeviz_smoothpars', smoothparsbase, /no_block

endif
end

;-----------------------------------------------------------------------
pro kubeviz_smoothpars_event, ev
; handle events generated by the zoom spectrum window
common kubeviz_state
common kubeviz_pars

widget_control, ev.id, get_uvalue=value
if (n_elements(value) eq 0) then value = ''
name = strmid(tag_names(ev, /structure_name), 7, 4)
case (name) of
  "BUTT": begin
     if (value eq "CANCEL") then widget_control, ev.top, /destroy
     if (value eq "APPLY")  then begin
	if pars.spatsmooth ge 1 and pars.specsmooth ge 1 then begin
	  widget_control, ev.top, /destroy
          widget_control, /hourglass
	  state.smooth = pars.spatsmooth
          state.specsmooth = pars.specsmooth
	
	  ; apply the new smoothing and replace bad data with zeros:
	  kubeviz_smooth, /datacube, /montecarlo
	  
	  kubeviz_linefit_init
	  kubeviz_medsum_image_update
	  kubeviz_plotspax
    	  kubeviz_plotinfo
    	  kubeviz_plotspec
    	  kubeviz_plotspeczoom
    	  kubeviz_linefit_update
	endif else printf, state.log_lun, '[WARNING] The smooth length must be >= 1.'  
     endif
  end
  "TEXT": begin
     if (value eq "SPATSMOOTH") then begin
        widget_control, ev.id, get_value=value
	pars.spatsmooth = value
     endif	
     if (value eq "SPECSMOOTH") then begin
        widget_control, ev.id, get_value=value
	pars.specsmooth = value
     endif	
  end
  else:
endcase  

end

;----------------------------------------------------------------------
pro kubeviz_fitall_range
common kubeviz_pars
common kubeviz_state
common kubeviz_widgetids

if (not (xregistered('kubeviz_fitallrange', /noshow))) then begin

pars.fitallrange = state.fitallrange
pars.Nloopadj = state.Nloopadj

fitallrangebase = widget_base(title='FITALL range',group_leader = widgetids.base1,  /col, /base_align_center, $
                 scr_xsize = 350, scr_ysize=175, xoffset=250, yoffset=55)
		 
dummy = widget_label(fitallrangebase, value='Parameters (pixel coordinates) for the FITALL button')
row1 = widget_base(fitallrangebase, /row)
dummy = widget_label(row1, xsize=150, /align_left, value='X range (min/max): ')
dummy = widget_text(row1, /editable, xsize=5, ysize=1, uvalue='FITXMIN', value=kubeviz_str(pars.fitallrange[0], format='(I4)'))
dummy = widget_text(row1, /editable, xsize=5, ysize=1, uvalue='FITXMAX', value=kubeviz_str(pars.fitallrange[1], format='(I4)'))
row2 = widget_base(fitallrangebase, /row)
dummy = widget_label(row2, xsize=150, /align_left, value='Y range (min/max): ')
dummy = widget_text(row2, /editable, xsize=5, ysize=1, uvalue='FITYMIN', value=kubeviz_str(pars.fitallrange[2], format='(I4)'))
dummy = widget_text(row2, /editable, xsize=5, ysize=1, uvalue='FITYMAX', value=kubeviz_str(pars.fitallrange[3], format='(I4)'))
row3 = widget_base(fitallrangebase, /row)
dummy = widget_label(row3, xsize=150, /align_left, value='Niter FITADJ: ')
dummy = widget_text(row3, /editable, xsize=5, ysize=1, uvalue='NLOOPADJ', value=kubeviz_str(pars.Nloopadj, format='(I4)'))
row4 = widget_base(fitallrangebase, /row)
dummy = widget_button(row4,  value=' Apply ', uvalue='APPLY')
dummy = widget_label(row4, xsize=50, /align_left, value='')
dummy = widget_button(row4,  value=' Defaults ', uvalue='DEFAULT')
dummy = widget_label(row4, xsize=50, /align_left, value='')
dummy = widget_button(row4,  value=' Cancel ', uvalue='CANCEL')

widget_control, fitallrangebase, /realize
xmanager, 'kubeviz_fitallrange', fitallrangebase, /no_block

endif
end

;-----------------------------------------------------------------------
pro kubeviz_fitallrange_event, ev
common kubeviz_state
common kubeviz_pars

widget_control, ev.id, get_uvalue=value
if (n_elements(value) eq 0) then value = ''
name = strmid(tag_names(ev, /structure_name), 7, 4)
case (name) of
  "BUTT": begin
     if (value eq "CANCEL") then Widget_control, ev.top, /destroy
     if (value eq "DEFAULT") then begin
       state.fitallrange = [0, state.Ncol-1, 0, state.Nrow-1]
       state.Nloopadj = 0
       Widget_control, ev.top, /destroy
     endif
     if (value eq "APPLY") then begin
	if pars.fitallrange[0] ge 0 and  pars.fitallrange[2] ge 0 and $
	   pars.fitallrange[1] lt state.Ncol and pars.fitallrange[3] lt state.Nrow and pars.Nloopadj ge 0 then begin
	  Widget_control, ev.top, /destroy
          state.fitallrange = pars.fitallrange
	  state.Nloopadj = pars.Nloopadj
	endif else printf, state.log_lun, '[WARNING] At least one value is out of limits'
     endif
  end
  "TEXT": begin
     widget_control, ev.id, get_value=numval
     case value of
      "NLOOPADJ": pars.Nloopadj = numval
      "FITXMIN": pars.fitallrange[0] = numval
      "FITXMAX": pars.fitallrange[1] = numval
      "FITYMIN": pars.fitallrange[2] = numval
      "FITYMAX": pars.fitallrange[3] = numval
      else:
     endcase	
  end
  else:
endcase  

end


;----------------------------------------------------------------------
pro kubeviz_gotomask
common kubeviz_pars
common kubeviz_state
common kubeviz_widgetids

if (not (xregistered('kubeviz_gotomask', /noshow))) then begin

pars.gotomask = state.imask

gotomaskbase = widget_base(title='Go To Mask',group_leader = widgetids.base1,  /col, /base_align_center, $
                 scr_xsize = 200, scr_ysize=100, xoffset=250, yoffset=55)
row1 = widget_base(gotomaskbase, /row)
dummy = widget_label(row1, xsize=150, /align_left, value='Go To Mask: ')
dummy = widget_text(row1, /editable, xsize=5, ysize=1, uvalue='GOTOMASK', value=kubeviz_str(pars.gotomask, format='(I3)'))
row3 = widget_base(gotomaskbase, /row)
dummy = widget_button(row3,  value=' Apply ', uvalue='APPLY')
dummy = widget_label(row3, xsize=55, /align_left, value='')
dummy = widget_button(row3,  value=' Cancel ', uvalue='CANCEL')

widget_control, gotomaskbase, /realize
xmanager, 'kubeviz_gotomask', gotomaskbase, /no_block

endif
end
;-----------------------------------------------------------------------
pro kubeviz_gotomask_event, ev
common kubeviz_state
common kubeviz_pars

widget_control, ev.id, get_uvalue=value
if (n_elements(value) eq 0) then value = ''
name = strmid(tag_names(ev, /structure_name), 7, 4)
case (name) of
  "BUTT": begin
     if (value eq "CANCEL") then Widget_control, ev.top, /destroy
     if (value eq "APPLY") then begin
	if pars.gotomask le state.Nmask  and pars.gotomask gt 0 then begin
	  Widget_control, ev.top, /destroy
          state.imask = pars.gotomask
	  kubeviz_medianspec
	  kubeviz_linefit_update
          kubeviz_plotspax 
          kubeviz_plotinfo
          kubeviz_plotspec
          kubeviz_plotspeczoom
	endif else printf, state.log_lun, '[WARNING] This mask index does not exist.'
     endif
  end
  "TEXT": begin
     if (value eq "GOTOMASK") then begin
        widget_control, ev.id, get_value=value
	pars.gotomask = value
     endif	
  end
  else:
endcase  

end


;----------------------------------------------------------------------
pro kubeviz_maskpars
common kubeviz_pars
common kubeviz_state
common kubeviz_widgetids

if (not (xregistered('kubeviz_maskpars', /noshow))) then begin

pars.maskradius = state.maskradius
pars.maskmode = state.maskmode

maskparsbase = widget_base(title='Mask Parameters',group_leader = widgetids.base1,  /col, /base_align_center, $
                 scr_xsize = 250, scr_ysize=110, xoffset=250, yoffset=55)
row1 = widget_base(maskparsbase, /row, /exclusive)
b0 = widget_button(row1, value='SINGLE',   uvalue='SINGLE')
b1 = widget_button(row1, value='CIRCULAR', uvalue='CIRCULAR')
b2 = widget_button(row1, value='SQUARE',   uvalue='SQUARE')
row2 = widget_base(maskparsbase, /row)
dummy = widget_label(row2, xsize=150, /align_left, value='Mask Radius (Size): ')
dummy = widget_text(row2, /editable, xsize=5, ysize=1, uvalue='maskpars', value=kubeviz_str(pars.maskradius, format='(I3)'))
row3 = widget_base(maskparsbase, /row)
dummy = widget_button(row3,  value=' Apply ', uvalue='APPLY')
dummy = widget_label(row3, xsize=55, /align_left, value='')
dummy = widget_button(row3,  value=' Cancel ', uvalue='CANCEL')

case pars.maskmode of 
  0: widget_control, b0, set_button=1
  1: widget_control, b1, set_button=1
  2: widget_control, b2, set_button=1
  else:
endcase

widget_control, maskparsbase, /realize
xmanager, 'kubeviz_maskpars', maskparsbase, /no_block

endif
end


;-----------------------------------------------------------------------
pro kubeviz_maskpars_event, ev
common kubeviz_state
common kubeviz_pars

widget_control, ev.id, get_uvalue=value
if (n_elements(value) eq 0) then value = ''
name = strmid(tag_names(ev, /structure_name), 7, 4)
case (name) of
  "BUTT": begin
     if (value eq "SINGLE")   then pars.maskmode = 0
     if (value eq "CIRCULAR") then pars.maskmode = 1
     if (value eq "SQUARE")   then pars.maskmode = 2
     if (value eq "CANCEL")   then Widget_control, ev.top, /destroy
     if (value eq "APPLY")    then begin
	if pars.maskradius ge 1 then begin
	  Widget_control, ev.top, /destroy
          state.maskmode   = pars.maskmode
	  state.maskradius = pars.maskradius
	endif else printf, state.log_lun, '[WARNING] The mask radius must be >= 1.'  
     endif
  end
  "TEXT": begin
     if (value eq "maskpars") then begin
        widget_control, ev.id, get_value=value
	pars.maskradius = value
     endif	
  end
  else:
endcase  

end

;----------------------------------------------------------------------
pro kubeviz_zcut_parameters
common kubeviz_state
common kubeviz_pars
common kubeviz_widgetids

if (not (xregistered('kubeviz_zcutpars', /noshow))) then begin

pars.zmin_ima = state.zmin_ima
pars.zmax_ima = state.zmax_ima

zcutparsbase = widget_base(title='Scale Parameters (zcut)',group_leader = widgetids.base1,  /col, /base_align_center, $
                 scr_xsize = 400, scr_ysize=300)

zcuthist = widget_draw(zcutparsbase, xsize=380, ysize=200 )

dummy = widget_base(zcutparsbase, ysize=3,  /row)
row1 = widget_base(zcutparsbase, /row)
dummy = widget_label(row1, xsize=100, /align_right, value='Limits     Low:')
dummy = widget_text(row1, /editable, xsize=7, ysize=1, uvalue='ZMIN', value=kubeviz_str(pars.zmin_ima, format='(F7.2)'))
dummy = widget_label(row1, xsize=40, /align_left, value=' High:')
dummy = widget_text(row1, /editable, xsize=7, ysize=1, uvalue='ZMAX', value=kubeviz_str(pars.zmax_ima, format='(F7.2)'))

dummy = widget_base(zcutparsbase, ysize=3,  /row)
row2 = widget_base(zcutparsbase, /row)
dummy = widget_button(row2,  value=' Apply ', uvalue='APPLY')
dummy = widget_label(row2, xsize=75, /align_left, value='')
dummy = widget_button(row2,  value=' Cancel ', uvalue='CANCEL')
dummy = widget_base(zcutparsbase, ysize=3,  /row)

widget_control, zcutparsbase, /realize
widget_control, zcuthist, get_value=wid5
widgetids.wid5 = wid5
kubeviz_plotzcuthist

xmanager, 'kubeviz_zcutpars', zcutparsbase, /no_block

endif
end

;-----------------------------------------------------------------------
pro kubeviz_zcutpars_event, ev
common kubeviz_state
common kubeviz_pars

widget_control, ev.id, get_uvalue=value
if (n_elements(value) eq 0) then value = ''
name = strmid(tag_names(ev, /structure_name), 7, 4)
case (name) of
  "BUTT": begin
     if (value eq "CANCEL") then Widget_control, ev.top, /destroy
     if (value eq "APPLY") then begin
	if pars.zmin_ima lt pars.zmax_ima then begin
	  Widget_control, ev.top, /destroy
          state.zmin_ima = pars.zmin_ima
          state.zmax_ima = pars.zmax_ima
	  state.zcuts = 1
	  kubeviz_plotspax
	endif else printf, state.log_lun, '[WARNING] Min value should be smaller than max value.'
     endif
  end
  "TEXT": begin
     if (value eq "ZMIN") then begin
        widget_control, ev.id, get_value=value
	pars.zmin_ima = value
	kubeviz_plotzcuthist 
     endif	
     if (value eq "ZMAX") then begin
        widget_control, ev.id, get_value=value
	pars.zmax_ima = value
	kubeviz_plotzcuthist 
     endif	
  end
  else:
endcase  

end

;-----------------------------------------------------------------------
pro kubeviz_plotzcuthist
common kubeviz_state
common kubeviz_widgetids
common kubeviz_pars

wset, widgetids.wid5
erase, 255
bin = 3*(max((*state.curr_ima), /NaN)-min((*state.curr_ima), /NaN))/((state.Ncol)*(state.Nrow))
;plothist, (*state.curr_ima), /NaN, bin=bin, /ylog, /fill, fcolor=254, color=254, background=255, /xstyle, /ystyle, $
;          title='Pixel Distribution', /device, position=[40, 20, 360, 180] , charsize=1.0

; overplot marker to highlight selected wavelength
oplot, [pars.zmin_ima, pars.zmin_ima], [1E-10,1E10], color=253, linestyle=0, thick=2.0
oplot, [pars.zmax_ima, pars.zmax_ima], [1E-10,1E10], color=252, linestyle=0, thick=2.0


end

;----------------------------------------------------------------------
pro kubeviz_selecthead
common kubeviz_pars
common kubeviz_state
common kubeviz_widgetids
common kubeviz_winsize

if (not (xregistered('kubeviz_selecthead', /noshow))) then begin

selectheadbase = widget_base(title='Sel. Head',group_leader = widgetids.base1,  /col, /base_align_center, $
                 scr_xsize = 180, scr_ysize=130, xoffset=(winsize.xwin1/2)-80, yoffset=winsize.ywin1/2.-30)

if ptr_valid(state.inprihead) then begin
 row0 = widget_base(selectheadbase, /row)
 dummy = widget_button(row0, xsize=150, /align_center, value='Primary Header', uvalue='PRIMARY')
endif
row1 = widget_base(selectheadbase, /row)
dummy = widget_button(row1, xsize=150, /align_center, value='Data Header', uvalue='DATA')
row2 = widget_base(selectheadbase, /row)
dummy = widget_button(row2, xsize=150, /align_center, value='Noise Header', uvalue='NOISE')
row3 = widget_base(selectheadbase, /row)
dummy = widget_button(row3,  value=' Cancel ', /align_center, uvalue='CANCEL')

widget_control, selectheadbase, /realize
xmanager, 'kubeviz_selecthead', selectheadbase, /no_block

endif
end

;------------------------------------------
pro kubeviz_selecthead_event, ev
common kubeviz_state
common kubeviz_pars
common kubeviz_widgetids

widget_control, ev.id, get_uvalue=value
if (n_elements(value) eq 0) then value = ''
name = strmid(tag_names(ev, /structure_name), 7, 4)

case (name) of
  "BUTT": begin
     if (value eq "CANCEL") then begin
       Widget_control, ev.top, /destroy
       return
     endif  
     if (value eq "PRIMARY") then begin
       pars.selhead = 0
       Widget_control, ev.top, /destroy
     endif  
     if (value eq "DATA")    then begin
       pars.selhead = 1
       Widget_control, ev.top, /destroy
     endif  
     if (value eq "NOISE")   then begin
       pars.selhead = 2
       Widget_control, ev.top, /destroy
     endif  
   end
  else:
endcase  

if (not (xregistered('kubeviz_data_header', /noshow))) then begin

headerbase = widget_base(group_leader = widgetids.base1, /column, $
                         /base_align_right, uvalue = 'header_base')

case pars.selhead of
  0: selhead = *state.inprihead
  1: selhead = *state.indatahead
  2: selhead = *state.innoisehead
endcase  

headlen = size(selhead)
headwid = max(strlen(selhead))

dummy = widget_text(headerbase, value=selhead, xsize=headwid, ysize=min([50,headlen[1]]), /scroll)

widget_control, headerbase, /realize
xmanager, 'kubeviz_data_header', headerbase, /no_block

endif

end

;-----------------------------------------------------------------------
pro kubeviz_help_whatisnew
common kubeviz_state
common kubeviz_widgetids

if (not (xregistered('kubeviz_help_whatisnew', /noshow))) then begin

helpnewbase = widget_base(group_leader = widgetids.base1, /column, $
                        /base_align_right, uvalue = 'help_instrbase')

findpro, 'kubeviz', dirlist=kubeviz_dir, /noprint 
whatisnew_file = kubeviz_dir[0]+'doc/whatsnew.txt'

h = 'Kubeviz version: ' + state.version
line = ''
openr, lun, whatisnew_file, /get_lun
while not EOF(lun) do begin
  readf, lun, line 
  h = [h, line] 
endwhile
close, lun

if size(h, /N_elements) lt 50 $
      then dummy = widget_text(helpnewbase, value=h, xsize=85, ysize=size(h, /N_elements)+1) $
      else dummy = widget_text(helpnewbase, value=h, xsize=85, ysize=50, /scroll)


widget_control, helpnewbase, /realize
xmanager, 'kubeviz_help_whatisnew', helpnewbase, /no_block

endif

end

;-----------------------------------------------------------------------
pro kubeviz_help_instructions
common kubeviz_state
common kubeviz_widgetids

if (not (xregistered('kubeviz_help_instructions', /noshow))) then begin

helpinstrbase = widget_base(group_leader = widgetids.base1, /column, $
                        /base_align_right, uvalue = 'help_instrbase')

findpro, 'kubeviz', dirlist=kubeviz_dir, /noprint 
instr_file = kubeviz_dir[0]+'doc/instructions.txt'

h = 'Kubeviz version: ' + state.version
line = ''
openr, lun, instr_file, /get_lun
while not EOF(lun) do begin
  readf, lun, line 
  h = [h, line] 
endwhile
close, lun

if size(h, /N_elements) lt 50 $
      then dummy = widget_text(helpinstrbase, value=h, xsize=95, ysize=size(h, /N_elements)+1) $
      else dummy = widget_text(helpinstrbase, value=h, xsize=95, ysize=50, /scroll)


widget_control, helpinstrbase, /realize
xmanager, 'kubeviz_help_instructions', helpinstrbase, /no_block

endif

end

;-----------------------------------------------------------------------
pro kubeviz_help_shortcuts
common kubeviz_state
common kubeviz_widgetids

if (not (xregistered('kubeviz_help_shortcuts', /noshow))) then begin

helpkeybase = widget_base(group_leader = widgetids.base1, /column, /base_align_right, uvalue = 'helpkey_base')

h = strarr(50)

 h[0] = 'Kubeviz version: ' + state.version
 h[1] = ''
 h[2] = 'Keyboard shortcuts:'
 h[3] = ''
 h[4] = 'cursor keys : step through spatial plane'
 h[5] = '<  >        : step through wavelength planes ' 
 h[6] = 'm           : select spaxel' 
 h[7] = 'n           : deselect spaxel' 
 h[8] = 'r           : clear all selected spaxels' 
 h[9] = 'k           : show previous spaxel map'
h[10] = 'l           : show next spaxel map'
h[11] = 'c           : fit a 2d gaussian to the current image'
h[12] = 'b           : toggle the flag for the selected spaxel' 
h[13] = 'f           : Fit using the linefit framework (set parameters first!)'
h[14] = 'a           : Fit using linefit and the initial guess'
h[15] = '              from adjacent spaxels results'
h[16] = 'q           : exit'
h[17] = ''
h[18] = 'In spectrum window only:'
h[19] = ''
h[20] = 's	     : Set wavelength (in combination with Sel Range button)'
h[21] = 'z	     : Reset the redshift assuming the current wavelength' 
h[22] = '	       to be that of the mainline'
h[23] = ''
h[24] = 'In spectral zoom window only:'
h[25]=	''
h[26]=	'g	     : Simple Gauss fit at current wavelength'



if size(h, /N_elements) lt 50 $
      then dummy = widget_text(helpkeybase, value=h, xsize=85, ysize=size(h, /N_elements)+1) $
      else dummy = widget_text(helpkeybase, value=h, xsize=85, ysize=50, /scroll)

widget_control, helpkeybase, /realize
xmanager, 'kubeviz_help_shortcuts', helpkeybase, /no_block

endif

end

;------------------------------------------------------------------------
pro kubeviz_batchmode, redshift=redshift, fit_all_lines=fit_all_lines

common kubeviz_state

if N_elements(redshift) eq 0 then begin
      print, '[ ERROR ] Batch mode does not work if the redshift keyword is not set. '
      kubeviz_destroy
      return
endif 

if total(state.lines) eq 0 then begin
      print, '[ ERROR ] No lines identified in the spectrum for the given redshift and lineset. '
      kubeviz_destroy
      return
endif 

if N_elements(fit_all_lines) gt 0 then begin
   Nlines = N_elements(kubeviz_chooselines((*state.linesets), state.lines, lineset=lineset))
   printf, state.log_lun, '[KUBEVIZ] Setting fit parameters for all lines in set.'
   for iline=0, Nlines-1 do begin
     state.pndofit[iline] = 1
     state.pcdofit[iline] = 1
   endfor			 
endif else begin
   state.pndofit[0] = 1
   state.pcdofit[0] = 1
endelse
case state.domontecarlo of
   0: str_errormethod = "Noise Cube"
   1: str_errormethod = "Bootstrap"
   2: str_errormethod = "MonteCarlo 1"
   3: str_errormethod = "MonteCarlo 2"
   4: str_errormethod = "MonteCarlo 3"
endcase
case state.linefit_type of
   0: str_fitmethod = 'Gauss'
   1: str_fitmethod = 'Moments'
endcase

if state.linefit_mode eq 0 then begin 
   printf, state.log_lun, '[KUBEVIZ] Fitting spaxels with method: '+str_fitmethod
   printf, state.log_lun, '[KUBEVIZ] Fitting spaxels with error method: '+str_errormethod
   str_resmode = 'spax'
   t_start = systime(/seconds)
   step_count = 0L
   abs_step = 1L*(state.Ncol * state.Nrow * state.percent_step)/100
   cur_col = state.col
   cur_row = state.row
   for col=0,state.Ncol-1 do begin
       state.col = col
       for row=0,state.Nrow-1 do begin
	   state.row = row
	   kubeviz_statusline, step_count, abs_step, t_start
	   kubeviz_linefit_dofit 
       endfor
   endfor
   kubeviz_statusline, /close
   state.col  = state.Ncol/2
   state.row  = state.Nrow/2
   state.wpix = state.Nwpix/2
   kubeviz_linefit_autoflag
 endif else begin
   printf, state.log_lun, '[KUBEVIZ] Fitting loaded masks with method: '+str_fitmethod
   printf, state.log_lun, '[KUBEVIZ] Fitting loaded masks with error method: '+str_errormethod
   str_resmode = 'mask'
   cur_mask = state.imask
   for mask = 1L, state.Nmask do begin
       state.imask = mask
       kubeviz_linefit_dofit
   endfor
   state.imask = cur_mask
  endelse
   
kubeviz_linefit_saveres, fname=state.outdir+strmid(state.filename,0,strlen(state.filename)-5)+'_res_'+strlowcase(strmid(str_fitmethod,0,3))+'_'+str_resmode+'.fits'
kubeviz_savesession, fname=state.outdir+strmid(state.filename,0,strlen(state.filename)-5)+'_'+strlowcase(strmid(str_fitmethod,0,3))+'_'+str_resmode+'.sav'
kubeviz_destroy

end

;------------------------------------------------------------------------
pro kubeviz_create, fname, scroll=scroll
;; create the graphical windows

common kubeviz_state
common kubeviz_widgetids
common kubeviz_winsize, winsize

;Set X environment
set_plot, 'X'                     ; initialize X-server plotting device
device, decomposed=0, retain=2    ; 
kubeviz_setcolour, 0              ; start with grey scale colour table
!x.thick = 1                      ; set default x-axis line thickness
!y.thick = 1                      ; set default y-axis line thickness

dimensions = GET_SCREEN_SIZE()

; define a zoomfac
if state.zoomfac eq 0 then begin
    case 1 of
      state.Ncol lt 100:                       state.zoomfac = fix(500./state.Ncol)
      state.Ncol ge 100 and state.Ncol lt 200: state.zoomfac = fix(600./state.Ncol)
      state.Ncol ge 200 and state.Ncol lt 800: state.zoomfac = fix(800./state.Ncol)
      state.Ncol ge 800 and state.Nrow lt dimensions[1]-200: state.zoomfac = 1
      else: begin
        print, '[ ERROR ] The spatial dimensions of the cube are too big for the current display.'
	print, '[ ERROR ] Increase the screen resolution or use the trim keyword to load a portion of the cube.'
	kubeviz_destroy, /abort
      end
    endcase
    ; If the window is too large in the y direction reduce zoomfac. 200 is the size of the fixed rows (menubar+info)
    if (state.Nrow*state.zoomfac) gt (dimensions[1]-200) then state.zoomfac = fix((dimensions[1]-200)/state.Nrow)
endif

if keyword_set(scroll) then state.scroll = 1

; Initialize parameters
state.col  = state.Ncol/2
state.row  = state.Nrow/2
state.wpix = state.Nwpix/2

state.imgmode = 0     ; set image mode to show individual slices
state.marker =  1     ; toggle spectrum marker line on/off 
state.scale  = -1     ; toggle user-scaling on/off

; define the size of the spaxel window 
winsize.xwin1 = max([500,state.Ncol * state.zoomfac])
winsize.ywin1 = state.Nrow * state.zoomfac

; define arrays for wavelength summed/median image
state.img1 = ptr_new(fltarr(state.Ncol, state.Nrow))
state.img2 = ptr_new(fltarr(state.Ncol, state.Nrow))
; and noise and badpix images:
state.nimg1 = ptr_new(fltarr(state.Ncol, state.Nrow))
state.nimg2 = ptr_new(fltarr(state.Ncol, state.Nrow))
state.bimg1 = ptr_new(fltarr(state.Ncol, state.Nrow))
state.bimg2 = ptr_new(fltarr(state.Ncol, state.Nrow))

; define arrays for the current image
state.curr_ima = ptr_new(fltarr(state.Ncol, state.Nrow))
state.curr_byteima = ptr_new(fltarr(state.Ncol, state.Nrow))
state.curr_byteima_grid = ptr_new(fltarr(state.Ncol*state.zoomfac, state.Nrow*state.zoomfac))

; set up the lay out of each window
kubeviz_setup_viewers, menubar
kubeviz_setup_linefit
if state.linefit_type eq 1 then kubeviz_linefit_typeswitch

; Display the starting image and spectra
kubeviz_plotspax
kubeviz_plotinfo
kubeviz_plotspec
kubeviz_plotspeczoom

; Call the window manager
xmanager, 'kubeviz_spax',       widgetids.base1, /no_block
xmanager, 'kubeviz_spec',       widgetids.base2, /no_block
xmanager, 'kubeviz_speczoom',   widgetids.base3, /no_block
xmanager, 'kubeviz_linefit',    widgetids.base5, /no_block

end

;-------------------------------------------------------------
pro kubeviz_libtest

;Check consistency of IDL setup. Warnings and errors are printed 
;to the terminal even if log_lun is set.

skip_routines = ['LIST', 'DELLOG', 'SETLOG', 'TRNLOG', 'DL_MAC', 'DL_DOS', 'BOOLEAN']
resolve_all, /continue_on_error, /quiet,  unresolved=unresolved, skip_routines=skip_routines ; compile every useful routine
if total(unresolved eq skip_routines) lt n_elements(unresolved) then begin
   match, unresolved, skip_routines, unresmatch, dummy
   remove, unresmatch, unresolved
   print, ''
   print, '[WARNING] ################################################################'
   print, '[WARNING] ## The following routines are missing in your IDL setup:      ##'
   print, '[WARNING] ##',unresolved
   print, '[WARNING] ## This can cause to code to behave improperly.               ##' 
   print, '[WARNING] ##               Continue at your own risk!                   ##'
   print, '[WARNING] ################################################################'
endif


;Check IDL version. The warning is printed to the terminal even if log_lun is set.
if float(!version.release) lt 6.1 then begin
  print, ''
  print, '[WARNING] ###############################################################'
  print, '[WARNING] ##    This code has been developed for IDL v6.1 or greater   ##' 
  print, '[WARNING] ##           Some features might now work correctly!         ##'
  print, '[WARNING] ###############################################################' 
endif

if !version.memory_bits eq 32 then begin
  print, ''
  print, '[WARNING] ###############################################################'
  print, '[WARNING] ##               This is a 32 bit IDL version.               ##' 
  print, '[WARNING] ##           Some features might now work correctly!         ##'
  print, '[WARNING] ###############################################################'
endif

;Check that MPFIT is in the idl path.
if file_which('mpfit.pro') eq '' then begin
  print, ''
  print, '[ ERROR ] ###############################################################'
  print, '[ ERROR ] ## It seems you do not have MPFIT in your IDL_PATH           ##'
  print, '[ ERROR ] ## Please download the most recent version from:             ##'
  print, '[ ERROR ] ##   https://www.physics.wisc.edu/~craigm/idl/fitting.html   ##'
  print, '[ ERROR ] ###############################################################'
  retall
endif

;Now check MPFIT version
mpfit_ok = 0 
mpfitver = '<unknown>'
CATCH, CATCHERR
IF CATCHERR EQ 0 THEN mpfit_ok = mpfit(/query, VERSION=mpfitver, MIN_VERSION='1.70')
CATCH, /CANCEL

if not mpfit_ok then begin
  print, ''
  print, '[ ERROR ] ###############################################################'
  print, '[ ERROR ] ## You must have MPFIT v1.70 or higher in your IDL_PATH      ##'
  print, '[ ERROR ] ## Version found: '+string(mpfitver,format='(A-43)')+'##'
  print, '[ ERROR ] ## Please download the most recent version from:             ##'
  print, '[ ERROR ] ##   https://www.physics.wisc.edu/~craigm/idl/fitting.html   ##'
  print, '[ ERROR ] ###############################################################'
  retall
endif  

;First check that astrolib exists. If mrdfits.pro is not present then print error
if file_which('mrdfits.pro') eq '' then begin
  print, ''
  print, '[ ERROR ] ###############################################################'
  print, '[ ERROR ] ## It seems you have no ASTROLIB in your IDL_PATH            ##'
  print, '[ ERROR ] ## Please download the most recent version from:             ##'
  print, '[ ERROR ] ##        http://idlastro.gsfc.nasa.gov/contents.html        ##'
  print, '[ ERROR ] ###############################################################'
  retall
endif

;Then attempt to check astrolib version by testing the airtovac routine that fails
;if astrolib is older than 2011.
astrolib_test = 0
astrolib_dir = file_dirname(file_which('airtovac.pro'))
CATCH, astrolib_test
IF astrolib_test EQ 0 THEN airtovac, 5000, output
CATCH, /cancel


if astrolib_test  NE 0 then begin
  print, ''
  print, '[ ERROR ] ###############################################################'
  print, '[ ERROR ] ## It seems your ASTROLIB installed at:                      ##'
  print, '[ ERROR ] ## '+string(astrolib_dir,format='(A-58)')+'##'
  print, '[ ERROR ] ## has not been updated from 2011.                           ##'
  print, '[ ERROR ] ## Please download the most recent version from:             ##'
  print, '[ ERROR ] ##        http://idlastro.gsfc.nasa.gov/contents.html        ##'
  print, '[ ERROR ] ###############################################################'
  retall
endif 

end



;-----------------------------------------------------------------------
; main program
;-----------------------------------------------------------------------

pro kubeviz, datafile, noisefile=noisefile, ext=ext, noise_ext=noise_ext, $ 
             trim=trim, fluxfac=fluxfac, smooth=smooth, $
	     specsmooth=specsmooth, transpose=transpose, help=help, version=version, $
             vacuum=vacuum, logarithmic=logarithmic, waveunit=waveunit, $
             redshift=redshift, do_mc_errors=do_mc_errors, bootstrap=bootstrap, $
             nmontecarlo=nmontecarlo, plot_mc_pdf=plot_mc_pdf, save_mc_pdf=save_mc_pdf, $ 
	     use_mc_noise=use_mc_noise, scroll=scroll, lineset=lineset, debug=debug, $
             mask_sn_thresh=mask_sn_thresh, mask_maxvelerr=mask_maxvelerr, $
	     mom_thresh=mom_thresh, instr=instr, band=band, $
             batch=batch, fit_all_lines=fit_all_lines, fix_ratios=fix_ratios, $
	     resfile=resfile, fittype=fittype, spmask=spmask, outdir=outdir, $
	     logunit=logunit, fitpars=fitpars, momwindowfile=momwindowfile, $
	     gaussinitfile=gaussinitfile, scalenoiseerr=scalenoiseerr


compile_opt idl2, hidden
!except = 0 ; Remove arithmetic underflow/overflow error messages

kubeviz_libtest

;Initialize common blocks and state structures
common kubeviz_state
common kubeviz_winsize
kubeviz_common 

cd, current = cwdir
state.cwdir = cwdir+'/'

if keyword_set(outdir) then state.outdir = outdir else state.outdir = state.cwdir
if file_test(state.outdir) eq 0 then file_mkdir, state.outdir

if keyword_set(logunit) then state.log_lun = logunit

if keyword_set(help) then begin
   doc_library, 'kubeviz'
   return
endif

if n_elements(version) ne 0 then begin
   printf, state.log_lun, 'Version: ', state.version
   return
endif

;-----------------------------------------------------------------------
printf, state.log_lun, ''
printf, state.log_lun, '                              *** KUBEVIZ ***'
printf, state.log_lun, ''
printf, state.log_lun, 'version: ', state.version
printf, state.log_lun, ''
;-----------------------------------------------------------------------

if n_elements(datafile) eq 0 then begin
   kubeviz_selectcube, dataname, dir, noisename, noisedir
   if n_elements(dataname) eq 0 then begin
      print, '[ ERROR ] no datacube file selected'
      return
   endif
   if n_elements(noisename) eq 0 then begin
      print, '[ ERROR ] no noisecube file selected'
      return
   endif
endif else begin
   kubeviz_splitpath, datafile, dataname, dir
   if file_test(dir+dataname) eq 0 then begin
      print, '[ ERROR ] Unable to locate file: '+dir+dataname
      return
   endif
   if strmid(dataname,strlen(dataname)-4) eq '.sav'  then begin
      kubeviz_loadsession, dataname, dir, spmask=spmask
      return
   endif
   if n_elements(noisefile) gt 0 then begin
      kubeviz_splitpath, noisefile, noisename, noisedir
      if strlen(noisedir) eq 0 then noisedir=dir
      if file_test(noisedir+noisename) eq 0 then begin
        print, '[ ERROR ] Unable to locate file: '+noisedir+noisename
        return
      endif
   endif else begin
      kubeviz_selectcube, dataname, dir, noisename, noisedir, /noiseonly
      if n_elements(noisename) eq 0 then begin
        print, '[ ERROR ] no noisecube file selected'
        return
      endif
   endelse 
endelse

state.filename = dataname
state.indir = dir

; Command line selection of the error method to be used

if (n_elements(bootstrap) ne 0) then begin 
     do_mc_errors = 1
     bootstrap_file = bootstrap
endif     

if (n_elements(do_mc_errors) ne 0) then begin
    case do_mc_errors of 
    0: state.domontecarlo = 0
    1: begin
	 kubeviz_splitpath, bootstrap_file, bootstrap_fname, bootstrap_path
         if strlen(bootstrap_fname) eq 0 then begin
	    printf, state.log_lun,  '[WARNING] invalid bootstrap filename. Use Noise Cube errors instead.' 
	 endif else begin 
	   state.domontecarlo = 1
	   ;If the bootstrap path is not given in the bootstrap keyword it is assumed to be the same as the datacube
	   if strlen(bootstrap_fname)-strlen(bootstrap) eq 0 then bootstrap_path = state.indir
	 endelse  
       end 
    2: begin
         state.domontecarlo = 2
         if n_elements(nmontecarlo) ne 0 then state.Nmc1 = nmontecarlo else state.Nmc1 = 100
       end 
    3: begin
         state.domontecarlo = 3  
         if n_elements(nmontecarlo) ne 0 then state.Nmc2 = nmontecarlo else state.Nmc2 = 100
       end
     4: begin
         state.domontecarlo = 4  
         if n_elements(nmontecarlo) ne 0 then state.Nmc3 = nmontecarlo else state.Nmc3 = 100
       end       
    endcase  
endif

if N_elements(resfile) gt 0 then begin
    kubeviz_splitpath, resfile, res_name, res_path
    if strlen(res_name) eq 0 then printf, state.log_lun, '[WARNING] Invalid results filename. ' else begin 
    ;If the results file path is not given it is assumed to be the same as the datacube
    validres = 1
    if strlen(res_path) eq 0 then res_path = state.indir
    endelse
endif      

if n_elements(plot_mc_pdf)    eq 0 then state.plotMonteCarlodistrib=0 else state.plotMonteCarlodistrib=1
if n_elements(save_mc_pdf)    eq 0 then state.saveMonteCarlodistrib=0 else state.saveMonteCarlodistrib=1
if n_elements(use_mc_noise)   eq 0 then state.useMonteCarlonoise=0    else state.useMonteCarlonoise=1
if n_elements(instr)          ne 0 then state.instr=instr
if n_elements(band)           ne 0 then state.band=band
if n_elements(fluxfac)        ne 0 then state.fluxfac=fluxfac
if n_elements(scroll)         eq 0 then scroll=0
if n_elements(lineset)        eq 0 then state.selected_lineset=0 else state.selected_lineset=lineset
if n_elements(smooth)         eq 0 then state.smooth=1 else state.smooth=smooth
if n_elements(specsmooth)     eq 0 then state.specsmooth=1 else state.specsmooth=specsmooth
if n_elements(transpose)      ne 0 then state.transpose = 1
if n_elements(vacuum)         ne 0 then state.vacuum = 1
if n_elements(debug)          ne 0 then state.debug= 1
if n_elements(redshift)       gt 0 then state.redshift = redshift
if n_elements(mask_sn_thresh) gt 0 then state.mask_sn_thresh = mask_sn_thresh
if n_elements(mask_maxvelerr) gt 0 then state.mask_maxvelerr = mask_maxvelerr
if n_elements(mom_thresh)     gt 0 then state.mom_thresh = mom_thresh
if n_elements(fix_ratios)     ne 0 then state.fitfixratios=fix_ratios
if n_elements(scalenoiseerr)  ne 0 then state.scalenoiseerrors=scalenoiseerr
if n_elements(fittype)        gt 0 then begin
  case fittype of
    'gauss' : state.linefit_type = 0
    'moments': state.linefit_type = 1
    else: printf, state.log_lun, '[WARNING]: fittype keyword value not recognized. Use "gauss" instead. '
    endcase
endif  
if n_elements(fitpars) gt 0 then kubeviz_set_user_fitpars, fitpars

state.zcuts = 4                  ; start with hist eq selected zcut
state.zoomrange = 32             ; start half-range of spectral zoom window

;Create the file where to store fit err msg
errmsg_file_info = file_info(state.cwdir+'fit_errmsg_'+strmid(state.filename,0,strlen(state.filename)-5)+'.txt')
if errmsg_file_info.exists eq 0 then spawn, 'touch '+state.cwdir+'fit_errmsg_'+strmid(state.filename,0,strlen(state.filename)-5)+'.txt'

;open the data cube
kubeviz_getdata, state.indir, state.filename, noisedir=noisedir, noisefile=noisefile, ext=ext, noise_ext=noise_ext, $
		 trim=trim, logarithmic=logarithmic, waveunit=waveunit, bootstrap_fname=bootstrap_fname, momwindowfile=momwindowfile, $
		 gaussinitfile=gaussinitfile
;initialize linefit:
kubeviz_linefit_init

;If a montecarlo method is set then create the corresponding cubes
case state.domontecarlo of
    1: kubeviz_readbootstrapcubes, bootstrap_fname, dir=bootstrap_path, momwindowfile=momwindowfile
    2: kubeviz_createmc1cubes
    3: kubeviz_createmc2cubes
    4: kubeviz_createmc3cubes
    else:
endcase

; mask bad pixels in bootstrap cubes ( & setup percentile arrays):
if state.domontecarlo gt 0 then kubeviz_setupmontecarlocubes

;read spaxels mask if provided
if n_elements(spmask) gt 0 then kubeviz_load_mask, fname=spmask

; Set the kubeviz main mode (0= interactive, 1=automatic)
if N_elements(batch) eq 0 then begin 
  ; create the windows
  if N_elements(validres) gt 0 then  kubeviz_linefit_loadres, resfile, dir=res_path
  kubeviz_create, fname, scroll=scroll
endif else  kubeviz_batchmode, redshift=redshift, fit_all_lines= fit_all_lines

end

