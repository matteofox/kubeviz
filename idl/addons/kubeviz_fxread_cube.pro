FUNCTION KUBEVIZ_FXREAD_CUBE, FILENAME, HEADER, P1, P2, P3, P4, P5, P6,        $
		NANVALUE=NANVALUE, PROMPT=PROMPT, 	$
		NOSCALE=NOSCALE, NOUPDATE=NOUPDATE,	$
		ERRMSG=ERRMSG, NODATA=NODATA, COMPRESS = COMPRESS,      $
		EXTENSION=EXTENSION0, SILENT=SILENT
;+
; NAME: 
;	KUBEVIZ_FXREAD_CUBE
; Purpose     : 
;	Read basic FITS files. Optimized to read datacubes with 3 dimesions
; Explanation : 
;	Read an image array from a disk FITS file.  Optionally allows the
;	user to read in only a subarray.
; Use         : 
;	FXREAD_CUBE, FILENAME, [, HEADER  [, I1, I2  [, J1, J2 ] [, K1, K2 ]]  
; Inputs      : 
;	FILENAME = String containing the name of the file to be read.
; Opt. Inputs : 
;	I1,I2	 = Data range to read in the first dimension.  If passed, then
;		   HEADER must also be passed.  If not passed, or set to -1,-1,
;		   then the entire range is read.
;	J1,J2	 = Data range to read in the second dimension.  If passed, then
;		   HEADER and I1,J2 must also be passed.  If not passed, or set
;		   to -1,-1, then the entire range is read.
;	K1,K2	 = Data range to read in the second dimension.  If passed, then
;		   HEADER and I1,J2 must also be passed.  If not passed, or set
;		   to -1,-1, then the entire range is read.
; Outputs     : 
;	DATA	 = Data array to be read from the file.
; Opt. Outputs: 
;	HEADER	 = String array containing the header for the FITS file.
; Keywords    : 
;       /COMPRESS - If this keyword is set and non-zero, then then treat
;                the file as gzip compressed.    By default FXREAD assumes
;                the file is gzip compressed if it ends in ".gz"
;	NANVALUE = Value signalling data dropout.  All points corresponding to
;		   IEEE NaN (not-a-number) are set to this value.  Ignored
;		   unless DATA is of type float or double-precision.
;       EXTENSION = FITS extension.  It can be a scalar integer,
;                indicating the extension number (extension number 0
;                is the primary HDU).  It can also be a scalar string,
;                indicating the extension name (EXTNAME keyword).
;                Default: 0 (primary HDU)
;	PROMPT	 = If set, then the optional parameters are prompted for at the
;		   keyboard.
;	NOSCALE	 = If set, then the output data will not be scaled using the
;		   optional BSCALE and BZERO keywords in the FITS header.
;		   Default is to scale, if and only if BSCALE and BZERO are
;		   present and nontrivial.
;	NOUPDATE = If set, then the optional BSCALE and BZERO keywords in the
;		   optional HEADER array will not be changed.  The default is
;		   to reset these keywords to BSCALE=1, BZERO=0.  Ignored if
;		   NOSCALE is set.
;	ERRMSG   = If defined and passed, then any error messages will be
;		   returned to the user in this parameter rather than
;		   depending on the MESSAGE routine in IDL.  If no errors are
;		   encountered, then a null string is returned.  
;       NODATA   = If set, then the array is not read in, but the
;                  primary header is read.
;       SILENT   = suppress error messages and informational messages
;
; Calls       : 
;	GET_DATE, FXADDPAR, FXHREAD, FXPAR, WHERENAN
; Common      : 
;	None.
; Restrictions: 
;	Groups are not supported.
;
;	The optional parameters I1, I2, only work with one two or there
;	dimensional arrays.  J1 and J2 only work with two- or three-dimensional
;	arrays. K1 and K2 only work with three-dimensional arrays
;
;
; Side effects: 
;	If the keywords BSCALE and BZERO are present in the FITS header, and
;	have non-trivial values, then the retureed array DATA is formed by the
;	equation
;
;			DATA = BSCALE*original + BZERO
;
;	However, this behavior can overridden by using the /NOSCALE keyword.
;
;	If the data is scaled, then the optional HEADER array is changed so
;	that BSCALE=1 and BZERO=0.  This is so that these scaling parameters
;	are not applied to the data a second time by another routine.  Also,
;	history records are added storing the original values of these
;	constants.  Note that only the returned array is modified--the header
;	in the FITS file itself is untouched.
;
;	If the /NOUPDATE keyword is set, however, then the BSCALE and BZERO
;	keywords are not changed.  It is then the user's responsibility to
;	ensure that these parameters are not reapplied to the data.  In
;	particular, these keywords should not be present in any header when
;	writing another FITS file, unless the user wants their values to be
;	applied when the file is read back in.  Otherwise, FITS readers will
;	read in the wrong values for the data array.
;	
; Category    : 
;	Data Handling, I/O, FITS, Generic.
; Prev. Hist. : 
;	W. Thompson, May 1992, based in part on READFITS by W. Landsman, and
;			       STSUB by M. Greason and K. Venkatakrishna.
;	W. Thompson, Jun 1992, added code to interpret BSCALE and BZERO
;			       records, and added NOSCALE and NOUPDATE
;			       keywords.
;	W. Thompson, Aug 1992, changed to call FXHREAD, and to add history
;			       records for BZERO, BSCALE.
; Minimium IDL Version:
;       V6.0 (uses V6.0 notation) 
; Written     : 
;	William Thompson, GSFC, May 1992.
; Modified    : 
;	Version 1, William Thompson, GSFC, 12 April 1993.
;		Incorporated into CDS library.
;	Version 2, William Thompson, GSFC, 17 November 1993.
;		Corrected bug with AVERAGE keyword on non-IEEE compatible
;		machines.
;		Corrected bug with subsampling on VAX machines.
;	Version 3, William Thompson, GSFC, 31 May 1994
;		Added ERRMSG keyword.
;       Version 4, William Thompson, GSFC, 23 June 1994
;               Modified so that ERRMSG is not touched if not defined.
;       Version 5, Zarro (SAC/GSFC), 14 Feb 1997 
;               Added I/O error checking
;       Version 6, 20-May-1998, David Schlegel/W. Thompson
;               Allow a single pixel to be read in.
;               Change the signal to read in the entire array to be -1
;       Version 7 C. Markwardt 22 Sep 2003
;               If the image is empty (NAXIS EQ 0), or NODATA is set, then
;               return only the header.  
;       Version 8 W. Landsman  29 June 2004
;               Added COMPRESS keyword, check for .gz extension  
;       Version 9, William Thompson, 19-Aug-2004
;               Make sure COMPRESS is treated as a scalar
;       Version 10, Craig Markwardt, 01 Mar 2004
;               Add EXTENSION keyword and ability to read different
;               extensions than the primary one.
;       Version 11,  W. Landsman   September 2006 
;               Assume since V5.5, remove VMS support
;       Version 11.1,  W. Landsman   November 2007
;               Allow for possibility number of bytes requires 64 bit integer
;       Version 12, William Thompson, 18-Jun-2010, update BLANK value.
;       Version 13, W. Landsman  Remove IEEE_TO_HOST, V6.0 notation
;       Version 14, M.Fossati, 13 Jun 2014  Removed AVERAGE and STEP, extend the 
;               subarray read to 3D datacubes. Converted into a function.
;               Introduced SILENT keyword.
;-
;
	ERRMSG=''
	ON_ERROR, 2
;    
;  This parameter will be used later.
;
	ALREADY_CONVERTED = 0
        READ_OK=0
;
;  Parse the input parameters.
;       
	CASE N_PARAMS() OF
		1:  BEGIN & I1=-1 & I2=-1 & J1=-1 & J2=-1 & K1=-1 & K2=-1 & END
		2:  BEGIN & I1=-1 & I2=-1 & J1=-1 & J2=-1 & K1=-1 & K2=-1 & END
		3:  BEGIN & I1=-1 & I2=-1 & J1=-1 & J2=-1 & K1=-1 & K2=-1 & END
		4:  BEGIN & I1=P1 & I2=P2 & J1=-1 & J2=-1 & K1=-1 & K2=-1 & END
		5:  BEGIN & I1=P1 & I2=P2 & J1=-1 & J2=-1 & K1=-1 & K2=-1 & END
		6:  BEGIN & I1=P1 & I2=P2 & J1=P3 & J2=P4 & K1=-1 & K2=-1 & END
		7:  BEGIN & I1=P1 & I2=P2 & J1=P3 & J2=P4 & K1=-1 & K2=-1 & END
		8:  BEGIN & I1=P1 & I2=P2 & J1=P3 & J2=P4 & K1=P5 & K2=P6 & END
	        9:  BEGIN & I1=P1 & I2=P2 & J1=P3 & J2=P4 & K1=P5 & K2=P6 & END
		ELSE:  BEGIN
			ERRMSG = 'Syntax:  FXREAD_CUBE(FILENAME, ' + $
				'[, HEADER [, I1, I2 [, J1, J2 ] [, K1, K2 ] ])'
			MESSAGE, ERRMSG, /con
		        return, -1
			END
	ENDCASE

	;; Extension number	
	IF N_ELEMENTS(EXTENSION0) EQ 0 THEN EXTENSION = 0L $
	ELSE EXTENSION = EXTENSION0[0]

	SZ = SIZE(EXTENSION)
	ETYPE = SZ[SZ[0]+1]
	IF ETYPE EQ 8 THEN BEGIN
		ERRMSG = 'EXTENSION must not be a structure'
		MESSAGE, ERRMSG, /con
		return, -1
	ENDIF


;
;  Determine if file is compressed, get the UNIT number, and open the file.
;
        IF NOT KEYWORD_SET(COMPRESS) THEN $
         COMPRESS = STRLOWCASE( STRMID(FILENAME, STRLEN(FILENAME)-3,3)) EQ '.gz'
	OPENR, UNIT, FILENAME, /GET_LUN, ERROR=ERROR,COMPRESS=COMPRESS[0], /SWAP_IF_LITTLE_ENDIAN
        IF ERROR NE 0 THEN BEGIN
	    ERRMSG='Error opening '+FILENAME
	    MESSAGE, ERRMSG, /con
            return, -1
        ENDIF
;
;  Read in the FITS header.
;

	;; Starting extension number is zero
	I_EXT = 0L
	FOUND_EXT = 0

        WHILE NOT FOUND_EXT DO BEGIN
            FXHREAD,UNIT,HEADER,STATUS
            IF STATUS NE 0 THEN BEGIN
               FREE_LUN,UNIT
                ERRMSG = 'Unable to read requested FITS header extension'
                MESSAGE, ERRMSG, /con
		return, -1
            ENDIF
;
;  Extract the keywords BITPIX, NAXIS, NAXIS1, ...
;
            START = 0L
            BITPIX = FXPAR(HEADER,'BITPIX', START=START)
            NAXIS = FXPAR(HEADER,'NAXIS', START=START)
            GCOUNT = FXPAR(HEADER,'GCOUNT', START=START)
            IF GCOUNT EQ 0 THEN GCOUNT = 1
            PCOUNT = FXPAR(HEADER,'PCOUNT', START=START)
            IF NAXIS GT 0 THEN BEGIN 
                DIMS = LONG64(FXPAR(HEADER,'NAXIS*')) ;Read dimensions
                NDATA = DIMS[0]
                IF NAXIS GT 1 THEN FOR I=2,NAXIS DO NDATA = NDATA*DIMS[I-1]
            ENDIF ELSE NDATA = 0
            NBYTES = LONG64(ABS(BITPIX) / 8) * GCOUNT * (PCOUNT + NDATA)
            NREC = (NBYTES + 2879L) / 2880L
	    
            IF ETYPE EQ 7 THEN BEGIN
                EXTNAME = STRTRIM(STRUPCASE(FXPAR(HEADER,'EXTNAME', $
                                                  START=START)),2)
                IF EXTNAME EQ EXTENSION THEN FOUND_EXT = 1
            END ELSE IF I_EXT EQ EXTENSION THEN FOUND_EXT = 1

            IF NOT FOUND_EXT THEN BEGIN
                ;; Check to be sure there are extensions
                IF I_EXT EQ 0 THEN BEGIN
                    IF NOT FXPAR(HEADER,'EXTEND', START=START) THEN BEGIN
		        FREE_LUN,UNIT
                        ERRMSG = 'Requested extension not found, and file ' + $
                          FILENAME + ' does not contain extensions'
                        MESSAGE, ERRMSG, /con
		        return, -1
                    ENDIF
                ENDIF
                
	        POINT_LUN, -UNIT, POINTLUN		;Current position
                MHEAD0 = POINTLUN + NREC*2880L
	        POINT_LUN, UNIT, MHEAD0			;Next FITS extension

                I_EXT++
            ENDIF
        ENDWHILE

        ;;
        ;; If we got here, then we have arrived at the requested
        ;; extension.  We still need to be sure that it is an image
        ;; and not a table (for extensions beyond the primary one,
        ;; that is).
        ;;
        IF I_EXT GT 0 THEN BEGIN
            XTENSION = STRTRIM(STRUPCASE(FXPAR(HEADER,'XTENSION', START=START)),2)
            IF (XTENSION NE 'IMAGE') THEN BEGIN
		FREE_LUN,UNIT
                ERRMSG = 'Extension ' + STRTRIM(EXTENSION,2) +		$
                  ' is not an image'
                MESSAGE, ERRMSG, /con
		return, -1
            ENDIF
        ENDIF
            
        ;; Handle case of empty image, or no data requested
        IF NAXIS EQ 0 OR KEYWORD_SET(NODATA) THEN BEGIN
            ;; Make DATA an undefined variable, reflecting no data
            ERRMSG = ''
            FREE_LUN,UNIT
            return, -1
        ENDIF

	DIMS = FXPAR(HEADER,'NAXIS*')
	N1 = DIMS[0]
	CASE NAXIS OF
	   2: BEGIN & N2 = DIMS[1] & N3 = 1       & END
	   3: BEGIN & N2 = DIMS[1] & N3 = DIMS[2] & END
	   ELSE: BEGIN & N2 = 1 & N3 = 1 & END
	ENDCASE   
;
;  Determine the array type from the keyword BITPIX.
;
	CASE BITPIX OF
		  8:	IDLTYPE = 1	; Byte
		 16:	IDLTYPE = 2	; Integer*2
		 32:	IDLTYPE = 3	; Integer*4
		-32:	IDLTYPE = 4	; Real*4
		-64:	IDLTYPE = 5	; Real*8
	ENDCASE
;
;  Set the default values for the optional parameters.
;
	IF (I1 EQ -1) && (I2 EQ -1) THEN BEGIN
           I1 = 0
           I2 = N1-1
        ENDIF
	IF (J1 EQ -1) && (J2 EQ -1) THEN BEGIN
           J1 = 0
           J2 = N2-1
        ENDIF
	IF (K1 EQ -1) && (K2 EQ -1) THEN BEGIN
           K1 = 0
           K2 = N3-1
        ENDIF
;
;  If the prompt keyword was set, the prompt for the parameters.
;
	IF KEYWORD_SET(PROMPT) THEN BEGIN
		ANSWER = ''
		READ,'Enter lower limit for X ['+STRTRIM(I1,2)+']: ', ANSWER
		IF ANSWER NE '' THEN I1 = (ANSWER)
;
		ANSWER = ''
		READ,'Enter upper limit for X ['+STRTRIM(I2,2)+']: ', ANSWER
		IF ANSWER NE '' THEN I2 = LONG(ANSWER)
;
		ANSWER = ''
		READ,'Enter lower limit for Y ['+STRTRIM(J1,2)+']: ', ANSWER
		IF ANSWER NE '' THEN J1 = LONG(ANSWER)
;
		ANSWER = ''
		READ,'Enter upper limit for Y ['+STRTRIM(J2,2)+']: ', ANSWER
		IF ANSWER NE '' THEN J2 = LONG(ANSWER)
		
		ANSWER = ''
		READ,'Enter lower limit for Z ['+STRTRIM(K1,2)+']: ', ANSWER
		IF ANSWER NE '' THEN K1 = LONG(ANSWER)
;
		ANSWER = ''
		READ,'Enter upper limit for Z ['+STRTRIM(K2,2)+']: ', ANSWER
		IF ANSWER NE '' THEN K2 = LONG(ANSWER)

	ENDIF
	
;
;  If any of the optional parameters were passed, then update the dimensions
;  accordingly.  First check I1 and I2.
;
	IF (I1 NE 0) || (I2 NE N1-1) THEN BEGIN
		IF NAXIS GT 3 THEN BEGIN
			FREE_LUN,UNIT
			ERRMSG = 'Range parameters can only be set for ' + $
				'one, two or three dimensional arrays'
			MESSAGE, ERRMSG, /con
		        return, -1
		ENDIF
		IF (MIN([I1,I2]) LT 0) OR (MAX([I1,I2]) GE DIMS[0]) THEN BEGIN
			FREE_LUN,UNIT
			ERRMSG = 'I1,I2 must be in the range 0 to ' +	$
				STRTRIM(DIMS[0]-1,2)
			MESSAGE, ERRMSG, /con
		        return, -1
		END ELSE IF I1 GT I2 THEN BEGIN
			ERRMSG = 'I2 must be >= I1'
			MESSAGE, ERRMSG, /con
		        return, -1
		ENDIF
		DIMS[0] = I2 - I1 + 1
	ENDIF
;
;  Next, check J1 and J2.
;
	IF (J1 NE 0) || (J2 NE N2-1) THEN BEGIN
		IF NAXIS lt 2 THEN BEGIN
			FREE_LUN,UNIT
			ERRMSG = 'J1, J2 can only be set for ' +	$
				'two or three-dimensional arrays'
			MESSAGE, ERRMSG, /con
		        return, -1
		ENDIF
		IF (MIN([J1,J2]) LT 0) OR (MAX([J1,J2]) GE DIMS[1]) THEN BEGIN
			FREE_LUN,UNIT
			ERRMSG = 'J1,J2 must be in the range 0 to ' +	$
				STRTRIM(DIMS[1]-1,2)
			MESSAGE, ERRMSG, /con
		        return, -1
		END ELSE IF J1 GT J2 THEN BEGIN
			ERRMSG = 'J2 must be >= J1'
			MESSAGE, ERRMSG, /con
		        return, -1
		ENDIF
		DIMS[1] = J2 - J1 + 1
	ENDIF

;
;  Next, check K1 and K2.
;
	IF (K1 NE 0) || (K2 NE N3-1) THEN BEGIN
		IF NAXIS NE 3 THEN BEGIN
			FREE_LUN,UNIT
			ERRMSG = 'K1, K2 can only be set for ' +	$
				'three-dimensional arrays'
			MESSAGE, ERRMSG, /con
		        return, -1
		ENDIF
		IF (MIN([K1,K2]) LT 0) OR (MAX([K1,K2]) GE DIMS[2]) THEN BEGIN
			FREE_LUN,UNIT
			ERRMSG = 'K1,K2 must be in the range 0 to ' +	$
				STRTRIM(DIMS[2]-1,2)
			MESSAGE, ERRMSG, /con
		        return, -1
		END ELSE IF K1 GT K2 THEN BEGIN
			ERRMSG = 'K2 must be >= K1'
			MESSAGE, ERRMSG, /con
		        return, -1
		ENDIF
		DIMS[2] = K2 - K1 + 1
	ENDIF	
		
;
;  Make the array.
;
	DATA = MAKE_ARRAY(DIMENSION=DIMS,TYPE=IDLTYPE,/NOZERO)
;
;  Find the start of the data to be read in.
;       

	POINT_LUN,-UNIT,OFFSET		;Current position
	DELTA = LONG64(N1)*N2*ABS(BITPIX)/8
	IF K1 NE 0 THEN BEGIN
		OFFSET = OFFSET + K1*DELTA
		POINT_LUN,UNIT,OFFSET
	ENDIF

;
;  If the I, J range is non-trivial, then read in the file plane by plane.  
;
        ON_IOERROR,QUIT
	IF (DIMS[0] NE N1) || (DIMS[1] ne N2) THEN BEGIN
	    
	    if N_elements(SILENT) eq 0 THEN MESSAGE,'Reading the FITS extension plane by plane', /INFO 
	     
	    CASE NAXIS OF
	     1: NK=1
	     2: NK=1 
	     3: NK=DIMS[2] 
	    ENDCASE
	     
	    FOR K = 0,NK-1 DO BEGIN

		   PLANE = MAKE_ARRAY(N1,N2,TYPE=IDLTYPE,/NOZERO)
	           READU,UNIT,PLANE
;
;  If I1,I2 do not match the array size, then extract the relevant subarray.
;
	           IF (I1 NE 0) || (I2 NE N1-1) || (J1 NE 0) || (J2 NE N2-1) THEN PLANE = PLANE[I1:I2,J1:J2]

;  Now store the PLANE in the data array.

	           DATA[*,*,K] = PLANE
	    ENDFOR
;
;  Otherwise, if the file doesn't have to be read in plane by plane, then just
;  read the data array.
;
	END ELSE BEGIN
	    if N_elements(SILENT) eq 0 THEN MESSAGE,'Reading the FITS extension in one go.', /INFO 
	    READU,UNIT,DATA
	ENDELSE    

;
;  Convert the data from IEEE to host format, keeping track of any IEEE NaN
;  values.  Don't do this if the conversion has already taken place.

       IF ~ALREADY_CONVERTED THEN BEGIN
	       IF (N_ELEMENTS(NANVALUE) EQ 1) && (IDLTYPE GE 4) &&     $
		       (IDLTYPE LE 6) THEN W = WHERENAN(DATA,COUNT) ELSE $
		       COUNT = 0
	       ;SWAP_ENDIAN_INPLACE,DATA, /SWAP_IF_LITTLE
       END ELSE COUNT = 0


;  If the parameters BZERO and BSCALE are non-trivial, then adjust the array by
;  these values.  Also update the BLANK keyword, if present.

       IF ~KEYWORD_SET(NOSCALE) THEN BEGIN
	       BZERO  = FXPAR(HEADER,'BZERO')
	       BSCALE = FXPAR(HEADER,'BSCALE')
	       BLANK  = FXPAR(HEADER,'BLANK',COUNT=NBLANK)
	       GET_DATE,DTE
	       IF (BSCALE NE 0) && (BSCALE NE 1) THEN BEGIN
		       DATA *= BSCALE
		       IF ~KEYWORD_SET(NOUPDATE) THEN BEGIN
			   FXADDPAR,HEADER,'BSCALE',1.
			   FXADDPAR,HEADER,'HISTORY',DTE +	       $
			     ' applied BSCALE = '+ STRTRIM(BSCALE,2)
			   IF NBLANK EQ 1 THEN BEGIN
			       print, bscale, blank
			       BLANK *= BSCALE
			       FXADDPAR,HEADER,'BLANK',BLANK
			   ENDIF
		       ENDIF
	       ENDIF
	       IF BZERO NE 0 THEN BEGIN
		       DATA += BZERO
		       IF ~KEYWORD_SET(NOUPDATE) THEN BEGIN
			   FXADDPAR,HEADER,'BZERO',0.
			   FXADDPAR,HEADER,'HISTORY',DTE +	       $
			     ' applied BZERO = '+ STRTRIM(BZERO,2)
			   IF NBLANK EQ 1 THEN BEGIN
			       BLANK += BZERO
			       FXADDPAR,HEADER,'BLANK',BLANK
			   ENDIF
		       ENDIF
	       ENDIF
       ENDIF
;
;  Store NANVALUE everywhere where the data corresponded to IEE NaN.
;
	IF COUNT GT 0 THEN DATA[W] = NANVALUE
;
;  Close the file and return
;
        READ_OK=1

QUIT:   ON_IOERROR,NULL
	FREE_LUN, UNIT
        IF NOT READ_OK THEN BEGIN
	    ERRMSG='Error reading file '+FILENAME
	    MESSAGE, ERRMSG, /con
            return, -1
	ENDIF
	IF N_ELEMENTS(ERRMSG) NE 0 THEN ERRMSG = ''
	return, DATA
	END
