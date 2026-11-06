; $Id: xpdmenu.pro,v 1.10 2005/02/01 20:26:16 scottm Exp $
;
; Copyright (c) 1991-2005, Research Systems, Inc.  All rights reserved.
;	Unauthorized reproduction prohibited.
;+
; NAME:
;	XPDMENU
;
; PURPOSE:
;	This procedure implifies setting up widget pulldown menus. XPDMENU
;	reads a description of the menu to be generated, and calls
;	the appropriate widget creation functions to generate it.
;
; CALLING SEQUENCE:
;	XPDMENU, Desc, Parent
;
; INPUTS:
;	DESC:	Either the name of a file that contains the description of the
;		pulldown menu to be generated, or a string array that
;		contains the description.  The rules for a pull-down menu
;		description are as follows:
;
;		Leading and trailing whitespace is ignored.  Lines starting
;		with the '#' character or blank lines are ignored.  All other
;		lines contain 2 fields, a button label and a value.  The label
;		should be quoted with any desired delimiter, usually single
;		or double quotes.  The value can be omitted, in which case the
;		label is used as the value.  To make a menu choice reveal
;		another pull-down menu, place a '{' character in its value
;		field.  Such a pulldown is terminated by a line containing
;		a '}' in the label field.
;
;		Example:
;			"Colors" {
;			    "Red"
;			    "Green"
;			    "Blue"	{
;				"Light"
;				"Medium"
;				"Dark"
;				"Navy"
;				"Royal"
;			    }
;			    "Cyan"
;			    "Magenta"
;			}
;			"Quit"		DONE
;
;		This example builds a menu bar with 2 buttons,
;		named "Colors" and "Quit". "Colors" is a pulldown
;		containing "Red", "Green", "Blue", "Cyan", and "Magenta".
;		"Blue" is a sub-pulldown containing shades of blue.
;		Such sub-menus can be nested to any desired level.
;		Most of the lines don't specify an explicit value. The
;		exception is "Quit", which has the value "DONE". It can
;		be instructive to run the following small program:
;
;			a = WIDGET_BASE()
;			XPDMENU, a, 'test'	; Test contains the above
;			widget_control, /REALIZE, a
;			uvalue=''
;			repeat begin
;			  event = widget_event(a)
;			  WIDGET_CONTROL, get_uvalue=uvalue, event.id
;			  print, uvalue
;			end until uvalue eq "EXIT"
;			WIDGET_CONTROL, /destroy, a
;			end
;
;		Note that if you choose to make DESC be a string array,
;		the arrays contents must be exactly the same as the file
;		would be (including the quotes around the fields). Each
;		element of the array corresponds to one line of a file.
;
;	PARENT:	Widget ID of the parent base widget for the pulldown menu.
;		If this argument is omitted, the menu base is a top-level base.
;
; KEYWORDS:
;	BASE:	A named variable to recieve the widget ID of the created base.
;
;	COLUMN:	If set, the buttons will be arranged in a column.  If unset,
;		the buttons will be arranged in a row.
;
;	FRAME:	The width, in pixels of the frame drawn around the base.  The
;		default is no frame.
;
;	TITLE:	If PARENT is not supplied, TITLE can be set a string to be
;		used as the title for the widget base.
;
;	FONT:	A string that contains the name of the font to use for the
;		menu buttons.
;
; OUTPUTS:
;	None.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	A pulldown menu widget heirarchy is created, but not realized.
;	Each button has the label specified by the first field of the
;	corresponding pulldown menu description line.  Each button has a
;	user value (uvalue) specified by the second field.
;
; RESTRICTIONS:
;	Very little syntax checking is done on the description file.
;	Incorrectly formated input can lead to unexpected results.
;
; EXAMPLE:
;	For an example of using XPDMENU, see the "Pull-Down Menu" example
;	in the "Simple Widget Examples".  To create the simple widget examples
;	main menu, enter WEXMASTER from the IDL prompt.
;
; MODIFICATION HISTORY:
;	4 October 1990, AB, RSI.
;	16 January 1991, AB	Added the option of DESC being a string
;				array containing the description.
;-

function mkpull_getline, unit, data, idx, n, label, value
; Reads the next non-comment line from the description. If Unit is
; non-zero, it represents an open file from which the description
; is read. Otherwise, the description comes from data(idx) and idx is
; incremented. In this case, EOF is defined as idx being equal to n.
;
; LABEL is set to the label part and VALUE to the value part.
; On EOF, the file unit, if any, is closed. The return value of
; the function is:
;	-1 - End of file was seen
;	0 - Line was single button
;	1 - Line is a pulldown button
;

ret = -1
value = ''

if (unit eq 0) then not_eof = (idx lt n) else not_eof = (not eof(unit))
while ((not_eof) and (ret eq -1)) do begin
  if (unit eq 0) then begin
    value = data[idx]
    idx = idx + 1
  endif else begin
    readf, unit, value
  endelse
  value = strtrim(value, 2)		; Leading/trailing whitespace
  delim = strmid(value, 0, 1)
  case delim of
    "" :
    "#" :
    "}" : ret = 2
    else: begin
      value = strmid(value, 1, 100000)
      pos = strpos(value, delim)
      if (pos eq -1) then begin
	message, "Bad delimiter in line: " + delim + value, /INFORM
      endif else begin
        label = strmid(value, 0, pos)
        value = strtrim(strmid(value, pos+1, 100000), 2)
        if (value eq "{") then begin
	  value = ""
	  ret = 1
        endif else begin
          if (strlen(value) eq 0) then value = label
	  ret = 0
        endelse
      endelse
    end
  endcase
if (unit eq 0) then not_eof = (idx lt n) else not_eof = (not eof(unit))
endwhile

if ((unit ne 0) and (ret eq -1)) then begin free_lun, unit & unit = 0 & end

return, ret

end







pro mkpull_pulldown, parent, unit, data, idx, n, font
;
; unit - A file LUN or 0.
; data - If Unit is 0, data is a string array containing the menu
;	description.
; idx - If Unit is 0, idx is an integer giving the current index into data.
; n - If Unit is 0, n is an integer giving the # of elements in data.

while 1 do begin
  ret = mkpull_getline(unit, data, idx, n, label, value)

  case ret of
    -1 : return
    0 : begin
	if font ne '' then begin
		but = WIDGET_BUTTON(parent,value=label,uvalue=value,font=font)
	endif else begin
		but = WIDGET_BUTTON(parent, value=label, uvalue=value)
	endelse
	end
    1 : begin
	    if font ne '' then begin
		but = WIDGET_BUTTON(parent, value=label, MENU = 2, font=font)
	    endif else begin
		but = WIDGET_BUTTON(parent, value=label, MENU = 2)
	    endelse
	    mkpull_pulldown, but, unit, data, idx, n, font
	end
    2 : return
  endcase
endwhile

end


pro XPDMENU, DESC, PARENT, BASE=BASE, FRAME=FRAME, TITLE=TITLE,	$
	COLUMN=COLUMN, FONT=FONT

  s = size(parent)
  if (s[s[0] + 1] eq 0) then begin
    ; No parent is specified.
    parent = 0
    if (not keyword_set(TITLE)) then TITLE='Menu'
  endif else begin
    if (s[0] ne 0) then message, 'PARENT must be a scalar value.'
    if (s[1] ne 3) then message, 'PARENT must be a long integer.'
  endelse

  s = size(desc)
  if (s[s[0]+1] ne 7) then $
	message,'Description argument must be of type string."
  if (s[0] eq 0) then begin
    openr, unit, desc, /GET_LUN
    n = 0
 endif else begin
    if (s[0] ne 1) then message, 'String array must be 1-D."
    unit = 0
    n = s[1]
  endelse

  if (not keyword_set(frame)) then frame = 0
  if (not keyword_set(font)) then font = ''

  if (parent eq 0) then $
      IF(KEYWORD_SET(COLUMN)) THEN $
      	  base = WIDGET_BASE(/COLUMN, TITLE=TITLE, FRAME=FRAME) $
      ELSE $
	  base = WIDGET_BASE(/ROW, TITLE=TITLE, FRAME=FRAME) $
  else $
      IF(KEYWORD_SET(COLUMN)) THEN $
      	  base = WIDGET_BASE(parent, /COLUMN, FRAME=FRAME) $
      ELSE $
	  base = WIDGET_BASE(parent, /ROW, FRAME=FRAME)

  mkpull_pulldown, base, unit, desc, 0, n, font

end
