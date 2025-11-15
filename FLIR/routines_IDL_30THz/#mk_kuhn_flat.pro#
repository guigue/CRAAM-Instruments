;
;+
;
; Make a flat field image using a set of spatially displaced
; data images.
;

PRO mk_kuhn_flat, images, thresh, flat, offsets=offsets, niter=niter

;
; INPUT PARMETERS:
;	images = cube of images.  Temporal variations must be removed.
;	thresh = threshold for selecting valid parts of images.  If
;		thresh > 0, valid data > thresh.  If thresh < 0,
;		valid data < |thresh|.
;
; INPUT KEYWORDS:
;       offsets = vector of approximate positional offsets of the images,
;               giving x,y shifts for all images. Dimensioned (2,nimages).
;               Offset is defined to be the location of a common place in the
;               images, in CCD pixels.  Offsets should be given to 
;		fractional pixels, to get shifts between images correct.
;		If the IDL Astronomy correlation procedures are used
;		to measure the shifts between images, the offsets
;		to this routine should be the negatives of the shifts.
;
; OUTPUT PARAMETERS:
;	flat = flat field image, normalized.  Zero values => no data.
;
; OPTIONAL KEYWORDS:
;	niter = Number of iterations to run.  Default = 10.
;
; PROCEDURE:
;	Uses the method of Kuhn, Lin, and Loranz, Pub. Astron.
;	Soc. Pacific, 103:1097, 1991.
;
; HISTORY:
;	Written October 1993  Barry LaBonte
;
;-

asize = SIZE(images)
ncols = asize(1)
nrows = asize(2)
nimages = asize(3)

; Iterations
IF( KEYWORD_SET(niter) LE 0) THEN niter = 10

; Construct the K array and the number array
masks = REPLICATE( BYTE(0), ncols*nrows, nimages)
FOR i=0,nimages-1 DO BEGIN
        IF( thresh GT 0) THEN BEGIN
                valid = WHERE( images(*,*,i) GT thresh )
                ENDIF ELSE BEGIN
                valid = WHERE( images(*,*,i) LT (-thresh) )
        ENDELSE
	masks(valid,i) = BYTE(1)
ENDFOR
masks = REFORM(masks, ncols, nrows, nimages)

karray = FLTARR(ncols, nrows)
numarray = INTARR(ncols, nrows)
PRINT, ' Building basic image.'
FOR i=0,nimages-2 DO BEGIN
	tmp = images(*,*,i)
	valid = WHERE( masks(*,*,i) GT 0 )
	logi = REPLICATE(0.,ncols,nrows)
        logi(valid) = ALOG(tmp(valid))
	
	FOR j=i+1, nimages-1 DO BEGIN
		dx = offsets(0,i) - offsets(0,j)
		IF( dx NE 0. ) THEN dx = FIX(dx/ABS(dx)) * FIX(ABS(dx)+0.5)
		dy = offsets(1,i) - offsets(1,j)
                IF( dy NE 0. ) THEN dy = FIX(dy/ABS(dy)) * FIX(ABS(dy)+0.5)
		tmp = images(*,*,j)
		valid = WHERE( masks(*,*,j) GT 0 )
		logj = REPLICATE(0.,ncols,nrows)
                logj(valid) = ALOG(tmp(valid))

; Identify pixel combinations with valid data in both images.
; Allow for pixels that are thresholded out.
; First, the unshifted I image and the shifted J image.
		sum  = logi( (dx>0):ncols-1+(0<dx), (dy>0):nrows-1+(0<dy) )
		mask = masks( (dx>0):ncols-1+(0<dx), (dy>0):nrows-1+(0<dy), i )
		sum  = sum - logj( (-dx>0):ncols-1-(dx>0), (-dy>0):nrows-1-(dy>0) )
		mask = mask < masks( (-dx>0):ncols-1-(dx>0), (-dy>0):nrows-1-(dy>0), j )
		sum = sum * mask
		karray( (dx>0):ncols-1+(0<dx), (dy>0):nrows-1+(0<dy) ) = karray( (dx>0):ncols-1+(0<dx), (dy>0):nrows-1+(0<dy) ) + sum
		numarray( (dx>0):ncols-1+(0<dx), (dy>0):nrows-1+(0<dy) ) = numarray( (dx>0):ncols-1+(0<dx), (dy>0):nrows-1+(0<dy) ) + mask

; Now the unshifted J image and the shifted I image.
		dx = -dx
		dy = -dy
                sum  = logj( (dx>0):ncols-1+(0<dx), (dy>0):nrows-1+(0<dy) )
                mask = masks( (dx>0):ncols-1+(0<dx), (dy>0):nrows-1+(0<dy), j )
                sum  = sum - logi( (-dx>0):ncols-1-(dx>0), (-dy>0):nrows-1-(dy>0) )
                mask = mask < masks( (-dx>0):ncols-1-(dx>0), (-dy>0):nrows-1-(dy>0), i )
                sum = sum * mask
		karray( (dx>0):ncols-1+(0<dx), (dy>0):nrows-1+(0<dy) ) = karray( (dx>0):ncols-1+(0<dx), (dy>0):nrows-1+(0<dy) ) + sum
		numarray( (dx>0):ncols-1+(0<dx), (dy>0):nrows-1+(0<dy) ) = numarray( (dx>0):ncols-1+(0<dx), (dy>0):nrows-1+(0<dy) ) + mask

	ENDFOR
ENDFOR
logi = 0
logj = 0
tmp = 0


valid = WHERE(numarray GT 0)
karray(valid) = karray(valid)/numarray(valid)
TVSCL,karray

; Now iterate to a flat 
change = 1.

; Smooth to get a start on the large scale structure
flat = karray
PRINT, 'Iteration, Change'
FOR iter=1,niter DO BEGIN
correction = REPLICATE(0., ncols, nrows)
; Sum up the corrections
	FOR i=0,nimages-2 DO BEGIN
		FOR j=i+1,nimages-1 DO BEGIN
                	dx = offsets(0,i) - offsets(0,j)
			IF( dx NE 0. ) THEN dx = FIX(dx/ABS(dx)) * FIX(ABS(dx)+0.5)
                	dy = offsets(1,i) - offsets(1,j)
			IF( dy NE 0. ) THEN dy = FIX(dy/ABS(dy)) * FIX(ABS(dy)+0.5)
			mask = masks( (dx>0):ncols-1+(0<dx), (dy>0):nrows-1+(0<dy), i ) < masks( (-dx>0):ncols-1-(dx>0), (-dy>0):nrows-1-(dy>0), j )
			sum = flat( (-dx>0):ncols-1-(dx>0), (-dy>0):nrows-1-(dy>0) )
			bad = WHERE( mask LE 0, bcount )
			IF( bcount GT 0 ) THEN sum( bad ) = 0.
			correction( (dx>0):ncols-1+(0<dx), (dy>0):nrows-1+(0<dy) ) = correction( (dx>0):ncols-1+(0<dx), (dy>0):nrows-1+(0<dy) ) + sum
	                dx = -dx
                	dy = -dy
			mask = masks( (dx>0):ncols-1+(0<dx), (dy>0):nrows-1+(0<dy), j ) < masks( (-dx>0):ncols-1-(dx>0), (-dy>0):nrows-1-(dy>0), i )
			sum = flat( (-dx>0):ncols-1-(dx>0), (-dy>0):nrows-1-(dy>0) )
			bad = WHERE( mask LE 0, bcount )
                        IF( bcount GT 0 ) THEN sum( bad ) = 0.
                	correction( (dx>0):ncols-1+(0<dx), (dy>0):nrows-1+(0<dy) ) = correction( (dx>0):ncols-1+(0<dx), (dy>0):nrows-1+(0<dy) ) + sum

		ENDFOR
	ENDFOR

	newflat = karray
	newflat(valid) = newflat(valid) + correction(valid) / numarray(valid)

	change = STDEV(newflat - flat)
	flat = newflat
	TVSCL,flat

	PRINT, iter, change
ENDFOR


flat(valid) = EXP(flat(valid))

END
