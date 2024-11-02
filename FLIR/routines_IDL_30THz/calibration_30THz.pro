
	pro calibration_30THz, files, index, flatdir=flatdir, outdir=outdir, flat_on=flat_on, center_on=center_on, rotate_on=rotate_on, $
	limbdark_on=limbdark_on  

	;common share1, level, year, month, day, hour, minute, second

	;-------------------------------------------------- Information ----------------------------------------------------------------------
	; name: "calibration_30THz.pro"
	; call: calibration_30THz, files, index
	;
	; Objective: This routine apply the calibration techniques developed for the Mid-IR 30 THz camera (10 um) 
	; The routine also update the header information. 
	;
	; ########## IMPORTANT: "The input images must be the output fits images of the routine "transform_fpf_to_fits_30THz.pro" ##################
	;
	; Procedures in order of realization:
	; 1) Normalization by Flat Field (routine: "correcting_flat_field_30THz.pro")
	; 2) Center the Sun's center with the center of the images (necessary for the limb darkening correction, routine: "center_30THz.pro")
	; 3) Rotation of the images to get Solar North UP -- (routine: "rotate_30THz.pro")
	;	4) Optional: Limb darkening correction (routine "correcting_limb_darkening_30THz.pro") -- IMPORTANT: LD REQUIRES CENTER THE IMAGES!!!!
	; 
	; Input:   
	; Data obtained by the 30 THz camera from Sao Paulo (telescope: BR30T) and El Leoncito (telescope: AR30T) are admitted    
	; The routine requires:
	; - files = the list of files to be calibrated. ".fits" images are accepted. 
	; - index = a vector containing the indexes of the images to calibrate (-1 results in the calibration of all images) 
	;
	; Keywords: 
	; flat_dir = 'path' = the path where the Flat Field image is located (not required if keyword /flat_on is not used)
	; outdir='path'     = The path to save the calibrated images (if not indicated, it is created in same folder than the input images)
	; /flat_on          = for applying the FLat field correction
	; /center_on        = for align the center of the Sun with the center of the image 
	; /rotate_on        = for rotating the images to solar north up (apply the camera angle and P0 corrections)
	; /limbdark_on      = for applying the limb darkening correction
	; 
	; Levels adopted for the 30 THz Data:
	; Level 0.0: Raw data in .fpf format as obtained by the 30 THz Camera (with no modifications) 
	; Level 0.2: Raw data (in fpf or fits format) with any correction applied (but not including the pre-processing techniques) 
	; Level 0.5: Raw data + Flat fielding correction + center of the Sun's center
	; Level 1.0: Raw data + Flat fielding + center of the Sun in the center of the image + rotations (Solar North up)
	; Level 1.5: Raw data + Flat fielding + centered images + rotations (Solar North up) + limb darkening corrections
	; Level 2.0: Level 1.5 images + Post-processing techniques (e.g. correction by jitter + elimination of noise) 
	; 
	; Output: 
	; The resulting filename is indicated as: 
	; "AR30T_YYYY-MM-DDTHH:MM:SS.SSS_levelxx.fits" ; "BR30T_YYYY-MM-DDTHH:MM:SS.SSS_levelxx.fits"  
	; --- AR30T: for images obtained at El Leoncito and BR30T for images obtained at CRAAM-Sau Paulo
	;
	; History:
	; written by Fernando M. Lopez (CRAAM-Mackenzie) --- October 2019
	; The routine is based in the work done by Carlos Francile and Franco Manini (OAFA-UNSJ, San Juan, Argentina)
	;
	; ---------------------------------------------------------------------------------------------------------------------------
	; Start estimating the time of the process
	t0 = systime(1)

	; Case of insufficient inputs:
	if n_params() lt 2 then begin
		err_mess = ' Input error: Include the list of files and indexes to calibrate.  Returning.'
		print, ' Input error: Include the list of files and indexes to calibrate.  Returning.'
		return
	endif

	; input mode: list 
	if size(files, /tname) eq 'STRING' then input_mode = 'file_list'

	ss_infil = index
		if ss_infil[0] eq -1 then ss_infil = indgen(n_elements(files))

	nimages=n_elements(ss_infil)

	; -------------------------------------------------------------------------------------------------------------------
	; check that all images are of the same format (.fits) and same level of processing

	dummy=intarr(nimages)
	dummy_level=intarr(nimages)
	for j=0l, nimages -1 do begin
		tipo=strmid(files[j],3,/reverse)
		level0=strmid(files[0],5,1,/reverse)+'.'+strmid(files[0],6,1,/reverse)

		dummy_level[j]=level0
		if tipo eq 'fits' or tipo eq '.fts' then dummy[j]=1 else if tipo eq '.fpf' then $
			message, '--- Error: Only fits images are allowed. Stopping Calibration --- '
	endfor

	re1=where(dummy eq 0, countfpf,ncomplement=countfits)
	suma_dummy=total(dummy_level)
	if suma_dummy ne nimages*level0 then message, '****** Error: The images have different processing levels. Stopping Calibration *******'

	if countfpf gt 0 and countfits gt 0 then message, 'Error: Different formats (fits and fpf) were found. Stopping Calibration'

	print, '****************************************************************'
	print, ' Preparing calibration for ', nimages, ' images of format:', tipo
	print, ' The images have a calibration level', level0
	print, '****************************************************************'

	level=level0

	; -----------------------------------------------------------------------------------------------------------------------------
	; Search for the Flat field image
	; extract the date of the observation 

		if tipo eq 'fits' or tipo eq '.fts' then begin
				Image_dummy0=readfits(files[0],header_dummy0)
				DATE_OBS=sxpar(header_dummy0, 'DATE_OBS')		
				date_flat = strmid(DATE_OBS, 0,4) +'-'+ strmid(DATE_OBS, 5,2) +'-'+ strmid(DATE_OBS, 8,2)
		endif
		
	if ~keyword_set(flat_on) then print, '****** Flat field correction will not be applied *******'

	; From here, this is only done if correction by Flat field will be applied
	if keyword_set(flat_on) then begin
		; read the Flat field image
		if ~keyword_set(flatdir) then begin 
					flatdir=''
					read, flatdir,prompt=' ***** Enter the path where the Flat field is located, [Enter] ***** : ' 
							
	; must search for a Flat of the respective date	
		endif					
	
		file_flat=file_search(flatdir+'*'+date_flat+'*FLAT.fits',count=nflats)
		if nflats eq 0 then begin
				print, '****************************************************************'
				print, 'NOT FLAT IMAGES FOUND FOR THE DATE OF THE OBSERVATIONS'
				print, '****************************************************************'
			  nameflat=''
			  READ, nameflat, PROMPT=' ***** Enter the Flat field name (no extension), [Enter] ***** : '
				file_flat=file_search(flatdir+nameflat+'*.fits',count=nflats)
		endif
		if nflats gt 1 then begin
				print, '***************************************************************'				
				print, 'MULTIPLES FLAT IMAGES FOUND FOR THE DATE OF THE OBSERVATIONS'
				print, '***************************************************************'
			  nameflat=''
			  READ, nameflat, PROMPT=' ***** Enter the Flat field name (no extension), [Enter] ***** : '
				file_flat=file_search(flatdir+nameflat+'*FLAT'+tipo,count=nflats)
		endif

		Iflat=readfits(file_flat,headerflat)
		
	endif

	; if outdir is not supplied then it creates a folder
	if ~keyword_set(outdir) then begin
		pos_barra=strpos(files[0],'/',/reverse_search)
		outdir=strmid(files[0],0,pos_barra+1) + 'calibrated/'
		file_mkdir, outdir
		
	endif
		print, '**********************************************************************************'
		print, 'SAVING CALIBRATED IMAGES IN:'
		print, outdir
		print, '**********************************************************************************'
	

	window,4,xs=512,ys=512, title='Image centered'
	window,1, xs=512,ys=512, title='Final Image'	
	loadct,8

	;*****************************************************************************************************************************
	; start calibration process for all images 
	for i=0l, nimages-1 do begin

	; read the files
		Imgin=readfits(files[i],header) 
		telescope=strmid(sxpar(header, 'TELESCOP'),0,5)
		; obtain the date and time for ephemerides!!!!
					date_obs=sxpar(header, 'DATE_OBS')
					year=strmid(date_obs,0,4)
					month=strmid(date_obs,5,2)
					day=strmid(date_obs,8,2)
					hour=strmid(date_obs,11,2)
					minute=strmid(date_obs,14,2)
					second=strmid(date_obs,17,6)
					 
			Imgoriginal=Imgin
;			Imgout=Imgoriginal
			corrections_list = '.'

	;*****************************************************************************************************************************
	; apply the correction by Flat Field: 
		
			If keyword_set(flat_on) then begin
					correcting_flat_field_30THz, Imgin, header, Iflat, headerflat, Imgout 
					Imgin=Imgout
					Img05=Imgout
					level='05'
;					corrections_list = corrections_list+' Flat Fielding Correction Applied'
			endif
	; Now Imgout is the image corrected by Flat field
	;*****************************************************************************************************************************
	; locate the center of the Sun and put it in the center of the image
			if keyword_set(center_on) then begin
					center_30THz, Imgin, header, Imgout
					Imgin=Imgout
					Img05=Imgout
					level='05'
			endif
	; Now Img05 are the images corrected by Flat Field and centered. The dimension of the new images level 0.5 are 800x800
	;*****************************************************************************************************************************
	; Apply the rotation to the images: 

		If keyword_set(rotate_on) then begin
				Imgin = Imgout
				Img10 = Imgout
				rotating_30THz, Imgin, header, Imgout
				level='10'				
		endif
	; Now Img10 is the image corrected by Flat field, centered and rotated 
	; The dimension of level 1.0 images is 800x800 pixels^2
	;*****************************************************************************************************************************
	; Apply the Limb darkening correction: 
		If keyword_set(limbdark_on) then begin
				correcting_limb_darkening_30THz, Imgin, header, Imgout
				Imgin = Imgout
				Img15 = Imgout
				level='15'
		endif
	;*****************************************************************************************************************************
	; save the images
		wset, 1 & tv, bytscl (congrid(Imgout,512,512), min=max(Imgout)-stdev(Imgout)/4,max=max(Imgout))
	;	sxaddpar,header, 'HISTORY', corrections_list
		
		filename = telescope+'_'+year+'-'+month+'-'+day+'T'+hour+':'+minute+':'+second+'_level'+level+'.fits'
		sxaddpar, header, 'FILENAME', filename		

		writefits, outdir+filename, Imgout, header
		
	endfor  

	t1 = systime(1)
	deltatime=(t1-t0)/60.
	print, '*****************************************************************************'
	print, 'Calibration procedure finished in: ', strtrim(string(deltatime),1), ' minutes'
	print, ' All calibrated images saved in ', outdir
	print, '*****************************************************************************'
	end


	;###############################################################################################################################################

	pro correcting_flat_field_30THz, Imgin, header, Iflat, headerflat, Imgout 
	;common share1, level, year, month, day, hour, minute, second

	; ----------------------- Explanation --------------------------------------------------------------------
	; This process performes the Flat fielding correction (uniform of the gain of the detector)  
	; It is very important to indicate the correct flat field image
	; The input is a level 0 raw image, with no pre-processing techniques previously applied
	; The output is a level 0.5 image (level 0.5 = level 0 + Flat field correction)
	;
	; Input: 
	; Iflat = the flat field image
	; Imgin = the raw image level 0.0 to be calibrated
	; header = the header of the image to be calibrated
  ; hederflat = the header of the flat field image
	; 
	; Output:
	; Imgout = image corrected by flat field (level 0.5) 
	;
	; written by Fernando Lopez -- Octobre - 2019  
	; --------------------------------------------------------------------------------------------------------
	 
		  ;Elinate cero values in the flat image
		  ceros=where(Iflat eq 0., count)
		  if count ne 0 then Iflat(ceros)=1.0
	
	telescope_img=sxpar(header,'TELESCOP')
	telescope_flat=sxpar(headerflat,'TELESCOP')
			
	if sxpar(header,'TELESCOP') ne sxpar(headerflat,'TELESCOP') then $
				message, 'THE FLAT CORRESPONDS TO A DIFFERENT TELESCOPE. STOPPING CALIBRATION'

	; size of the flat image:
	sflat=size(Iflat)
	simage=size(Imgin)

	; gives error message if the dimenson of the flat and of the image are different
	if sflat[1] ne simage[1] then message, 'The Flat and the Image have different dimensions',returning
	if sflat[2] ne simage[2] then message, 'The Flat and the Image have different dimensions',returning

	 ;--------------------
		 ;Apply Flat correction
		 Imgout=Imgin/Iflat

	; img05 is the final image corrected by Flat field (level 0.5)
	
	sxaddpar,header,'LVL_NUM', '0.5'
	sxaddpar,header, 'HISTORY', 'Flat Fielding Correction Applied '
	
	return
	end


	;###############################################################################################################################################
	pro center_30THz, Imgin, header, Imgout

	; ------------------------------------------- Explanation ------------------------------------------------------ 
	; Name: "center_30THz.pro"
	; Project: Part of the package of IDL routines destinated to calibrate the 30 THz (10um) IR telescope
	;
	; Objective: The routine determines the pixels location of the Sun's center and move it to the center of the image. 
	; 					 Important to correct for movements during the tracking of the telescope. 
	; 
	; Input: 
	; Imgin:  The images level 0.5, with Flat field previously corrected. 
	; header: the header of the image
	;
	; Output: 
	; Imgout: The level 0.5 image corrected by flat field and centered (imgages are rebinned to 800x800 pix^2, but the Sun is not rebinned.
	; 				Then the pixel scale (cdelt) does not change)
	; header: The modified header of the centered image (with new pixels of the center of the Sun)
	;
	; History:
	; written by Fernando M. Lopez (CRAAM-Mackenzie) --- October 2019
	; The routine is based in the work done by Carlos Francile and Franco Manini (OAFA-UNSJ, San Juan, Argentina)
	;------------------------------------------------------------------------------------------------------------------------------------- 
	; 

	; using get_sun this below is the correct format: 
	; result=get_sun('2019-05-15 19:00')  ; epheme2 is equal to get_sun

	; The final dimension of the images centered. This values can be changed if desired. 
	dimxf=800
	dimyf=800

	; obtain the date and time for ephemerides!!!!
	date_obs=sxpar(header, 'DATE_OBS')
	year=strmid(date_obs,0,4)
	month=strmid(date_obs,5,2)
	day=strmid(date_obs,8,2)
	hour=strmid(date_obs,11,2)
	minute=strmid(date_obs,14,2)
	second=strmid(date_obs,17,2)
	 
	SUN_DIST=sxpar(header,'SUN_D')

	; dimension of the images
	ss=size(Imgin)
	dimx=ss[1] & dimy=ss[2]

	; The center of the image 800x800 is located at:
	Imgcenter=[(dimx/2.),(dimy/2.)]

	; the image to be centered
	Idummy = Imgin
	Idummy = Idummy - min(Idummy)
	Idummy = Idummy/max(Idummy)


	Imgout = fltarr(dimxf,dimyf)+min(Imgin)

	Imgout[((dimxf-dimx)/2.): dimxf-((dimxf-dimx)/2.)-1, ((dimyf-dimy)/2.):dimyf-((dimyf-dimy)/2.)-1] = Imgin
	;if dimxf gt dimy then Imgout[((dimxf-dimx)/2.): dimxf-((dimxf-dimx)/2.)-1, ((dimyf-dimy)*.7):dimyf-((dimyf-dimy)*.3)-1] = Imgin


	; define the limb as I/I0=0.5
	intcenter = mean(Idummy[dimx/2.-20:dimx/2+20,dimy/2.-20:dimy/2.+20])
	resultcenter=WHERE(Idummy lt intcenter/2., countador2, COMPLEMENT=notresultcenter)
		  if countador2 eq 0 then print, '----- FAILED TO FIND OUT THE CENTER ----'

		  ;--------------------------------------------------------------------------------------------------------
		  ;Generate in Icalcu and image where the limb is determined by the method of Roberts for edge detection 
		  Icalcu=FLTARR(dimx,dimy)+1.
		  Icalcu(resultcenter)=0.
		  ;Icalcu(notresultcenter)=1.
		  Icalcu=ROBERTS(Icalcu)

		  resultcentered=where(Icalcu eq 1)
		  resultcentered3 = ARRAY_INDICES(Icalcu, resultcentered)

		  ;-------------------------------------------------------------------------------------------------
		  ;interpolate a circle to the position of the Limb of Icalcu, using the package mpfit (L-Markwardt)
		  start_parms=[367.,367.,Imgcenter[0]-1.,Imgcenter[1]-1,0.]
		  parms = MPFITELLIPSE(FLOAT(REFORM(resultcentered3(0,*))), FLOAT(REFORM(resultcentered3(1,*))), start_parms,/CIRCULAR, QUIET=1)
		  ;      P[0]   Ellipse semi axis 1
		  ;      P[1]   Ellipse semi axis 2   ( = P[0] if CIRCLE keyword set)
		  ;      P[2]   Ellipse center - x value
		  ;      P[3]   Ellipse center - y value
		  ;      P[4]   Ellipse rotation angle (radians) if TILT keyword set

		  ;---------------------------------------------------------------------------------------------------
		  ;check for validity of the value of the Sun's center obtained. 
		  if parms(2) gt (dimx/2.)-0.5*(dimx/2.) and parms(2) lt (dimx/2.)+0.5*(dimx/2.) and $
			parms(3) gt (dimy/2.)-0.5*(dimy/2.) and parms(3) lt (dimy/2.)+0.5*(dimy/2.) then begin
			;if parms(2) gt 0 and parms(2) lt dimx and parms(3) gt 0 and parms(3) lt dimy then begin
					
		      x0=parms(2)        ; pixel x of the solar center
		      y0=parms(3)        ; pixel y of the solar center 
		      r0=round(parms(1)) ; solar radius in pixels
		  endif else message, 'Could not find correctly the center of the Sun.. check for the image',/continue
	
			if n_elements(x0) eq 0 then x0=dimxf/2 
			if n_elements(y0) eq 0 then y0=dimyf/2  
			if n_elements(r0) eq 0 then r0=1			
			print, '******** Center and Radius of the Sun Found *********', x0,y0,r0
			
			; Sun's center found:
			Suncenter=[x0,y0]

					;CDELT1, CDELT2 are the resolution of the pixels in x and y directions
					;Solar radius has a value 6.957 10^8 m  (de Stix) o  r = 6,96 10^8 m = 0,00465247 AU
					;AU = 149.597.870 +-2 KM   (Stix),   149.597.870.700 m (IAU)
					;r0 is the solar radius in pixels
					;SUNDIST is the distance Sun-Earth for the date_obs
					;Solar radius in arcseconds = arctg (solar radius / distance Sun-Earth)
					angle= ATAN(6.957D8 / (SUN_DIST*149597870700.D)) * (180.D/!DPI) * 60.D * 60.D ;en segundos de arco
					resolution=angle/r0


	; solar center coordinates in the expanded image
	newx0 = x0 + ((dimxf-dimx)/2.)
	newy0 = y0 + ((dimyf-dimy)/2.) 

	; new desplacelent 
	deltax_new =  (dimxf/2.) - newx0
	deltay_new = (dimyf/2.) - newy0

	; vector to desplacement of the Sun deltar = vector center image  -  vector center sun
	; Now have to move all the Sun to coincidite with the center of the image

	deltar = Imgcenter - [x0,y0]        ; vector with distance in x,y between center of the image and Sun's center
	deltar_mod = sqrt((deltar[0])^2. + (deltar[1])^2.)  ; module of the distance en pixels

	deltax = (Imgcenter[0] - x0) 
	deltay = (Imgcenter[1] - y0)

	;if abs(deltax) gt (dimxf-dimx) or abs(deltay) gt (dimyf-dimy) then $
	;		message, 'Message: Large displacement of the Sun: Require increase the dimension of the images'

	;Imgout = shift(Imgout,deltax_new,deltay_new)
	Imgout = shift(Imgout,deltax,deltay)

	Img45=Imgout

	;imagen3=shift(imagen3,640-(round(x)+320),480-(round(y)+240))

	Img45[dimxf/2-1:dimxf/2+1,*]=10000
	Img45[*,dimyf/2-1:dimyf/2+1]=10000

	Img45[newx0-3:newx0+3,newy0-3:newy0+3]=0
   wset, 4 & tv, bytscl (congrid(Img45,512,512), min=max(Imgout)-stdev(Imgout)/4,max=max(Imgout))

	; #######################################################################
	; update the header with the new modifications to the image
	
		sxaddpar,header,'CRPIX1', (dimxf+1)/2
		sxaddpar,header,'CRPIX2', (dimyf+1)/2
		sxaddpar,header,'RSUN_LF', r0, ' Limb fit solar radius in pixels',before='SUN_D',format='f10.5'
		sxaddpar,header,'CDELT1', resolution, after='CRVAL2',format='f8.5'
		sxaddpar,header,'CDELT2', resolution, after='CDELT1',format='f8.5'
		sxaddpar,header,'CTYPE1', 'arcsec' , after='CDELT2'
		sxaddpar,header,'CTYPE2', 'arcsec' , after='CTYPE1'			
		sxaddpar,header,'LVL_NUM', '0.5'
		sxaddpar,header, 'HISTORY', 'Centered'
		if r0 eq 1 then sxaddpar,header, 'HISTORY', 'NOT POSSIBLE TO DETERMINE LIMB POSITION IN THE IMAGE'
	; #######################################################################


	return
	end

	;################################################################################################################################################
pro rotating_30THz, Imgin, header, Imgout
	; ------------------------------------------- Explanation ------------------------------------------------------ 
	; Name: "rotating_30THz.pro"
	; Project: Part of the package of IDL routines destinated to calibrate the 30 THz (10um) IR telescope
	;
	; Objective: The routine rotates the images to obtain the Solar North UP. Necessary to compare with other observations
	; Corrections applied: 
	;  - Angle desviation of the vertical axis camera with the perpendicular direction to the Ecliptic (derived by Manini et al. 2017, BAAA) 
	;  - Correction of the P0 angle (put the solar axis perpendicular to Ecliptic - The value of P0 is obtained by ephemerides)
  ;  - reverse to obtain correct position of E-W 
	; 
	; Input: 
	; Imgin:  The images level 0.5, with Flat field correction and Sun centered. 
	; header: the header of the image
	;
	; Output: 
	; Imgout: The level 1.0 image rotated 
	;					(imgages are rebinned to 800x800 pix^2, but the Sun is not rebinned. Then the pixel scale (cdelt) does not change)
	; header: The modified header of the rotated image
	;
	; History:
	; written by Fernando M. Lopez (CRAAM-Mackenzie) --- October 2019
	; The routine is based in the work done by Carlos Francile and Franco Manini (OAFA-UNSJ, San Juan, Argentina)
	;------------------------------------------------------------------------------------------------------------------------------------- 

		; correct the E-W positions 		
		Imgout = reverse(Imgin,1)
		x=sxpar(header,'X0 ')
    y=sxpar(header,'Y0 ')
    P0=sxpar(header,'P0 ')
		desviation = sxpar(header,'CROTA')
		nullvalues= min(Imgin)
		x0=sxpar(header,'CRPIX1')
		y0=sxpar(header,'CRPIX2')
		Imgout = rot(Imgout, (P0+desviation), 1.0, x0, y0, CUBIC=-0.5, /PIVOT, MISSING=nullvalues)
	
		sxaddpar,header,'LVL_NUM', '1.0'
		sxaddpar,header, 'HISTORY', 'Image reversed and rotated (clockwise): '+ strtrim(string(P0+desviation),1) +' deg -- Solar North Up'
		sxaddpar, header, 'CROTA', strtrim(string((P0+desviation)),1), ' Rot angle (deg), clockwise, from Y-axis.'
return
end


	;################################################################################################################################################
	pro correcting_limb_darkening_30THz, Imgin, header, Imgout
	;common share1, level, year, month, day, hour, minute, second

	;----------------------------- Explanation ------------------------------------------------------------------------
	; This program apply correction of limb darkening to the 30THz images.
	; The input is a raw (fpf or fits) image, with Flat field correction previously applied
	; The output image is the image with the correction of Flat field and limb darkening!!!!
	; The required step is to perform the normalization by Flat field (uniform the gain of the detector)
	;
	; Input: 
	; 	- Imgin = the image previously corrected by flat field (level0.5)
	; 	- centro = [xc,yc,r0]; array containing: the pixels of the Sun's center (xc,yc) and solar radius (r0)
	; Output:
	;		- Imgout = the image corrected by limb darkening (level 1.0 if corrected previously by Flat Field)
	;-------------------------------------------------------------------------------------------------------------------

	centro=fltarr(3)
	; p0 the parameters obtained in Manini et al. (2018) - BAAA
	p0=[0.188731, -0.0626018] 
	centro=[sxpar(header,'CRPIX1'), sxpar(header,'CRPIX2'), sxpar(header,'RSUN_LF')]	

	; calculate distances to solar centre

	isize=size(Imgin)
	xi=isize(1)
	yi=isize(2)

	mascara=fltarr(xi,yi)+1.

		  for i=0, xi-1 do begin
		      for j=0, yi-1 do begin
		         dist=SQRT((i-centro(0))^2. + (j-centro(1))^2.)/centro(2)
		         if dist lt 1. then Mascara(i,j)=(1. - p0(0) - p0(1) + p0(0)*cos(asin(dist)) + p0(1)*cos(asin(dist))^2.)
		     endfor
		  endfor

	;correct imput image for limb darkening effects
	;lo que hace es restar el minimo, corregir, y luego aplicar de nuevo el minimo.
	;el programa manda temperatura en base a lo que calcula, entonces hay un gap entre la base de la curva
	;y el cero del eje, donde la curva de oscurecimiento se aplica, entonces debo bajar la curva
	;aplicar la funcion de osc y luego volver a subir

	minimo=min(Imgin)
	Imgin=Imgin-minimo
	Imgout=Imgin/mascara
	Imgout=Imgout+minimo
	sxaddpar, header,'LVL_NUM', '1.5'
	sxaddpar,header, 'HISTORY', 'Limb Darkening Correction Applied'
	; Img10 is the final image corrected by Flat field (level 0.5) and Limb darkening (level 1.0)
	;sxaddpar, header, 'HISTORY', 'Limb darkening corrected'
	return
	end



