pro transform_fpf_to_fits_30THz, telescope=telescope, indir=indir, level=level


;*****************************************Explanation *********************************************************************************
; call:  transform_fpf_to_fits_30THz, telescope='telescope', indir='path'
;
; This routine reads the 30 THz fpf images and then save them in agreement of their Date_obs field value. 
; The routine also creates a header for the fits images
; The routine reads the raw .fpf data and save the new images in fits format 
; Creates a new forder named "fits" located in the imput folder. 
; The name for the new images are given by: "AR30T_YYYY_MM_DDTHH:MM:SS.SSS_levelXX.fits"
; the output images of this routine are the input images for the calibration routines ("calibration_30THz.pro")
; The routine also creates a header for the fits images
;
; Contains two function to create the headers. One for the AR30T and BR30T
; Keywords: 
; indir = 'the path where the raw fpf images are located'
; telescope='AR30T' or 'BR30T'
; level= 'levelXX' ; the level of the images(eg. level00 = raw) - optional: if not indicated it is assumed as raw 
;
; 
; Output:
; header= the header of the .fits image
;
; History: Written by Fernando M. Lopez (CRAAM-Mackenzie) --- October 2019
; The routine is based in the work done by Carlos Francile and Franco Manini (OAFA-UNSJ, San Juan, Argentina)
;**************************************************************************************************************************************

	; If not indicated must enter the path where the original images are located
	if not keyword_set(indir) then begin 
			indir=''
			READ, indir, PROMPT='Indicate the path of the input images [Enter]: '
	endif
	if not keyword_set(telescope) then begin
			telescope=''
			READ, telescope, PROMPT='Indicate the telescope (AR30T or BR30T) [Enter]: '
	endif
	
	; check that indir and telescope keywords are string
	if size(indir,/type) ne 7 then indir=string(indir)
	if size(telescope,/type) ne 7 then telescope=strtrim(string(telescope),1)

	if not keyword_set(level) then level='level00'

	if level eq 'level00' then print, ' ----- The images are RAW level 0.0 ------ ' else print, ' ----- The images are PROCESSED  -------'

	; check last character of indir be a /
	total_indir=strlen(indir)
	pos_bar=strpos(indir,'/',/reverse_search)
	if strmid(indir,pos_bar,1) ne '/' then indir=indir+'/' 
		
	
	; create the output folder
	file_mkdir, indir + 'fits/'

	outdir = indir+'fits/'

	; read the images
	files = file_search(indir+'*.fpf',count=n)

	if n eq 0 then begin
		print, ' Not .fpf images found, looking for .fits images'
		files = file_search(indir+'*.fits',count=n)
		if n eq 0 then message, 'Error : not images found--check for the path'
	endif

	; the format of the images (should be .fpf)
	formato = strmid(STRTRIM(string(files(0)),1),3,/reverse)

	if formato eq '.fpf' then begin

		for i=0,n-1 do begin
			
			fpf = read_fpf(files[i]) 
			image = fpf.data
			pos_barra=strpos(files[i],'/',/reverse_search)
			total_strings=strlen(files[i])
		
			original_name = strmid(files[i],pos_barra+1,total_strings-pos_barra)

		
		; ----------------------------- create the header for the fits file to save -------------------------------------------------------
		; make a header for the telescope: AR30T or BR30T 	
			if telescope eq 'AR30T' then header = makeheader_AR30T(fpf, level)
			if telescope eq 'BR30T' then header = makeheader_BR30T(fpf, level)
					
						sxaddpar, header, 'ORIGIFIL', original_name,' Original fpf filename', before='DATE_OBS'
	; -----------------------------------------------------------------------------------------------------------------------------------
			filename=sxpar(header,'FILENAME')
			file = outdir+filename+'.fits'
			print, 'Processing image ', i , ' from a total of: ', n , ' for telescope: ', telescope 	
			writefits, file, image, header

		endfor
	endif else message, '-------- Error:Input images are not .fpf --------' 
	print, '-------------------------------- Process finishded --------------------------------'
	end

; ############################################################################################################################################
function makeheader_AR30T, fpf, level

   
    ;This function generates the header for the images
		image=fpf.data
		SS=size(fpf.data)
		datatype=4
		NAXISX = SS[1]
		NAXISY=SS[2]
		naxisarr=[NAXISX,NAXISY] 
		level2 = strmid(level,5,1)+'.'+strmid(level,6,1)
			mkhdr, header, datatype, naxisarr
    
		year = STRTRIM(string(fpf.DATE_YEAR),1)
    month = STRTRIM(string(fpf.DATE_MONTH),1) 
    day = STRTRIM(string(fpf.DATE_DAY),1)
    hour = STRTRIM(string(fpf.DATE_HOUR),1)
    minute = STRTRIM(string(fpf.DATE_MINUTE),1)
    second = STRTRIM(string(fpf.DATE_SECOND),1) + '.' + STRTRIM(string(fpf.DATE_MILLISECOND),1)

  			if strlen(month) eq 1 then month='0'+month
        if strlen(day) eq 1 then day='0'+ day
        if strlen(hour) eq 1 then hour = '0'+ hour
        if strlen(minute) eq 1 then minute = '0' + minute
        if strlen(STRTRIM(string(fpf.DATE_SECOND),1)) eq 1 then second = '0' + second


     ;generate variables dateobs and timeobs
    dateobs=year+'/'+month+'/'+day
    date_obs=year+'-'+month+'-'+day+'T'+hour+':'+minute+':'+second+'Z'
    time_obs=hour+':'+minute+':'+second

   	filename = 'AR30T_'+year+'-'+month+'-'+day+'T'+hour+':'+minute+':'+second+'_'+level	

;**************************************************************************************************************************************
		; add the keywords to the header
    sxaddpar, header, 'FILENAME', filename, before='COMMENT'
    sxaddpar, header, 'DATE_OBS', date_obs, before='COMMENT'
    sxaddpar, header, 'T_OBS', time_obs, before='COMMENT'
    sxaddpar, header, 'DATE-OBS', dateobs, before='COMMENT'
    sxaddpar, header, 'TELESCOP', 'AR30T', before='COMMENT', FORMAT="A5"
    sxaddpar, header, 'INSTRUME', '30 THz camera', before='COMMENT'
    sxaddpar, header, 'DETECTOR', 'FLIR A645sc', before='COMMENT'
    sxaddpar, header, 'EXPTYPE', 'EXPOSURE', before='COMMENT'
    sxaddpar, header, 'WAVELNTH', 10 , ' in micrometers',before='COMMENT'
; 		sxaddpar, header, 'WAVEUNIT', 'Micrometers' ,before='COMMENT'
    sxaddpar, header, 'OBSERVAT', ' Obs. Astronomico Felix Aguilar', before='COMMENT'
    sxaddpar, header, 'PLACE', ' El Leoncito - San Juan - ARGENTINA', before='COMMENT'
    sxaddpar, header, 'LONGITUD', ' -69 19.8', before='COMMENT'
    sxaddpar, header, 'LATITUDE', ' -31 48.1', before='COMMENT'


		sxaddpar, header, 'LVL_NUM', level2, before='COMMENT'
		sxaddpar, header, 'CRPIX1', (NAXISX+1)/2 , before='COMMENT'
    sxaddpar, header, 'CRPIX2', (NAXISY+1)/2 , before='COMMENT'
    sxaddpar, header, 'CRVAL1', '0' , before='COMMENT'
    sxaddpar, header, 'CRVAL2', '0' , before='COMMENT'
    sxaddpar, header, 'CROTA', '-15.774' , ' Rot angle (deg), CCW, from Y.axis.', before='COMMENT'

; calculate the ephemerides
		sun, year, month, day, hour+minute/60.+second/3600., DIST=DIST, SD=SD, TRUE_LONG = TRUE_LONG, TRUE_LAT = TRUE_LAT, $
	APP_LONG = APP_LONG, TRUE_RA = TRUE_RA, TRUE_DEC = TRUE_DEC, APP_RA = APP_RA, APP_DEC=APP_DEC,  LAT0 = LAT0, LONG0=LONG0, $
	PA = PA, CARRINGTON=CARRINGTON

; KEYWORD PARAMETERS:
;       Keywords:
;         /LIST displays values on screen.
;         DIST = distance in AU.
;         SD = semidiameter of disk in arc seconds.
;         TRUE_LONG = true longitude (deg).
;         TRUE_LAT = 0 always.
;         APP_LONG = apparent longitude (deg).
;         APP_LAT = 0 always.
;         TRUE_RA = true RA (hours).
;         TRUE_DEC = true Dec (deg).
;         APP_RA = apparent RA (hours).
;         APP_DEC = apparent Dec (deg).
;         LAT0 = latitude at center of disk (deg).
;         LONG0 = longitude at center of disk (deg).
;         PA = position angle of rotation axis (deg).
;	 			  CARRINGTON = Carrington rotation number.


;**************************************************************************************************************************************
   ; Parameters obtained from ephemerides
;    sxaddpar, header, 'JULDAT', par05, before='COMMENT'
    sxaddpar, header, 'L0', LONG0, before='COMMENT', FORMAT="F10.5"
    sxaddpar, header, 'B0', LAT0, before='COMMENT', FORMAT="F10.5"
    sxaddpar, header, 'P0', PA, before='COMMENT', FORMAT="F10.5"
    sxaddpar, header, 'POS_A', PA, before='COMMENT', FORMAT="F10.5"
    sxaddpar, header, 'SUN_D', DIST,' Distance of the Earth to the Sun in AU', before='COMMENT', FORMAT="F11.8"
;    sxaddpar, header, 'SUN_REF', SD,' Semidiameter of disk in arc seconds', before='COMMENT', FORMAT="F12.5"
    sxaddpar, header, 'SUN_LE', TRUE_LONG, before='COMMENT', FORMAT="F10.5"
    sxaddpar, header, 'SUN_DE', APP_DEC, ' The apparent declination of the Sun in deg',before='COMMENT', FORMAT="F10.5"
    sxaddpar, header, 'SUN_RA',  APP_RA,' The apparent right ascension of the Sun in deg', before='COMMENT', FORMAT="F10.5"
    sxaddpar, header, 'UNITS',' Kelvin', ' Physical units', before='COMMENT'
    sxaddpar, header, 'ORIENT', 0, before='COMMENT'
		sxaddpar, header, 'CARRINGTON', CARRINGTON, ' Carrington rotation number', before='COMMENT', FORMAT="F10.5"

;**************************************************************************************************************************************
; Parameters imported from the original fpf image
;    sxaddpar, header, 'IM01', ,      'image_fpf_id', before='COMMENT'
    sxaddpar, header, 'IM02', fpf.image_version,     ' version', before='COMMENT'
    sxaddpar, header, 'IM03', fpf.image_pixel_offset,' image_pixel_offset', before='COMMENT'
    sxaddpar, header, 'IM04', fpf.image_type,        ' image_image_type', before='COMMENT'
    sxaddpar, header, 'IM05', fpf.image_pixel_format,' image_pixel_format', before='COMMENT'
    sxaddpar, header, 'IM06', fpf.image_trig_count,  ' image_trig_count', before='COMMENT'
    sxaddpar, header, 'IM07', fpf.image_frame_count, ' image_frame_count', before='COMMENT'

    sxaddpar, header, 'CAM01', fpf.camera_partn,          ' camera_partn', before='COMMENT'
    sxaddpar, header, 'CAM02', fpf.camera_sn,             ' camera_sn',before='COMMENT'
    sxaddpar, header, 'CAM03', fpf.camera_range_t_min,    ' camera_range_t_min',before='COMMENT'
    sxaddpar, header, 'CAM04', fpf.camera_range_t_max,    ' camera_range_t_max',before='COMMENT'
    sxaddpar, header, 'CAM05', fpf.camera_lens_name,      ' camera_lens_name',before='COMMENT'
    sxaddpar, header, 'CAM06', fpf.camera_lens_partn,     ' camera_lens_partn',before='COMMENT'
    sxaddpar, header, 'CAM07', fpf.camera_lens_sn,        ' camera_lens_sn',before='COMMENT'
    sxaddpar, header, 'CAM08', fpf.camera_filter_name,    ' filter_name',before='COMMENT'
    sxaddpar, header, 'CAM09', fpf.camera_filter_part_n,  ' filter_part_n',before='COMMENT'
    sxaddpar, header, 'CAM10', fpf.camera_filter_part_sn, ' filter_part_sn',before='COMMENT'

    sxaddpar, header, 'OBJ01', fpf.object_emissivity,    ' emissivity', before='COMMENT'
    sxaddpar, header, 'OBJ02', fpf.object_distance,      ' object_distance', before='COMMENT'
    sxaddpar, header, 'OBJ03', fpf.object_amb_temp,      ' amb_temp', before='COMMENT'
    sxaddpar, header, 'OBJ04', fpf.object_atm_temp,      ' atm_temp',before='COMMENT'
    sxaddpar, header, 'OBJ05', fpf.object_compu_tao,     ' compu_tao',before='COMMENT'
    sxaddpar, header, 'OBJ06', fpf.object_estim_tao,     ' estim_tao',before='COMMENT'
    sxaddpar, header, 'OBJ07', fpf.object_ref_temp,      ' ref_temp',before='COMMENT'
    sxaddpar, header, 'OBJ08', fpf.object_ext_opt_temp,  ' ext_opt_temp',before='COMMENT'
    sxaddpar, header, 'OBJ09', fpf.object_ext_opt_trans, ' ext_opt_trans',before='COMMENT'

    sxaddpar, header, 'SCALING1', fpf.scaling_t_min_cam,  ' t_min_cam',before='COMMENT'
    sxaddpar, header, 'SCALING2', fpf.scaling_t_max_cam,  ' t_max_cam',before='COMMENT'
    sxaddpar, header, 'SCALING3', fpf.scaling_t_min_calc, ' t_min_calc',before='COMMENT'
    sxaddpar, header, 'SCALING4', fpf.scaling_t_max_calc, ' t_max_calc',before='COMMENT'
    sxaddpar, header, 'SCALING5', fpf.scaling_t_min_scale,' t_min_scale',before='COMMENT'
    sxaddpar, header, 'SCALING6', fpf.scaling_t_max_scale,' t_max_scale',before='COMMENT'

   ; sxaddpar, header, 'HISTORY', '' ,before='COMMENT'

 return, header
end


; ###############################################################################################################################################
function makeheader_BR30T, fpf, level
   
    ;This function generates the header for the images
		image=fpf.data
		SS=size(fpf.data)
		datatype=4
		NAXISX = SS[1]
		NAXISY=SS[2]
		naxisarr=[NAXISX,NAXISY] 
		level2 = strmid(level,5,1)+'.'+strmid(level,6,1)
			mkhdr, header, datatype, naxisarr
    
		year = STRTRIM(string(fpf.DATE_YEAR),1)
    month = STRTRIM(string(fpf.DATE_MONTH),1) 
    day = STRTRIM(string(fpf.DATE_DAY),1)
    hour = STRTRIM(string(fpf.DATE_HOUR),1)
    minute = STRTRIM(string(fpf.DATE_MINUTE),1)
    second = STRTRIM(string(fpf.DATE_SECOND),1) + '.' + STRTRIM(string(fpf.DATE_MILLISECOND),1)

  			if strlen(month) eq 1 then month='0'+month
        if strlen(day) eq 1 then day='0'+ day
        if strlen(hour) eq 1 then hour = '0'+ hour
        if strlen(minute) eq 1 then minute = '0' + minute
        if strlen(STRTRIM(string(fpf.DATE_SECOND),1)) eq 1 then second = '0' + second


     ;generate variables dateobs and timeobs
    dateobs=year+'/'+month+'/'+day
    date_obs=year+'-'+month+'-'+day+'T'+hour+':'+minute+':'+second+'Z'
    time_obs=hour+':'+minute+':'+second

   	filename = 'BR30T_'+year+'-'+month+'-'+day+'T'+hour+':'+minute+':'+second+'_'+level

;**************************************************************************************************************************************
		; add the keywords to the header
    sxaddpar, header, 'FILENAME', filename, before='COMMENT'
    sxaddpar, header, 'DATE_OBS', date_obs, before='COMMENT'
    sxaddpar, header, 'T_OBS', time_obs, before='COMMENT'
    sxaddpar, header, 'DATE-OBS', dateobs, before='COMMENT'
    sxaddpar, header, 'TELESCOP', 'BR30T', before='COMMENT', FORMAT="A5"
    sxaddpar, header, 'INSTRUME', '30 THz camera', before='COMMENT'
    sxaddpar, header, 'DETECTOR', 'FLIR A20', before='COMMENT'
    sxaddpar, header, 'EXPTYPE', 'EXPOSURE', before='COMMENT'
    sxaddpar, header, 'WAVELNTH', 10 , ' in micrometers',before='COMMENT'
; 		sxaddpar, header, 'WAVEUNIT', 'Micrometers' ,before='COMMENT'
    sxaddpar, header, 'OBSERVAT', 'CRAAM', before='COMMENT'
    sxaddpar, header, 'PLACE', 'São Paulo - BRAZIL', before='COMMENT'
    sxaddpar, header, 'LONGITUD', '-46 39.1', before='COMMENT'
    sxaddpar, header, 'LATITUDE', '-23 32.8', before='COMMENT'
		sxaddpar, header, 'LVL_NUM', level2, before='COMMENT'
		sxaddpar, header, 'CRPIX1', (NAXISX+1)/2 , ' The location of disk center in x direction', before='COMMENT'
    sxaddpar, header, 'CRPIX2', (NAXISY+1)/2 , ' The location of disk center in y direction', before='COMMENT'
    sxaddpar, header, 'CRVAL1', 0.0 , before='COMMENT'
    sxaddpar, header, 'CRVAL2', 0.0 , before='COMMENT'
;    sxaddpar, header, 'CROTA', '-15.774' , 'Rot angle (deg), CCW, from Y.axis.', before='COMMENT'

; calculate the ephemerides
		sun, year, month, day, hour+minute/60.+second/3600., DIST=DIST, SD=SD, TRUE_LONG = TRUE_LONG, TRUE_LAT = TRUE_LAT, $
	APP_LONG = APP_LONG, TRUE_RA = TRUE_RA, TRUE_DEC = TRUE_DEC, APP_RA = APP_RA, APP_DEC=APP_DEC,  LAT0 = LAT0, LONG0=LONG0, $
	PA = PA, CARRINGTON=CARRINGTON

; KEYWORD PARAMETERS:
;       Keywords:
;         /LIST displays values on screen.
;         DIST = distance in AU.
;         SD = semidiameter of disk in arc seconds.
;         TRUE_LONG = true longitude (deg).
;         TRUE_LAT = 0 always.
;         APP_LONG = apparent longitude (deg).
;         APP_LAT = 0 always.
;         TRUE_RA = true RA (hours).
;         TRUE_DEC = true Dec (deg).
;         APP_RA = apparent RA (hours).
;         APP_DEC = apparent Dec (deg).
;         LAT0 = latitude at center of disk (deg).
;         LONG0 = longitude at center of disk (deg).
;         PA = position angle of rotation axis (deg).
;	 			  CARRINGTON = Carrington rotation number.

;**************************************************************************************************************************************
  ; Parameters obtained from ephemerides
;    sxaddpar, header, 'JULDAT', par05, before='COMMENT'
    sxaddpar, header, 'L0', LONG0, before='COMMENT', FORMAT="F10.5"
    sxaddpar, header, 'B0', LAT0, before='COMMENT', FORMAT="F10.5"
    sxaddpar, header, 'P0', PA, before='COMMENT', FORMAT="F10.5"
    sxaddpar, header, 'POS_A', PA, before='COMMENT', FORMAT="F10.5"
    sxaddpar, header, 'SUN_D', DIST,' Distance of the Earth to the Sun in AU', before='COMMENT', FORMAT="F11.8"
;    sxaddpar, header, 'SUN_REF', SD,'Semidiameter of disk in arc seconds', before='COMMENT', FORMAT="F12.5"
    sxaddpar, header, 'SUN_LE', TRUE_LONG, before='COMMENT', FORMAT="F10.5"
    sxaddpar, header, 'SUN_DE', APP_DEC, ' The apparent declination of the Sun in deg',before='COMMENT', FORMAT="F10.5"
    sxaddpar, header, 'SUN_RA',  APP_RA,' The apparent right ascension of the Sun in deg', before='COMMENT', FORMAT="F10.5"
    sxaddpar, header, 'UNITS',' Celsius', ' Physical units', before='COMMENT'
    sxaddpar, header, 'ORIENT', 0, before='COMMENT'
		sxaddpar, header, 'CARRINGTON', CARRINGTON, ' Carrington rotation number', before='COMMENT', FORMAT="F10.5"

;**************************************************************************************************************************************
; Parameters imported from the original fpf image
;    sxaddpar, header, 'IM01', fpf.image_fpf_id,      'image_fpf_id', before='COMMENT'
    sxaddpar, header, 'IM02', fpf.image_version,     ' version', before='COMMENT'
    sxaddpar, header, 'IM03', fpf.image_pixel_offset,' image_pixel_offset', before='COMMENT'
    sxaddpar, header, 'IM04', fpf.image_type,        ' image_image_type', before='COMMENT'
    sxaddpar, header, 'IM05', fpf.image_pixel_format,' image_pixel_format', before='COMMENT'
    sxaddpar, header, 'IM06', fpf.image_trig_count,  ' image_trig_count', before='COMMENT'
    sxaddpar, header, 'IM07', fpf.image_frame_count, ' image_frame_count', before='COMMENT'

    sxaddpar, header, 'CAM01', fpf.camera_partn,          ' camera_partn', before='COMMENT'
    sxaddpar, header, 'CAM02', fpf.camera_sn,             ' camera_sn',before='COMMENT'
    sxaddpar, header, 'CAM03', fpf.camera_range_t_min,    ' camera_range_t_min',before='COMMENT'
    sxaddpar, header, 'CAM04', fpf.camera_range_t_max,    ' camera_range_t_max',before='COMMENT'
    sxaddpar, header, 'CAM05', fpf.camera_lens_name,      ' camera_lens_name',before='COMMENT'
    sxaddpar, header, 'CAM06', fpf.camera_lens_partn,     ' camera_lens_partn',before='COMMENT'
    sxaddpar, header, 'CAM07', fpf.camera_lens_sn,        ' camera_lens_sn',before='COMMENT'
    sxaddpar, header, 'CAM08', fpf.camera_filter_name,    ' filter_name',before='COMMENT'
    sxaddpar, header, 'CAM09', fpf.camera_filter_part_n,  ' filter_part_n',before='COMMENT'
    sxaddpar, header, 'CAM10', fpf.camera_filter_part_sn, ' filter_part_sn',before='COMMENT'

    sxaddpar, header, 'OBJ01', fpf.object_emissivity,    ' emissivity', before='COMMENT'
    sxaddpar, header, 'OBJ02', fpf.object_distance,      ' object_distance', before='COMMENT'
    sxaddpar, header, 'OBJ03', fpf.object_amb_temp,      ' amb_temp', before='COMMENT'
    sxaddpar, header, 'OBJ04', fpf.object_atm_temp,      ' atm_temp',before='COMMENT'
    sxaddpar, header, 'OBJ05', fpf.object_compu_tao,     ' compu_tao',before='COMMENT'
    sxaddpar, header, 'OBJ06', fpf.object_estim_tao,     ' estim_tao',before='COMMENT'
    sxaddpar, header, 'OBJ07', fpf.object_ref_temp,      ' ref_temp',before='COMMENT'
    sxaddpar, header, 'OBJ08', fpf.object_ext_opt_temp,  ' ext_opt_temp',before='COMMENT'
    sxaddpar, header, 'OBJ09', fpf.object_ext_opt_trans, ' ext_opt_trans',before='COMMENT'

    sxaddpar, header, 'SCALING1', fpf.scaling_t_min_cam,  ' t_min_cam',before='COMMENT'
    sxaddpar, header, 'SCALING2', fpf.scaling_t_max_cam,  ' t_max_cam',before='COMMENT'
    sxaddpar, header, 'SCALING3', fpf.scaling_t_min_calc, ' t_min_calc',before='COMMENT'
    sxaddpar, header, 'SCALING4', fpf.scaling_t_max_calc, ' t_max_calc',before='COMMENT'
    sxaddpar, header, 'SCALING5', fpf.scaling_t_min_scale,' t_min_scale',before='COMMENT'
    sxaddpar, header, 'SCALING6', fpf.scaling_t_max_scale,' t_max_scale',before='COMMENT'

 ;   sxaddpar, header, 'HISTORY', 'Convertion to .fits by "routine transform_fpf_to_fits_30THz.pro"' ,before='COMMENT'

 return, header
end







