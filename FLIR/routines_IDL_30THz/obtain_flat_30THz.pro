PRO  obtain_flat_30THz, indir=indir, telescope=telescope


   ; ****************************************************************************************************************************************
    ;Este proceso obtiene una mascara de Flat para la camara 30THz utilizando el metodo de Kuhn, Lin & Loranz 1991
    ;Utiliza la rutina mk_kuhn_flat del solar soft
    ;Lee las imagenes de entrada en formato fpf (Flir Public Format) o fits del directorio dirimag.
    ;Escibe la imagen de Flat en formato fits y jpg en directorio de salida, con nombre la fecha de las imagenes
    ;Para aplicar la corrección de Flat hay que hacer Imagen=Imagen/Imagenflat

    ;INPUT 
    ;Define input directory where the N flat images are located
    ; indir = 'path of images' (eg: indir = '\home\Fernando\flatimages\')
    ; telescope='AR30T' or 'BR30T'
    ;
    ;
    ; History:
    ; written by Fernando M. Lopez (CRAAM-Mackenzie) --- October 2019
    ; The routine is based in the work done by Carlos Francile and Franco Manini (OAFA-UNSJ, San Juan, Argentina)
    ; ****************************************************************************************************************************************


    loadct,0

   
		if not keyword_set(indir) then message, '*** MUST INDICATE AS INPUT KEYWORD THE PATH OF THE IMAGES, STOPPING PROCESS ***'
		if not keyword_set(telescope) then message, '*** MUST INDICATE AS INPUT KEYWORD THE TELESCOPE, STOPPING PROCESS ***'
		if telescope ne 'AR30T' then if telescope ne 'BR30T' then $
						message, '***** INCORRECT TELESCOPE (AR30T - BR30T), STOPPING CALIBRATION *****' 

  	 outdir=indir+'final/'

		; create the output folder
		file_mkdir, outdir

    ;Define variables de trabajo
    cantimag=100 	;nro imagenes para calcular del set total, forma un cubo de (640,480,cantimag) puede saturar memoria 
    radiomascara=330.   ;radio en pixeles a considerar en las imagenes de trabajo para eliminar el efecto del limbo que da problema
    threshold= 0.1       ;nivel para la rutina mk_kuhn_flat sobre el cual considera validos los datos de imagenes
                        ; este valor de threshold deberíamos obtenerlo del histograma!!! tomar que el corte sea en 4sigmas para evitar los sectores brillantes!!
    iteraciones= 10      ;numero de iteraciones que utilizará la rutina mk_kuhn_flat (se probaron mas iteraciones y no mejora, sino parece empeorar

    seed = 1000L         ;Numero arbitrario para la secuencia aleatoria de eleccion de imagenes del set total

    ;------------------------------------------------------------
    ;Define matrices y variables generales
    window,0, xsize=640, Ysize=480

    Idummy=FLTARR(640,480)
    Icalcu=FLTARR(640,480)
    Icalcu2=FLTARR(640,480)
    Iflat=FLTARR(640,480)
    desplazamiento=FLTARR(2,cantimag)
    Icubo=FLTARR(640,480,cantimag)
    indice=0
    indice2=0
    HeaderFile=STRARR(200)

    ;------------------------------------------------------------
    ;Crea una mascara circular con unos y ceros que se utiliza para eliminar los bordes de imágenes
    ;(dejar fuera el limbo)
    mascara=FLTARR(2560,1920)

    for i=0,2559 do begin
        for j=0,1919 do begin
           if SQRT((i-1280.)^2 + (j-960.)^2) lt radiomascara then mascara(i,j) = 1.0
        endfor
    endfor

    ;------------------------------------------------------------
    ;Busca imagenes de trabajo presentes en el directorio
    estaimag=file_search(indir+'*.fts', COUNT=count)
    if count ne 0 then begin
        tipo='fits'
        estaimag=FINDFILE(dirimag+'*.fts', COUNT=count)
        estaimag=estaimag(SORT(estaimag))
        Idummy=FLOAT(READFITS(estaimag(0), HeaderFile, /Silent))
        dateobs= sxpar(HeaderFile, 'DATE-OBS')

        ;Define nombre del archivo Flat de salida
        nombre=STRMID(dateobs,0,4)+STRMID(dateobs,5,2)+STRMID(dateobs,8,2)
    endif

    estaimag=file_search(indir+'*.fpf', COUNT=count)
    if count ne 0 then begin
        tipo='fpf'
        estaimag=file_search(indir+'*.fpf', COUNT=count)
        estaimag=estaimag(SORT(estaimag))
        fpf = read_fpf(estaimag(0))
        año=STRTRIM(STRING(fpf.date_year),2)
        mes=STRTRIM(STRING(fpf.date_month),2)
        if strlen(mes)eq 1 then mes='0'+ mes
        dia=STRTRIM(STRING(fpf.date_day),2)
        if strlen(dia)eq 1 then dia='0'+ dia

        ;Define nombre del archivo Flat de salida
        nombre=telescope+'_'+año+'-'+mes+'-'+dia+'_FLAT'

    endif

    if tipo eq 'nada' then begin
        print, 'ERROR, no se encontraron imagenes'
        goto, final
    endif

;-----------------------------------------------------------------------------
;PASO 1
;-----------------------------------------------------------------------------
    ;------------------------------------------------------------
    ;Genera un indice aleatorio de las imagenes del set total
    aleatorio=FIX(randomu(seed,count) * count)

    ;------------------------------------------------------------
    ;Comienza iteración
    WHILE indice lt cantimag DO BEGIN

       ;---------------------------------------------
       ;lee una imagen del indice aleatorio
       if tipo ne 'fits' then begin
            flagformatofpf=1
            fpf = read_fpf(estaimag(aleatorio(indice2)))
            Idummy=fpf.data

       endif else begin
            flagformatofpf=0
            Idummy=FLOAT(READFITS(estaimag(aleatorio(indice2)), HeaderFile, /Silent))
            Idummy= Idummy-MIN(Idummy)
            Idummy= Idummy/MAX(Idummy)
       endelse
       ;---------------------------------------------
        ;Calcula ofsett del centro de la imagen

        maximo=MAX(Idummy)

        result=WHERE(Idummy lt (maximo-maximo/5.), count2, COMPLEMENT=notresult)
        if count2 eq 0 then begin
            indice=indice+1
            goto, siguiente
       endif
       ;---------------------------------------------
       ;Genera en Icalcu un arreglo donde esta definido el limbo del disco solar
       ;utilizando el filtro de Roberts
        Icalcu(result)=0.
        Icalcu(notresult)=1.
        Icalcu=ROBERTS(Icalcu)

        result2=where(Icalcu eq 1)
        Result3 = ARRAY_INDICES(Icalcu, result2)

       ;---------------------------------------------
       ;interpola un circulo a la posición del limbo de Icalcu, utilizando el paquete mpfit (L-Markwardt)
        start_parms=[320.,320.,240.0,240.0,300.]
        parms = MPFITELLIPSE(FLOAT(REFORM(Result3(0,*))), FLOAT(REFORM(Result3(1,*))), start_parms,/CIRCULAR, QUIET=1);, /CIRCULAR)
        ;      P[0]   Ellipse semi axis 1
        ;      P[1]   Ellipse semi axis 2   ( = P[0] if CIRCLE keyword set)
        ;      P[2]   Ellipse center - x value
        ;      P[3]   Ellipse center - y value
        ;      P[4]   Ellipse rotation angle (radians) if TILT keyword set
       ;---------------------------------------------
       ;si el valor del centro del disco solar obtenido esta dentro de ciertos límites, lo considera valido
        if parms(2) gt 100 and parms(2) lt 539 and parms(3) gt 100 and parms(3) lt 379 then begin
            ;---------------------------------------------
            ;se guarda en desplazamiento del centro de la imagen considerada
            desplazamiento(*,indice)=[parms(2), parms(3)]


            ;---------------------------------------------
            ;Guarda en el cubo el sector central de la imagen considerada utilizando la mascara, prescindiendo de los bordes.
            Idummy=Idummy*mascara(1280-320-FIX(parms(2)-320.):1280+319-FIX(parms(2)-320.),960-240-FIX(parms(3)-240.):960+239-FIX(parms(3)-240.))
            rrr=where(Idummy eq 0.)
            Idummy(rrr)=0.001     ;elimina valores cero de la mascara para evitar divisiones por cero
            Icubo(*,*,indice)=Idummy

            ;---------------------------------------------
            ;muestra la imagen
            tvscl, Idummy


         ;---------------------------------------------
            ;actualiza indices
            indice=indice+1
        endif

        siguiente:
       ;---------------------------------------------
        ;actualiza indices
        indice2=indice2+1
        if indice2 ge cantimag then indice2=0
       wait, 0.1
    ENDWHILE


    ;-----------------------------------------------------------------------------
    ;Definido el set de imagenes de trabajo con sus centros y acondicionadas estas, pasa a calcular el Flat
    thresh=threshold
    niter=iteraciones

    mk_kuhn_flat, Icubo, thresh, Iflat,  offsets=desplazamiento,  niter=niter

    ceros=where(Iflat eq 0, Cuenta)
    if Cuenta ne 0 then Iflat(ceros)=1.0

    Iflat2=Iflat

    print, 'Paso 2, recalcula sobre imagenes corregidas por el flat anterior'
;-----------------------------------------------------------------------------------------
;PASO 2;
;-----------------------------------------------------------------------------------------
    Idummy(*,*)=0;=FLTARR(640,480)
    Icalcu(*,*)=0;=FLTARR(640,480)
    Icalcu2(*,*)=0;=FLTARR(640,480)
    Iflat(*,*)=0;=FLTARR(640,480)
    desplazamiento(*,*)=0;=FLTARR(2,cantimag)
    Icubo(*,*,*)=0;=FLTARR(640,480,cantimag)


    seed = 1001L  ;elige otra secuencia aleatoria
    indice=0
    indice2=0

    ;------------------------------------------------------------
    ;Genera un indice aleatorio de las imagenes del set total
    aleatorio=FIX(randomu(seed,count) * count)


    ;------------------------------------------------------------
    ;Comienza iteración
    WHILE indice lt cantimag DO BEGIN

       ;---------------------------------------------
       ;lee una imagen del indice aleatorio
       if tipo ne 'fits' then begin
            fpf = read_fpf(estaimag(aleatorio(indice2)))
            Idummy=fpf.data
            Idummy=Idummy/Iflat2

       endif else begin
            Idummy=FLOAT(READFITS(estaimag(aleatorio(indice2)), HeaderFile, /Silent))
            Idummy= Idummy-MIN(Idummy)
            Idummy= Idummy/MAX(Idummy)
       endelse
       ;---------------------------------------------
        ;Calcula ofsett del centro de la imagen

        maximo=MAX(Idummy)

        result=WHERE(Idummy lt (maximo-maximo/5.), count2, COMPLEMENT=notresult)
        if count2 eq 0 then begin
            indice=indice+1
            goto, siguiente2
       endif
       ;---------------------------------------------
       ;Genera en Icalcu un arreglo donde esta definido el limbo del disco solar
       ;utilizando el filtro de Roberts
        Icalcu(result)=0.
        Icalcu(notresult)=1.
        Icalcu=ROBERTS(Icalcu)

        result2=where (Icalcu eq 1)
        Result3 = ARRAY_INDICES(Icalcu, result2)

       ;---------------------------------------------
       ;interpola un circulo a la posición del limbo de Icalcu, utilizando el paquete mpfit (L-Markwardt)
        start_parms=[320.,320.,240.0,240.0,300.]
        parms = MPFITELLIPSE(FLOAT(REFORM(Result3(0,*))), FLOAT(REFORM(Result3(1,*))), start_parms,/CIRCULAR, QUIET=1);, /CIRCULAR)
        ;      P[0]   Ellipse semi axis 1
        ;      P[1]   Ellipse semi axis 2   ( = P[0] if CIRCLE keyword set)
        ;      P[2]   Ellipse center - x value
        ;      P[3]   Ellipse center - y value
        ;      P[4]   Ellipse rotation angle (radians) if TILT keyword set
       ;---------------------------------------------
       ;si el valor del centro del disco solar obtenido esta dentro de ciertos límites, lo considera valido
        if parms(2) gt 100 and parms(2) lt 539 and parms(3) gt 100 and parms(3) lt 379 then begin

            ;---------------------------------------------
            ;se guarda en desplazamiento del centro de la imagen considerada
            desplazamiento(*,indice)=[parms(2), parms(3)]


            ;---------------------------------------------
            ;Guarda en el cubo el sector central de la imagen considerada utilizando la mascara, prescindiendo de los bordes.
            Idummy=Idummy*mascara(1280-320-FIX(parms(2)-320.):1280+319-FIX(parms(2)-320.),960-240-FIX(parms(3)-240.):960+239-FIX(parms(3)-240.))
            rrr=where(Idummy eq 0.)
            Idummy(rrr)=0.001     ;elimina valores cero de la mascara para evitar divisiones por cero
            Icubo(*,*,indice)=Idummy

            ;---------------------------------------------
            ;muestra la imagen
            ;tvscl, Idummy*Iflat2*Iflat


            ;---------------------------------------------
            ;actualiza indices
            indice=indice+1
        endif

        siguiente2:
       ;---------------------------------------------
        ;actualiza indices
        indice2=indice2+1
        if indice2 ge cantimag then indice2=0
    wait, 0.1
    ENDWHILE


    ;-----------------------------------------------------------------------------
    ;Definido el set de imagenes de trabajo con sus centros y acondicionadas estas, pasa a calcular el Flat
    thresh=threshold
    niter=iteraciones

    mk_kuhn_flat, Icubo, thresh, Iflat,  offsets=desplazamiento,  niter=niter

    ceros=where(Iflat eq 0, Cuenta)
    if Cuenta ne 0 then Iflat(ceros)=1.0


    ;Finalizado el Paso 2, el Flat final será el producto de ambos
    Iflat=Iflat*Iflat2
    ;-------------------------------------------------------
    ;Graba la imagen de flat
    filename = nombre+'.fits'
		dateobs=año+'-'+mes+'-'+dia
; ##########################################################################################################################
    MKHDR, headerfile, Iflat
		; add few keywords to the header
		if telescope eq 'AR30T' then begin
				sxaddpar, headerfile, 'FILENAME', filename, before='COMMENT'
				sxaddpar, headerfile, 'DATE-OBS', dateobs, before='COMMENT'
				sxaddpar, headerfile, 'TELESCOP', 'AR30T', before='COMMENT'
				sxaddpar, headerfile, 'INSTRUME', '30 THz camera', before='COMMENT'
				sxaddpar, headerfile, 'DETECTOR', 'FLIR A645sc', before='COMMENT'
				sxaddpar, headerfile, 'EXPTYPE', 'FLAT', before='COMMENT'
				sxaddpar, headerfile, 'WAVELNTH', 10 , ' in micrometers',before='COMMENT'
				sxaddpar, headerfile, 'OBSERVAT', ' Obs. Astronomico Felix Aguilar', before='COMMENT'
				sxaddpar, headerfile, 'PLACE', ' El Leoncito - San Juan - ARGENTINA', before='COMMENT'
				sxaddpar, headerfile, 'LONGITUD', ' -69 19.8', before='COMMENT'
				sxaddpar, headerfile, 'LATITUDE', ' -31 48.1', before='COMMENT'		
		endif

		if telescope eq 'BR30T' then begin
				sxaddpar, headerfile, 'FILENAME', filename, before='COMMENT'
				sxaddpar, headerfile, 'DATE-OBS', dateobs, before='COMMENT'
				sxaddpar, headerfile, 'TELESCOP', 'BR30T', before='COMMENT'
				sxaddpar, headerfile, 'INSTRUME', '30 THz camera', before='COMMENT'
				sxaddpar, headerfile, 'DETECTOR', 'FLIR A20', before='COMMENT'
				sxaddpar, headerfile, 'EXPTYPE', 'FLAT', before='COMMENT'
				sxaddpar, headerfile, 'WAVELNTH', 10 , ' in micrometers',before='COMMENT'
				sxaddpar, headerfile, 'OBSERVAT', 'CRAAM', before='COMMENT'
				sxaddpar, headerfile, 'PLACE', 'São Paulo - BRAZIL', before='COMMENT'
				sxaddpar, headerfile, 'LONGITUD', '-46 39.1', before='COMMENT'
				sxaddpar, headerfile, 'LATITUDE', '-23 32.8', before='COMMENT'
		endif
	
    writefits, outdir+filename, Iflat, headerfile


    ;-------------------------------------------------------
    ;Muestra la imagen de flat y la graba en formato jpg
    tvscl, Iflat*Iflat2
    TVLCT, R, G, B, /GET
    Imagejpg = TVRD(  TRUE=1)
    WRITE_JPEG , outdir+nombre+'.jpg', Imagejpg, QUALITY=100, TRUE=1


final:
print, 'Termino el proceso, grabo flat '+nombre
END
