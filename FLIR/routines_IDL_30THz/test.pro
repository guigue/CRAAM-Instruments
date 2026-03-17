PRO test
  
    ;Define variables de trabajo
    cantimag=100 	;nro imagenes para calcular del set total, forma un cubo de (640,480,cantimag) puede saturar memoria 
    radiomascara=330.   ;radio en pixeles a considerar en las imagenes de trabajo para eliminar el efecto del limbo que da problema
    threshold= 0.1       ;nivel para la rutina mk_kuhn_flat sobre el cual considera validos los datos de imagenes
                        ; este valor de threshold deberíamos obtenerlo del histograma!!! tomar que el corte sea en 4sigmas para evitar los sectores brillantes!!
    iteraciones= 10      ;numero de iteraciones que utilizará la rutina mk_kuhn_flat (se probaron mas iteraciones y no mejora, sino parece empeorar

    seed = 1000L                ;Numero arbitrario para la secuencia aleatoria de eleccion de imagenes del set total
    INDIR = './flats/'
  
    ;------------------------------------------------------------
    ;Define matrices y variables generales
    window,0, xsize=640, Ysize=480

    Idummy=FLTARR(640,480)
    Icalcu=FLTARR(640,480)
    Icalcu2=FLTARR(640,480)
    Iflat=FLTARR(640,480)
    INDICES=[0,13,26,28,30]
    cantimag=(size(indices))(1)
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
    estaimag=file_search(indir+'*.fits', COUNT=count)
    
    print,count
    if count ne 0 then begin
        tipo='fits'
        estaimag=FINDFILE(indir+'*.fits', COUNT=count)
        estaimag=estaimag(SORT(estaimag))
        Idummy=FLOAT(READFITS(estaimag(0), HeaderFile, /Silent))
        dateobs= sxpar(HeaderFile, 'DATE-OBS')

        ;Define nombre del archivo Flat de salida
        nombre=STRMID(dateobs,0,4)+STRMID(dateobs,5,2)+STRMID(dateobs,8,2)
    endif


    WHILE indice lt cantimag DO BEGIN
       flagformatofpf=0
       fitsname = estaimag[indices[indice2]]
       Idummy=FLOAT(READFITS(fitsname, HeaderFile, /Silent))
       Idummy= Idummy-MIN(Idummy)
       Idummy= Idummy/MAX(Idummy)

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
       ;---------------------------------------------
       ;si el valor del centro del disco solar obtenido esta dentro de ciertos límites, lo considera valido
       print,"  Fits file: ",fitsname,",  Radius = ",parms[0]," Center = (",parms[2],parms[3],")"
        
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
            ;actualiza indices
            indice=indice+1
        endif

        siguiente:
       ;---------------------------------------------
        ;actualiza indices
        indice2=indice2+1
        if indice2 ge cantimag then indice2=0
        
    ENDWHILE


    ;-----------------------------------------------------------------------------
    ;Definido el set de imagenes de trabajo con sus centros y acondicionadas estas, pasa a calcular el Flat
    thresh=threshold
    niter=iteraciones

    mk_kuhn_flat, Icubo, thresh, Iflat,  offsets=desplazamiento,  niter=niter

    ceros=where(Iflat eq 0, Cuenta)
    if Cuenta ne 0 then Iflat(ceros)=1.0

;    Iflat2=Iflat
    save,iflat,filename='flat.save'

    RETURN
    END
