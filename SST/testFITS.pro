;+----------------------------------------------------------------------
;
; testFITS
;
; This is the second part (IDL) of the test of the conversion procedures
; used to convert SST raw data in FITS. The first part is testRBD.py
; Every FITS created corresponds to one of the different formats used
; by SST in its history, since 1999-05-01.  The IDL checks that the
; convertion is perfect, meaning with this, that the difference
; between the original Raw Binary Data has difference 0 with the
; created FITS file.
;
; Since some other functionalities were introduced in the python
; module oRBD, like concat() and reduced(), these are tested as well.
;
; If the IDL pro does not stop means the test was succesfull. It also
; has now an MD5 checksum included
;
; Use:
;    testFITS
;
; Written by Guigue@gcastro.net - 2017-08-09
;
;-------------------------------------------------------------------

pro testFITS

  journal,'testFITS.out'        ; log the outputs
  
  RBD_bifiles=['bi1010822','bi1021019','bi1021202','bi1021221']

  RBD_rsfiles=['rs1020715.1300','rs1021205.2200',$
                           'rs1061206.2100','rs990909.1700']

  FITS_bifiles = ['sst_auxiliary_2001-08-22T14:04:40.034-23:45:04.714_level0.fits',$
                  'sst_auxiliary_2002-10-19T11:24:39.704-21:37:04.235_level0.fits',$
                  'sst_auxiliary_2002-12-02T02:01:34.154-23:18:05.244_level0.fits',$
                  'sst_auxiliary_2002-12-21T11:15:57.732-22:26:07.212_level0.fits']
  
  FITS_rsfiles = ['sst_integration_2002-07-15T13:24:19.605-13:59:59.592_level0.fits',$
                  'sst_integration_2002-12-05T21:59:58.728-22:08:31.215_level0.fits',$
                  'sst_integration_2006-12-06T20:59:58.423-21:54:05.746_level0.fits',$
                  'sst_integration_1999-09-09T16:59:52.056-17:57:40.420_level0.fits']

  for i=0, n_elements(RBD_bifiles)-1 do begin
     read_sst,rbd,RBD_bifiles[i],/mon,/close,recr=1000000
     trbd=tag_names(rbd)
     
     fits=mrdfits(FITS_bifiles[i],1,h,/unsigned)
     tfits=tag_names(fits)

     if (n_elements(trbd) ne n_elements(tfits)) then begin
        print,''
        print,'Structures are different for files '+RBD_bifiles[i]+' and '+FITS_bifiles[i]
        print,'Aborting....'
        print,''
        stop
     endif

     print,RBD_bifiles[i]+' and '+FITS_bifiles[i]
     for itag=0,n_elements(trbd)-1 do begin
        command='m=moment(rbd.'+trbd[itag]+'-fits.'+tfits[itag]+')'
        r=execute(command)
        if (abs(m[0]) gt 0) then begin
           print,'Tag = '+trbd(itag)+' , '+tfits(itag)
           print,'       Difference = '+string(m[0])+'  sigma = '+string(sqrt(m[1]))
           stop
        endif
        
     endfor
     
  endfor

  for i=0, n_elements(RBD_rsfiles)-1 do begin
     read_sst,rbd,RBD_rsfiles[i],/close,recr=1000000
     trbd=tag_names(rbd)
     
     fits=mrdfits(FITS_rsfiles[i],1,h,/unsigned)
     tfits=tag_names(fits)

     if (n_elements(trbd) ne n_elements(tfits)) then begin
        print,''
        print,'Structures are different for files '+RBD_rsfiles[i]+' and '+FITS_rsfiles[i]
        print,'Aborting....'
        print,''
        stop
     endif

     print,RBD_rsfiles[i]+' and '+FITS_rsfiles[i]
     for itag=0,n_elements(trbd)-1 do begin
        command='m=moment(rbd.'+trbd[itag]+'-fits.'+tfits[itag]+')'
        r=execute(command)
        if (abs(m[0]) gt 0) then begin
           print,'Tag = '+trbd(itag)+' , '+tfits(itag)
           print, 'Difference = '+string(m[0])+'  sigma = '+string(sqrt(m[1]))
           stop
        endif
        
     endfor
     
  endfor
  
  fits=mrdfits('rs1150621.17-18-concat.fits',1,h,/unsigned)
  read_sst,d1,'rs1150621.1700',/close,recr=1000000
  read_sst,d2,'rs1150621.1800',/close,recr=1000000
  rbd=[d1,d2]
  print, 'rs1150621.1700   rs1150621.1800 rs1150621.17-18-concat.fits'
  
  trbd=tag_names(rbd)
  tfits=tag_names(fits)

  if (n_elements(trbd) ne n_elements(tfits)) then begin
     print,''
     print,'Structures are different for files rs1150621.17-18-concat.fits and '+$
           'rs1150621.1700 + rs1150621.1800'
     print,'Aborting....'
     print,''
     stop
  endif

  for itag=0,n_elements(trbd)-1 do begin
     command='m=moment(rbd.'+trbd[itag]+'-fits.'+tfits[itag]+')'
     r=execute(command)
     if (abs(m[0]) gt 0) then begin
        print,'Tag = '+trbd(itag)+' , '+tfits(itag)
        print, 'Difference = '+string(m[0])+'  sigma = '+string(sqrt(m[1]))
        stop
     endif

  endfor
  
  fits=mrdfits('rs1150621.17-18-red-concat.fits',1,h,/unsigned)
  tfits=tag_names(fits)

  for itag=0,n_elements(tfits)-1 do begin
     itrbd = (where(trbd eq tfits[itag]))(0)
     command='m=moment(rbd.'+trbd[itrbd]+'-fits.'+tfits[itag]+')'
     r=execute(command)
     if (abs(m[0]) gt 0) then begin
        print,'Tag = '+trbd(itrbd)+' , '+tfits(itag)
        print, 'Difference = '+string(m[0])+'  sigma = '+string(sqrt(m[1]))
        stop
     endif

  endfor

  journal                       ; stop logging
  MD5 = '00B6121EF75BFA12FCDAA6FD50A736D0'xULL
  spawn,"sed '1,5d' testFITS.out > t; mv -f t testFITS.out" 
  spawn,'md5sum testFITS.out',r
  r=strsplit(r,/extract)
  comando="md5sum='"+strupcase(r[0])+"'xULL"
  c=execute(comando)
  if (md5sum eq MD5) then begin
     print, ''
     print,''
     print,'Tests run Ok!'
     print,''
     print,''
  endif else begin
     print,''
     print,''
     print,'MD5 checksum differs'
     print,'Test does not pass'
     print,'Original = ', MD5
     print,'Obtained = ',md5sum
     print,''
     print,''
  endelse
  
end
