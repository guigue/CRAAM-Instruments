;+----------------------------------------------------------------------
;
; testFITS
;
; This is a simple IDL routine to check the FPF -> FITS procedure.
; The name(s) of the files can be changed (added) in the array
;'files'.
;
; The routine takes the difference between the image data for every
; pixel and compares against the ZERO variable, case it is greater,
; it is considered that the files are different, the program shows
; a message and ends.
;
; Use:
;    testFITS
;
; Written by Guigue@gcastro.net - 2017-10-24
;
;-------------------------------------------------------------------

pro testFITS
  
  !PATH=!PATH+':~/solar/IR'
  
  journal,'testFITS.out'        ; log the outputs
  files=['img0323']

  ZERO = 1.0D-08                ; Definition of Zero

  for n=0,n_elements(files)-1 do begin
     
     read_fpf,files[n]+'.fpf',fpf
     fit=mrdfits(files[n]+'.fits',0,hfit)
     
     Xsize=fpf.header.image_data.xsize
     Ysize=fpf.header.image_data.ysize

     for i=0, xsize - 1 do begin
        for j=0, ysize - 1 do begin
           if ( (fpf.data[i,j]-fit[i,j]) gt ZERO ) then begin
              print, files[n]+'.fpf and '+files[n]+'.fit are different '
              print, 'Element ('+string(i)+','+string(j)+') = ', fpf.data[i,j], fit[i,j]
              print, 'Exiting....'
              return
           endif
        endfor
     endfor

     print, ' '
     print, 'Congrats!!  '+files[n]+'.fpf and '+files[n]+'.fit are identical '
     print, ' '
  endfor
  return
end

  
  
