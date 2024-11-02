function getsize,struc

n=n_tags(struc)

sum=0
for i=0,n-1 do begin
    nelements=size(struc.(i),/n_elements)
    tipo=size(struc.(i),/type)
    case tipo of
        1:    nbytes=1 ; byte
        2:    nbytes=2 ; integer
        3:    nbytes=4 ; long integer
        4:    nbytes=4 ; floating 
        5:    nbytes=8 ; double precision floating
        6:    nbytes=8 ; complex single precision
        9:    nbytes=16; complex double precision
        12:   nbytes=2 ; unsigned integer
        13:   nbytes=4 ; long unsigned integer
        else: nbytes=0
    endcase
    sum = sum + nelements * nbytes

endfor
return,sum
end

;+
;
; NAME: sst_struc
;
; PURPOSE: create structure of data to readu SST/ITA/12GHZ binary archives
;
; INPUT: (You must pass the keyword: n_struct)
;
; OUPUT: 
;          raw_struc: a named structure, the name is composed by the
;                     identification file and the date.  v.g. RS1021225, 'slow'
;                     (integrated SST dat) of 2002/12/25
;
; KEYWORD: odata: observational data.  If set the structure will be of type
;                 struct_data of buffers.h  
;          monitor: Monitoring Data.  If set the structure will be of type
;                   struct_monitor in buffers.h
;          display: Display Data.  If set the structure will be of type
;                   struct_send_data1 in buffers.h
;          n_struct: name of the structure
;
; COMMENTS: This new version of sst_struc replaces the old one.  In this
; new version, structures are declared here and the system take into
; account the date to know the different kinds of structures.  Every time,
; buffers.h is changed, this file should reflect the changes.
;
; HISTORY
;         Guigue - Dec 2002, Updated 2017-07-29 (see comments at the end)
;
;-

pro sst_struc, raw_struc, raw_size, odata=odata, monitor=monit, display=display, mdata=mdata, $
                    n_struc=n_struc,def_file=def_file

fecha = long(strmid(n_struc,2,strlen(n_struc)-1))
kind  = strmid(n_struc,0,2) 

if keyword_set(monit) then begin

   if (fecha gt 1021213 ) then begin
 
        nchan = 6
        nradiom = 6
        cmd='raw_struc={'+n_struc+',time:0l, azipos:0.0, elepos:0.0, azierr:0.0, eleerr:0.0,'
        cmd=cmd+'adc:uintarr(nchan), sigma:fltarr(nchan), gps_status:0,acq_gain:0,'
        cmd=cmd+'target:0b, opmode:0b, off:intarr(nradiom), hot_temp:0.0, amb_temp:0.0,'
        cmd=cmd+'opt_temp:0.0,if_board_temp:0.0,radome_temp:0.0,humidity:0.0,temperature:0.0,opac_210:0.0,'
        cmd=cmd+'opac_405:0.0,elevation:0.0, pressure :0.0,burst :0b,errors:0l}'

    endif


    if (fecha gt 1021123 and fecha le  1021213 ) then begin
 
        nchan = 6
        nradiom = 6
        cmd='raw_struc={'+n_struc+',time:0l, azipos:0.0, elepos:0.0, azierr:0.0, eleerr:0.0,'
        cmd=cmd+'adc:lonarr(nchan), sigma:fltarr(nchan), gps_status:0, daq_status:0, acq_gain:0,'
        cmd=cmd+'target:0b, opmode:0b, off:intarr(nradiom), hot_temp:0.0, amb_temp:0.0,'
        cmd=cmd+'opt_temp:0.0,if_board_temp:0.0,radome_temp:0.0,humidity:0.0,temperature:0.0,opac_210:0.0,'
        cmd=cmd+'opac_405:0.0,elevation:0.0, pressure :0.0,burst :0b,errors:0l}'

    endif


    if (fecha gt 1020915 and fecha le 1021123) then begin
 
        nchan   = 8
        nradiom = 6
        cmd='raw_struc={'+n_struc+',time:0l, azipos:0.0, elepos:0.0, azierr:0.0, eleerr:0.0,'
        cmd=cmd+'adc:intarr(nchan), sigma:fltarr(nchan), gps_status:0, daq_status:0, acq_gain:0,'
        cmd=cmd+'target:0b, opmode:0b, att:intarr(nradiom),off:intarr(nradiom),  mix_volt:fltarr(nradiom),'
        cmd=cmd+'mix_curr:fltarr(nradiom),hot_temp:0.0, amb_temp:0.0,opt_temp:0.0,'
        cmd=cmd+'if_board_temp:0.0,radome_temp:0.0,humidity:0.0,temperature:0.0,opac_210:0.0,'
        cmd=cmd+'opac_405:0.0,elevation:0.0, pressure :0.0,burst :0b,errors:0l}'

    endif

    if (fecha le 1020915 ) then begin

        nradiom = 6
        nchan   = 8

        cmd='raw_struc={'+n_struc+',time:0l, azipos:0.0, elepos:0.0, azierr:0.0, eleerr:0.0,'
        cmd=cmd+'adc:intarr(nchan), sigma:fltarr(nchan), gps_status:0, daq_status:0, acq_gain:0,'
        cmd=cmd+'target:0b, opmode:0b,att:intarr(nradiom),off:intarr(nradiom),mix_volt:fltarr(nradiom),'
        cmd=cmd+'mix_curr:fltarr(nradiom),hot_temp:0.0, amb_temp:0.0,opt_temp:0.0,'
        cmd=cmd+'if_board_temp:0.0,radome_temp:0.0,humidity:0.0,wind:0.0,opac_210:0.0,'
        cmd=cmd+'opac_405:0.0,elevation:0.0, seeing:0.0,burst :0b,errors:0l}'

    endif

    res=execute(cmd)
    for i=0,n_tags(raw_struc)-1 do raw_struc.(i)=0    
    raw_size=getsize(raw_struc)
    return

endif

if keyword_set(display) then begin

    nchan   = 6
    nradiom = 6

    cmd = 'raw_struc={'+n_struc+'date:0L,time:0L,adcval:uintarr(nchan),adc_min:uintarr(nchan),'
    cmd = cmd + 'adc_max:uintarr(nchan),sigma:fltarr(nchan),azipos:0.0,elepos:0.0,azierr:0.0,'
    cmd = cmd + 'eleerr:0.0,x_off:0.0,y_off:0.0,off:intarr(nradiom),call_mirr_pos:0L,'
    cmd = cmd + 'gps_status:0,hot_temp:0.0,amb_temp:0.0,opt_temp:0.0,'
    cmd = cmd + 'if_board_temp:0.0,radome_temp:0.0,humidity:0.0,temperature:0.0,opac_210:0.0,'
    cmd = cmd + 'opac_405:0.0,pressure:0.0, burst:0B,target:0B,opmode:0B,dummy:0B,errors:0L}'

    res=execute(cmd)
    for i=0,n_tags(raw_struc)-1 do raw_struc.(i)=0    
    raw_size=getsize(raw_struc)
    return

endif

;---------------- observational data ----------------------------------------
if keyword_set(odata) then begin

;
;  Itapetinga Data 14m dish
;
    if (kind eq 'is' or kind eq 'if' or kind eq 'ir') then begin

        nchan   = 8
        nradiom = 5
        cmd = 'raw_struc={'+n_struc+',time:0l, adcval:intarr(nchan), pos_time:0l,'
        cmd = cmd + 'azipos:0l, elepos:0l, azivel:0, elevel:0, azierr:0l,'
        cmd = cmd + 'eleerr:0l, x_off:0, y_off:0, att:intarr(nradiom), '
        cmd = cmd + 'target:0b, opmode:0b, gps_status:0, daq_status:0, recnum:0l}'
;        if (kind eq 'ir') then begin
;           nradiom = 6
;           cmd = 'raw_struc={'+n_struc+',time:0l, adcval:intarr(nchan), pos_time:0l,'
;           cmd = cmd + 'azipos:0l, elepos:0l, azivel:0, elevel:0, azierr:0l,'
;           cmd = cmd + 'eleerr:0l, x_off:0, y_off:0, att:intarr(nradiom), '
;           cmd = cmd + 'target:0b, opmode:0b, gps_status:0, recnum:0l}'
;        endif

    endif else begin

;
; SST Data
;
        
        nradiom = 6

        if (fecha gt 1021213) then begin

            nchan   = 6
            cmd ='raw_struc={'+n_struc+',time:0l, adcval:uintarr(nchan), pos_time:0l,'
            cmd = cmd + 'azipos:0l, elepos:0l, pm_daz:0, pm_del:0, azierr:0l,'
            cmd = cmd + 'eleerr:0l, x_off:0, y_off:0, off:intarr(nradiom),'
            cmd = cmd + 'target:0b, opmode:0b, gps_status:0, recnum:0l}'
            
        endif

        if (fecha gt 1021203 and fecha le 1021213) then begin

            nchan   = 6
            cmd = 'raw_struc={'+n_struc+',time:0l, adcval:lonarr(nchan), pos_time:0l,'
            cmd = cmd + 'azipos:0l, elepos:0l, pm_daz:0, pm_del:0, azierr:0l,'
            cmd = cmd + 'eleerr:0l, x_off:0, y_off:0, off:intarr(nradiom),'
            cmd = cmd + 'target:0b, opmode:0b, gps_status:0, recnum:0l}'
            
        endif

        if (fecha gt 1020520 and fecha le 1021203) then begin

            nchan   = 8
            
            cmd = 'raw_struc={'+n_struc+',time:0l, adcval:intarr(nchan), pos_time:0l,'
            cmd = cmd + 'azipos:0l, elepos:0l, pm_daz:0, pm_del:0, azierr:0l,'
            cmd = cmd + 'eleerr:0l, x_off:0, y_off:0, att:intarr(nradiom),'
            cmd = cmd + 'target:0b, opmode:0b, gps_status:0, daq_status:0, recnum:0l}'

        endif

        if (fecha gt 990501 and fecha le 1020520) then begin

            nchan   = 8
            cmd = 'raw_struc={'+n_struc+',time:0l, adcval:intarr(nchan), pos_time:0l,'
            cmd = cmd + 'azipos:0l, elepos:0l, azivel:0, elevel:0, azierr:0l,'
            cmd = cmd + 'eleerr:0l, x_off:0, y_off:0, att:intarr(nradiom), '
            cmd = cmd + 'target:0b, opmode:0b, gps_status:0, daq_status:0, recnum:0l}'
            
        endif

        if (fecha le 990501) then begin

            nchan   = 8
            cmd = 'raw_struc={'+n_struc+',time:0l, adcval:intarr(nchan), pos_time:0l,'
            cmd = cmd + 'azipos:0l, elepos:0l, azivel:0, elevel:0, azierr:0l,'
            cmd = cmd + 'eleerr:0l, att:intarr(nradiom), target:0b, opmode:0b,'
            cmd = cmd + 'gps_status:0, daq_status:0, recnum:0l}'
       
        endif

    endelse

    res=execute(cmd)
    for i=0,n_tags(raw_struc)-1 do raw_struc.(i)=0    
    raw_size=getsize(raw_struc)

endif


return

;---------------------------------------------------------------------------------------
;
; Comments to the 2012-07-29 version
;
; In the process to transform SST data to FITS files I discovered an error in 
; the auxiliary (monitoring) data staructure that affected the reading for the
; timespan 2002-09-16 to 2002-11-23. In the 2002-12 version off was set as
; a float and mix_curr as an integer.  In fact off is an integer and mix_curr is 
; a float. 
;
;----------------------------------------------------------------------------------------

end
