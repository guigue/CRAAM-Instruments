# External methods
import string, os, struct, glob
import numpy as np
import datetime as dt
from enum import IntEnum

# Our methods
import astronomical_methods as am
import oRBD

class SST(object):

###############################################################################################
#
#    SST: A python class to use with SST data. It is the OOP evolution of the original
#         sst_methods.py 
#
#    Author:  Guigue@Sampa
#             guigue@craam.mackenzie.br  , guigue@gcastro.net
#             http://www.guigue.gcastro.net
#
#    Change Record :  
#         - 2017-06-15 : First written
#         - 2017-06-17 : time now is python datetime format compatible.
#
###############################################################################################

    def str2datetime(self,isotime):
        try: 
            ymd,time = isotime.strip().split(" ",1)
        except:
            ymd = isotime.strip()
            time = ''
    
        year  = int(ymd[0:4])
        month = int(ymd[5:7])
        day   = int(ymd[8:10])
  
        hms = time.strip().split(":")

        ss_int = 0
        usec   = 0
        mm     = 0
        hh     = 0

        if len(hms) == 3 :
            ss = float(hms[2])
            ss_int = int(ss)
            ss_frac=ss-ss_int
            usec = int(float(ss_frac)*1e6)
            mm = int(hms[1])
            hh = int(hms[0])  
        if len(hms) == 2 :
            mm = int(hms[1])
            hh = int(hms[0])  
        if len(hms) == 1 :
            hh = int(hms[0])  
            
        return dt.datetime(year,month,day,hh,mm,ss_int,usec)


    def base_name(self,oTime):
        """
        base_name

        Aim: 
             

        Input: oTime an object of time datetime.datetime
 
        Output: a list with 
                [0] a string in the format YYYMMDD
                [1] an integer with the time in hundred of microseconds since 0 UTC
                [2] a string in the format YYYY-MM-DD

        Examples:
        >>> import oSST, datetime
        >>> d=oSST.SST()
        >>> date=datetime.datetime(2014,10,11,14,25,33)
        >>> ['1141011', 519330000L, '2014-10-11'

        Change record:
            First written by Guigue on a sunny day of October 2014, in Sao Paulo
            2017-06-17 : Adapted to use datetime python classes (Guigue, Sao Paulo)

        """
        bn  = str(int(oTime.strftime("%Y"))-1900) + oTime.strftime("%m") + oTime.strftime("%d")
        ymd = oTime.strftime("%Y")+ '-' + oTime.strftime("%m") + '-' + oTime.strftime("%d")
        hus =  int((float(oTime.strftime("%H")) * 3600 + \
                    float(oTime.strftime("%M")) * 60   + \
                    float(oTime.strftime("%S")) + \
                    float(oTime.strftime("%f")) / 1e6 ) * 10000 )

        return [bn,hus,ymd]

    def define_fmt(self,sst_file_type,sst_file_date):
        """
        Aim: to return the format of the binary record to be read. 
             This is maybe the most sensitive procedure. The strings that
             represent the binary format depend on the file archive and date.
             Don't change unless you know very well what you ar doing.
             In python language the following applys:
              i : long signed integer , 4 bytes
              h : short signed integer , 2 bytes
              H : short unsigned integer , 2 bytes
              B : unsigned byte , 1 byte
              f : float, 4 bytes

        Input:  
        sst_file_type : scalar string == 'bi' | 'rs' | 'rf'
        sst_file_date : scalar string with an ISO date format of the day

        Output:
        a scalar string with the binary format
        """

        if (sst_file_type == 'rs') or (sst_file_type == 'rf'):
            if (sst_file_date > '2002-12-13'):
                sst_bin_fmt='=iHHHHHHiiihhiihhhhhhhhBBhi'
            elif (sst_file_date > '2002-12-03' and sst_file_date <= '2002-12-13'):
                sst_bin_fmt='=iiiiiiiiiihhiihhhhhhhhBBhi'
            elif (sst_file_date > '1999-05-01' and sst_file_date <= '2002-12-03'):
                sst_bin_fmt='=ihhhhhhhhiiihhiihhhhhhhhBBhhi'
            else: 
                sst_bin_fmt=[0,'']
        elif (sst_file_type == 'bi'):
            if (sst_file_date > '2002-12-13'):
                sst_bin_fmt='=iffffHHHHHHffffffhhBBhhhhhhfffffffffffBi'
            elif (sst_file_date > '2002-11-23') and (sst_file_date <= '2002-12-13'):
                sst_bin_fmt='=iffffhhhhhhffffffhhhBBhhhhhhfffffffffffBi'
            elif (sst_file_date > '2002-09-15') and (sst_file_date <= '2002-11-23'):
                sst_bin_fmt='=iffffhhhhhhhhffffffffhhhBBhhhhhhhhhhhhfffffffffffffffffffffffBi'
            elif (sst_file_date <= '2002-09-15') : 
                sst_bin_fmt='=iffffhhhhhhhhffffffffhhhBBhhhhhhhhhhhhfffffffffffffffffffffffBi'
            else:
                sst_bin_fmt=''
        else:
            sst_bin_fmt=''

        return sst_bin_fmt      


    def read_one_record(self,sst_fd,sst_fmt):
        sst_record=os.read(sst_fd,struct.calcsize(sst_fmt))
        return struct.unpack(sst_fmt,sst_record)


    def unpack_one_record(self,ur,sst_file_date,sst_file_type):

        year  = int(self.initial_time.strftime("%Y"))
        month = int(self.initial_time.strftime("%m"))
        day   = int(self.initial_time.strftime("%d"))

        # The time stamp is in python datetime format
        hours        = ur[0] / 36000000L
        minutes      = (ur[0] % 36000000L)/600000L
        seconds      = ((ur[0] % 36000000L)%600000L)/1.0E+04
        seconds_int  = int(seconds)
        seconds_frac = seconds - int(seconds)
        useconds     = int(seconds_frac * 1e6)

        time_stamp = dt.datetime(year,month,day,hours,minutes,seconds_int,useconds)

        jd = self.am.julian_date(time_stamp.isoformat(" "))
        #####
        #
        #  RS & RF files
        #
        #####
        if (sst_file_type == 'rs') or (sst_file_type == 'rf'):
            if (sst_file_date > '2002-12-13'):
                sst_unpacked_record={'hus_time':ur[0], \
                                     'adc':ur[1:7], \
                                     'ant_temp':[0.0E0,0.0E0,0.0E0,0.0E0,0.0E0,0.0E0], \
                                     'pos_time':ur[7], \
                                     'ant_coord':[ur[8]/1000.0,ur[9]/1000.0], \
                                     'pm': [ur[10]/1000.0,ur[11]/1000.0], \
                                     'ant_coord_err': [ur[12]/1000.0,ur[13]/1000.0], \
                                     'scan_off':  [ur[14]/1000.0,ur[15]/1000.0], \
                                     'rec_offset':ur[16:22] , \
                                     'target':ur[22] % 32, \
                                     'mirror_pos':ur[22] /32 , \
                                     'opmode':ur[23] , \
                                     'julday': jd , \
                                     'time':time_stamp }  
            elif (sst_file_date > '2002-12-03') and (sst_file_date <= '2002-12-13'):
                sst_unpacked_record={'hus_time':ur[0], \
                                     'adc':ur[1:7], \
                                     'ant_temp':[0.0E0,0.0E0,0.0E0,0.0E0,0.0E0,0.0E0], \
                                     'pos_time':ur[7], \
                                     'ant_coord':[ur[8]/1000.0,ur[9]/1000.0], \
                                     'pm': [ur[10]/1000.0,ur[11]/1000.0], \
                                     'ant_coord_err': [ur[12]/1000.0,ur[13]/1000.0], \
                                     'scan_off':  [ur[14]/1000.0,ur[15]/1000.0], \
                                     'rec_offset':ur[16:22] , \
                                     'target':ur[23] % 32, \
                                     'mirror_pos':ur[23] /32 , \
                                     'opmode':ur[24] , \
                                     'julday': jd , \
                                     'time':time_stamp }  
            elif (sst_file_date > '2002-05-20' and sst_file_date <= '2002-12-03'):
                sst_unpacked_record={'hus_time':ur[0], \
                                     'adc':ur[1:9], \
                                     'ant_temp':[0.0E0,0.0E0,0.0E0,0.0E0,0.0E0,0.0E0], \
                                     'pos_time':ur[9], \
                                     'ant_coord':[ur[10]/1000.0,ur[11]/1000.0], \
                                     'pm': [ur[12]/1000.0,ur[13]/1000.0], \
                                     'ant_coord_err': [ur[14]/1000.0,ur[15]/1000.0], \
                                     'scan_off':  [ur[16]/1000.0,ur[17]/1000.0], \
                                     'rec_att':ur[18:24] , \
                                     'target':ur[24] % 32, \
                                     'mirror_pos':ur[24] /32 , \
                                     'opmode':ur[25] , \
                                     'julday': jd , \
                                     'time':time_stamp }  
            elif (sst_file_date > '1999-05-01' and sst_file_date <= '2002-05-20'):
                sst_unpacked_record={'hus_time':ur[0], \
                                     'adc':ur[1:9], \
                                     'ant_temp':[0.0E0,0.0E0,0.0E0,0.0E0,0.0E0,0.0E0,0.0E0,0.0E0] , \
                                     'pos_time':ur[9], \
                                     'ant_coord':[ur[10]/1000.0,ur[11]/1000.0], \
                                     'ant_vel': [ur[12]/1000.0,ur[13]/1000.0], \
                                     'ant_coord_err': [ur[14]/1000.0,ur[15]/1000.0], \
                                     'scan_off':  [ur[16]/1000.0,ur[17]/1000.0], \
                                     'rec_att':ur[18:24] , \
                                     'target':ur[24] % 32, \
                                     'mirror_pos':ur[24] / 32 , \
                                     'opmode':ur[25] , \
                                     'julday': jd , \
                                     'time':time_stamp }  
        ####
        #
        # BI files
        #
        ####
        elif (sst_file_type == 'bi'):
            if (sst_file_date > '2002-12-13'):
                sst_unpacked_record={'hus_time':ur[0], \
                                     'ant_coord':[ur[1],ur[2]], \
                                     'ant_coord_err': [ur[3],ur[4]], \
                                     'adc':ur[5:11], \
                                     'ant_temp':[0.0E0,0.0E0,0.0E0,0.0E0,0.0E0,0.0E0], \
                                     'adc_sigma':ur[11:17], \
                                     'gps_status': ur[17], \
                                     'acq_gain': ur[18] , \
                                     'target':ur[19] % 32, \
                                     'mirror_pos':ur[19] /32 , \
                                     'opmode' : ur[20] , \
                                     'adc_offset': ur[21:27] , \
                                     'hot_temp': ur[27] , \
                                     'amb_temp': ur[28] , \
                                     'opt_temp': ur[29] , \
                                     'if_board_temp': ur[30] , \
                                     'radome_temp': ur[31] , \
                                     'humidity': ur[32] , \
                                     'temperature': ur[33] , \
                                     'opac_210': ur[34] , \
                                     'opac_405': ur[35] , \
                                     'elevation': ur[36] , \
                                     'pressure' : ur[37] , \
                                     'burst' : ur[38] , \
                                     'errors' : ur[39] , \
                                     'julday': jd , \
                                     'time':time_stamp }  
            elif (sst_file_date > '2002-11-23' and sst_file_date <= '2002-12-13'):
                sst_unpacked_record={'hus_time':ur[0], \
                                     'ant_coord':[ur[1],ur[2]], \
                                     'ant_coord_err': [ur[3],ur[4]], \
                                     'adc':ur[5:11], \
                                     'ant_temp':[0.0E0,0.0E0,0.0E0,0.0E0,0.0E0,0.0E0], \
                                     'adc_sigma':ur[11:17], \
                                     'gps_status': ur[17], \
                                     'acq_gain': ur[18] , \
                                     'target':ur[19] % 32, \
                                     'mirror_pos':ur[19] /32 , \
                                     'opmode' : ur[20] , \
                                     'adc_offset': ur[21:27] , \
                                     'hot_temp': ur[27] , \
                                     'amb_temp': ur[28] , \
                                     'opt_temp': ur[29] , \
                                     'if_board_temp': ur[30] , \
                                     'radome_temp': ur[31] , \
                                     'humidity': ur[32] , \
                                     'temperature': ur[33] , \
                                     'opac_210': ur[34] , \
                                     'opac_405': ur[35] , \
                                     'elevation': ur[36] , \
                                     'pressure' : ur[37] , \
                                     'burst' : ur[38] , \
                                     'errors' : ur[39] , \
                                     'julday': jd , \
                                     'time':time_stamp }  
            elif (sst_file_date > '2002-09-15' and sst_file_date <= '2002-11-23'):
                sst_unpacked_record={'hus_time':ur[0], \
                                     'ant_coord':[ur[1],ur[2]], \
                                     'ant_coord_err': [ur[3],ur[4]], \
                                     'adc':ur[5:13], \
                                     'ant_temp':[0.0E0,0.0E0,0.0E0,0.0E0,0.0E0,0.0E0,0.0E0,0.0E0], \
                                     'adc_sigma':ur[13:21], \
                                     'gps_status': ur[21], \
                                     'daq_status': ur[22], \
                                     'acq_gain': ur[23] , \
                                     'target':ur[24] % 32, \
                                     'mirror_pos':ur[24] /32 , \
                                     'opmode' : ur[25] , \
                                     'adc_attenuators': ur[26:32] , \
                                     'adc_offset': ur[32:38] , \
                                     'mix_voltage': ur[38:44] , \
                                     'mix_current': ur[44:50] , \
                                     'hot_temp': ur[50] , \
                                     'amb_temp': ur[51] , \
                                     'opt_temp': ur[52] , \
                                     'if_board_temp': ur[53] , \
                                     'radome_temp': ur[54] , \
                                     'humidity': ur[55] , \
                                     'temperature': ur[56] , \
                                     'opac_210': ur[57] , \
                                     'opac_405': ur[58] , \
                                     'elevation': ur[59] , \
                                     'pressure' : ur[60] , \
                                     'burst' : ur[61] , \
                                     'errors' : ur[62] , \
                                     'julday': jd , \
                                     'time':time_stamp }  
            elif (sst_file_date <= '2002-09-15'):
                sst_unpacked_record={'hus_time':ur[0], \
                                     'ant_coord':[ur[1],ur[2]], \
                                     'ant_coord_err': [ur[3],ur[4]], \
                                     'adc':ur[5:13], \
                                     'ant_temp':[0.0E0,0.0E0,0.0E0,0.0E0,0.0E0,0.0E0,0.0E0,0.0E0], \
                                     'adc_sigma':ur[13:21], \
                                     'gps_status': ur[21], \
                                     'daq_status': ur[22], \
                                     'acq_gain': ur[23] , \
                                     'target':ur[24] % 32, \
                                     'mirror_pos':ur[24] /32 , \
                                     'opmode' : ur[25] , \
                                     'adc_attenuators': ur[26:32] , \
                                     'adc_offset': ur[32:38] , \
                                     'mix_voltage': ur[38:44] , \
                                     'mix_current': ur[44:50] , \
                                     'hot_temp': ur[50] , \
                                     'amb_temp': ur[51] , \
                                     'opt_temp': ur[52] , \
                                     'if_board_temp': ur[53] , \
                                     'radome_temp': ur[54] , \
                                     'humidity': ur[55] , \
                                     'wind': ur[56] , \
                                     'opac_210': ur[57] , \
                                     'opac_405': ur[58] , \
                                     'elevation': ur[59] , \
                                     'seeing' : ur[60] , \
                                     'burst' : ur[61] , \
                                     'errors' : ur[62] , \
                                     'julday': jd , \
                                     'time':time_stamp }  

        return sst_unpacked_record


    def files(self):
        """
        files: It looks the files referent to the interval to be read.
               The list is filled in self.data_files

        use:  self.files()

        Inputs: none

        Outputs: none
              
        Requisites: files() needs
               self.initial_time : a datetime time of the interval beginning.
               self.final_time : a datetime time of the interval end.
               self.data_path : a string where the data files are.
               self.data_type : a string with rs || rf || bi

        Change Record:
               2017-06-17: First Written for the OOP programm - Guigue@Sampa

        """

        # Empty the list before a new search starts
        self.data_files = []

        # Get the base name, i.e. YYYMMDD of SST files
        t1=self.base_name(self.initial_time)
        t2=self.base_name(self.final_time)
        # Base names should be the same. No chance to read data from different days!
        if (t1[0] != t2[0]): return []  

        # Get the hour of the interval
        h1 = int(self.initial_time.strftime("%H"))
        h2 = int(self.final_time.strftime("%H"))
        hours = range(h1,h2+1) 

        # Look for the files
        if (self.data_type.lower() == 'rs'):
            for i in hours : self.data_files.append(self.data_path+'/rs'+t1[0]+'.'+'{0:=02d}'.format(i)+'00')
        elif (self.data_type.lower() == 'rf'):
            for i in hours :
                sst_pattern=self.data_path+'/rf'+t1[0]+'.'+'{0:=02d}'.format(i)+'*' 
                self.data_files = self.data_files + glob.glob(sst_pattern)
        elif (self.data_type.lower() == 'bi'):
            self.data_files = [self.data_path+'/bi'+t1[0]]

        # Retunr the ordered list
        self.data_files.sort()
        return 

    def read(self) :
        """
        read

        Aim: read SST data

        Input : None

        Output: a list of dictionaries. The structure depends on the data type and epoch
                (SST data changed over time....)

        Requisites: files() needs
               self.initial_time : a datetime time of the interval beginning.
               self.final_time : a datetime time of the interval end.
               self.data_path : a string where the data files are.
               self.data_type : a string with rs || rf || bi

        Examples:
        >>> import oSST
        >>> d=oSST.SST()
        >>> d.data_path='/path/to/the/data'
        >>> d.initial_time=d.str2datetime('YYYY-MM-DD HH:MM:SS.SSSSS')
        >>> d.final_time=d.str2datetime('YYYY-MM-DD HH:MM:SS.SSSSSS')
        >>> d.data_type='rf'
        >>> d.read()
        >>> import matplotlib.pyplot as plt
        >>> t=oSST.TimeAxis()
        >>> t.getTimeAxis(d,'dt')
        >>> tp=oSST.TotalPower()
        >>> tp.getTotalPower(d,2)
        >>> import matplotlib.pyplot as plt
        >>> plt.plot(t.time,tp.tp)
        >>> plt.show()
 
        Change record:
         2015-02-19: First written by Guigue@Sampa
         2015-04-01: (fool's day) Minor corrections, Guigue@Sampa
         2016-04-26: Fix the opmode reading in 'rs' and 'rf' files. Guigue@Sampa 
         2016-04-27: Big Change.  Now it reads RF files. Guigue@Sampa
         2017-06-17: First Written for the OOP programm - Guigue@Sampa
         2017-06-18: Enum classes created for Opmodes, targets, etc. 

        """
        self.files()
        self.data = []

        t1=self.base_name(self.initial_time)
        t2=self.base_name(self.final_time)

        for fname in self.data_files :
            if os.path.exists(fname) :
                fd     = os.open(fname,os.O_RDONLY)
                fmt    = self.define_fmt(self.data_type,t1[2])

                if ( len(fmt) < 1 ) : return

                nrec   = os.fstat(fd).st_size / struct.calcsize(fmt)
                for irec in range(nrec) :
                    ur = self.read_one_record(fd,fmt)
                    if (ur[0] >= t1[1] and ur[0] <= t2[1]):
                        sst_record = self.unpack_one_record(ur,t1[2],self.data_type)
                        self.data.append(sst_record)
                os.close(fd)
            else:
                print 'File '+fname+'  not found'

        return

    def __init__(self):
        self.data = []
        self.am = am
        self.initial_time=dt.datetime.now()
        self.final_time=dt.datetime.now()
        self.data_type='rs'
        self.data_path=''
        self.data_files=[]
        return 

####################################################################################################

class yAxis(object):
    """
    yAxis

       This class can be used to retrieve an N-Dimensional Array (ndarray)
       values stored in the list of dictionaries of the SST Data.  

       The SST.data object member is a list of dictionaries, that structures well the information
       but is  hard to manage with, v.g.,  matplotlib.  

       yAxis solves this problem extracting one file from the dictionary and creating an array.
       The class is supplemented with information about the extracted data.

       At the moment, not all of the dictionary fields cane be extracted. The list can be
       increased easily however.

       TotalPower class is a specialized yAxis.

    Input:
       yaxisname = a string with the name of the dictionary field

    Output:
       none

    Example:
       y=oSST.yAxis()
       y.getValues(d,'opmode')

       (where d is an oSST Object)

    Change Record:
       2017-06-19 : First created Guigue@Sampa

    """

    def getValues(self,SST,yaxisname='ant_coord'):

        ndata=len(SST.data)

        if yaxisname == 'ant_coord':
            self.AxisName=['Azimuth','Elevation']
            self.AxisUnits='degrees'
            try:
                self.ant_coord = np.zeros([ndata,2])
                for i in np.arange(ndata) : self.ant_coord[i] = SST.data[i]['ant_coord']
            except:
                self.ant_coord = []

        if yaxisname == 'ant_coord_err':
            self.AxisName=['Azimuth','Elevation']
            self.AxisUnits='degrees'
            try:
                self.ant_coord_err = np.zeros([ndata,2])
                for i in np.arange(ndata) : self.ant_coord_err[i] = SST.data[i]['ant_coord_err']
            except:
                self.ant_coord_err = []
                
        elif (yaxisname == 'adc') :
            self.AxisName =['Ch 1','Ch 2','Ch 3','Ch 4','Ch 5','Ch 6']
            self.AxisUnits = 'ADC Units'
            try:
                self.adc = np.zeros([ndata,6],dtype=np.int)
                for i in np.arange(ndata): self.adc[i] = SST.data[i]['adc'][0:6]
            except:
                self.adc = []

        elif (yaxisname == 'scan_off') and (SST.data_type == 'rs' or SST.data_type == 'rf') :
            self.AxisName= ['Azimuth Offset','Elevation Offset']
            self.AxisUnits='degrees'
            try:
                self.scan_off = np.zeros([ndata,2])
                for i in np.arange(ndata) : self.scan_off[i] = SST.data[i]['scan_off']
            except:
                self.scan_off = []
        elif (yaxisname == 'pnt_corr') and (SST.data_type == 'rs' or SST.data_type =='rf') :
            self.AxisName= ['Delta Azimuth','Delta Elevation']
            self.AxisUnits='degrees'
            try:
                self.pnt_corr = np.zeros([ndata,2])
                for i in np.arange(ndata) : self.pnt_corr[i] = SST.data[i]['pm']
            except:
                self.pnt_corr = []
        elif (yaxisname == 'target') :
            self.AxisName= 'Target'
            self.AxisUnits=''
            try:
                self.target = np.zeros(ndata,dtype=np.int)
                for i in np.arange(ndata) : self.target[i] = SST.data[i]['target']
            except:
                self.target = []
        elif (yaxisname == 'opmode') :
            self.AxisName= 'Opmode'
            self.AxisUnits=''
            try:
                self.opmode = np.zeros(ndata,dtype=np.int)
                for i in np.arange(ndata) : self.opmode[i] = SST.data[i]['opmode']
            except:
                self.opmode = []

        elif (yaxisname == 'mirror_pos') : 
            self.AxisName= 'Calibration Mirror Position'
            self.AxisUnits=''
            try:
                self.mirror_pos = np.zeros(ndata,dtype=np.int)
                for i in np.arange(ndata) : self.mirror_pos[i] = SST.data[i]['mirror_pos']
            except:
                self.mirror_pos = []
        elif (yaxisname == 'adc_sigma') and (SST.data_type == 'bi'):
            self.AxisName =['sigma Ch 1','sigma Ch 2','sigma Ch 3','sigma Ch 4','sigma Ch 5','sigma Ch 6']
            self.AxisUnits = 'ADC Units'
            try:
                self.adc_sigma = np.zeros([ndata,6],dtype=np.int)
                for i in np.arange(ndata): self.adc_sigma[i] = SST.data[i]['adc_sigma'][0:6]
            except:
                self.adc_sigma=[]
        elif (yaxisname == 'mix_voltage') and (SST.data_type == 'bi'):
            self.AxisName = ['Ch 1','Ch 2','Ch 3','Ch 4','Ch 5','Ch 6']
            self.AxisUnits = 'V'
            try:
                self.mix_voltage = np.zeros([ndata,6],dtype=np.int)
                for i in np.arange(ndata): self.mix_voltage[i] = SST.data[i]['mix_voltage']
            except:
                self.mix_voltage = []
        elif (yaxisname == 'mix_current') and (SST.data_type == 'bi'):
            self.AxisName = ['Ch 1','Ch 2','Ch 3','Ch 4','Ch 5','Ch 6']
            self.AxisUnits = 'mA'
            try:
                self.mix_current = np.zeros([ndata,6],dtype=np.int)
                for i in np.arange(ndata): self.mix_current[i] = SST.data[i]['mix_current']
            except:
                self.mix_current = []
        elif (yaxisname == 'hot_temp') and (SST.data_type == 'bi'):
            self.AxisName = 'Hot Load Temperature'
            self.AxisUnits = 'C'
            try:
                self.hot_temp = np.zeros(ndata,dtype=np.int)
                for i in np.arange(ndata): self.hot_temp[i] = SST.data[i]['hot_temp']
            except:
                self.hot_temp = []
        elif (yaxisname == 'amb_temp') and (SST.data_type == 'bi'):
            self.AxisName = 'Ambient Load Temperature'
            self.AxisUnits = 'C'
            try:
                self.amb_temp = np.zeros(ndata,dtype=np.int)
                for i in np.arange(ndata): self.amb_temp[i] = SST.data[i]['hot_temp']
            except:
                self.amb_data = []
        elif (yaxisname == 'opt_temp') and (SST.data_type == 'bi'):
            self.AxisName = 'Optical Box Temperature'
            self.AxisUnits = 'C'
            try:
                self.opt_temp = np.zeros(ndata,dtype=np.int)
                for i in np.arange(ndata): self.opt_temp[i] = SST.data[i]['opt_temp']
            except:
                self.opt_temp = []
        elif (yaxisname == 'if_temp') and (SST.data_type == 'bi'):
            self.AxisName = 'IF Board Temperature'
            self.AxisUnits = 'C'
            try:
                self.if_temp = np.zeros(ndata,dtype=np.int)
                for i in np.arange(ndata): self.if_temp[i] = SST.data[i]['if_board_temp']
            except:
                self.if_temp = []
        elif (yaxisname == 'radome_temp') and (SST.data_type == 'bi'):
            self.AxisName = 'Radome Temperature'
            self.AxisUnits = 'C'
            try:
                self.radome_temp = np.zeros(ndata,dtype=np.int)
                for i in np.arange(ndata): self.radome_temp[i] = SST.data[i]['radome_temp']
            except:
                self.radome_temp = []
        elif (yaxisname == 'temperature') and (SST.data_type == 'bi'):
            self.AxisName = 'Atmospheric Temperature'
            self.AxisUnits = 'C'
            try:
                self.temperature = np.zeros(ndata,dtype=np.int)
                for i in np.arange(ndata): self.temperature[i] = SST.data[i]['temperature']
            except:
                self.temperature = []
        elif (yaxisname == 'humidity') and (SST.data_type == 'bi'):
            self.AxisName = 'Atmospheric Humidity'
            self.AxisUnits = '%'
            try:
                self.humidity = np.zeros(ndata,dtype=np.int)
                for i in np.arange(ndata): self.humidity[i] = SST.data[i]['humidity']
            except:
                self.humidity = []
        elif (yaxisname == 'pressure') and (SST.data_type == 'bi'):
            self.AxisName = 'Atmospheric Pressure'
            self.AxisUnits = 'mmHg'
            try:
                self.pressure = np.zeros(ndata,dtype=np.int)
                for i in np.arange(ndata): self.pressure[i] = SST.data[i]['pressure']
            except:
                self.pressure = []

        return

    def __init__(self):
        self.AxisName = ['Azimuth','Elevation']
        self.AxisUnits= 'degrees'
        return

##############################################################################################

class TotalPower(object):

    def getTotalPower(self,SST,ch=1,adc='y'):

        """
        TotalPower
           Extracts either ADC values or Antenna Temperatures
           for a determined channel. 

        Inputs:
          ch: integer number 1...6 representing the receiver data to be read
          adc: string either 'y' (ADC values) or 'n' (Antenna Temperature)
               default is 'y'

        Output:
          A numpy ndarray

        Change Record:
          2017-06-15: First written by Guigue@Sampa

        """

        ndata=len(SST.data)
        if (ch < 1 or ch > 6) :
            print 'Wrong channel'
            tp = []
        else:
            tp = np.empty(ndata)
            if (adc == 'y') : 
                for i in np.arange(ndata) : tp[i]=SST.data[i]['adc'][ch]
            elif (adc == 'n') :
                for i in np.arange(ndata) : tp[i]=SST.data[i]['ant_temp'][ch]
            else:
                print 'Wrong choice adc=y||n'
                tp=[]

        self.initial_time = SST.initial_time 
        self.final_time = SST.final_time
        self.data_type = SST.data_type
        self.channel   = ch
        self.adc = adc
        self.tp = tp
        return

    def __init__(self):
        self.tp = np.zeros(0)
        self.am = am
        self.initial_time=dt.datetime.now()
        self.final_time=dt.datetime.now()
        self.data_type='rs'
        self.channel = 1
        self.adc = 'y'
        return 

####################################################################################################

class TimeAxis(object):

    def getTimeAxis(self,SST,ttype='ms'):
        """
        timeaxis

            It returns the time acis of SST data in three different formats: 
               1) milliseconds (ms) from the beginning of the day
               2) julian day (jd)
               3) ISO format (iso) : YYYY-MM-DD HH:MM:SS.SSSS
            Default is milieseconds

        Input: 
            ttype : string either 'ms', 'dt' or 'jd'

        Output:
            a ndarray (numpy) with the time axis

        Change Record:
            2017-06-15: First written by Guigue@Sampa
        """

        ndata=len(SST.data)
        if ttype == 'ms' :
            time = np.array(np.empty(len(SST.data)), dtype=np.uint32)
            for i in np.arange(ndata) : time[i]=SST.data[i]['hus_time']
        elif ttype == 'jd':
            time = np.array(np.empty(len(SST.data)), dtype=np.float64)
            for i in np.arange(ndata) : time[i]=SST.data[i]['julday']
        elif ttype == 'dt':
            #time = np.array(np.empty(len(SST.data)), dtype=np.string_)
            time = []
            for i in np.arange(ndata) : time.append(SST.data[i]['time'])
        else :
            print 'Unknown type format.'
            time = []

        self.initial_time = SST.initial_time 
        self.final_time = SST.final_time
        self.data_type = SST.data_type
        self.time = time 
                 
        return

    def __init__(self):
        self.time = np.zeros(0)
        self.am = am
        self.initial_time=dt.datetime.now()
        self.final_time=dt.datetime.now()
        self.data_type='rs'
        return 

######################################################################################

class Beams(IntEnum):
    """ 
    Beams

        An enum class to define the beams. We normally use numbers from 1 to 6
        while the array index goes from 0 to 5. 

    """
    one   = 0
    two   = 1
    three = 2
    four  = 3
    five  = 4
    six   = 5

class Opmodes(IntEnum):
    """
    Opmodes

      An Enum class with the Observing Modes of the SST
      Obtained from SST_OPMODE.pro

    """

    track       =  0
    map_radec   =  1
    map_azel    =  2
    map_radial1 =  3
    map_interm  =  4
    scan_az     =  5
    scan_el     =  6
    scan_ra     =  7
    scan_dec    =  8
    scan_interm =  9
    scan_tau    = 10
    scan_radec  = 11
    scan_azel   = 12
    fast_8scan  = 14
    onoff       = 15
    onon        = 16
    offpoint    = 17
    on_interm   = 18
    map_decra   = 21
    map_elaz    = 22
    map_radial2 = 23
    max_att     = 40
    stall       = 50
    ant_locked  = 55
    unknown_mode= 99

class MirrorPos(IntEnum):

    """
    MirrorPos

        An Enum class with the Calibration Mirror positions.
        antenna means: observing the sky (astronomical)
        amb_load, hot_load: the calibration loads
        cal_ref: I don't remember... 
        unknown : the mirror is moving

        Taken from SST_TARGET.pro

    """     
    antenna     = 0
    amb_load    = 1
    hot_load    = 2
    cal_ref     = 3
    unknown_cal = 7

class Target(IntEnum):
    """
    Target

       An Enum class with the target codification, created when Pluto
       was still a planet...

       Taken from SST_TARGET.pro

    """

    sky         =  0
    mercury     =  1
    venus       =  2
    earth       =  3
    mars        =  4
    jupiter     =  5
    saturn      =  6
    uranus      =  7
    neptune     =  8
    pluto       =  9
    moon        = 10
    sun         = 11
    ar          = 12
    star        = 13
    beacon      = 20
    service     = 21
    manual      = 22
    last_obj    = 30
    unknown_obj = 31
    

