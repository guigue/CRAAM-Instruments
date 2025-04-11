import os
import xml.etree.ElementTree as xmlet
import numpy as np
from astropy.io import fits
import collections
import warnings
import pdb
import datetime as dt
import pickle
from scipy.optimize import curve_fit
from astropy import units as u
from astropy import constants as c

############ Global Variables ##########
__version__       = "2025-04-11T1610ART"
__DATA_FILE__     = "hats_data_rbd.bin"
__HUSEC_FILE__    = "hats_husec.bin"
__RECORD_SIZE__   = 38
short_array       = collections.deque()
#######################################


###########################################################################################################################################
#
# HATS: A python class to read and deconvolve HATS RBD and AUX data.
#
# Usage:
#   > export HATSXMLTABLES=/my/directory/with/HATS/XMLTABLES
#   >>> import HATS
#
#   ### RBD data
#   >>> h=HATS.hats('2021-08-26 1800')      # Create the Object for ISO date 2021-08-26 16:00 UT 

#
#
##########################################################################################################################################
#
# Data Structures:
#
#   HATS
#   MetaData:   Dictionary with
#                               InputPath : string with directory where files are
#                               XMLDataDescriptionFile : XML file with RBD data description
#                               File_Name : string with the RBD name file
#                               Date      : ISO data date
#                               Hour      : Data Hour (UTC)
#                               FFTProgram: full path to FFT module (C)
#   rData: A numpy array with raw data
#                               sample    : sample number (np.uint32 )
#                               sec       : seconds from 1970-01-01, Unix time (np.uint32)
#                               ms        : milliseconds of the seconds (np.uint16)
#                               husec     : hundred of microsends since 0 UT, as in SST (np.uint64)
#                               golay     : Golay output in ADC units (np.uint32)
#                               chopper   : Power Supply in ADC units (np.uint32)
#                               temp_env
#                               temp_hics
#                               temp_golay
#
#   cData: A numpy array with converted data
#                               golay     : Golay output in ADC units (np.uint32)
#                               chopper   : Power Supply in ADC units (np.uint32)
#                               temp_env
#                               temp_hics
#                               temp_golay
#   Deconv: a Dictionary with
#                               amplitude : FFT amplitude in ADC units (np.float64)
#                               husec     : time of 'amplitude' in hundred of microsends since 0 UT, as in SST (np.uint64)
#                               time      : time of 'amplitude' in datetime()
#
#
#
#   AUX
#   MetaData:   Dictionary with
#                               InputPath             : string with directory where files are
#                               XMLAuxDescriptionFile : string with XML data description file
#                               File_Name             : string with the AUX name file
#                               Date                  : ISO data date
#                               Hour                  : Data Hour (UTC)
#   Data: a numpy array with auxliary data
#                               husec                 : hundred of microseconds (0.1 ms) since 0 UT
#                               jd                    : julian day
#                               sid                   : sidereal time
#                               elevation             : elevation (degrees)
#                               azimuth               : azimuth (degrees)
#                               right_ascension       : right ascension (degrees)
#                               declination           : declination (degrees)
#                               ra_rate               : right ascension speed (degrees / second)
#                               dec_rate              : declination speed (degrees / second)
#                               object                : pointed object ID (integer number)
#                               opmode                : operation mode ID (integer number)
#   Time: an array with time in datetime format 
#
#
#
#   ENV
#   Data: a numpy array with environmental data
#                                husec                : hundred of microseconds (0.1 ms) since 0 UT
#                                time                 : time in datetime format
#                                Tgolay               : Golay cell temperature
#                                Thics                : HICS internal temperature
#                                Tenv                 : External temperature
#
#   WS
#   
#
####################################################################################################################################
#
#   METHODS
#
#   getTimeAxis(husec,MetaData)  : returns a datetime ndarray
#                                  Examples : taux=h.getTimeAxis(h.aux.Data['husec'],h.MetaData)  # for auxiliary data
#                                             trbd=h.getTimeAxis(h.rbd.rData['husec'],h.MetaData) # for rbd data
#
#   check()                      : checks data integrity
#   getFFT()                     : Deconvolves RBD data, creating Deconv structure
#   extract(tr,save=False,pklname={file_name}): input a time range tr=[datetime(initial),datetime(end)] and returns a HATS object for the time range.
#   extrac_scans(stype='rigt_ascension' | 'declination')
#                                : Extract scans from a hats class object. Returns a list of dictionaries with the scans.
#   __add__                      : h=HATS.hats('2021-12-13T1500')
#                                  g=HATS.hats('2021-12-13T1600')
#                                  i=h+g
#   plot()                       : Simple plot of data
#                                  Returns the time axis
#   
####################################################################################################################################
#   
# Author: @guiguesp - São Paulo - returning to presencial classes - 2021-10-12
#                     Added more descriptions in XML file
#                     Simplified the readout of the Analog Devices variables
#                     Added the '+' method
#                     Added extract() method
#                     Many changes - 2023-02-01 From OAFA/Cesco
#                           - Deconv is now a structured ndarray
#                           - Fixed __add__ (problems with concatenaion)
#                           - Deletion of wrong records depends on numpy version.
#                           - Added env class that integrates environmental data in 5-minute periods
#                           - aux structure changed, now includes the time
#                           - extract_scans()
#                     2023-03-31 - Sampa
#                           - extract() corrected some bugs derived from the changes above.
#                     2023-04-14 - Sampa
#                           - seqLims() : corrected bug when only one sequence was present
#                           - smooth()  : updated to Deconv as an ndarray (before was a dictionary)
#                     2023-04-18 - Sampa
#                           _ __add__() : Fix the MetaData dictionary of the sum
#                     2023-05-29 - Sampa
#                            - np.int replaced by np.int32
#                     2024-09-11 - OAFA
#                            - added save a pickle file when extracted
#                     2024-10-10 - Sampa
#                            - weather station data included here
#                            - sky dip method included
#                     2025-04-08 - OAFA
#                            - Class ws computes the Precipitable Water Vapor Content.
#
####################################################################################################################################


########### Global Functions ##########
def contiguo(index):

    nu = len(index)
    if (nu < 2):
        return np.zeros([2,1],int)
    
    n = nu-1
    ad = 0
    
    dif = index[1:] - index[:-1]

#    pdb.set_trace()
    
    if ((dif[0] != 1) & (dif[1] == 1)):
        index = index[1:]
        dif  = index[1:] - index[:-1]
        ad   = 1
        nu   = nu - 1
        n    = n - 1
    
    if ((dif[n-1] != 1) & (dif[n-2] == 1)):
        index = index[0:n-1]
        dif   = index[1:] -index[:-1]
        nu    = nu-1
        n     = n-1
    
    tt, = np.where(dif == 1)
#    pdb.set_trace()    
    if (len(tt) == 0): 
        return np.zero([2,1],int)
        
    if (len(tt) == nu-1):
        pp = np.zeros([2,1],int)
        pp[0,0] = 0
        pp[1,0] = nu-1
        return pp+ad
        
    k, = np.where(dif != 1)
    dif[k] = 0
    dif2 = dif[1:]-dif[:-1]
        
    pi, = np.where(dif2 == -1)
    pf, = np.where(dif2 == 1)
        
    if (pi[0] < pf[0]):
        tmp = np.zeros(len(pi)+1,int)
        for i in range(len(pi)):
            tmp[i+1] = pi[i]
        pi = tmp
        
    if (pi[len(pi)-1] < pf[len(pf)-1]):
        tmp = np.zeros(len(pf)+1,int)
        for i in range(len(pf)):
            tmp[i] = pf[i]
        tmp[len(pf)] = nu-1
        pf = tmp        
            
    pi[1:]  = np.array(pi[1:]) + 2 
    pp      = np.zeros([2,len(pi)],int)
    pp[0,:] = pi
    pp[1,:] = pf
        
    return (pp+ad).T

##############################################


class hats(object):
    
    def __str__(self):
        return 'A Class representing HATS Data'
    
    def __init__(self,date='',PathToXML='',pklname=''):

        warnings.filterwarnings('ignore')

        if (len(pklname) > 0):
            try:
                with open(pklname,'rb') as file:
                    h=pickle.load(file)
                    self.MetaData = h.MetaData
                    self.aux = h.aux
                    self.env = h.env
                    self.rbd = h.rbd
                    return 
            except:
                print('\n\n File {0:s} not found\n\n'.format(pklname))

        # PathToXML should point to tha directory where the XML tables are copied
        # When not defined, look at the environment
        if (isinstance(PathToXML,str) and len(PathToXML) == 0):
            if ('HATSXMLPATH' in os.environ.keys()):
                PathToXML=os.environ['HATSXMLPATH']
            else:
                PathToXML='./'
        else:
            PathToXML=PathToXML

        if PathToXML[-1] != '/' :
            PathToXML=PathToXML+'/'
        
        if (os.path.exists(PathToXML)):
            
            if os.path.exists(PathToXML+"HATSAuxFormat.xml"):
                self.MetaData = {"XMLAuxDescriptionFile":PathToXML+"HATSAuxFormat.xml"}
            else:
                raise ValueError("No XML file found {}".format(PathToXML+"HATSAuxFormat.xml"))
                return
            
            if os.path.exists(PathToXML+"HATSDataFormat.xml"):            
                self.MetaData.update({"XMLDataDescriptionFile":PathToXML+"HATSDataFormat.xml"})
            else:
                raise ValueError("\n\n No XML file found {}\n\n".format(PathToXML+"HATSDataFormat.xml"))
                return
            
        else:
            raise ValueError("\n\n No XML directory {0:s} found \nExit\n\n".format(PathToXML))

        InputPath = os.getenv('HATS_DATA_InputPath')
        if InputPath == None :
            InputPath = './'
            
        if (InputPath[-1] != '/') :
            InputPath+='/'        
        if (os.path.exists(InputPath)):
            self.MetaData.update({'InputPath':InputPath})
        else:
            raise ValueError("\n\n No Input Path {} found \n\n".format(InputPath))

        self.rbd = rbd()
        self.aux = aux()
        self.env = env()
        
        if (len(date) >= 13):
            YMD = date[:10]
            H   = date[11:13]
            
            self.MetaData.update({'Date':YMD})
            self.MetaData.update({'Hour':H+'00'})

            fname = 'hats-'+YMD+'T'+H+'00.rbd'    # RBD file
            self.from_file(fname)

            fname = 'hats-'+YMD+'T'+H+'00.aux'    # AUX file
            self.from_file(fname)

            self.env.integrate(self.rbd.rData['husec'],self.rbd.cData)
            self.env.Data['time'] = self.getTimeAxis(self.env.Data['husec'],self.MetaData)

        return

    def from_file(self,fname):

        if (fname[-3:] == 'aux') & (fname[:4] == 'hats') :
            fullpathname = self.MetaData['InputPath']+'aux/'+fname
            if (not os.path.exists(str(fullpathname))):
                raise ValueError('\n\n File {0:s} not found \n\n'.format(fullpathname))
            Data = self.aux.from_file(fullpathname,self.MetaData)
            Time = self.getTimeAxis(Data['husec'],self.MetaData)
            eData = np.empty(Data.shape, dtype=Data.dtype.descr + ([('time', object)]))
            for field in Data.dtype.names:
                eData[field] = Data[field]
            eData['time'] = Time
            self.aux.Data = eData
            self.MetaData.update({'AUXFilename':fullpathname})

        elif (fname[-3:] == 'rbd') & (fname[:4] == 'hats'):
            fullpathname = self.MetaData['InputPath']+fname
            if (not os.path.exists(str(fullpathname))):
                raise ValueError('\n\n File {0:s} not found \n\n'.format(fullpathname))
            self.rbd.from_file(fullpathname,self.MetaData)
            self.getFFT()
            self.MetaData.update({'RBDFilename':fullpathname})

        else:
            print("\n\n  {0:s} : wrong file \n\n".format(fname))
            
        return
    
    def getTimeAxis(self,husec,MetaData):
        """

        getTimeAxis: Class method to convert the us time axis used in RBD files to a Python
                     datetime ndarray that can be used with matplotlib.pyplot.

        Change Record:  First written by Guigue @ Sampa
                        2017-11-04 St Charles Day !!!

        """

        N = husec.shape[0]
        if isinstance(MetaData['Date'], list):
            year  = np.zeros(N,dtype=np.int32) + int(MetaData['Date'][0][0:4])
            month = np.zeros(N,dtype=np.int32) + int(MetaData['Date'][0][5:7])
            day   = np.zeros(N,dtype=np.int32) + int(MetaData['Date'][0][8:10])
        else:
            year  = np.zeros(N,dtype=np.int32) + int(MetaData['Date'][0:4])
            month = np.zeros(N,dtype=np.int32) + int(MetaData['Date'][5:7])
            day   = np.zeros(N,dtype=np.int32) + int(MetaData['Date'][8:10])
            
        dt_time = list(map(self.husec2dt,husec,year,month,day))
        return np.asarray(dt_time)

    def husec2dt(self,husec,year,month,day):
        import datetime as dt

        ms           = husec
        hours        = int(ms // 36000000)
        minutes      = int((ms % 36000000) // 600000)
        seconds      = ((ms % 36000000) % 600000) / 1.0E+04
        seconds_int  = int(seconds)
        seconds_frac = seconds - int(seconds)
        useconds     = int(seconds_frac * 1e6)
        return dt.datetime(year,month,day,hours,minutes,seconds_int,useconds)

    def getFFT(self,steps=32,window_size=128,recoffset=0,nrecords=-1,calibrated=True):

        if (recoffset < 0):
            return
        
        Ntotal = self.rbd.rData['husec'].shape[0]
        if (nrecords == -1):
            s = slice(recoffset,Ntotal)
        elif (nrecords > -1):
            if (recoffset+nrecords > Ntotal):
                return
            s = slice(recoffset,recoffset+nrecords)

        self.rbd.rData['husec'][s].tofile(__HUSEC_FILE__)
        if calibrated:
            self.rbd.cData['golay'][s].tofile(__DATA_FILE__)
        else:
            self.rbd.rData['golay'][s].astype(np.float64).tofile(__DATA_FILE__)
            
        if os.path.exists(self.MetaData['FFTProgram']):
            os.spawnl(os.P_WAIT,self.MetaData['FFTProgram'],os.path.basename(self.MetaData['FFTProgram']),str(window_size),str(steps))
        else:
            print("\n\n No {0:s}  module found in {1:s}\nExit\n\n".format(os.path.basename(self.MetaData['FFTProgram']),
                                                                          os.path.dirname(self.MetaData['FFTProgram'])))
            return            
        
        husec  = np.fromfile(__HUSEC_FILE__ , dtype=np.uint64  , count=-1)
        fftAmp = np.fromfile(__DATA_FILE__  , dtype=np.float64 , count=-1)

        os.remove(__HUSEC_FILE__)
        os.remove(__DATA_FILE__)

        t = self.getTimeAxis(husec,self.MetaData)
        NRecords = husec.shape[0]
        self.rbd.Deconv = np.zeros(NRecords,dtype=[('time',object),('husec','uint64'),('amplitude','float')])
        self.rbd.Deconv['time']      = t
        self.rbd.Deconv['husec']     = husec
        self.rbd.Deconv['amplitude'] = fftAmp

        return
    
    def extract(self,time_range,save=False,pklname='hats_extract.pkl'):

        ######################################
        #
        # If save=True it saves a pickle file
        # the filename is given with pklname option in the command line
        #
        # Others file formats are on the way...
        #     @guiguesp 2024-09-11 (at OAFA after surviving a 24-hour Zonda wind)
        ##########################################################################
        
        _temp_ = hats()
        _temp_.MetaData.update({'XMLRBDDescriptionFile':self.MetaData['XMLDataDescriptionFile']})
        _temp_.MetaData.update({'XMLAuxDescriptionFile':self.MetaData['XMLAuxDescriptionFile']})
        _temp_.MetaData.update({'InoutPath':self.MetaData['InputPath']})
        _temp_.MetaData.update({'RBDFilename':self.MetaData['RBDFilename']})
        _temp_.MetaData.update({'Date': self.MetaData['Date']})
        _temp_.MetaData.update({'TimeRange':[time_range[0].isoformat(),time_range[1].isoformat()]})

        if isinstance(self.MetaData['Date'],list):
            ts0    = dt.datetime(int(self.MetaData['Date'][0][:4]),int(self.MetaData['Date'][0][5:7]),int(self.MetaData['Date'][0][8:]),0,0,0).timestamp()           
        else:
            ts0    = dt.datetime(int(self.MetaData['Date'][:4]),int(self.MetaData['Date'][5:7]),int(self.MetaData['Date'][8:]),0,0,0).timestamp()

            
        time_range_husec = [(time_range[0].timestamp() - ts0) * 10000, (time_range[1].timestamp() - ts0) * 10000]

        X = (self.rbd.rData['husec'] >= time_range_husec[0] ) & (self.rbd.rData['husec'] <= time_range_husec[1])
        _rData_  = np.copy(self.rbd.rData[X])
        _cData_  = np.copy(self.rbd.cData[X])

        X = (self.env.Data['husec'] >= time_range_husec[0] ) & (self.env.Data['husec'] <= time_range_husec[1])
        _eData_  = np.copy(self.env.Data[X])

        X = (self.rbd.Deconv['husec'] >= time_range_husec[0] ) & (self.rbd.Deconv['husec'] <= time_range_husec[1])
        _Deconv_ = np.copy(self.rbd.Deconv[X])
        
        X = (self.aux.Data['husec'] >= time_range_husec[0] ) & (self.aux.Data['husec'] <= time_range_husec[1])
        _Data_ = self.aux.Data[X]

        _temp_.rbd.rData  = _rData_
        _temp_.rbd.cData  = _cData_
        _temp_.rbd.Deconv = _Deconv_
        _temp_.env.Data   = _eData_
        _temp_.aux.Data   = _Data_

        if save:
            with open(pklname,'wb') as file:
                pickle.dump(_temp_,file,pickle.HIGHEST_PROTOCOL)                                            

        return _temp_
    
    def __add__(self,h):

        _temp_ = hats()

        # Build the new MetaData Dictionary
        
        if isinstance(self.MetaData['Date'],list):
            Date = self.MetaData['Date']
        else:
            Date = [self.MetaData['Date']]
        if isinstance(h.MetaData['Date'],list):
            Date = Date + h.MetaData['Hour']
        else:
            Date.append(h.MetaData['Date'])

        if isinstance(self.MetaData['Hour'],list):    
            Hour = self.MetaData['Hour']
        else:
            Hour = [self.MetaData['Hour']]
        if isinstance(h.MetaData['Hour'],list):    
            Hour = Hour + h.MetaData['Hour']
        else:
            Hour.append(h.MetaData['Hour'])
            
        if isinstance(self.MetaData['RBDFilename'],list):    
            FNam = self.MetaData['RBDFilename']
        else:
            FNam = [self.MetaData['RBDFilename']]
        if isinstance(h.MetaData['RBDFilename'],list):    
            FNam = Hour + h.MetaData['RBDFilename']
        else:
            FNam.append(h.MetaData['RBDFilename'])
        
        if isinstance(self.MetaData['XMLDataDescriptionFile'],list):    
            XMLRBD = self.MetaData['XMLDataDescriptionFile']
        else:
            XMLRBD = [self.MetaData['XMLDataDescriptionFile']]
        if isinstance(h.MetaData['XMLDataDescriptionFile'],list):    
            XMLRBD = Hour + h.MetaData['XMLDataDescriptionFile']
        else:
            XMLRBD.append(h.MetaData['XMLDataDescriptionFile'])

        if isinstance(self.MetaData['XMLAuxDescriptionFile'],list):    
            XMLAUX = self.MetaData['XMLAuxDescriptionFile']
        else:
            XMLAUX = [self.MetaData['XMLAuxDescriptionFile']]
        if isinstance(h.MetaData['XMLAuxDescriptionFile'],list):    
            XMLAUX = Hour + h.MetaData['XMLAuxDescriptionFile']
        else:
            XMLAUX.append(h.MetaData['XMLAuxDescriptionFile'])
            
        if isinstance(self.MetaData['FFTProgram'],list):    
            FFT = self.MetaData['FFTProgram']
        else:
            FFT = [self.MetaData['FFTProgram']]
        if isinstance(h.MetaData['FFTProgram'],list):    
            FFT = Hour + h.MetaData['FFTProgram']
        else:
            FFT.append(h.MetaData['FFTProgram'])

        if isinstance(self.MetaData['InputPath'],list):    
            PATH = self.MetaData['InputPath']
        else:
            PATH = [self.MetaData['InputPath']]
        if isinstance(h.MetaData['InputPath'],list):    
            PATH = Hour + h.MetaData['InputPath']
        else:
            PATH.append(h.MetaData['InputPath'])
            
        _temp_.MetaData = {'XMLDataDescriptionFile': XMLRBD,
                           'XMLAuxDescriptionFile': XMLAUX,
                           'InputPath': PATH,
                           'RBDFilename': FNam,
                           'Date': Date,
                           'Hour': Hour,
                           'FFTProgram': FFT}

        _temp_.rbd.rData = np.concatenate((self.rbd.rData,h.rbd.rData))
        
        _temp_.rbd.cData = {}        
        _temp_.rbd.cData = np.concatenate((self.rbd.cData,h.rbd.cData))

        _temp_.env.Data = {}        
        _temp_.env.Data = np.concatenate((self.env.Data,h.env.Data))
        _temp_.rbd.Deconv = np.concatenate((self.rbd.Deconv,h.rbd.Deconv))
        _temp_.aux.Data = np.concatenate((self.aux.Data,h.aux.Data))
                        
        return _temp_
    
    def slice_mean(self,item):
        global short_array
        old = short_array.popleft()            # pop oldest from left
        short_array.append(item)               # push newest in from right
        try :
            m = np.mean(short_array)
            return m
        except:
            print(short_array)
            return 0

    def smooth(self,M,truncate=False,wrap=False,zero=False,mirror=False):
        
    ############################################################################
    # smooth: Finds the mean for the points in a sliding window (fixed size)
    #         as it is moved from left to right by one point at a time.
    # Inputs:
    #         M: number of items in sliding window
    # Options: The options follow the IDL 8.7 smooth options. Only one is valid.
    #         Boolean variables:
    #            truncate
    #            wrap
    #            zero
    #            mirror
    #  Otput:
    #      self.rbd.Deconv['smoothed']: ndarray the same size of y
    #
    #  Author: Adapted from a public code.
    #        Guiguesp @ Sao Paulo - 2018-03-01
    #        Introduced the map() function to speed the process
    #############################################################################

        y = np.copy(self.rbd.Deconv['amplitude'])
        if not (M % 2) :
            M+=1        # must be odd

        N = y.shape[0]
        edge=True

        M2 = M//2

        if wrap:
            y=np.concatenate((y[-M2+1:],y,y[0:M2]))
        elif zero:
            y=np.concatenate((np.zeros(M2),y,np.zeros(M2)))
        elif mirror:
            y=np.concatenate((y[M2:0:-1],y,y[N-2:N-2-M2:-1]))
        elif truncate:
            y=np.concatenate((np.zeros(M2)+y[0],y,np.zeros(M2)+y[-1]))
        else:
            edge=False

        N = y.shape[0]

        # Load deque (d) with first window of seq
        global short_array
        short_array = collections.deque(y[0:M])
        means = [np.mean(short_array)]             # contains mean of first window

        meansF = list((map(self.slice_mean,y[M:])))
        means  = means + meansF
        means  = np.asarray(means)

        if (not edge) :
            means = np.concatenate((y[0:M2+1],means,y[-M2+1:]))

        NRecords = means.shape[0]
        _Deconv_              = np.zeros(NRecords,dtype=[('time',object),('husec','uint64'),('amplitude','float'),('smooth','float')])
        _Deconv_['time']      = self.rbd.Deconv['time']
        _Deconv_['husec']     = self.rbd.Deconv['husec']
        _Deconv_['amplitude'] = self.rbd.Deconv['amplitude']
        _Deconv_['smooth']    = means

        self.rbd.Deconv = _Deconv_
        
        if 'Smooth_Window' in self.rbd.MetaData:
            self.rbd.MetaData['Smooth_Window'] = M
        else:
            self.rbd.MetaData.update({'Smooth_Window':M})

        return                    # return an ndarray
    
    def seqLims(self,array,condition):

        index, = np.where(array==condition)
        nu = index.shape[0]
        if (nu < 2):
            return np.zeros([2,1],int)
    
        n = nu-1
        ad = 0
    
        dif = index[1:] - index[:-1]

        if ((dif[0] != 1) & (dif[1] == 1)):
            index = index[1:]
            dif  = index[1:] - index[:-1]
            ad   = 1
            nu   = nu - 1
            n    = n - 1
    
        if ((dif[n-1] != 1) & (dif[n-2] == 1)):
            index = index[:n-1]
            dif   = index[1:] -index[:-1]
            nu    = nu-1
            n     = n-1
    
        tt, = np.where(dif == 1)

        if (tt.shape[0] == 0): 
            return np.zero([2,1],int)
        
        if (len(tt) == nu-1):
            pp = np.zeros([1,2],int)
            pp[0,0] = 0
            pp[0,1] = nu-1
            return pp
        
        k, = np.where(dif != 1)
        dif[k] = 0
        dif2 = dif[1:]-dif[:-1]
        
        pi, = np.where(dif2 == -1)
        pf, = np.where(dif2 == 1)
        
        if (pi[0] < pf[0]):
            tmp = np.zeros(pi.shape[0]+1,int)
            for i in range(pi.shape[0]):
                tmp[i+1] = pi[i]
            pi = tmp
        
        if (pi[len(pi)-1] < pf[len(pf)-1]):
            tmp = np.zeros(pf.shape[0]+1,int)
            for i in range(pf.shape[0]):
                tmp[i] = pf[i]
            tmp[len(pf)] = nu-1
            pf = tmp        
            
        pi[1:]  = np.array(pi[1:]) + 2 
        pp      = np.zeros([2,pi.shape[0]],int)
        pp[0,:] = pi
        pp[1,:] = pf
        
        return (pp+ad).T

    def extract_scans(self,stype='right_ascension'):

        if stype=='right_ascension':
            opmode = 8
        elif stype=='declination':
            opmode = 7
        else:
            opmode = 8

        wscan, = np.where(self.aux.Data['opmode']==opmode)
        icon  = self.seqLims(self.aux.Data['opmode'],opmode)

        NScans = icon.shape[0]
        
        scans = {}
        for iscan in np.arange(0,NScans):
            itr = (self.rbd.Deconv['time']>= self.aux.Data['time'][wscan[icon[iscan,0]]]) &  \
                (self.rbd.Deconv['time']<= self.aux.Data['time'][wscan[icon[iscan,1]]])
            _husec_ = self.rbd.Deconv['husec'][itr]
            if _husec_.shape[0] != 0:
                _amplitude_,      = self.rbd.Deconv['amplitude'][itr],
                _right_ascension_ = np.interp(_husec_,self.aux.Data['husec'],self.aux.Data['right_ascension'])
                _declination_     = np.interp(_husec_,self.aux.Data['husec'],self.aux.Data['declination'])
                _elevation_       = np.interp(_husec_,self.aux.Data['husec'],self.aux.Data['elevation'])
                _azimuth_         = np.interp(_husec_,self.aux.Data['husec'],self.aux.Data['azimuth'])
            scans.update({
                'time'            : self.rbd.Deconv['time'][itr],
                'husec'           : _husec_,
                'amplitude'       : _amplitude_,
                'right_ascension' : _right_ascension_,
                'declination'     : _declination_,
                'elevation'       : _elevation_,
                'azimuth'         : _azimuth_
            })
            
        return scans

    def skyModel(self,x,Toff,Ts,tau):
        
        return Toff + Ts * ( 1 - np.exp(-tau/np.sin(np.radians(x))))

    def SkyDip(self):

        x, = np.where((self.aux.Data['opmode']==10))
        c  = contiguo.contiguo(x)
        time_interval = [self.rbd.Deconv['time'][x[0]],self.rbd.Deconv['time'][x[-1]]]

        elevation = []
        mVolts    = []
    
        for i in np.arange(len(c)):
            el = self.aux.Data['elevation'][x[c[i,0]]+2:x[c[i,1]]-2].mean()
            xx = (self.rbd.Deconv['husec'] >= self.aux.Data['husec'][x[c[i,0]]]) & (self.rbd.Deconv['husec'] <= self.aux.Data['husec'][x[c[i,1]]])
            mV = self.rbd.Deconv['amplitude'][xx].mean()
            elevation.append(el)
            mVolts.append(mV)

        par0     = [1.0e+02,1.0E+02,0.5]
        RERR     = 1.0
        TOLL     = 0.15
                
        while (RERR > TOLL):
            try:
                par, cov = curve_fit(self.skyModel, elevation, mVolts, p0 = par0)
                dfit     = self.skyModel(elevation,par[0],par[1],par[2])
                RERR     = np.sqrt(np.diag(cov))[2]/par[2]
                i = 0
                for val in mVolts:
                    if (np.abs(dfit[i]-val)/val > TOLL):
                        elevation.pop(i)
                        mVolts.pop(i)
                    i+=1
            except:
                self.skydip = -99
                return

        if par[2] < 0:
            self.skydip = -99
            return
        
        perr     = (np.sqrt(np.diag(cov)))
        dfit     = self.skyModel(elevation,par[0],par[1],par[2])

        self.skydip = {'tau':par[2],'stau':perr[2],
                       'Toff':par[0],'sToff': perr[0],
                       'Ts':par[1],'sTs': perr[1],
                       'Elevation':np.asarray(elevation),'mVolts':np.asarray(mVolts),
                       'Model':dfit,'Data':mVolts,'Trange':time_interval}

        return

    def plot(self):

        from matplotlib import pyplot as plt
        from matplotlib import dates

        thats = self.getTimeAxis(self.rbd.Deconv['husec'],self.MetaData)
        
        fig=plt.figure()
        fig.set_size_inches(((25*u.cm).to(u.imperial.inch)).value,
                            ((15*u.cm).to(u.imperial.inch)).value)
        ax=fig.add_subplot(1,1,1)
        pos=[0.1,0.1,0.85,0.85]
        ax.set_position(pos)
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
        ax.plot(thats,self.rbd.Deconv['amplitude'],'-r')
        ax.set_xlabel('UT')
        ax.set_ylabel('mV')
        if isinstance(self.MetaData['Date'],list):
            dia=self.MetaData['Date'][0]
        else:
            dia=self.MetaData['Date']
            
        ax.set_title(dia)
                         
        return thats


class env(object):

    def __init__(self):

        self.Data = np.empty((0))
        return
    
    def __str__(self):
        return 'A Class representing HATS Environment Data'

    def version(self):
        self.version = __version__
        return self.version

    def integrate(self,husec,cData):

        IntegrationInterval = 5000
        frac = husec.shape[0] // IntegrationInterval
        modulus = husec.shape[0] % IntegrationInterval
        NRecords = frac
        if modulus > 0 :
            NRecords+=1
        self.Data = np.zeros(NRecords,dtype=[('time',object),('husec','uint64'),('Tgolay','float'),
                                            ('Thics','float'),('Tenv','float')]
                            )
        for i in np.arange(0,NRecords):
            first = i*IntegrationInterval
            last  = (i+1)*IntegrationInterval
            if (last > husec.shape[0]):
                last = husec.shape[0]-1

            self.Data['husec'][i]  = husec[(first+last)//2]
            self.Data['Tgolay'][i] = cData['temp_golay'][first:last].mean()
            self.Data['Thics'][i]  = cData['temp_hics'][first:last].mean()
            self.Data['Tenv'][i]   = cData['temp_env'][first:last].mean()

        return
        
class aux(object):
    
    def __init__(self):

        self.Data = np.empty((0))
        self.MetaData = {}
        return
    
    def __str__(self):
        return 'A Class representing HATS Auxiliary Data'

    def version(self):
        self.version = __version__
        return self.version
    
    def from_file(self, fname, MetaData, offset=0, nrecords=-1):

        if 'Date' not in MetaData:
            MetaData.update({'Date':fname[-19:-9]})
            MetaData.update({'Hour':fname[-8:-4]})
        
        xml = xmlet.parse(MetaData['XMLAuxDescriptionFile']).getroot()
        self._tblheader = dict()
        self._tblheader = collections.OrderedDict()

        for child in xml:
            var_name = child[0].text
            var_dim  = int(child[1].text)
            var_type = child[2].text
            var_unit = child[3].text

            if var_type == "xs:int":
                np_type = np.int32
            if var_type == "xs:unsignedInt":
                np_type = np.uint32                
            if var_type == "xs:unsignedShort":
                np_type = np.uint16            
            if var_type == "xs:unsignedLong":
                np_type = np.uint64
            if var_type == "xs:double":
                np_type = np.float64

            self._tblheader.update({var_name:[var_dim, np_type, var_unit]})

        dt_list = list()
        for key, value in self._tblheader.items():
            dt_list.append((key, value[1], value[0]))

        fhandler = open(fname,'rb')
        whole_data = np.fromfile(fhandler, dt_list,count=nrecords,offset=offset)     # raw Data, no calibration

        flag = whole_data['husec'] < int(MetaData['Hour'][:2])*36000000              # catch the interval with wrong husecs
        if np.version.full_version > '1.19':
            Data = np.delete(whole_data,flag,0)                                     # delete wrong records            
        else:
            Data = np.delete(whole_data,np.where(flag),0)                                     # delete wrong records            
        
        deleted_records = whole_data.shape[0] - self.Data.shape[0]
        self.MetaData.update({'RBDRecDeleted':deleted_records})        
        self.MetaData.update({'N_Records_Read':whole_data.shape[0]})
        return Data
    
class rbd(object):

    def __init__(self,InputPath='./'):

        self.rData = np.empty((0))
        self.cData = np.empty((0))
        self.Deconv = np.empty((0))
        self.MetaData = {}
        self.Units = {'raw':{}, 'calibrated':{}}
        
        return

    def __str__(self):
        return 'A Class representing HATS Raw Binary Data'

    def version(self):
        self.version = __version__
        return self.version
    
    def from_file(self, fname, MetaData, recoffset=0, nrecords=-1):

        if 'Date' not in MetaData:
            MetaData.update({'Date':fname[-19:-9]})
            MetaData.update({'Hour':fname[-8:-4]})
            
        offset = recoffset * __RECORD_SIZE__
        xml = xmlet.parse(MetaData['XMLDataDescriptionFile']).getroot()
        FFTProgram = os.getenv('HATS_FFTProgram')
        if FFTProgram == None:
            MetaData.update({'FFTProgram':'./'})
        else:
            MetaData.update({'FFTProgram':FFTProgram})
        MetaData.update({'RBDFilename':fname})
        
        self._tblheader = collections.OrderedDict()

        for child in xml:
            var_name  = child[0].text
            var_dim   = int(child[1].text)
            var_type  = child[2].text
            var_unit  = child[3].text
            var_orig  = child[4].text
            var_conv  = child[5].text
            var_ord   = float(child[6].text)
            var_slop  = float(child[7].text)
            var_cunit = child[8].text

            if var_type == "xs:int":
                np_type = np.int32
            if var_type == "xs:unsignedInt":
                np_type = np.uint32                
            if var_type == "xs:unsignedShort":
                np_type = np.uint16            
            if var_type == "xs:unsignedLong":
                np_type = np.uint64

            self._tblheader.update({var_name:[var_dim, np_type, var_unit,var_orig,var_conv,var_ord,var_slop]})
            self.Units['raw'].update({var_name:var_unit})
            self.Units['calibrated'].update({var_name:var_cunit})

        dt_list = list()
        for key, value in self._tblheader.items():
            dt_list.append((key, value[1], value[0]))

        fhandler = open(fname,'rb')
        whole_data = np.fromfile(fhandler, dt_list,offset=offset,count=nrecords)     # raw Data, no calibration
        flag = whole_data['husec'] < int(MetaData['Hour'][:2])*36000000              # catch the interval with wrong husecs
        
        if np.version.full_version > '1.19':
            self.rData = np.delete(whole_data,flag,0)                                     # delete wrong records
        else:
            self.rData = np.delete(whole_data,np.where(flag),0)                                     # delete wrong records

        self.MetaData.update({'N_Records_Deleted':whole_data.shape[0] - self.rData.shape[0]})        
        
        # Number convertion for data coming from the Analog Devices A/D board.
        # They come in a 4 bytes integer, however, only the last 3 are significant.
        # The sign is in the third byte, the first bit.
        for key,val in self._tblheader.items():
            if (val[3] == "ad7770"):
                adcu = self.rData[key]&0x00FFFFFF             # Get the last 2 bytes
                negative = (self.rData[key]&0x0800000)>0      # Get the sign bit
                adcu[negative] = adcu[negative]-0x1000000    # Subtract the 'excess' to have negative numbers.
                self.rData[key] = adcu

        # Calibrated Data. 
        self.cData = {}
        self.CalFac = {}
        dt_list = list()
        conv_list = dict()
        for key, value in self._tblheader.items():
            if (value[4] == 'yes'):
                dt_list.append((key, np.float64,value[0]))
                conv_list.update({key:[value[6],value[5]]})

                
        self.cData = np.ndarray(self.rData.shape[0],dtype=dt_list)
        for key in dt_list:
            self.cData[key[0]] = self.rData[key[0]] * conv_list[key[0]][0] + conv_list[key[0]][1]
            self.CalFac.update({key:{'Ordinate':conv_list[key[0]][1],'Slope':conv_list[key[0]][0]}})


        if 'N_Records_Read' in self.MetaData:
            h.rbd.MetaData['N_Records_Read'] = whole_data.shape[0]
        else:
            self.MetaData.update({'N_Records_Read':whole_data.shape[0]})

        if 'Read_Offset_Record' in self.MetaData:
            h.rbd.MetaData['Read_Offset_Record'] = recoffset
        else:
            self.MetaData.update({'Read_Offset_Record':recoffset})

        return

    def check(self):

        t_initial = self.husec2dt(self.rData['husec'][0])
        t_final   = self.husec2dt(self.rData['husec'][-1])

        ds        = self.rData['sample'][1:]-self.rData['sample'][:-1]
        leap_ds   = (ds != 1)
        list_sample_leaps = np.where(leap_ds)

        dhusec     = self.rData['husec'][1:]-self.rData['husec'][:-1]
        leap_husec = (dhusec != 10)
        list_husec_leaps = np.where(leap_husec)

        sec      = np.empty(self.rData['sec'].shape[0],dtype=np.float64)
        sec      = self.rData['sec']+self.rData['ms']/1000
        dsec     = sec[1:]-sec[:-1]
        leap_sec = (np.abs(dsec) < 0.0009) | (np.abs(dsec) > 0.00101) 
        list_sec_leaps = np.where(leap_sec)
        
        # Print Stats
        print("\n\n")
        print(" - - - STATS - - - \n\n")
        print("Total Records : {0:d}".format(self.rData['sample'].shape[0]))
        print("Initial Date  : {0:s} \nFinal Date    : {1:s} \nTotal Duration: {2:d} seconds \n\n".format(t_initial.strftime("%c"),
                                                                                         t_final.strftime("%c"),
                                                                                         (t_final-t_initial).seconds ))
        print("Mean Values")
        print("          Golay : {0:f} mV".format(self.cData['golay'].mean()))
        print("   Temp HICS    : {0:f} C".format(self.cData['temp_hics'].mean()))
        print("   Temp Ambient : {0:f} C".format(self.cData['temp_env'].mean()))
        print("   Temp Golay   : {0:f} C".format(self.cData['temp_golay'].mean()))
        print("\n")

        if (list_sample_leaps[0].shape[0] > 0):
            print("    Sample leaps > 1\n")
            i=1
            for y in list_sample_leaps:
                for x in y:
                    init = self.rData['sample'][x]
                    final = self.rData['sample'][x+1]
                    if (final >= init):
                        leap = final-init
                    else:
                        leap = -int(init-final)
                        
                    print("      # : {0:d}    Index: {1:6d}  Init = {2:6d} Next = {3:6d} Leap = {4:6d}".format(i, x,init,final,leap))
                    i+=1

        if (list_husec_leaps[0].shape[0] > 0):
            i=1
            print("    Husec leaps > 10\n")
            for y in list_husec_leaps:
                for x in y:
                    init = self.rData['husec'][x]
                    final = self.rData['husec'][x+1]
                    if (final >= init):
                        leap = final-init
                    else:
                        leap = -int(init-final)
                    
                    print("      # : {0:d}    Index: {1:6d}  Init = {2:6d} husec Next = {3:d} husec , Leap = {4:f} seconds".format(i, x, init, final,leap*1e-4))
                    i+=1

        if (list_sec_leaps[0].shape[0] > 0):
            i=1
            print("    Sec leaps > 0.001\n")
            for y in list_sec_leaps:
                for x in y:
                    print("      # : {0:d}    Index: {1:6d}  , Leap = {2:f} seconds".format(i, x, sec[x+1]-sec[x]))
                    i+=1
                    
        print(" - - - - - - - - - - -  \n\n")
        return list_husec_leaps,list_sample_leaps
                                  
class ws(object):
    
    def __str__(self):
        return 'A Class representing HATS Weather Station Data'
    
    def __init__(self,date=''):
        
        warnings.filterwarnings('ignore')
        InputPath=os.getenv('HATS_WS_InputPath')
        if InputPath == None:
            InputPath = './'
            
        if (InputPath[-1] != '/') :
            InputPath+='/'        
        if (os.path.exists(InputPath)):
            self.MetaData={"InputPath":InputPath}
        else:
            raise ValueError("\n\n No Input Path {} found \n\n".format(InputPath))

        self.MetaData = {'InputPath':InputPath}

        if (len(date) == 10):
            YMD = date[:10]
        else:
            date = dt.datetime.now()
            YMD="{0:4d}-{1:02d}-{2:02d}".format(date.year,date.month,date.day)

        self.MetaData.update({'Date':YMD})
            
        fname = 'hats-'+YMD+'.ws'    # ASCII file
        self.MetaData.update({'Filename':fname})
        self.data = {}
        self.from_file(fname)
        self.Units = {'Temperature':'°C','Humidity':'%','Pressure':'HPa'}
        
        return

    def from_file(self,fname):
        
        fullpathname = self.MetaData['InputPath']+fname
        if (not os.path.exists(str(fullpathname))):
            raise ValueError('\n\n File {0:s} not found \n\n'.format(fullpathname))

        DateTime = []
        Temperature = []
        Humidity = []
        Pressure = []
        
        f = open(fullpathname,'r',errors='ignore')
        for line in f:
            s = line.split(sep=',')
            if (len(s) == 5):
                try:
                    DateTime.append(dt.datetime.fromisoformat(s[0]))
                    Temperature.append(float(s[2].split(sep='=')[1][:-1]))
                    Humidity.append(float(s[3].split(sep='=')[1][:-1]))
                    temp = s[4].split(sep='=')[1]
                    Pressure.append(float(temp[:temp.find('H')]))
                except:
#                    pdb.set_trace()
                    pass
                    
        f.close()
        self.data.update({'time':np.asarray(DateTime),
                          'temperature':np.asarray(Temperature),
                          'humidity':np.asarray(Humidity),
                          'pressure':np.asarray(Pressure)})    
        return

    def to_csv(self):
        
        csv_fname=self.MetaData['Filename'][:-3]+'.csv'
        f=open(csv_fname,'w')
        f.write("time,temperature [C],Humidity [%],Pressure [Hpa]\n")
        for i in np.arange(self.data['time'].shape[0]):
            line = self.data['time'][i].isoformat()            + ','
            line+='{0:6.1f}'.format(self.data['temperature'][i]) + ','
            line+='{0:6.1f}'.format(self.data['humidity'][i])    + ','
            line+='{0:6.1f}'.format(self.data['pressure'][i])    + '\n'
            f.write(line)

        f.close()

        return
    
    def dew(self,T,r):
    ############################
    #
    # Dew Point Temperature
    #
    #     Inputs:  T, dry-bulb temperature (Celsius)
    #              r, relative humidity in fraction (0...1)
    #
    #     Output:  Dew point temperature in Celsius
    #
    #
    #     Source: https://en.wikipedia.org/wiki/Dew_point (2021-07-31)
    #
    ####################################################################
    
        a = 6.1121 * u.mbar
        b = 18.678 * u.Unit('')
        c = 257.14 * u.Celsius
        d = 234.50 * u.Celsius

        Psm = a * np.exp( (b-T/d) * (T/(c+T)) )
        gm  = np.log(r * np.exp( (b-T/d) * (T/(c+T))) )
        return   c * gm / (b - gm)

    def P0(self,T,h):
    ##############################
    #
    # Water Vapor Partial Pressure
    #
    #    Inputs:  T, dry-bulb temperature in Celsius
    #             h, relative humidity in %
    #
    #    Output:  Water Vapor partial pressure in mbar
    #
    #    Source: ALMA memo 237 - Butler, B (1998)
    #
    ######################################################
    
        D    = self.dew(T, h/100 )
        expo = 1.81E+00 + (1.727E+01 * D /(D + 237.3 * u.Celsius))
        return  np.exp(expo) * u.mbar
    
    def pwv(self,h2o_scale=2.0):

    #################################
    #
    # Precipitable Water Vapor
    #
    #    Parameter inputs:
    #           temperature = dry-bulb temperature in Celsius
    #           humidity    = relative humidity in %
    #           Hh2o        = Water scale height in km (typically 1.5 - 2 )
    #
    #    Output: Precipitable Water Vapor in mm
    #
    #    Source: ALMA memo 237 - Butler, B (1998)
    #
    ########################################################

        temperature = self.data['temperature'] * u.Celsius
        humidity    = self.data['humidity'] * u.Unit('')
        Hh2o        = h2o_scale*u.km
        
        amu      = 1.66053906660e-24 * u.g  # Atomic Mass Unit in g
        mw       = 18 * amu                 # Water mole mass
        rhol     = 1 * u.g / u.cm**3        # water density
        constant = mw / ( rhol * c.k_B.cgs)
        
        num = constant * self.P0(temperature,humidity).cgs * Hh2o.to('cm')
        den = temperature.to(u.Kelvin, equivalencies=u.temperature())

        self.data.update({'PWV':(num/den).to('mm')})
        self.MetaData.update({'H2O_scale_height':Hh2o})
        
        return  
