import os
from pathlib import Path
import xml.etree.ElementTree as xmlet
import numpy as np
from astropy.io import fits
import collections
import pdb

__version__       = "2021-11-09T2045BRT"
__SIGNAL_LENGTH__ =  8192
__DATA_FILE__     = "hats_data_rbd.bin"
__HUSEC_FILE__    = "hats_husec.bin"

##############################################
#
# HATS: A python class to read and deconvolve HATS RBD data.
#
# Usage:
#   >>> import HATS
#   >>> h=HATS.rbd(PathToXML='./XMLTables')        # Create the Object. PathToXML should point where XML shema files are
#   >>> h.from_file('hats-2021-08-26T1800.rbd')    # Import data from RBD file
#   >>> h.getFFT()                                 # Deconvolve the data.
#
# Data Structures:
#   MetaData:   Dictionary with
#                               InputPath : string with directory where files are
#                               PathToXML : string with XML schema files
#                               File_Name : string with the RBD name file
#                               Date      : ISO data date
#                               Hour      : Data Hour (UTC)
#  Data: A numpy array with
#                               sample    : sample number (np.uint32 )
#                               sec       : seconds from 1970-01-01, Unix time (np.uint32)
#                               ms        : milliseconds of the seconds (np.uint16)
#                               husec     : hundred of microsends since 0 UT, as in SST (np.uint64)
#                               adcu      : Golay output in ADC units (np.uint32)
#                               power_supp: Power Supply in ADC units (np.uint32)
#
# Deconv: a Dictionary with
#                               amplitude : FFT amplitude in ADC units (np.float64)
#                               husec     : time of 'amplitude' in hundred of microsends since 0 UT, as in SST (np.uint64)
#                               time      : time of 'amplitude' in datetime()
#
# Author: @guiguesp - São Paulo - returning to presencial classes - 2021-10-12
#
############################################################################

class rbd(object):

    def __init__(self,FirmwareVersion='1.10',PathToXML='',InputPath='./'):
        self.MetaData = {}
        self.Data     = np.empty((0))

        self.MetaData.update({'InputPath':InputPath})
        
        # PathToXML should point to tha directory where the XML tables are copied
        # When not defined, look at the environment
        if (isinstance(PathToXML,str) and len(PathToXML) == 0):
            if ('HATSXMLPATH' in os.environ.keys()):
                self.MetaData.update({'PathToXML' :os.environ['HATSXMLPATH']})
            else:
                self.MetaData.update({'PathToXML':'./'})
        else:
            self.MetaData.update({'PathToXML':PathToXML})

        if self.MetaData['PathToXML'][-1] != '/' :
            self.MetaData['PathToXML'] = self.MetaData['PathToXML'] + '/'
            
        self.MetaData['Firmware_Version'] = FirmwareVersion

    def __str__(self):
        return 'A Class representing HATS Raw Binary Data'

    def version(self):
        self.version = __version__
        return self.version
    
    def from_file(self, fname):

        fullpathname = self.MetaData['InputPath']+'/'+fname
        if (not os.path.exists(str(fullpathname))):
            print('\n\n File {0:s} not found \n\n'.format(fullpathname))
            return

        self.MetaData.update({'File_Name':fname})
        if ( (self.MetaData['File_Name'][:4] != "hats") & (self.MetaData['File_Name'][-4:] != "rbd") ) :
            raise ValueError("Invalid file type {}".format(self.filename))
        else:
            self.MetaData.update({'Date':fname[5:15]})
            self.MetaData.update({'Hour':fname[16:18]})

        xml = xmlet.parse(self.MetaData['PathToXML'] / Path("HATSDataFormat.xml")).getroot()
        self._tblheader = dict()
        self._tblheader = collections.OrderedDict()
        for child in xml:
            var_name = child[0].text
            var_dim  = int(child[1].text)
            var_type = child[2].text
            var_unit = child[3].text
            var_orig = child[4].text
            var_conv = child[5].text
            var_ord  = float(child[6].text)
            var_slop = float(child[7].text)

            if var_type == "xs:int":
                np_type = np.int32
            if var_type == "xs:unsignedInt":
                np_type = np.uint32                
            if var_type == "xs:unsignedShort":
                np_type = np.uint16            
            if var_type == "xs:unsignedLong":
                np_type = np.uint64

            self._tblheader.update({var_name:[var_dim, np_type, var_unit,var_orig,var_conv,var_ord,var_slop]})

        dt_list = list()
        for key, value in self._tblheader.items():
            dt_list.append((key, value[1], value[0]))
        
        fhandler = open(fullpathname,'rb')
        self.rData = np.fromfile(fhandler, dt_list)     # raw Data, no calibration

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
        for key,val in self._tblheader.items():
            if (val[4] == "yes"):                                            # Data must be calibrated
                self.cData.update({key:np.zeros(self.rData[key].shape[0])})  # Create room for the data
                self.cData[key] = self.rData[key] * val[6] + val[5]          # apply a linear conversion

                                   
        return

    def husec2dt(self,husec):
        import datetime as dt

        ms           = husec
        hours        = int(ms // 36000000)
        minutes      = int((ms % 36000000) // 600000)
        seconds      = ((ms % 36000000) % 600000) / 1.0E+04
        seconds_int  = int(seconds)
        seconds_frac = seconds - int(seconds)
        useconds     = int(seconds_frac * 1e6)
        year         = int(self.MetaData['Date'][0:4])
        month        = int(self.MetaData['Date'][5:7])
        day          = int(self.MetaData['Date'][8:10])
        return dt.datetime(year,month,day,hours,minutes,seconds_int,useconds)

    def getTimeAxis(self):
        """

        getTimeAxis: Class method to convert the us time axis used in RBD files to a Python
                     datetime ndarray that can be used with matplotlib.pyplot.

        Change Record:  First written by Guigue @ Sampa
                        2017-11-04 St Charles Day !!!

        """

        dt_time = list(map(self.husec2dt,self.rData['husec']))
        return np.asarray(dt_time)

    def getFFT(self):

        Nchunks = self.rData['husec'].shape[0] // __SIGNAL_LENGTH__
        Lchunk = self.rData['husec'].shape[0] % __SIGNAL_LENGTH__
        amp = []
        hus = []
        
        self.rData['husec'].tofile(__HUSEC_FILE__)
        self.rData['adcu'].tofile(__DATA_FILE__)
        os.spawnl(os.P_WAIT,'./HATS_fft','HATS_fft')
        
        husec  = np.fromfile(__HUSEC_FILE__ , dtype=np.uint64  , count=-1)
        fftAmp = np.fromfile(__DATA_FILE__  , dtype=np.float64 , count=-1)

        os.remove(__HUSEC_FILE__)
        os.remove(__DATA_FILE__)
        dt_time = list(map(self.husec2dt,husec))
        self.Deconv = {'amplitude': fftAmp, 'husec':husec, 'time':np.asarray(dt_time)}

        return 
