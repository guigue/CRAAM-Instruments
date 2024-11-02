# External methods
import sys, string, os, struct, glob
import numpy as np
import xml.etree.ElementTree as xmlet
from astropy.io import fits
import CASLEO
import oRBD
import pdb


################################################################
#                                                              #
#  Version Number. Change everytime the code is changed.       #
#                                                              #
_Version = '2024-01-23T11:37BRT'                               #
#                                                              #
################################################################
class RBD(object):

###############################################################################################
#
#    RBD: A python class to use with SST data, its main objective is to create SST level-0 files.
#         The class reads one file (the full path can be introduced with the name) and writes
#         one fits file in the same directory where the programm is running. A python script can be
#         written in order to automate the conversion process.  The fits files are an exact copy
#         of the original SST files, although I believe some changes are introduced in the conversion
#         process because of the binary representation of the numbers. python reads the binary files
#         and stores them in internal memories with its own binary format. When it stores the data
#         in the fits file, it 'casts' the numbers to the fits binary representation.
#
#         Some Glossary:
#             RBD : Raw Binary Data. The original SST files in the formats 'slow' (RS), 'fast' (RF)
#                   and 'monitoring' (BI)
#             FITS level-0: fits file produced from an RBD file.  It should be an identical copy
#                  of the RBD, the only difference are the headers with a description of the data
#                  including the units
#
#         The fits level-0 file has two headers. The primary one has general information about the
#         data, including frequencies, geographical coordinates of the observatory, etc. The second
#         header includes th description of the data represented as a binary table. The description
#         includes the data units.
#
#         The description of the data is included in the XML files. These files should be in the same
#         directory of the python program.
#
#         New methods were added to this version.
#           writeFITSwithName() : allows to give a non-standard name to the fits file
#           concat() : allows to concatenate two or more RBD objects
#           reduced() : creates a reduced version of the RBD object, leaving the most used variables.
#
#        Also the calling sequence was changed. Now when creating an object only some basic parameters
#        are defined. For instance the input and output directories. You may changed them afterwards
#        >>> d.InputPath = /path/to/files
#        >>> d.OutputPath = /path/to/fits
#
#    Tests:
#
#         In order to assess the accuracy of the conversion to fits, a test was carried out using the
#         following files
#             bi1010822 --> sst_auxiliary_2001-08-22T14:04:40.034-23:45:04.714_level0.fits
#             bi1021019 --> sst_auxiliary_2002-10-19T11:24:39.704-21:37:04.235_level0.fits
#             bi1021202 --> sst_auxiliary_2002-12-02T02:01:34.154-23:18:05.244_level0.fits
#             bi1021221 --> sst_auxiliary_2002-12-21T11:15:57.732-22:26:07.212_level0.fits
#             rs1020715.1300 --> sst_integration_2002-07-15T13:24:19.605-13:59:59.592_level0.fits
#             rs1021205.2200 --> sst_integration_2002-12-05T21:59:58.728-22:08:31.215_level0.fits
#             rs1061206.2100 --> sst_integration_2006-12-06T20:59:58.423-21:54:05.746_level0.fits
#             rs990909.1700 --> sst_integration_1999-09-09T16:59:52.056-17:57:40.420_level0.fits
#
#         Every file was later opened using SSW's mrdfits() and SST read_sst procedures.  The data
#         was compared doing a substraction, and obtaining the mean and standard deviation. Were the
#         mean different from zero, the program stopped and issue a message.  We found many tiny errors
#         in the python procedures oRBD, that were corrected, repeating the test until mean = 0.
#         See testRBD.py and testFITS.pro
#
#    Use:
#            It is recommended to set the path where the XML files are. For instance
#            $ export RBDXMLPATH=/path/to/xmlfiles
#            > setenv RBDXMLPATH /path/to/xmlfiles
#
#            >>> import oRBD
#            >>> d=oRBD.RBD(InputPath=[path],OutputPath=[path])
#            >>> d.readRBDinDictionary([RBD filename])
#            >>> d.writeFITS()
#            If you want to put a non-standar name use
#            >>> d.writeFITSwithName([FITS filename])
#
#            Check
#            >>> from astropy.io import fits
#            >>> h=fits.open([fits filename])
#            >>> print h.info()
#            >>> print(repr(h[0].header))             # Primary Header
#            >>> print(repr(h[1].header))             # Table header
#            >>> h[1].Data                            # Actual Data
#
#             Another check can be done with Solar Software
#             IDL> r=mrdfits([fits filename],1,h,/unsigned) ; VERY important the unsigned parameter!!
#             IDL> read_sst,d,[RBD filname],/close,/recr=1000000 [,/mon]
#             Then r and d are indentical structures.
#
#    Author:  Guigue@Sampa
#             guigue@craam.mackenzie.br  , guigue@gcastro.net
#             http://www.guigue.gcastro.net
#
#
#    Change Record :
#         - 2017-06-15 : First written
#         - 2017-08-19 : xml files implementation
#         - 2017-08-29 : writeFITS implementation
#         - 2017-11-02 : Check that PathToXML points to the XML tables repositories.
#                        readRBDinDictionary() and writeFITS() methods return True or False.
#
#         - 2018-10-09 : adapted for python 3
#                        print --> print()
#                        integer division now is //
#
#         - 2018-10-12 : adapted for python 3
#                        print --> print()
#                        integer division now is //
#                        viewkeys() -> keys()
#                        reduce() method was altered
#         - 2020-04-26 : added two new methods to RBD
#                        __add__ : overload of '+' it allows to do >> d = d1+d2, with d1, d2 two RBD objects
#                        extract : extracts an interval of time example
#                             >>> time_interval = [datetime.datetime(2020,4,26,20,30,15,350), datetime.datetime(2020,4,26,20,40,28,550)]
#                             >>> d1 = d.extract(time_interval)
#         - 2020-05-29 : now CASLEO information is in separate function.
#         - 2020-06-18 : stupid mistake corrected !
#
###############################################################################################

    @property
    def version(self):
        return self._version

    def __init__(self,PathToXML='',InputPath='./',OutputPath='./'):

    # PathToXML should point to tha directory where the XML tables are copied
    # When not defined, look at the environment
        if (isinstance(PathToXML,str) and len(PathToXML) == 0):
            if ('RBDXMLPATH' in os.environ.keys()):
                self.PathToXML = os.environ['RBDXMLPATH']
            else:
                self.PathToXML = './'
        else:
            self.PathToXML=PathToXML

        if self.PathToXML[-1] != '/' :
            self.PathToXML = self.PathToXML+'/'

    # Check the existence of the XML tables
        if not self.CheckXMLTables() :
            return

        self.OutputPath = OutputPath
        self.InputPath  = InputPath

        self.Data     = {}
        self.MetaData = {}
        self.History  = []
        self._version = _Version

        return

    def __str__(self):
        return 'A Class representing SST Raw Binary Data'

    def extract(self,time_interval):
        _temp_ = oRBD.RBD()
        _temp_.header = self.header
        _temp_.MetaData = self.MetaData
        _temp_.History.append('Extracted Interval : ' + str(time_interval))

        t0 = int( 1.0E+04 * (time_interval[0].hour * 3600 + time_interval[0].minute * 60 + (time_interval[0].second + time_interval[0].microsecond/1.0E+06)))
        t1 = int( 1.0E+04 * (time_interval[1].hour * 3600 + time_interval[1].minute * 60 + (time_interval[1].second + time_interval[1].microsecond/1.0E+06)))

        TagList = list(self.Data.keys())
        for iTag in TagList:
            _temp_.Data[iTag] = self.Data[iTag][ (self.Data['time']>= t0 ) & ( self.Data['time'] <= t1)]

        return _temp_

    def __add__(self,d):

        _temp_ = RBD()

        _temp_.Data      = {}
        _temp_.Metadata  = {}
        _temp_.header    = self.header
        _temp_.History   = self.History

        if isinstance(self.MetaData['ISOTime'],str):
            ISOTime     = [self.MetaData['ISOTime']]
        else:
            ISOTime     = self.MetaData['ISOTime']
        ISOTime.append(d.MetaData['ISOTime'])
        
        if isinstance(self.MetaData['RBDFileName'],str):
            RBDFileName = [self.MetaData['RBDFileName']]
        else:
            RBDFileName = self.MetaData['RBDFileName']
            
        RBDFileName.append(d.MetaData['RBDFileName'])
        SSTType     = self.MetaData['SSTType']
        ISODate     = self.MetaData['ISODate']

        _temp_.MetaData = { 'ISODate': ISODate,
                            'ISOTime': ISOTime,
                            'RBDFileName': RBDFileName,
                            'SSTType':SSTType}

        TagList = list(self.Data.keys())
        for iTag in TagList:
            _temp_.Data[iTag] = np.concatenate((self.Data[iTag],d.Data[iTag]))

        _temp_.History.append('Concatenated Data')

        return _temp_

    def getTimeAxis(self):
        """

        getTimeAxis: Class method to convert the us time axis used in RBD files to a Python
                     datetime ndarray that can be used with matplotlib.pyplot.

        Change Record:  First written by Guigue @ Sampa
                        2017-11-04 St Charles Day !!!

        """

        import datetime as dt
        ndata = self.Data['time'].shape[0]

        ssttime = np.array(np.empty(ndata),dtype=dt.datetime)
        year  = int(self.MetaData['ISODate'][0:4])
        month = int(self.MetaData['ISODate'][5:7])
        day   = int(self.MetaData['ISODate'][8:])
        for i in np.arange(ndata):
            ms = self.Data['time'][i]
            hours =  ms // 36000000
            minutes = (ms % 36000000) // 600000
            seconds = ((ms % 36000000) % 600000) / 1.0E+04
            seconds_int  = int(seconds)
            seconds_frac = seconds - int(seconds)
            useconds     = int(seconds_frac * 1e6)
            ssttime[i] = dt.datetime(year,month,day,hours,minutes,seconds_int,useconds)

        return ssttime

    """-------------------------------------------------------------------------------------"""
    def CorrectAuxiliary(self):
        """
        CorectAuxiliary: For some unknown reason Auxiliary files (a.k.a. BI files) have a 0 record
                         when SST start_obs command is given. A 0 record means all fields are 0.
                         This method corrects this problem deleting the flawed records.
                         To verify the record does not correspond to O UT, we check that ADC are 0 too.

                         It only works for Auxiliary files.

        Change Record: First written by Guigue @ Sampa
                       2017-11-04 St Charles Day !!
                       2017-11-07 Completely changed, using python tools. Ethernal LOVE to Python!!

        """

        if (self.MetaData['SSTType'] != 'Auxiliary'):
            return

        Mask = (self.Data['time'] != 0 & np.all(self.Data['adc'] != 0) )
        TagList = self.Data.keys()
        for tag in TagList:
            v=self.Data[tag][Mask]
            self.Data[tag]=v

        return

    """------------------------------------------------------------------------------------ """

    def timeSpan(self):
        """
        timeSpan:
             looks for the first and last non-zero time in the time array and converts to
             ISO time (HH:MM:SS.SSS)

        Change Record:
             First written by Guigue @ Sampa - 2017-08-26

        """

        # Sometimes time=0 is registered. We do not normally observe at 0 UT, in the contrary
        # this false value is created by a wrong program. In order to know the starting and end
        # time of the data, we serch for values other than 0
        _nonzero_=self.Data['time'].nonzero()

        # Return the first and last _nonzero_ time converted to ISO standards
        return self.getISOTime(self.Data['time'][_nonzero_[0][0]]) ,self.getISOTime(self.Data['time'][_nonzero_[0][-1]])


    """-------------------------------------------------------------------------------------"""
    def base_name(self,fname):
        """
        base_name:
            Simple procedure to get the base name of a full described file name
            /path/to/file/filename --> filename

        Change Record:
            Guigue @ Sampa - 2017-08-26

        """
        _s_ = fname.strip().split('/')
        return _s_[-1]

    """------------------------------------------------------------------------------------ """

    def getISOTime(self,hustime):
        """
        getISOTime:
           Convert from Hus (hundred of microseconds) to HH:MM:SS.SSSS

        inputs:
           hustime : a 4bytes integer with the time in hundred of microseconds since 0 UT.
                     hustime is the SST time format.

        Change Record:
           First written Guigue @ Sampa - 2017-08-26

        """
        _hours_ = hustime // 36000000
        _minutes_ = (hustime % 36000000) // 600000
        _secs_ = (hustime - (_hours_ * 36000000 + _minutes_ * 600000)) // 10000.0
        return '{0:=02d}'.format(_hours_)+':'+ '{0:=02d}'.format(_minutes_) +':'+'{0:=06.3f}'.format(_secs_)

    """ -------------------------------------------------------------------------------------"""

    def getISODate(self):
        """
        getISODate:
        Converts an SST filename in a ISO Date format. e.g.: RS1070812.1500 -> 2007-08-12
        It also sets the data type to Data (for RS and RF files) and Auxiliary (for BI)

        Input:
        the SST filename.

        Output:
        A list [Type,ISODate]

        Change Record:
        First written by Guigue @ Sampa on 2017-08-19
        """
        _bname_=self.base_name(self.RBDfname)                                # Get the base name (remove /path/to/file )
        try:
            _name1_,_name2_=_bname_.strip().split(".")                       # Does it have hours?
        except:
            _name1_ = _bname_.strip()
            _name2_ = ''

        _SSTprefix_ = _name1_[0:2].upper()
        _SSTtype_ = ''

        if (_SSTprefix_ == 'RS'): SSTtype='Integration'
        if (_SSTprefix_ == 'RF'): SSTtype='Subintegration'
        if (_SSTprefix_ == 'BI'): SSTtype='Auxiliary'

        # The SST RBD file names are formed of three parts:
        #   2 characters to describe the type of data
        #   8 or 9 characters to describe year, month and day
        #   The separator dot .
        #   For Data files (RS,RF) the time in the format HHMM
        # Therefore a data file opened on 2017-08-31 at 17:32 hs
        # is named RF1170831.1732. When the file correspond to
        # Integrated Data the minutes are removed: RS1170831.1700
        # For Auxiliary data, there is no indication of the time: BI1170831
        #
        # SST started to observe in 1999. The mechanism to create their filenames
        # is wrong, since it takes the number of years since 1900.  Therefore
        # the first year files started with 99, having the day description
        # with 8 characters while, since 2000, files start with 100 and the
        # day description has 9 characters.
        if (len(_name1_) == 8):
            date=str(int(_name1_[2:4])+1900) + '-' + _name1_[4:6] + '-' + _name1_[6:8]
        elif (len(_name1_) == 9):
            date=str(int(_name1_[2:5])+1900) + '-' + _name1_[5:7] + '-' + _name1_[7:9]
        else:
            print (self.RBDfname+ '  is a wrong RBD Filename. Aborting...')
            return False

        # From the time description of the RBD file name we get
        if (len(_name2_) == 4) :
            time=_name2_[0:2]+':'+_name2_[2:4]
        else:
            time='00:00'

        self.MetaData.update({'RBDFileName' : _bname_})
        self.MetaData.update({'ISODate'     : date   })
        self.MetaData.update({'ISOTime'     : time   })
        self.MetaData.update({'SSTType'     : SSTtype})

        return True

    """-----------------------------------------------------------------------------"""

    def define_fmt(self):
        """

        define_fmt:
        It defines three elements to be used for reading  the raw binary data (RBD) files.
        1) fmt: The binary format in pythons rules. It is a string. With this string python
                can read one record and put it in a list (vector) of 'unpacked' data
        2) ranges : it is a list of 2-members lists. Every member points to the begginning and the end
                    of a variable in the 'unpacked' list (vector)
        3) fieldnames : it is a list of field names in the binary file. It will be used to create the FITS
                        header

        Change Record:
        First created by Guigue @ Sampa 2017.08.20 : checked against the xml files only.

        """
        # Internal variables destroyed at the end of procedure
        _fmt_        = '='
        _fieldnames_ = []
        _ranges_     = []
        _fielditer_  = 0
        _TotalDim_   = 0

        for child in self.header:

            # xml table. Children have three fields
            _VarName_ = child[0].text                                          # Variable Name
            _VarDim_  = int(child[1].text)                                     # Variable Dimension
            _TotalDim_ += _VarDim_
            _VarType_ = child[2].text                                          # Variable type

            _fieldnames_.append(_VarName_)
            _ranges_.append([_fielditer_,_fielditer_+_VarDim_-1])
            _fielditer_+=_VarDim_

            for i in range(_VarDim_):
                # A short table that converts from xml types to python bynary formats
                if ( _VarType_ == 'xs:int')           : _fmt_ = _fmt_ + 'i'  # a 4 bytes integer
                if ( _VarType_ == 'xs:unsignedShort') : _fmt_ = _fmt_ + 'H'  # a 2 bytes unsigned integer
                if ( _VarType_ == 'xs:short')         : _fmt_ = _fmt_ + 'h'  # a 2 bytes integer
                if ( _VarType_ == 'xs:byte')          : _fmt_ = _fmt_ + 'B'  # a byte
                if ( _VarType_ == 'xs:float')         : _fmt_ = _fmt_ + 'f'  # a 4 bytes float

            bin_header = {'names':_fieldnames_, 'ranges':_ranges_,'fmt':_fmt_,'TotalDim':_TotalDim_}
            self.bin_header=bin_header
        return

    """-----------------------------------------------------------------------------------------"""

    def read_xml_header(self):
        if not self.getISODate():
            return False
        _tt_        = oRBD.DataTimeSpan(self.PathToXML)
        _hfname_    = _tt_.findHeaderFile(SSTType=self.MetaData['SSTType'],SSTDate=self.MetaData['ISODate'])
        _xmlheader_ = xmlet.parse(self.PathToXML + _hfname_)
        self.hfname = _hfname_
        self.header = _xmlheader_.getroot()
        return True

    """-----------------------------------------------------------------------------"""

    def readRBDinDictionary(self,RBDfname='rs990501'):

        """
        readRBDinDictionary
           It is the class method used to read a SST Raw Binary Data (RBD). The data is
           represented with a dictionary, each element represents one SST variable, and data
           is stored in a numpy ndarray. Every ndarray has the python dtype corresponding to
           the original SST data.

        Output:
           It returns True or False

        Change Record:
           First Written by Guigue @ Sampa - 2017-08-26

        """

        self.RBDfname = RBDfname
        if not self.read_xml_header():
            return

        self.define_fmt()
        self.CleanPaths()
        self.History = []

        _fmt_    = self.bin_header['fmt']
        _header_ = self.bin_header['names']
        _ranges_ = self.bin_header['ranges']
        _Nfields_= len(_header_)

        if os.path.exists(self.InputPath+self.RBDfname) :
            _fd_         = os.open(self.InputPath+self.RBDfname,os.O_RDONLY)
            _nrec_       = os.fstat(_fd_).st_size // struct.calcsize(_fmt_)

            for child in self.header:

                # xml table. Children have three fields
                _VarName_ = child[0].text                                          # Variable Name
                _VarDim_  = int(child[1].text)                                     # Variable Dimension
                _VarType_ = child[2].text                                          # Variable type

                 # A short table that converts from xml types to python bynary formats
                if ( _VarType_ == 'xs:int') :
                     _nptype_ = np.int32   # a 4 bytes integer
                elif ( _VarType_ == 'xs:unsignedShort') :
                     _nptype_ = np.uint16  # a 2 bytes unsigned integer
                elif ( _VarType_ == 'xs:short')         :
                     _nptype_ = np.int16   # a 2 bytes integer
                elif ( _VarType_ == 'xs:byte')          :
                     _nptype_ = np.byte    # a byte
                elif ( _VarType_ == 'xs:float')         :
                     _nptype_ = np.float32 # a 4 bytes float
                else:
                     _nptype_ = np.float32

                if (_VarDim_ == 1):
                    self.Data.update( {_VarName_ : np.array(np.empty([_nrec_],_nptype_)) })
                else:
                    self.Data.update( {_VarName_ : np.array(np.empty([_nrec_,_VarDim_],_nptype_))})

            for _irec_ in range(_nrec_) :
                _one_record_  = os.read(_fd_,struct.calcsize(_fmt_))
                _ur_          = struct.unpack(_fmt_,_one_record_)

                for _field_ in range(_Nfields_):
                    if (_ranges_[_field_][0] == _ranges_[_field_][1]) :
                        self.Data[_header_[_field_]][_irec_]= _ur_[_ranges_[_field_][0]]
                    else:
                        self.Data[_header_[_field_]][_irec_,...] = _ur_[_ranges_[_field_][0]:_ranges_[_field_][1]+1]

            os.close(_fd_)
        else:
            print ('File '+self.InputPath+self.RBDfname+'  not found. Aborting...')
            return False

        return True

    """-----------------------------------------------------------------------------"""
    def CleanPaths(self):

        """
        It seems that fits.writeto does not understand the meaning of '~/'
        We chhange for $HOME/
        """
        path = self.OutputPath.strip().split('/')
        if (path[0] == '~') : path[0]=os.environ['HOME']
        newpath=''
        for ipath in path:
            newpath+=ipath+'/'
        if newpath[-2]=='/' : newpath=newpath[0:-1]
        self.OutputPath=newpath

        path = self.InputPath.strip().split('/')
        if (path[0] == '~') : path[0]=os.environ['HOME']
        newpath=''
        for ipath in path:
            newpath+=ipath+'/'
        if newpath[-2]=='/' : newpath=newpath[0:-1]
        self.InputPath=newpath

        return

    """-----------------------------------------------------------------------------"""

    def writeFITS(self) :

        _hhmmss_ = self.timeSpan()
        self.History.append('Converted to FITS level-0 with oRBD.py version '+self.version)
        return self.writeFITSwithName('sst_'  +
                               self.MetaData['SSTType'].lower() + '_' +
                               self.MetaData['ISODate'] + 'T' +
                               _hhmmss_[0]+'-' + _hhmmss_[1] +
                               '_level0.fits')


    """-----------------------------------------------------------------------------"""

    def writeFITSwithName(self,FITSfname):

        """
        writeFITS:
             A method to write the SST data as a FITS file. The fits format used is a binary table.
             Then the file has a primary and a table or secondary header. It also defines the name of the fits file as

             sst_[integration | subintegration | auxiliary]_YYYY-MM-DDTHH:MM:SS.SSS-HH:MM:SS.SSS_level0.fits

             integration is a RS file, subintegration is a RF file and auxiliary is a BI file.
             If the file already exist, the system breaks.

             The system implements two headers. The primary header has general information, while the secondary
             header is specific for the table, including the units of the columns.

        Output:
             It returns True or False on success or failure respectively.

        Change Record:
             First written by Guigue @ Sampa - 2017-08-26
             Return value added on 2017-11-02

        """

        cg = CASLEO.Observatory_Coordinates()

        self.CleanPaths()
        self.MetaData.update({'FITSfname':FITSfname})

        _isodate_ = self.MetaData['ISODate']
        _hhmmss_ = self.timeSpan()

        _hdu_ = fits.PrimaryHDU()

        #
        # This is the Primary (global) header. It gives information about the instrument, and the data.
        # It probably would be better to include this information in the XML decriptive files.
        # Anyway, oRBD.py and the XML files are a linked together.
        #
        _hdu_.header.append(('origin','CRAAM/Universidade Presbiteriana Mackenzie',''))
        _hdu_.header.append(('telescop','Solar Submillimeter Telescope',''))
        _hdu_.header.append(('observat','CASLEO',''))

        _hdu_.header.append(('station','Lat = {0:12.8f}, Lon = {1:12.8f}, Height = {2:5.3f} km'.format(cg.lat.value, cg.lon.value, cg.height.value/1000),''))
        _hdu_.header.append(('tz','GMT-3',''))

        _hdu_.header.append(('date-obs',_isodate_,''))
        _hdu_.header.append(('t_start',_isodate_+'T'+ _hhmmss_[0],''))
        _hdu_.header.append(('t_end',_isodate_+'T'+ _hhmmss_[1],''))
        _hdu_.header.append(('data_typ',self.MetaData['SSTType'],''))
        if isinstance(self.MetaData['RBDFileName'],list) :
            for iRBD in self.MetaData['RBDFileName']:_hdu_.header.append(('origfile',iRBD,'SST Raw Binary Data file'))
        else:
            _hdu_.header.append(('origfile',self.MetaData['RBDFileName'],'SST Raw Binary Data file'))

        _hdu_.header.append(('frequen','212 GHz ch=1,2,3,4; 405 GHz ch=5,6',''))

        # About the Copyright
        _hdu_.header.append(('comment','COPYRIGHT. Grant of use.',''))
        _hdu_.header.append(('comment','These data are property of Universidade Presbiteriana Mackenzie.'))
        _hdu_.header.append(('comment','The Centro de Radio Astronomia e Astrofisica Mackenzie is reponsible'))
        _hdu_.header.append(('comment','for their distribution. Grant of use permission is given for Academic '))
        _hdu_.header.append(('comment','purposes only.'))

        for i in range(len(self.History)):
            _hdu_.header.append(('history',self.History[i]))

        _fits_cols_ = []
        for _child_ in self.header:

            # xml table. Children have four fields
            _ttype_   = _child_[0].text      # Name
            _tform_   = _child_[1].text      # (Dimension) Format
            _tunit_   = _child_[3].text      # Unit
            _tzero_   = 0                    # Effective 0 (to mimmic an unsigned integer)
            _tscal_   = 1.0                  # Data Scaling Factor
            _VarType_ = _child_[2].text      # Variable type

            if ( _VarType_ == 'xs:int') :
                _tform_ += 'J'
                _np_form_ = np.dtype('i4')

            if ( _VarType_ == 'xs:unsignedShort') :
                # FITS does not have unsignedShort representation
                # The way to represent them correctly is to use
                # 4 bytes integers and add a 32768 offset
                #
                # Important! Use /unsigned with mrdfits()
                _tform_ += 'I'
                _tzero_ = 32768
                _np_form_ = np.dtype('u2')

            if ( _VarType_ == 'xs:short'):
                _tform_ += 'I'
                _np_form_ = np.dtype('i2')

            if ( _VarType_ == 'xs:byte') :
                _tform_ += 'B'
                _np_form_=np.dtype('b')

            if ( _VarType_ == 'xs:float') :
                _tform_ += 'E'
                _np_form_=np.dtype('f4')

            _c_ = fits.Column(name   = _ttype_  ,
                              format = _tform_  ,
                              unit   = _tunit_  ,
                              bscale = _tscal_  ,
                              bzero  = _tzero_  ,
                              array  = self.Data[_ttype_])
            _fits_cols_.append(_c_)

        _coldefs_ = fits.ColDefs(_fits_cols_)
        pdb.set_trace()
        _tbhdu_   = fits.BinTableHDU.from_columns(_coldefs_)
        # About the units
        _tbhdu_.header.append(('comment','Time is in hundred of microseconds (Hus) since 0 UT',''))
        _tbhdu_.header.append(('comment','ADCu = Analog to Digital Conversion units. Proportional to Voltage',''))
        _tbhdu_.header.append(('comment','mDeg = milli degree',''))
        _tbhdu_.header.append(('comment','Temperatures are in Celsius',''))

        _hduList_ = fits.HDUList([_hdu_,_tbhdu_])

        self.CleanPaths()

        try:
            _hduList_.writeto(self.OutputPath+self.MetaData['FITSfname'])
            return True
        except OSError as err:
            print("\n\n Write Error\n\n: {0}".format(err))
            return False

    """------------------------------------------------------------------------------------ """

    def concat(self,RBDlist) :
        """
        concat:
               A method to concatenate two or more RBD objects.
               The procedures doesn't check for the similarity of the objects.
               The user must provide two objects with the same Data structure.

               The procedure doesn't check neither the time sequence. The user
               must take care of it.

        use:
              import oRBD
              d1=oRBD.RBD()
              d1.readRBDinDictionary([RBD filename1])
              d2=oRBD.RBD()
              d2.readRBDinDictionary([RBD filename2])
              d=oRBD.RBD()
              d.concat([d1,d2])

        Change Record:
              First written by Guigue @ Sampa - 2017-09-08

        """


        self.RBDfname =[]
        self.Data = {}
        self.Metadata = {}
        self.header = RBDlist[0].header
        self.History = RBDlist[0].History

        ISOTime = []
        RBDFileName = []
        SSTType = RBDlist[0].MetaData['SSTType']
        ISODate = RBDlist[0].MetaData['ISODate']

        primDimension= 0
        TagList = list(RBDlist[0].MetaData.keys())
        for iRBD in RBDlist :
            primDimension += iRBD.Data['time'].shape[0]
            ISOTime.append(iRBD.MetaData['ISOTime'])
            RBDFileName.append(iRBD.MetaData['RBDFileName'])

        self.MetaData = { 'ISODate': ISODate,
                          'ISOTime': ISOTime,
                          'RBDFileName': RBDFileName,
                          'SSTType':SSTType}


        TagList = list(RBDlist[0].Data.keys())
        for iTag in TagList:
            secDimension=0
            if len(RBDlist[0].Data[iTag].shape) > 1 :
                secDimension = RBDlist[0].Data[iTag].shape[1]
            if secDimension == 0 :
                self.Data.update( { iTag: np.array(np.empty(primDimension,dtype=RBDlist[0].Data[iTag].dtype)) } )
            else:
                self.Data.update( { iTag: np.array(np.empty([primDimension,secDimension],dtype=RBDlist[0].Data[iTag].dtype)) } )

        iFirst=0
        for iRBD in RBDlist:
            iLast=iFirst+iRBD.Data['time'].shape[0]
            for iTag in TagList: self.Data[iTag][iFirst:iLast]=iRBD.Data[iTag]
            iFirst=iLast

        self.History.append('Concatenated Data')

        return


    """------------------------------------------------------------------------------------ """

    def reduced(self):
        """
        reduced:
             It produces a reduced form of the original data. It saves only the following columns in the
             table:

             time    : time in Hus
             azipos  : encoder's azimuth
             elepos  : encoder's elevation
             adc or adcval : receiver's output
             opmode  : oberving mode
             target  : target observed
             x_off   : scan offset in azimuth
             y_off   : scan offset in elevation

        use:
            import oRBD
            d=oRBD.RBD([RBD filename])
            d.readRBDinDictionary()
            d.reduced()

        Change Record:

            First written by Guigue @ Sampa - 2017-08-30

        """

        ListToPreserve = ['time','adc','adcval','elepos','azipos','opmode','target','x_off','y_off']
        ListToDelete = ['pm_daz', 'azierr', 'gps_status', 'pm_del', 'recnum', 'eleerr', 'pos_time', 'off']

        for iKey in ListToDelete:
            try:
                del self.Data[iKey]
            except:
                continue

        ChildToRemove=[]
        for iChild in self.header:
            if ListToPreserve.count(iChild.attrib['VarName']) == 0 : ChildToRemove.append(iChild)

        for iChildToRemove in ChildToRemove : self.header.remove(iChildToRemove)
        self.History.append('Reduced Data File. Selected Variables saved')

        del self.bin_header

        return

    """------------------------------------------------------------------------------------ """
    def CheckXMLTables(self):

        if  not os.path.exists(self.PathToXML+'SSTDataFormatTimeSpanTable.xml')  :
            print ('  ')
            print ('File : '+ self.PathToXML+'SSTDataFormatTimeSpanTable.xml not found')
            print ('Exiting...')
            return False

        if  not os.path.exists(self.PathToXML+'DataFormat-2002-12-14_to_2100-01-01.xml')  :
            print ('  ')
            print ('File : '+ self.PathToXML+'DataFormat-2002-12-14_to_2100-01-01.xml not found')
            print ('Exiting...')
            return False

        if  not os.path.exists(self.PathToXML+'DataFormat-2002-12-04_to_2002-12-13.xml')  :
            print ('  ')
            print ('File : '+ self.PathToXML+'DataFormat-2002-12-04_to_2002-12-13.xml not found')
            print ('Exiting...')
            return False

        if  not os.path.exists(self.PathToXML+'DataFormat-1999-05-02_to_2002-05-20.xml')  :
            print ('  ')
            print ('File : '+ self.PathToXML+'DataFormat-1999-05-02_to_2002-05-20.xml not found')
            print ('Exiting...')
            return False

        if  not os.path.exists(self.PathToXML+'DataFormat-1900-01-01_to_1999-05-01.xml')  :
            print ('  ')
            print ('File : '+ self.PathToXML+'DataFormat-1900-01-01_to_1999-05-01.xml not found')
            print ('Exiting...')
            return False

        if  not os.path.exists(self.PathToXML+'AuxiliaryDataFormat-2002-12-14_to_2100-01-01.xml')  :
            print ('  ')
            print ('File : '+ self.PathToXML+'AuxiliaryDataFormat-2002-12-14_to_2100-01-01.xml not found')
            print ('Exiting...')
            return False

        if  not os.path.exists(self.PathToXML+'AuxiliaryDataFormat-2002-11-24_to_2002-12-13.xml')  :
            print ('  ')
            print ('File : '+ self.PathToXML+'AuxiliaryDataFormat-2002-11-24_to_2002-12-13.xml not found')
            print ('Exiting...')
            return False

        if  not os.path.exists(self.PathToXML+'AuxiliaryDataFormat-2002-09-16_to_2002-11-23.xml')  :
            print ('  ')
            print ('File : '+ self.PathToXML+'AuxiliaryDataFormat-2002-09-16_to_2002-11-23.xml not found')
            print ('Exiting...')
            return False

        if  not os.path.exists(self.PathToXML+'AuxiliaryDataFormat-1900-01-01_to_2002-09-15.xml')  :
            print ('  ')
            print ('File : '+ self.PathToXML+'AuxiliaryDataFormat-1900-01-01_to_2002-09-15.xml not found')
            print ('Exiting...')
            return False

        return True

    """------------------------------------------------------------------------------------ """


######################################################################################

class DataTimeSpan(object):
    """

    DataTimeSpan

    This class has the unique effect of reading the SSTDataFormatTimeSpanTable
    and find the right RBD xml description file.

    Example of use:
       import oRBD
       tt=oRBD.DataTimeSpan()
       type='Auxiliary'
       date='2017-08-19'
       print (tt.findHeaderFile(SSTType=type,SSTDate=date))

    Change Record:

        Firts created by Guigue @ Sampa on 2017.08.19 (very cold indeed)

    """

    def __str__(self):
        return 'A class to read the SST Data Format Table and find the right xml description file.'

    def findHeaderFile(self,SSTType="Integration",SSTDate="1899-12-31"):
        """
        findHeaderFile:
        Given the DataType and the Date, look for the RBD xml description file.
        I guess there is a better and more efficient way to search in a b-tree, but
        I'm tired to find it ;-)

        Input:
            SSTType = "Subintegration" | "Integration" | "Auxiliary"
            SSTDate = YYYY-MM-DD

        Output:
            a string with the name of the RBD xml decription file.

        Change Record
            First created by Guigue @ Sampa on 2017-08-19
            Corrected on 2017-08-24

        """
        DataDescriptionFileName=''
        if (SSTType == 'Integration') or (SSTType == 'Subintegration') :
            filetype='Data'
        else:
            filetype='Auxiliary'

        for child in self.table:
            if (child[0].text == filetype) and (child[1].text <= SSTDate) and (child[2].text >= SSTDate) :
                DataDescriptionFileName=child[3].text
        return DataDescriptionFileName

    """------------------------------------------------------------------------"""

    def __init__(self,PathToXML):

        _tt_ = xmlet.parse(PathToXML+'SSTDataFormatTimeSpanTable.xml')
        self.table = _tt_.getroot()
        return
