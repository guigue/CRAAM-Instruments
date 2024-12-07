import os
import datetime as dt
from pathlib import Path
import xml.etree.ElementTree as xmlet
import numpy as np
import collections
import astropy.constants as ctes
from astropy.io import fits
import pdb

#################################################
#
# POEMAS: A class to deal with the radiotelescope
# "POlarization Emission of Millimeter Activity at the Sun" (POEMAS)
#
#
# Use:
# > import poemas
# > p=poemas.Poemas()
# > p.read([poemas file name],subintegration=[True|False])
#      subintegration=True is the default, 10 ms data is read
#      subintegration=False integrates data inside a second 
#
# Methods:
#  to_fits()
#  extract()
#  getTimeAxis()
#
# Author:
#    @Guiguesp - from a version written by Julia (?)
#
# Version: 
__Version__ = '2024-12-07T17:30BST'
#
# Notes:
#    Presently only TRK files can be read.
#
####################################################

class Poemas(object):

    def __init__(self):
        self.filename = []
        self.type = ""
        self.data = {}
        self.history = []
        self.level = 'L0' 

    def __str__(self):
        return "A class for POEMAS data"
        
    @property
    def version(self):
        return __Version__
    
    @property
    def headercolumns(self):
        """Returns the names of the header columns in a tuple."""

        return self.headerdata.dtype.names

    @property
    def columns(self):
        """Returns the names of the columns in a tuple."""

        return self.data.dtype.names

    def read(self,path, subintegration=True, path_to_xml=None):

        if not path_to_xml:
        # __file__ is a python variable that stores where the module is located.
            path_to_xml = Path(__file__).parent / Path("XMLTables/")
        else:
            path_to_xml = Path(path_to_xml)

        if not isinstance(path, bytes):
            path = Path(path).expanduser()
            if not path.exists():
                raise FileNotFoundError("File not found: {}".format(path))
            name = path.name
        
        if not path_to_xml.exists():
            raise ValueError("Invalid path to XML: {}".format(path_to_xml))

        return self.from_file(path, name, path_to_xml, subintegration)

    def __find_header(self, path_to_xml, xml_type):

        """
        Method for finding the correct description file.
        Returns a dict representing the description found,
        the key is the variable name and the value is a list
        containing the var dimension, type and unit respectively.
        """

        if xml_type == "head":
            xml = xmlet.parse(path_to_xml / Path("POEMASDataFormatHead.xml")).getroot()
        elif xml_type == "tbl":
            xml = xmlet.parse(path_to_xml / Path("POEMASDataFormat.xml")).getroot()
        else:
            raise ValueError("Invalid xml type: {}".format(xml_type))

        header = dict()
        header = collections.OrderedDict()
        for child in xml:
            var_name = child[0].text
            var_dim = int(child[1].text)
            var_type = child[2].text
            var_unit = child[3].text

            if var_type == "xs:int":
                np_type = np.int32
            elif var_type == "xs:float":
                np_type = np.float32

            header.update({var_name:[var_dim, np_type, var_unit]})

        return header

    def tb2sfu(self,tb):
        # Convert from brightness temperature to flux
        
        d45 = 4.5E-01
        A45 = (d45/2.0)**2 * np.pi
        Aeff45 = 0.45
        Ae45 = A45 * Aeff45
          
        d90 = 1.65E-01
        A90 = (d90/2.0)**2 * np.pi
        Aeff90 = 0.71
        Ae90 = A90 * Aeff90

        W2sfu = 1.0E+22

        return (tb[:,0]+tb[:,1]) * ctes.k_B.value / Ae45 * W2sfu,(tb[:,2]+tb[:,3]) * ctes.k_B.value / Ae90 * W2sfu
    
    def from_file(self, path, name, path_to_xml, subintegration):

        """
        Read data from a file.
        Parameters
        ----------
            path : pathlib.Path
                Location of the POEMAS file in the file system.
            name : str
                Name of the POEMAS file.
            path_to_xml : Path, optional
                Location of the POEMAS xml description files in the file system.
        Raises
        ------
        ValueError
            If the filename is invalid.

        ------
            TRK only
        
        """
        
        if (not os.path.exists(str(path))) :
            print('\n\n File not found \n\n')
            return
        
        if name[-4:] != ".TRK":
            raise ValueError("Invalid file type {}".format(self.filename))
        else:
            self.type = "TRK"

        self.filename.append(name)
        self._header = self.__find_header(str(path_to_xml),"head")

        hdt_list = list()
        for key, value in self._header.items():
            hdt_list.append((key, value[1], value[0]))

        self._tblheader = self.__find_header(path_to_xml,"tbl")

        dt_list = list()
        for key, value in self._tblheader.items():
            dt_list.append((key, value[1], value[0]))

        if isinstance(path, bytes):
            self.headerdata = np.frombuffer(path, hdt_list, count = 1)
            fhandler = open(str(path),'rb')
            fhandler.seek(28,0)
            data = np.fromfile(fhandler, dt_list)
        else:
            self.headerdata = np.fromfile(str(path), hdt_list, count = 1)
            fhandler = open(str(path),"rb")
            fhandler.seek(28,0)
            data = np.fromfile(fhandler, dt_list)

        self.t_init  = (dt.timedelta(seconds=float(data['sec'][0]))  + dt.datetime(2001,1,1)).isoformat()
        self.t_end   = (dt.timedelta(seconds=float(data['sec'][-1])) + dt.datetime(2001,1,1)).isoformat()

        if subintegration:
            tb  = data['TB'].reshape(data.shape[0] * 100,4)
            ms  = np.empty(0)
            sms = [] 
            for s in data['sec']:
                t   = dt.timedelta(seconds=float(s)) + dt.datetime(2001,1,1)
                ms1 = int(t.second *1000 + t.minute * 60000 + t.hour * 3600000)
                ms  = np.concatenate((ms, ms1 + np.arange(100)*10))
                sms.append(ms1)
            sms = np.asarray(sms)
            ele = np.interp(ms,sms,data['ele_ang'])
            azi = np.interp(ms,sms,data['zi_ang'])
            flux45, flux90 = self.tb2sfu(tb)
            
        else:
            ms  = []
            ele = []
            azi = []
            tb  = np.empty([data['sec'].shape[0],4])
            TB  = data['TB'].reshape(data.shape[0] * 100,4)
            for i in np.arange(data['sec'].shape[0]):
                t   = dt.timedelta(seconds=float(data['sec'][i])) + dt.datetime(2001,1,1)
                ms.append(int(t.second *1000 + t.minute * 60000 + t.hour * 3600000))
                tb[i,0] = TB[i*100:(i+1)*100,0].mean()
                tb[i,1] = TB[i*100:(i+1)*100,1].mean()
                tb[i,2] = TB[i*100:(i+1)*100,2].mean()
                tb[i,3] = TB[i*100:(i+1)*100,3].mean()
                ele.append(data['ele_ang'][i])
                azi.append(data['zi_ang'][i])

            tb  = np.asarray(tb)
            ele = np.asarray(ele)
            azi = np.asarray(azi)
            ms  = np.asarray(ms)
            flux45, flux90 = self.tb2sfu(tb)
            
        self.data = {'ms':ms,
                     'tbl45':tb[:,0],
                     'tbr45':tb[:,1],
                     'tbl90':tb[:,2],
                     'tbr90':tb[:,3],
                     'flux45':flux45,
                     'flux90':flux90,
                     'azi':azi,
                     'ele':ele}

        return self
    
    def getTimeAxis(self):
        """

        getTimeAxis: Class method to convert the ms time axis used to Python
                     datetime ndarray that can be used with matplotlib.pyplot.

        Change Record:  Taken from SST class oRBD 
                        First written by Guigue @ Sampa
                        2017-11-04 St Charles Day !!!

        """
        ptime = []
        year  = int(self.t_init[:4])
        month = int(self.t_init[5:7])
        day   = int(self.t_init[8:10])
        
        for ms in self.data['ms']:
            hours        =  int(ms // 3600000)
            minutes      = int((ms % 3600000) // 60000)
            seconds      = ((ms % 3600000) % 60000) / 1.0E+03
            seconds_int  = int(seconds)
            seconds_frac = seconds - int(seconds)
            useconds     = int(seconds_frac * 1e6)
            ptime.append(dt.datetime(year,month,day,hours,minutes,seconds_int,useconds))

        return np.asarray(ptime)

    def extract(self,time_interval):
        _temp_ = Poemas()
        _temp_.type = self.type
        t0 = int( 1.0E+03 * (time_interval[0].hour * 3600 + time_interval[0].minute * 60 + (time_interval[0].second + time_interval[0].microsecond/1.0E+06)))
        t1 = int( 1.0E+03 * (time_interval[1].hour * 3600 + time_interval[1].minute * 60 + (time_interval[1].second + time_interval[1].microsecond/1.0E+06)))

        TagList = list(self.data.keys())
        for iTag in TagList:
            _temp_.data[iTag] = self.data[iTag][ (self.data['ms']>= t0 ) & ( self.data['ms'] <= t1)]

        _temp_.t_init  = time_interval[0].isoformat()
        _temp_.t_end   = time_interval[1].isoformat()
        _temp_.history.append('Extracted interval: {} - {}'.format(time_interval[0].isoformat(),time_interval[1].isoformat()))
        
        return _temp_

    def __add__(self,p):

        if (p.type != 'TRK'):
            return
        
        _temp_ = Poemas()
        _temp_.type = 'TRK'
        _temp_.t_init = self.t_init
        _temp_.t_end  = p.t_end
        _temp_.filename.append(p.filename[0])

        TagList = list(self.data.keys())
        for iTag in TagList:
            _temp_.data[iTag] = np.concatenate((self.data[iTag],p.data[iTag]))

        _temp_.history.append('Concatenated Data')

        return _temp_
        
    def add_history(self,history):

        self.history.append(history)
        return

    def change_level(self,level):

        self.level=level
        return
  
    def to_fits(self):

        """
        Write POEMAS data to a FITS file.
        By default the name of the fits file is defined as:
        poemas_[TRK |  | ]_YYYY-MM-DD_HH:MM:SS.SSSTHH:MM:SS.SSS_LLL.fits
        The file has two HDUs. The primary containing just a header with general
        information such as the origin, telescope, time zone. The second is a BinaryTable
        containing the data and a header with data specific information.

        TRK files only

        Parameters
        ----------
        name : str, optional
            Name of the fits file.
        
        output_path : str, pathlib.Path, optional
            Output path of the fits file. By default
            is where the script is being called from.
        
        Raises
        ------
        FileExistsError
            If a file with the same name already exists
            in the output path.
        """

        name = "poemas_{}_{}-{}_{}.fits".format(self.type.lower(),
                                                self.t_init.replace(':','_'),
                                                self.t_end[11:].replace(':','_'),
                                                self.level)
        name = Path(name)

        if (name).exists():
            raise FileExistsError("File {} already exists.".format(str(name)))

        hdu = fits.PrimaryHDU()
        hdu.header.append(('origin', 'CRAAM/Universidade Presbiteriana Mackenzie', ''))
        hdu.header.append(('telescop', 'POEMAS - POlarization Emission of Millimeter Activity at the Sun', ''))
        hdu.header.append(('observat', 'CASLEO', '')) 
        hdu.header.append(('station', 'Lat = -31.79897222, Lon = -69.29669444, Height = 2.491 km', ''))
        hdu.header.append(('tz', 'GMT-3', ''))

        hdu.header.append(('date-obs', self.t_init[:10], ''))
        hdu.header.append(('t_start', self.t_init,''))
        hdu.header.append(('t_end', self.t_end, ''))
        hdu.header.append(('data_typ', self.type, ''))
        if isinstance(self.filename, list) :
            for fname in self.filename: hdu.header.append(('origfile', fname, 'POEMAS Raw Binary Data file')) 
        else:
            hdu.header.append(('origfile',self.filename, 'POEMAS Raw Binary Data file')) 
            
        hdu.header.append(('frequen', '45 GHz ch=R,L; 90 GHz ch=R,L', ''))

        # About the Copyright
        hdu.header.append(('comment', 'COPYRIGHT. Grant of use.', ''))
        hdu.header.append(('comment', 'These data are property of Universidade Presbiteriana Mackenzie.'))
        hdu.header.append(('comment', 'The Centro de Radio Astronomia e Astrofisica Mackenzie is reponsible'))
        hdu.header.append(('comment', 'for their distribution. Grant of use permission is given for Academic ')) 
        hdu.header.append(('comment', 'purposes only.'))

        #History
        hdu.header.append(("history", "Flux density obtained from brightness temperatures")) #modificar        
        hdu.header.append(("history", "Converted to FITS level = {} with poemas.py".format(self.level))) #modificar

        for hist in self.history:
            hdu.header.append(("history", hist))

        dscal = 1.0
        fits_cols = list()
        offset = 0
 
        fits_cols.append(fits.Column(name='MS',format='1K',unit='ms',bscale=dscal,bzero=offset,array=self.data['ms']))
        fits_cols.append(fits.Column(name='TBL45',format='1E',unit='K',bscale=dscal,bzero=offset,array=self.data['tbl45']))
        fits_cols.append(fits.Column(name='TBR45',format='1E',unit='K',bscale=dscal,bzero=offset,array=self.data['tbr45']))
        fits_cols.append(fits.Column(name='TBL90',format='1E',unit='K',bscale=dscal,bzero=offset,array=self.data['tbl90']))
        fits_cols.append(fits.Column(name='TBR90',format='1E',unit='K',bscale=dscal,bzero=offset,array=self.data['tbr90']))
        fits_cols.append(fits.Column(name='FLUX45',format='1E',unit='K',bscale=dscal,bzero=offset,array=self.data['flux45']))
        fits_cols.append(fits.Column(name='FLUX90',format='1E',unit='K',bscale=dscal,bzero=offset,array=self.data['flux90']))
        fits_cols.append(fits.Column(name='AZI',format='1E',unit='Deg',bscale=dscal,bzero=offset,array=self.data['azi']))
        fits_cols.append(fits.Column(name='ELE',format='1E',unit='Deg',bscale=dscal,bzero=offset,array=self.data['ele']))
        
        tbhdu   = fits.BinTableHDU.from_columns(fits.ColDefs(fits_cols))
        # About the units
        tbhdu.header.append(('comment','Time is ms since 0 UT',''))
        hduList = fits.HDUList([hdu,tbhdu])

#        self.CleanPaths()

        try:
            hduList.writeto(name)
            return True
        except OSError as err:
            print("\n\n Write Error\n\n: {0}".format(err))
            return False

        return
    

