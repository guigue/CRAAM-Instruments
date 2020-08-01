from pathlib import Path
import xml.etree.ElementTree as xmlet
import numpy as np
from astropy.io import fits
from utils import julday
import collections


def openP(path, name=None, path_to_xml=None, ms=None):

    """
    Function to open a POEMAS file and return an `POEMAS` object.
    Parameters
    ----------
    path : str, pathlib.Path, buffer
        File to be opened.
    name : str, optional
        Name of the POEMAS file. Only needed if path
        is a buffer.
    path_to_xml : str, pathlib.Path, optional
        Location of the POEMAS xml description files in the file system.
        If not defined it is assumed that the path is XMLTbles
        within the module's own directory.

    Raises
    ------
    FileNotFoundError
        If the POEMAS file was not found.

    ValueError
        If the path to the xml files is invalid.
    """

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

    return POEMAS().from_file(path, name, path_to_xml, ms)


class POEMAS(object):

    def __init__(self):
        self.filename = ""
        self.poemas_type = ""
        self.type = ""
        self.date = ""
        self.time = ""
        self.records = 0
        self.headerdata = np.empty((0))
        self.data = np.empty((0))
        self.new_data = np.empty((0))
        self.file_info = {}
        self.data = np.empty((0))
        self.history = list()

    def get_time_span(self):

        """
        Returns a tuple containing the ISO time of the
        first and last record found in the data.
        """

        nonzero = self.data["sec"].nonzero()

        return (julday.time(int(self.data["sec"][nonzero[0][0]])),
        julday.time(int(self.data["sec"][nonzero[0][-1]])))

    @property
    def headercolumns(self):
        """Returns the names of the header columns in a tuple."""

        return self.headerdata.dtype.names

    @property
    def columns(self):
        """Returns the names of the columns in a tuple."""

        return self.data.dtype.names

    def get_date(self):

        """
        Returns a string containing the ISO date and time of the
        first record found in the data.
        """

        nonzero = self.data["sec"].nonzero()

        date = str(julday.date(int(self.data["sec"][nonzero[0][0]])) + " " +
        julday.time(int(self.data["sec"][nonzero[0][0]])))

        return date


    def to_fits(self, name=None, output_path=None):

        """Writes the POEMAS data to a FITS file.
        By default the name of the fits file is defined as:
        poemas_[TRK |  | ]_YYYY-MM-DD_HH:MM:SS.SSSTHH:MM:SS.SSS_level0.fits
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

        t_start, t_end = self.get_time_span()


        if not name:
            name = "poemas_{}_{}_{}T{}_level0.fits".format(self.type.lower(), self.date, t_start, t_end)

        else:
            if not name.endswith(".fits"):
                name += ".fits"


        name = Path(name)

        if not output_path:
            output_path = "."

        output_path = Path(output_path).expanduser()

        if (output_path / name).exists():
            raise FileExistsError("File {} already exists.".format(str(name)))

        hdu = fits.PrimaryHDU()
        hdu.header.append(('origin', 'CRAAM/Universidade Presbiteriana Mackenzie', ''))
        hdu.header.append(('telescop', 'POEMAS - POlarization Emission of Millimeter Activity at the Sun', ''))
        hdu.header.append(('observat', 'CASLEO', ''))
        hdu.header.append(('station', 'Lat = -31.79897222, Lon = -69.29669444, Height = 2.491 km', ''))
        hdu.header.append(('tz', 'GMT-3', ''))

        hdu.header.append(('date-obs', self.date, ''))
        hdu.header.append(('t_start', self.date + 'T' + t_start,''))
        hdu.header.append(('t_end', self.date + 'T' + t_end, ''))
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
        hdu.header.append(("history", "Converted to FITS level-0 with poemas.py"))

        for hist in self.history:
            hdu.header.append(("history", hist))

        dscal = 1.0
        fits_cols = list()

        for column, values in self._header.items():

            var_dim = str(values[0])

            offset = 0
            if values[1] == np.int32:
                var_dim += "J"
            else:
                var_dim += "E"

            if self.poemas_type == "Subintegration":
                header_array = np.repeat(self.headerdata, self.records * 100 , axis=0)
            else:
                header_array = np.repeat(self.headerdata, self.records, axis=0)
            fits_cols.append(fits.Column(name=column,
                                         format=var_dim,
                                         unit=values[2],
                                         bscale=dscal,
                                         bzero=offset,
                                         array= header_array[column]))

        id = 0
        for column, values in self._newtblheader.items():

            var_dim = str(values[0])

            offset = 0
            if values[1] == np.int32:
                var_dim += "J"
            else:
                var_dim += "E"

            if id == 0 and self.poemas_type == "Subintegration":
                _name = 'msec'
            else:
                _name = column


            fits_cols.append(fits.Column(name=_name,
                                        format=var_dim,
                                        unit=values[2],
                                        bscale=dscal,
                                        bzero=offset,
                                        array=self.new_data[id]))
            id += 1


        tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(fits_cols))


        hdulist = fits.HDUList([hdu, tbhdu])

        hdulist.writeto(output_path / name)


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
        elif xml_type == "sub":
            xml = xmlet.parse(path_to_xml / Path("POEMASSubintegrationDataFormat.xml")).getroot()
        elif xml_type == "int":
            xml = xmlet.parse(path_to_xml / Path("POEMASIntegrationDataFormat.xml")).getroot()
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

    def subintegration(self):

        dt_array = [ [], [], [], [], [], [], [] ]

        for i in range (0,self.records):

            for j in range (0,4):

                if(j == 0):

                    msec = julday.msec(int(self.data[i][j]))
                    msec_array = [n for n in range(msec, msec+(100*10),10)]
                    dt_array[j].extend(msec_array)

                elif(j>0 and j<=2):

                    dt_array[j].extend([self.data[i][j]]*100)

                else:

                    tb_array = self.data[i][j].reshape(100,4)
                    tb_array = tb_array.transpose()

                    cont = 3
                    for tb in tb_array:
                        dt_array[cont].extend(tb)

                        cont+=1

        self.new_data = np.array(dt_array)

        return self

    def integration(self):

        dt_array = [ [], [], [],    [], [], [], [],    [], [], [], []]

        for i in range (0,self.records):

            for j in range (0,4):

                if(j == 0):
                    sec = julday.msec(int(self.data[i][j]))
                    dt_array[j].append(sec)


                elif(j>0 and j<=2):
                    dt_array[j].extend([self.data[i][j]])

                else:

                    tb_array = self.data[i][j].reshape(100,4)
                    tb_array = tb_array.transpose()

                    cont = 3
                    for tb in tb_array:
                        mean = np.mean(tb)
                        dt_array[cont].append(mean)
                        cont+=1

                    for tb in tb_array:
                        std = np.std(tb)
                        dt_array[cont].append(std)
                        cont+=1

        self.new_data = np.array(dt_array)

        return self

    def from_file(self, path, name, path_to_xml, ms):

        """Loads data from a file and returns an `POEMAS` object.
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

        self.filename = name

        if self.filename[-4:] != ".TRK":
            raise ValueError("Invalid file type {}".format(self.filename))
        else:
            self.type = "TRK"

        self._header = self.__find_header(path_to_xml,"head")
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
            self.data = np.frombuffer(path, dt_list)
        else:
            self.headerdata = np.fromfile(str(path), hdt_list, count = 1)
            fhandler = open(str(path),"rb")
            fhandler.seek(28,0)
            self.data = np.fromfile(str(path), dt_list)


        self.date, self.time = self.get_date().split(" ")
        self.records = self.headerdata[0][1]

        if ms == 0 or ms == None:
            self.poemas_type = "Subintegration"
            self._newtblheader = self.__find_header(path_to_xml,"sub")
            self.subintegration()
        elif ms == 1:
            self.poemas_type = "Integration"
            self._newtblheader = self.__find_header(path_to_xml,"int")
            self.integration()

        t_start, t_end = self.get_time_span()

        self.file_info.update({"Filename": name ,
                               "Date": self.date,
                               "Initial Time": t_start,
                               "Final Time":  t_end,
                               "POEMAStype": self.poemas_type})

        return self
