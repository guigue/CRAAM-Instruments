# External methods
import sys, string, os, struct, glob
import numpy as np
import datetime as dt
import pdb

from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, EarthLocation
from astropy.coordinates import get_body_barycentric, get_body, SkyCoord
from astropy import units as u
from astropy.coordinates import AltAz
from astropy.utils import iers

import CASLEO

##################################
Version = '20200702T1847BRT'     #
##################################

#######################################################
#
# Flux
#
# This class should have all the methods to get Flux
# from the SST. The first and simplest one, namely, the
# 'analytical Gaussian all-212GHz-beams-equal' is implemented
# here, extracted from the (quite) old sstpos.pro
# I did not try to make it more efficient, I only tried to
# reproduce identical results with a very well tested program.
#
# The method has a drawback, that one needs to fulfill part
# of the Bpos structure (a dictionary in python), the rest is implemented as
# a method that reads the beam_pos.asc file. But, frankly, this is
# a problem with SST data as well, which is not integrated! (Our fault)
#
# So, an example of use is as follows
#
# >>> bpos = {'hpbw212':4,'hpbw405':2,'eta212':0.35,'eta405':0.2}
# >>> f=mb.Flux(d,bpos)
# >>> f.Analytic_Method(d)
#
# The first line defines some properties of the beams. The _init__()
# method reads the beam_pos.asc file (you can change the filename in the
# method invocation.
#
# There is a flag temperature=True to get temperature instead of SFU.
#
# Written by: @guiguesp - 2020-04-28 under quarrantine... :-(
#
##########################################################################

class Flux(object):

    def version(self):
        return Version

    def __init__(self,d,bpos,bposfile='beam_pos.asc'):

        self.MetaData = d.MetaData
        self.MetaData.update({'BPos_filename':bposfile})
        self.Bpos     = bpos
        self.get_Bpos()

        self.History  = d.History

        return

    def get_Bpos(self):

        self.Bpos.update({'off': np.zeros(6),
                          'el' : np.zeros(6) } )

        if os.path.exists(self.MetaData['BPos_filename']) :
            _fd_         = open(self.MetaData['BPos_filename'],'r')
            for i in np.arange(6):
                l=_fd_.readline().split()
                self.Bpos['off'][i] = float(l[0])
                self.Bpos['el'][i]  = float(l[1])
            _fd_.close()
        else:
            print ('File ',self.MetaData['BPos_filename'],' not found')
            print ('Exit...\n\n')

        return

    def att(self,x,y):

    # Beam determined in November 10, 2003

        xc  = self.Bpos['off'][4]
        yc  = self.Bpos['el'][4]
        sx  = 2 / np.sqrt(np.log(64))
        sy  = 4 / np.sqrt(np.log(64))
        ang = np.radians(19.15)  #  0.334226

        xx  = x - xc
        yy  = y - yc
        xx1 = xx * np.cos(ang) - yy * np.sin(ang)
        yy1 = xx * np.sin(ang) + yy * np.cos(ang)
        d2  = xx1**2 / (2 * sx**2)  + yy1**2 / (2 * sy**2)

        return np.exp(d2)

    def calcular_xy(self,t0,t1,t2,x43,y43,x42,y42,b_a,x_off,y_off,delta):

        x = np.log(t2/t0) * y43 - np.log(t2 / t1) * y42
        y = np.log(t2/t1) * x42 - np.log(t2 / t0) * x43

        x = (x/b_a + x_off) / (2 * delta)
        y = (y/b_a + y_off) / (2 * delta)

        return x,y

    def to_datetime(self,husecs):

        iso   = self.MetaData['ISODate'].split('-')
        year  = int(iso[0])
        month = int(iso[1])
        day   = int(iso[2])
        ms    = husecs
        hours =  ms // 36000000
        minutes = (ms % 36000000) // 600000
        seconds = ((ms % 36000000) % 600000) / 1.0E+04
        seconds_int  = int(seconds)
        seconds_frac = seconds - int(seconds)
        useconds     = int(seconds_frac * 1e6)

        return dt.datetime(year,month,day,hours,minutes,seconds_int,useconds)

        
    def sst2xy(self,eqOff,Off_hour,pntb=5,quiet=True):

        from sstMap import Sun
        sun = Sun()
        
        cg = CASLEO.Observatory_Coordinates()
        
        if (pntb < 1) | (pntb > 6):
            pntch = 4
        else:
            pntch=pntb-1
        if ~quiet:
            print("\n\n Pointing Beam = {0:3d} \n".format(pntb))

        time0=dt.datetime(self.MBSol['dtime'].year,self.MBSol['dtime'].month,self.MBSol['dtime'].day,Off_hour)
        
        Sradec0,Saltaz0,L0,B0,P,R_Sun,Carrington,North_Orientation = sun.Get_Sun_Coordinates(cg,time0)
        self.MetaData.update({'Sun_Coordinates': {'RaDec':radec,
                                                  'AltAz': altaz,
                                                  'L0': L0, 'B0': B0,
                                                  'P':P, 'R_Sun': R_Sun,
                                                  'Carrington': Carrington,
                                                  'North_Orientation':-North_Orientation},
                              'Observatory': cg})
        
        self.MetaData.update({'Parallactic_Angle':self.Compute_Parallactic_Angle()})
        daz = self.MBSol['off'] - self.Bpos['off'][pntch]
        d_el = self.MBSol['el']  - self.Bpos['el'][pntch]
        dalpha = daz * np.cos(self.MetaData['Parallactic_Angle']) - d_el * np.sin(self.MetaData['Parallactic_Angle'])
        ddec   = daz * np.sin(self.MetaData['Parallactic_Angle']) + d_el * np.cos(self.MetaData['Parallactic_Angle'])

        self.MBSol.update({'ra':self.MetaData['RaDec'][0]+eqOff[0]+dalpha})
        self.MBSol.update({'dec':self.MetaData['RaDec'][1]+eqOff[1]+ddec})
        
        return

    def Compute_Parallactic_Angle(self):
    ########################
    #  Compute the parallactic angle.
    #
    #    Formula taken from A. Guyonnet.
    #    eta = sin(ha) / ( cos(dec) * tan(latitude) - sin(dec) * cos(ha) )
    #
    #######################

        latitude = self.MetaData['Observatory'].lat.rad
        ha = np.radians(self.MetaData['Sun_Coordinates']['RaDec'].ra.hour)
        dec = self.MetaData['Sun_Coordinates']['RaDec'].dec.rad
        parallactic_angle = np.degrees(np.arctan(np.sin(ha) / (np.cos(dec) * np.tan(latitude) - np.sin(dec) * np.cos(ha))))
        return parallactic_angle

######################################################################################    
    def Analytic_Method(self, data, limite=1000, eliptic=False, temperature=False):

# Original comments from sstpos.pro follow
#
# SSTPOS
#
# Compute positions and flux assuming a point-like source
#
# CALLING SEQUENCE
#   sstpos,data,bpos,res,[limite=limite]
#
# INPUT
#   data:   structure with antenna temperatures (see sst_read)
#   bpos:   structure with antenna beams position
#
# OUTPUT
#   res:    structure with positions, flux and more...
#          time:  teh same as in data
#          flux212 [atemp212] : 212 GHz flux density in s.f.u. if /temperature set, in K
#          flux405 [atemp405] : 405 GHz flux density in s.f.u. if /temperature set, in K
#          off: azimuth source positions in arcmin of the antenna coordinates system
#          el : elevation source positions in arcmin of the antenna coordinates system
#
# OPTIONAL PARAMETERS
#   limite: A tricky parameter. It tries to bind the solutions
#           within this limit.  For example, limite=0.1 will not
#           allow position solutions to have a separation of more
#           than sqrt(2) * 0.1 arcmin
#   /eliptico: If set, consider an elliptical beam 5.
#   /temperature: If set, give the result in the same Units of d.adcval
#              structure
#
# HISTORY
#   First written by Guigue in the early days of the SST (1999.999999)
#   Last written by Guigue, of course, 2003  :-)
#   Added /temperature keyword.
#   Last version May/2011 (CASLEO)
#----------------------------------------------------------------------------------
#   Translated to python on 2020-04-26 at home because of the CoVid-19 quarrantine
#
#

# coefficient b, same for all the beams
        b_a212 = np.log(16) / self.Bpos['hpbw212']**2	#assumed antenna b
        b_a405 = np.log(16) / self.Bpos['hpbw405']**2	#assumed antenna b

# coeff_bolt=S(flux)/2k so with S in sfu (10^22 W/m-2..)
        Bcoef = 2 * 0.138
        tp    = data.Data['time'].shape[0]
        Output = 'FluxDensity'

# Computed values
        tp = data.Data['time'].shape[0]
        if (temperature):
            Output = 'Temperature'
            k212 = Bcoef
            k405 = Bcoef
            res  = {'off'        : np.array(np.empty([tp],float)),
                    'el'         : np.array(np.empty([tp],float)),
                    'AntTemp212' : np.array(np.empty([tp],float)),
                    'AntTemp405' : np.array(np.empty([tp],float)),
                    'uloff'      : np.array(np.empty([tp],float)),
                    'ulel'       : np.array(np.empty([tp],int))  ,
                    'toz'        : np.array(np.empty([tp],float)),
                    'time'       : data.Data['time']             ,
                    'dtime'      : np.array(np.empty([tp],dtype='datetime64[us]'))
            }

        else:
            k212 = (1.5/2)**2 * np.pi * self.Bpos['eta212']
            k405 = (1.5/2)**2 * np.pi * self.Bpos['eta405']
            res  = {'off'     : np.array(np.empty([tp],float)) ,
                    'el'      : np.array(np.empty([tp],float)) ,
                    'Flux212' : np.array(np.empty([tp],float)) ,
                    'Flux405' : np.array(np.empty([tp],float)) ,
                    'uloff'   : np.array(np.empty([tp],float)) ,
                    'ulel'    : np.array(np.empty([tp],int))   ,
                    'toz'     : np.array(np.empty([tp],float)) ,
                    'time'    : data.Data['time']              ,
                    'dtime'   : np.array(np.empty([tp],dtype='datetime64[us]'))
            }

        x42 = self.Bpos['off'][3] - self.Bpos['off'][1] # az(4) - az(2)
        x43 = self.Bpos['off'][3] - self.Bpos['off'][2] # az(4) - az(3)

        y42 = self.Bpos['el'][3] - self.Bpos['el'][1] #el(3) - el(1)
        y43 = self.Bpos['el'][3] - self.Bpos['el'][2] #el(3) - el(2)

        r42 = self.Bpos['off'][3]**2 - self.Bpos['off'][1]**2 + self.Bpos['el'][3]**2  - self.Bpos['el'][1]**2
        r43 = self.Bpos['off'][3]**2 - self.Bpos['off'][2]**2 + self.Bpos['el'][3]**2  - self.Bpos['el'][2]**2

        delta = x42 * y43 - x43 * y42
        x_off = r42 * y43 - r43 * y42
        y_off = r43 * x42 - r42 * x43

        for ii in np.arange(tp):

            t0 = data.Data['adcval'][ii,1]
            t1 = data.Data['adcval'][ii,2]
            t2 = data.Data['adcval'][ii,3]
            t405 = data.Data['adcval'][ii,4]

            if ((t0 <= 0.0) or (t1 <= 0.0) or (t2 <= 0.0) ):

                if (ii != 0):
                    res['off'][ii] = res['off'][ii-1]
                    res['el'][ii]  = res['el'][ii-1]
                else:
                    jj = ii
                    while (t0 <= 0.0):
                        jj = jj + 1
                        t0 = data.Data['adcval'][jj,1]
                        res['toz'][ii] = 0.0

                    jj = ii
                    while (t1 <= 0.0):
                        jj = jj + 1
                        t1 = data.Data['adcval'][jj,2]
                        res['toz'][ii] = 0.0
                    jj = ii
                    while (t2 <= 0.0):
                        jj = jj + 1
                        t2 = data.Data['adcval'][jj,3]
                        res['toz'][ii] = 0.0

            res['off'][ii], res['el'][ii] = self.calcular_xy(t0,t1,t2,x43,y43,x42,y42,b_a212,x_off,y_off,delta)

            if (ii != 0):

                if (np.abs(res['off'][ii] - res['off'][ii-1]) > limite):
                    res['off'][ii]  = res['off'][ii-1] + limite * np.sign(res['off'][ii]-res['off'][ii-1])
                    res[ii].uloff = 0.0

                if (np.abs(res['el'][ii] - res['el'][ii-1]) > limite):
                    res['el'][ii] = res['el'][ii-1] + limite * np.sign(res['el'][ii]-res['el'][ii-1])
                    res['ulel'][ii] = 0.0

            c212 = t0 * Bcoef / k212 * np.exp(b_a212 *( (res['off'][ii] - self.Bpos['off'][1])**2 + (res['el'][ii]  - self.Bpos['el'][1] )**2 ) )

            if (eliptic):
                c405 = t405 *  Bcoef / k405 * att(res['off'][ii],res['el'][ii])

            else:
                c405 = t405 *  Bcoef / k405 * np.exp(b_a405*( (res['off'][ii] - self.Bpos['off'][4])**2 + (res['el'][ii]  - self.Bpos['el'][4] )**2 ) )

            if (temperature):
                res['AntTemp212'][ii] = c212
                res['AntTemp405'][ii] = c405
            else:
                res['Flux212'][ii] = c212
                res['Flux405'][ii] = c405

            res['dtime'][ii] = self.to_datetime(res['time'][ii])

        self.MetaData.update({'Version': self.version(),
                              'Method' : 'Analytical',
                              'Output' : Output        }
        )

        self.History.append('Solved the Multibeam System')
        self.MBSol = res

        return

    def writeFITS(self,FITSfname):

        """
        writeFITS:
             A method to write the SST MB Solution in FITS format as a binary table.
             The file has a primary header and a table or secondary header.

             The system implements two headers. The primary header has general information,
             while the secondary header is specific for the table, including the units of the columns.

             Method taken from oRBD.py

        Output:
             It returns True or False on success or failure, respectively.

        Change Record:
             First written by Guigue @ Sampa - 2017-08-26
             Return value added on 2017-11-02
             Adapted to MB Solution on 2020-04-29

        """

        self.MetaData.update({'FITSfname':FITSfname})

        _isodate_ = self.MetaData['ISODate']
        _hhmmss_  = [self.MBSol['dtime'][0].astype('str'),self.MBSol['dtime'][-1].astype('str')]
        _hdu_     = fits.PrimaryHDU()

        #
        # This is the Primary (global) header. It gives information about the instrument, and the data.
        #
        _hdu_.header.append(('origin','CRAAM/Universidade Presbiteriana Mackenzie',''))
        _hdu_.header.append(('telescop','Solar Submillimeter Telescope',''))
        _hdu_.header.append(('observat','CASLEO',''))
        _hdu_.header.append(('station','Lat = -31.79897222, Lon = -69.29669444, Height = 2.491 km',''))
        _hdu_.header.append(('tz','GMT-3',''))

        _hdu_.header.append(('date-obs',_isodate_,''))
        _hdu_.header.append(('t_start',_hhmmss_[0],''))
        _hdu_.header.append(('t_end',_hhmmss_[1],''))
        _hdu_.header.append(('data_typ',self.MetaData['Output'],''))
        if isinstance(self.MetaData['RBDFileName'],list) :
            for iRBD in self.MetaData['RBDFileName']:_hdu_.header.append(('origfile',iRBD,'SST Raw Binary Data file'))
        else:
            _hdu_.header.append(('origfile',self.MetaData['RBDFileName'],'SST Raw Binary Data file'))

        _hdu_.header.append(('frequen','212 GHz; 405 GHz',''))

        # About the Copyright
        _hdu_.header.append(('comment','COPYRIGHT. Grant of use.',''))
        _hdu_.header.append(('comment','These data are property of Universidade Presbiteriana Mackenzie.'))
        _hdu_.header.append(('comment','The Centro de Radio Astronomia e Astrofisica Mackenzie is reponsible'))
        _hdu_.header.append(('comment','for their distribution. Grant of use permission is given for Academic '))
        _hdu_.header.append(('comment','purposes only.'))

        for i in range(len(self.History)):
            _hdu_.header.append(('history',self.History[i]))

        # Secondary Header (Table)
        _fits_cols_ = [  fits.Column(name   = 'OFF'   ,
                                     format = '1E'    ,
                                     unit   = 'arc minutes' ,
                                     bscale = 1.0     ,
                                     bzero  = 0       ,
                                     array  = self.MBSol['off']),
                         fits.Column(name   = 'EL'    ,
                                     format = '1E'    ,
                                     unit   = 'arc minutes' ,
                                     bscale = 1.0     ,
                                     bzero  = 0       ,
                                     array  = self.MBSol['off']),
                         fits.Column(name   = 'ULOFF' ,
                                     format = '1E'    ,
                                     unit   = ' '     ,
                                     bscale = 1.0     ,
                                     bzero  = 0       ,
                                     array  = self.MBSol['uloff']),
                         fits.Column(name   = 'ULEL'  ,
                                     format = '1E'    ,
                                     unit   = ' '     ,
                                     bscale = 1.0     ,
                                     bzero  = 0       ,
                                     array  = self.MBSol['ulel']),
                         fits.Column(name   = 'TOZ'   ,
                                     format = '1E'    ,
                                     unit   = ' '     ,
                                     bscale = 1.0     ,
                                     bzero  = 0       ,
                                     array  = self.MBSol['toz']),
                         fits.Column(name   = 'TIME'  ,
                                     format = '1J'    ,
                                     unit   = 'ms'    ,
                                     bscale = 1.0     ,
                                     bzero  = 0       ,
                                     array  = self.MBSol['time'] // 10   ) # Convert husecs to milliseconds
                         ]

        if (self.MetaData['Output'] == 'Temperature'):
            _fits_cols_.append( fits.Column(name    = 'TEMP212',
                                            format  = '1E'     ,
                                            unit    = ' K '    ,
                                            bscale  = 1.0      ,
                                            bzero   = 0.0      ,
                                            array   = self.MBSol['AntTemp212']))
            _fits_cols_.append( fits.Column(name    = 'TEMP405',
                                            format  = '1E'     ,
                                            unit    = ' K '    ,
                                            bscale  = 1.0      ,
                                            bzero   = 0.0      ,
                                            array   = self.MBSol['AntTemp405']))

        if (self.MetaData['Output'] == 'FluxDensity'):
            _fits_cols_.append( fits.Column(name    = 'FLUX212',
                                            format  = '1E'     ,
                                            unit    = 'SFU'    ,
                                            bscale  = 1.0      ,
                                            bzero   = 0.0      ,
                                            array   = self.MBSol['Flux212']))
            _fits_cols_.append( fits.Column(name    = 'FLUX405',
                                            format  = '1E'     ,
                                            unit    = 'SFU'    ,
                                            bscale  = 1.0      ,
                                            bzero   = 0.0      ,
                                            array   = self.MBSol['Flux405']))

        _coldefs_ = fits.ColDefs(_fits_cols_)
        _tbhdu_   = fits.BinTableHDU.from_columns(_coldefs_)
        # About the units
        _tbhdu_.header.append(('comment','Time is milliseconds from 0 UT',''))
        _hduList_ = fits.HDUList([_hdu_,_tbhdu_])

        if os.path.exists(self.MetaData['FITSfname']) :
            print ('File '+ self.MetaData['FITSfname']+ '  already exist. Aborting....')
            return False
        else:
            _hduList_.writeto(self.MetaData['FITSfname'])

        return True
