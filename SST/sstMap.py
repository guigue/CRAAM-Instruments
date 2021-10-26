import numpy as np
import matplotlib.pyplot as plt
import pdb
import scipy.interpolate as intl
from scipy.ndimage import rotate

import statistics as stat
import datetime as dt

from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, EarthLocation
from astropy.coordinates import get_body_barycentric, get_body
from astropy import units as u
from astropy.coordinates import AltAz
from astropy.utils import iers
from astropy.io import fits
from sunpy import coordinates as Sun_Coordinates
#
# This is needed because IERS data repository is not
# accepting downloads. The precision is lower, but
# it's enough for our computations
#
iers.conf.auto_download = False

# Local Packages
import CASLEO
from CraamTools import contiguo as c
from CraamTools.fit import Circle
from CraamTools.AstroTools import julian
from CraamTools.AstroTools import get_pb0
from CraamTools.filters.RunningMean import rm1d

_Version = '20200801T1654BRT'

#=======================================================================================
#
#  sstMap
#
#  Author:  @guiguesp - during the SARS-CoV-2 outbreak dark days of 2020
#                     - Major Changes on 2020-05-29
#
#  sstMap intends to be a python class to extract, calibrate, correct, and fit
#  sst maps.  It uses an oRBD object as input. The basic class is Scans
#  that extract every scan on a channel data and creates an object. Since the
#  scans have different length, Scans determines the "mode" lenght and then interpolates
#  all the scans to this new length. It does the same with the x_off.
#  Map is the other class, that makes an array with the scans and
#  interpolates the final image to a squared box.
#
#  More detailed information within the classes and methods descriptions.
#
#  sstMap now is more "robust" in terms of data security. It uses astropy and
#  sunpy (1.1) to get astronomical coordinates. It has a separated CASLEO set of functions.
#
#  Usage :
#    import oRBD
#    d = oRBD.RBD()
#    d.readRBDinDictionary(<SST RBD file name>)
#    from CraamTools.AstroTools import sstMap
#    m0=sstMap.sstMap(d,<channel number: 0...5>)
#    m0.Scans.Extract_Scans()
#    m0.Create_Image(image_size=<image size>)
#    m0.Rotate_North_Image()
#    m0.Calibrate_Map_in_T()
#    m0.Write_FITS(<fits file name>)
#
#  To correct the zig-zag in maps use
#    m0.Scans.Correct_Scans(eshift=<val1>, oshift=<val2>)
#
#=====================================================================================

class Scans(object):
    #
    #
    # Scans
    #
    #   This class has most of the original data from the RBD object.
    #   It has also properties derived from the data, like the
    #   Sun center and radius which may be used for further analysis.
    #
    #   Methods
    #       Get_Bpos
    #       Select_Data
    #       Extract_Scans
    #       Correct_Scans
    #       Limb_fit
    #       Get_Version
    #
    #_____________________________________________________________________________

    def __str__(self):
        return 'A Class to reduce raw scans.'

    @property
    def version(self):
        return self._version

    def Select_Data(self,direction='az'):
        #
        #
        # Select Data:
        #
        #   This procedure makes indices to locate the scans data.
        #   Since we implemented only AzEl maps, other types of maps
        #   are not extracted.
        #
        #   A drawback of the method is that if more than one map is
        #   present in the data, it cannot separate them. The user
        #   must supply an RBD object containing only one map.
        #
        #   @guiguesp - 2020-04-05
        #
        #_____________________________________________________________________________

        if (direction == 'az'):
            iMap, = np.where(self.d.Data['opmode'] == 2)
            cMap  = c.contiguo(iMap)
        else:
            print('\n Only Azimuth - Elevation Maps Implemented yet\n')
            return

        iInt, = np.where(self.d.Data['opmode'] == 4)
        cInt  = c.contiguo(iInt)

        self.direction = 'az'
        self.Indices = {'cMap':cMap,'iMap':iMap, 'iInt':iInt, 'cInt':cInt}
        return

    def Get_Bpos(self):
        #
        #
        # Get_Bpos
        #
        #   This method returns a dict with the beam positions.
        #   This implementeation is for 2017-09 only
        #   A general procedure must be implemented.
        #
        #   @guiguesp 2020-04-05
        #
        #_____________________________________________________________________________

        off = np.array([ 14.2,  16.27, 12.34,  14.47,  14.75, 14.08])
        ele = np.array([-6.03, -12.2, -12.04, -15.93, -13.69, -5.92])
        return {'OFF':off,'ELE':ele}

    def Correct_Scans(self,oshift=0,eshift=0):
        #
        #
        # Correct_Scans
        #
        #    This method "aligns" Scans. It was
        #    experimentally observed that scans are misalagned bewteen them,
        #    the reason is still unknow, however I suspect that the serial
        #    connection between computer and Orbit Mount is the cause.
        #    Anyway, following discussions with JFValle Silva, my idea
        #    is to plot the central scans over the Sun in function of the
        #    abs(x_off). Then, apply a shift (numpy.roll) until the limbs
        #    (left and right) coincides. This process must be done
        #    before applying this method, and inspected by eye. An automatic
        #    implementation can be done, in the future.
        #
        #    It was observed that "odd" scans have generally different shifts
        #    than "even" scans. And after applying the correction the center
        #    of the solar disc is displaced. Thereforebefore shifts we determine the
        #    original solar disc center, then apply the scans corrections
        #    and then move the solar disc center to the original position.
        #    For this reason, the method calls many times limb_fit()
        #
        #    @guiguesp - 2020-04-05
        #
        #_____________________________________________________________________________

        self.Limb_fit()
        x0_orig = self.limb['X0']

        for iscan in np.arange(self.Ny):
            if (np.mod(iscan,2)):
                self.T[iscan] = np.roll(self.T[iscan],oshift)
            else:
                self.T[iscan] = np.roll(self.T[iscan],eshift)

        self.Limb_fit()
        dx = (self.x_off[self.Ny//2,1:]-self.x_off[self.Ny//2,0:-1]).mean()

        delta = int((x0_orig - self.limb['X0']) / dx)
        for iscan in np.arange(self.Ny):
            self.T[iscan] = np.roll(self.T[iscan],delta)

        self.Limb_fit()
        return

    def Extract_Scans(self,direction='az'):
        #
        #
        # Extract_Scans
        #
        #   This method (1) separate the different scans, (2) reverse the order
        #   for the "odd" scans (3) subtracts a base line, making all scans
        #   to have a base=0 (4) interpolates all scans to have the same
        #   length (5) creates a matrix of uncalibrated temperatures, x_off and
        #   y_off.
        #
        #   The base line is obtained by fitting a linear polynomium to the
        #   mean values of the intermediate measurements. Intermediate readouts
        #   are coded with an opmode=4
        #
        #   The lenght of the final scans is the same of the mode of the lengths of the
        #   original scans. The positions and uncalibrated temperatures are interpolated.
        #
        #   If more than one map is found in the data, the system will extract all, mixing
        #   the maps. Watch up!
        #
        #   @guiguesp - 2020-04-05
        #
        #_______________________________________________________________________________

        _MD2DEG_ = 1.0e-03
        _MD2SEC_ = 3.6e+00

        scans = []
        x_offs = []
        y_offs = []
        azi    = []
        ele    = []

        scanL = []

        cMap = self.Indices['cMap']
        cInt = self.Indices['cInt']
        iMap = self.Indices['iMap']
        iInt = self.Indices['iInt']

        N = len(cMap)
        for iscan in np.arange(N):
            scan = self.d.Data['adcval'][iMap[cMap[iscan,0]]:iMap[cMap[iscan,1]],self.channel]
            xoff = self.d.Data['x_off'][iMap[cMap[iscan,0]]:iMap[cMap[iscan,1]]] * _MD2SEC_
            xoff = np.linspace(xoff[0],xoff[-1],xoff.shape[0])
            yoff = self.d.Data['y_off'][iMap[cMap[iscan,0]]:iMap[cMap[iscan,1]]] * _MD2SEC_
            yoff = np.linspace(yoff[0],yoff[-1],yoff.shape[0])
            sazi = self.d.Data['azipos'][iMap[cMap[iscan,0]]:iMap[cMap[iscan,1]]] * _MD2DEG_
            sele = self.d.Data['elepos'][iMap[cMap[iscan,0]]:iMap[cMap[iscan,1]]] * _MD2DEG_

            scan.astype(float)
            xoff.astype(float)
            yoff.astype(float)

            if (np.mod(iscan,2)):
                scan = np.flip(scan)
                xoff = np.flip(xoff)
                yoff = np.flip(yoff)
                sazi = np.flip(sazi)
                sele = np.flip(sele)

            # Background
            b0 = float(self.d.Data['adcval'][iInt[cInt[iscan,0]],self.channel])
            b1 = float(self.d.Data['adcval'][iInt[cInt[iscan,1]],self.channel])
            b  = (b1-b0)/(cMap[iscan,1]-cMap[iscan,0])
            a  = b0 - b * iMap[cMap[iscan,0]]
            x  = np.arange(iMap[cMap[iscan,0]],iMap[cMap[iscan,1]])
            y  = a + b * x

            scan = scan - y
            scanL.append(len(scan))
            x_offs.append(xoff)
            y_offs.append(yoff)
            azi.append(sazi)
            ele.append(sele)
            scans.append(scan)

        scanF   = stat.mode(scanL)
        xf      = np.arange(scanF)
        scansF  = []
        xoffF   = []
        yoffF   = []
        sazi    = []
        sele    = []

        for iscan in np.arange(N):
            x = np.arange(len(scans[iscan]))
            s = np.interp(xf,x,scans[iscan])
            scansF.append(s)

            xo = np.interp(xf,x,x_offs[iscan])
            xoffF.append(xo)

            yo = np.interp(xf,x,y_offs[iscan])
            yoffF.append(yo)

            s = np.interp(xf,x,azi[iscan])
            sazi.append(s)

            s = np.interp(xf,x,ele[iscan])
            sele.append(s)

        self.x_off  = np.asarray(xoffF)
        self.y_off  = np.asarray(yoffF)
        self.azi    = np.asarray(sazi)
        self.ele    = np.asarray(sele)
        self.T      = np.asarray(scansF)
        self.Nx     = scanF
        self.Ny     = N
        self.direction = direction

        return

    def Limb_fit(self):
        #
        #
        # limb_fit
        #
        # The function limb_fit fits the solar limb from the Scans.
        # It uses the x_off and y_off ccordinates (which are in arc seconds)
        # Afterward it finds the center in Alt Azimuth coordinates.
        #
        # This function is adpated from the original JERCosta limb_fit.pro
        # later used in my tpoint.pro. From this I obatined the definition
        # of "observed" coordinates. This definition was discussed with P. Wallace
        # wee back in 2006/11
        #
        # 2020-04-01 - Guigue
        #
        #_______________________________________________________________________________

        _SEC2DEG_ = 1.0 / 3.6e+03

        mo = np.zeros([self.Ny, self.Nx],float)
        for i in np.arange(self.Nx):
            for j in np.arange(self.Ny):
                mo[j,i] = self.T[j][i]

        smooth=0
        if self.channel > 3 :
            smooth=10
        x,y       = Circle.get_limb(self.x_off,self.y_off,mo,smooth=smooth)
        par,cov   = Circle.fit(x,y)
        self.limb = {'off':x,'el':y,'X0':par[0],'Y0':par[1],'R':par[2],'Cov':cov}
        self.image = mo
        return

    def Get_Version(self):
        return  '20200529T1210'

    def __init__(self,d,channel):
        self._version = _Version
        self.channel = channel
        self.d = d
        return

##########################################################################################################################

class Map(object):

    #
    # Map Class
    #
    #    This class has the method needed to create a square image from
    #    the scans extracted with the Scans methods.
    #
    #    Methods
    #      Create_Image
    #      Get_Version
    #      Rotate_North_Image
    #
    #_____________________________________________________________________________


    def __init__(self):
        self.version = _Version
        return

    def __str__(self):
        return 'A Class with fundamental Map operations.'

    @property
    def version(self):
        return self._version

    def Get_Version(self):
        return  '20200530T1819'

    def Get_Observed_Sun_Center(self):
        ##################
        # Get Alt Azi coordinated for the "observed" sun center
        #
        _SEC2DEG_ = 1.0 / 3.6e+03

        if not hasattr(self.Scans,'limb'):
            self.Scans.Limb_fit()

        obsazi = self.MetaData['Sun_Coordinates']['AltAz'].az.degree  - self.Scans.limb['X0'] * _SEC2DEG_ / np.cos(self.MetaData['Sun_Coordinates']['AltAz'].alt.radian)
        obsele = self.MetaData['Sun_Coordinates']['AltAz'].alt.degree - self.Scans.limb['Y0'] * _SEC2DEG_
        self.MetaData.update({'Observed_Sun_Center':{'Azi' : obsazi, 'Ele' : obsele}})

        return

    def Create_Image(self,image_size=512):

        #
        # Create_Image
        #
        #   This is the most important method of the class. It takes the
        #   Scans and converts them to a square matrix, with new axes.
        #   It determines the solar disc center in "matrix coordinates"
        #
        #   @guiguesp - 2020-04-05
        #             - 2020-05-30 : Changed the method to determine the Sun Center in Matrix Units
        #
        #______________________________________________________________________________

        _SEC2DEG_ = 1.0 / 3.6e+03

        self.Get_Observed_Sun_Center()

        mo = np.zeros([self.Scans.Ny, self.Scans.Nx],float)
        for i in np.arange(self.Scans.Nx):
            for j in np.arange(self.Scans.Ny):
                mo[j,i] = self.Scans.T[j][i]

        mf = intl.interp2d(np.arange(self.Scans.Nx),np.arange(self.Scans.Ny),mo)
        P  = image_size
        x  = np.arange(0, self.Scans.Nx, self.Scans.Nx / P)
        y  = np.arange(0, self.Scans.Ny, self.Scans.Ny / P)
        self.image = mf(x,y)

        x     = np.linspace(0,self.Scans.Nx,P)
        x_off = []
        y_off = []
        azi   = []
        ele   = []

        # This are mean x_off and y_off created from the original ones.
        pixel_size_x = np.abs(self.Scans.x_off[self.Scans.Ny//2,-1] - self.Scans.x_off[self.Scans.Ny//2,0])/P

        x_off_model = np.linspace(
            self.Scans.x_off[self.Scans.Ny//2,0],
            self.Scans.x_off[self.Scans.Ny//2,-1],
            P)

        pixel_size_y =  np.abs(self.Scans.y_off[-1,self.Scans.Nx//2] - self.Scans.y_off[0,self.Scans.Nx//2])/P
        y_off_model = np.linspace(
            self.Scans.y_off[0,self.Scans.Nx//2],
            self.Scans.y_off[-1,self.Scans.Nx//2],
            P)

        azi_model = np.linspace(self.MetaData['Observed_Sun_Center']['Azi'] - self.Scans.x_off[self.Scans.Ny//2,0] * _SEC2DEG_ ,
                                self.MetaData['Observed_Sun_Center']['Azi'] + self.Scans.x_off[self.Scans.Ny//2,-1] * _SEC2DEG_, P)

        ele_model = np.linspace(self.MetaData['Observed_Sun_Center']['Ele'] - self.Scans.x_off[0,self.Scans.Nx//2] * _SEC2DEG_ ,
                                self.MetaData['Observed_Sun_Center']['Ele'] + self.Scans.x_off[-1,self.Scans.Ny//2] * _SEC2DEG_, P)

        for iscan in np.arange(P):
            x_off.append(x_off_model)
            y_off.append(y_off_model)
            azi.append(azi_model)
            ele.append(ele_model)


        self.x_off = np.asarray(x_off)
        self.y_off = (np.asarray(y_off)).transpose()
        self.azi   = np.asarray(azi)
        self.ele   = (np.asarray(ele)).transpose()
        self.MetaData.update({'Pixel_Size':{'X': pixel_size_x, 'Y':pixel_size_y, 'Units':'arcsec'}})

        smooth=0
        if self.MetaData['Frequency'] == '405' :
            smooth=10
        x,y       = Circle.get_limb(np.arange(P),np.arange(P),self.image,smooth=smooth)
        par,cov   = Circle.fit(x,y)
        self.MetaData.update({'Sun_Center_Matrix':{'X0':int(par[0]),'Y0':int(par[1])}})

        return

    def Rotate_North_Image(self):
        #
        # Rotate_North_Image
        #
        #   SST maps are Alt Azimuth oriented. This method
        #   rotates the matrix to have the Sun North Pole pointing
        #   'up' in ordinates.  To do this it needs the Parallactic
        #   angle  and the Sun P angle both provided by sstMap methods
        #
        #  It also creates a sun_off_x/y coordinate system that makes the
        #  Sun center in (0,0) and has arcsec units.  It can be used to compare
        #  with other observations.
        #
        #  The method computes the new solar sun center, while it should be
        #  (0,0), the coarseness of the original data makes that interpolations
        #  displace a little bit the sun position.
        #
        #  @guiguesp - 2020-04-05
        #
        # A note on the Sun orientation. There is a funciton in sunpy.coordinates.sun
        # to get the "orientation angle": sunpy.coordinates.sun.orientation(location,time)
        # The orientation angle is the angle between the local zenith and the solar north
        # measured eastward from local zenith. This angle takes into account the
        # P angle and the parallactic angle. However, the value is very different.
        #
        # I leave both results in the MetaData. This version of sstMap rotates a map
        # in 180-orientation angle, something that seems to be in accord with
        # SDO/AIA images. But I have no clear explanation yet for what I'm doing.
        # It should be checked!
        #
        # @guiguesp - 2020-08-01
        #______________________________________________________________________________

        if self.FLAGS['Rot_North'] :
            # Done!!
            return

        dx = self.image.shape[1]//2-self.MetaData['Sun_Center_Matrix']['X0']
        dy = self.image.shape[1]//2-self.MetaData['Sun_Center_Matrix']['Y0']
        angle = self.MetaData['Sun_Coordinates']['North_Orientation'].value - 180

        im = np.roll(np.roll(self.image, dx,1),dy,0)
#        im = rotate(im,-self.MetaData['Parallactic_Angle'],reshape=False)
#        im = rotate(im,self.MetaData['Sun_Coordinates']['P'].value,reshape=False)
        im = rotate(im,angle,reshape=False)
        im = np.roll(np.roll(im, -dx, 1), -dy, 0)
        self.image = im

        xoff = np.roll(np.roll(self.x_off,dx,1),dy,0)
#        xoff = rotate(xoff,-self.MetaData['Parallactic_Angle'],reshape=False)
#        xoff = rotate(xoff,self.MetaData['Sun_Coordinates']['P'].value,reshape=False)
        xoff = rotate(xoff,angle,reshape=False)
        xoff = np.roll(np.roll(xoff,dx,1), dy,0)
        self.x_off = xoff

        yoff = np.roll(np.roll(self.y_off,dx,1),dy,0)
#        yoff = rotate(yoff,-self.MetaData['Parallactic_Angle'],reshape=False)
#        yoff = rotate(yoff,self.MetaData['Sun_Coordinates']['P'].value,reshape=False)
        yoff = rotate(yoff,angle,reshape=False)
        yoff = np.roll(np.roll(yoff,dx,1), dy,0)
        self.y_off = yoff

        azi = np.roll(np.roll(self.azi,dx,1),dy,0)
#        azi = rotate(azi,-self.MetaData['Parallactic_Angle'],reshape=False)
#        azi = rotate(azi,self.MetaData['Sun_Coordinates']['P'].value,reshape=False)
        azi = rotate(azi,angle,reshape=False)
        azi = np.roll(np.roll(azi,-dx,1), -dy,0)
        self.azi = azi

        ele = np.roll(np.roll(self.ele,dx,1),dy,0)
#        ele = rotate(ele,-self.MetaData['Parallactic_Angle'],reshape=False)
#        ele = rotate(ele,self.MetaData['Sun_Coordinates']['P'].value,reshape=False)
        ele = rotate(ele,angle,reshape=False)
        ele = np.roll(np.roll(ele,-dx,1), -dy,0)
        self.ele = ele

        dx = []
        dy = []

        dx_model = (np.arange(self.image.shape[1]) - self.MetaData['Sun_Center_Matrix']['X0']) * self.MetaData['Pixel_Size']['X']
        for iscan in np.arange(self.image.shape[0]):
            dx.append(dx_model)
        dy_model = (np.arange(self.image.shape[0]) - self.MetaData['Sun_Center_Matrix']['Y0']) * self.MetaData['Pixel_Size']['Y']
        for iscan in np.arange(self.image.shape[1]):
            dy.append(dy_model)
        dx = np.asarray(dx)
        dy = (np.asarray(dy)).transpose()

        x,y       = Circle.get_limb(dx,dy,self.image)
        par,cov   = Circle.fit(x,y)

        self.sun_off_x  = dx
        self.sun_off_y  = dy
        self.MetaData.update({'limb_off':{'X0':par[0],'Y0':par[1],'R':par[2]}})
        self.FLAGS['Rot_North'] = True

        return

class Sun(object):

    def __init__(self):
        self.quiet_sun_temp = {'212': 5.9E+03, '405':5.1E+03}

    def __str__(self):
        return 'A Class to compute Sun ephemeris for El Leoncito.'

    def Get_Sun_Coordinates(self,cg,ctime):
        #
        #  Get_Sun_Coordinates
        #
        #     This method obtains the Sun Coordinates both in (Ra,Dec) and (Alt, Az)
        #     While astropy may use very accurate calculations and corrections from
        #     IERS, for some reason it cannot download the data. So, I decided to
        #     work off-line with a lesser precision, although enoug for us.
        #
        #
        #  @guiguesp - 2020-04-05
        #            - 2020-05-30 : Now using sunpy 1.1
        #_____________________________________________________________________________

        t = Time(ctime)

        with solar_system_ephemeris.set('builtin'):
            sun_radec=get_body('sun',t,cg)

        aa = AltAz(location=cg,obstime=t)
        sun_altaz = sun_radec.transform_to(aa)

        L0 = Sun_Coordinates.sun.L0(ctime)
        B0 = Sun_Coordinates.sun.B0(ctime)
        P  = Sun_Coordinates.sun.P(ctime)
        R_Sun = Sun_Coordinates.sun.angular_radius(ctime)
        Carrington  = Sun_Coordinates.sun.carrington_rotation_number(ctime)
        North_Orientation = Sun_Coordinates.sun.orientation(cg,t)

        return sun_radec,sun_altaz,L0,B0,P,R_Sun,Carrington,North_Orientation

##########################################################################################################################

class sstMap(Map):

    #
    # sstMap class
    #
    #   This is the Main class. "It has"  two other Classes: Scans and Map
    #   It initializes the class, and computes the astronomical data used
    #   later. Some properties are the Sun ephemerides in RaDec and AltAz
    #   the parallactic angle and the convertion factor ("calibration")
    #   to brightness temperature. It also gathers the MetaData.
    #
    #   Other properties of interest are obsazi and obsele, derived from
    #   the original scans, can be used for the pointing model.
    #
    #   Methods
    #      Get_Version
    #      Get_Time
    #      Get_Sun_Coordinates
    #      Get_Observed_Sun_Center
    #      Compute_Parallactic_Angle
    #      Get_Calibration
    #      Calibrate_Map_in_T
    #
    #  @guiguesp - 2020-04-05
    #            - 2020-05-29 : Write_FITS() mtehod added
    #_____________________________________________________________________________

    def __str__(self):
        return 'The Class fo Solar maps with SST.'

    @property
    def version(self):
        return self._version

    def __init__(self,d,channel=0):

        if (d.MetaData['SSTType'] != 'Integration'):
            print('\n  Only Integration (RS) data accepted\n')
            return

        self._version = _Version
        self.FLAGS   = {'Cal':False, 'Rot_North': False}
        self.Sun     = Sun()

        self.image   = np.empty([256,256])
        frequency = '212'
        reference_Tb = self.Sun.quiet_sun_temp['212']
        if (channel > 3):
            frequency = '405'
            reference_Tb = self.Sun.quiet_sun_temp['405']
        self.MetaData = d.MetaData
        self.MetaData.update({'Frequency': frequency})
        self.MetaData.update({'Channel':channel})
        self.MetaData.update({'reference_Tb': reference_Tb})

        self.Scans = Scans(d,channel)
        self.Scans.Select_Data(direction='az')

        t_start = d.Data['time'][self.Scans.Indices['iMap'][0]]
        t_end   = d.Data['time'][self.Scans.Indices['iMap'][-1]]
        t_mean = (t_start +  t_end) // 2
        self.MetaData.update({'t_start':self.Get_Time(d.MetaData,t_start)})
        self.MetaData.update({'t_end':self.Get_Time(d.MetaData,t_end)})
        self.MetaData.update({'ctime':self.Get_Time(d.MetaData,t_mean)})
        cg = CASLEO.Observatory_Coordinates()
        self.MetaData.update({'Observatory':cg})
        radec,altaz,L0,B0,P,R_Sun,Carrington,North_Orientation = self.Sun.Get_Sun_Coordinates(cg,self.MetaData['ctime'])
        self.MetaData.update({'Sun_Coordinates': {'RaDec':radec,
                                                  'AltAz': altaz,
                                                  'L0': L0, 'B0': B0,
                                                  'P':P, 'R_Sun': R_Sun,
                                                  'Carrington': Carrington,
                                                  'North_Orientation':-North_Orientation},
                                'Observatory': cg})

        self.MetaData.update({'Parallactic_Angle':self.Compute_Parallactic_Angle()})
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

    def Get_Time(self,md,husec):
    #
    #
    #    Get_Time: Class method to convert us time a Python datetime.
    #
    #    Change Record:  First written by Guigue @ Sampa for oRBD
    #                    2017-11-04 St Charles Day !!!
    #                    2020-04-01 Fool's Day (adapted for sstMap)
    #
    #_____________________________________________________________________________

        year  = int(md['ISODate'][0:4])
        month = int(md['ISODate'][5:7])
        day   = int(md['ISODate'][8:])

        hours =  husec // 36000000
        minutes = (husec % 36000000) // 600000
        seconds = ((husec % 36000000) % 600000) / 1.0E+04
        seconds_int  = int(seconds)
        seconds_frac = seconds - int(seconds)
        useconds     = int(seconds_frac * 1e6)

        return dt.datetime(year,month,day,hours,minutes,seconds_int,useconds)

    def Get_Calibration(self):

    #
    #  Get_Calibration
    #
    #    This procedure determines the multiplicative factor to convert
    #    the background subtracted scans in ADCu to solar brightness temperature.
    #    From an histogram determines the maximum of the secondary peak (the primary
    #    peak is near zero an represents the sky) and attributes to it the
    #    quiet Sun brightness temperature determined by Adriana (2005) paper.
    #
    #  @guiguesp - 2020-04-05
    #
    #_____________________________________________________________________________

        nbins = 200
        h,xh  = np.histogram(np.ravel(self.Scans.T),bins=nbins)
        qs    = xh[nbins//2+np.argmax(h[nbins//2:])]
        #
        # Reference Quiet Sun Brightness Temperature from A. Silva (2005)
        #

        # Default is 212 GHz

        self.MetaData.update({'adc2T': self.MetaData['reference_Tb']/qs})
        return

    def Calibrate_Map_in_T(self):

    #
    #  Calibrate_Map_in_T
    #
    #    Applies the calibration to the data.
    #
    #_____________________________________________________________________________

        if self.FLAGS['Cal'] :
            # Done!!
            return

        self.Get_Calibration()
        self.image = self.MetaData['adc2T'] * self.image
        self.FLAGS['Cal'] = True

        return

    def Write_FITS(self,FITSname=''):
        if (not hasattr(self,'image')):
            return

        if len(FITSname) < 1:
            FITSname = 'sst_map_' + self.MetaData['Frequency'] + 'GHz_Channel_' + str(self.MetaData['Channel']) + '_' + self.MetaData['ctime'].isoformat()[0:16]+'.fits'

        self.MetaData.update({'FITSfname':FITSname})
        _isodate_ = self.MetaData['ISODate']
        _hdu_ = fits.PrimaryHDU(self.image)

                #
                # This is the Primary (global) header. It gives information about the instrument, and the data.
                # It probably would be better to include this information in the XML decriptive files.
                # Anyway, oRBD.py and the XML files are a linked together.
                #
        _hdu_.header.append(('origin','CRAAM/Universidade Presbiteriana Mackenzie',''))
        _hdu_.header.append(('telescop','Solar Submillimeter Telescope',''))
        _hdu_.header.append(('observat','CASLEO',''))

        _hdu_.header.append(('station',
        'Lat = ' + str(self.MetaData['Observatory'].lat.degree) + ', Lon = ' + \
         str(self.MetaData['Observatory'].lon.degree) + ', Height = ' + \
         str(self.MetaData['Observatory'].height.value/1000)[0:4] + ' km',''))

        _hdu_.header.append(('tz','GMT-3',''))

        _hdu_.header.append(('date-obs',_isodate_,''))
        _hdu_.header.append(('t_start',self.MetaData['t_start'].isoformat()[0:16],''))
        _hdu_.header.append(('t_end',self.MetaData['t_end'].isoformat()[0:16],''))
        if isinstance(self.MetaData['RBDFileName'],list) :
            for iRBD in self.MetaData['RBDFileName']:_hdu_.header.append(('origfile',iRBD,'SST Raw Binary Data file'))
        else:
            _hdu_.header.append(('origfile',self.MetaData['RBDFileName'],'SST Raw Binary Data file'))
        _hdu_.header.append(('frequen',self.MetaData['Frequency'],'GHz'))

        # The World Coordinate System
        _hdu_.header.append(('CTYPE1','HPLN-TAN',''))
        _hdu_.header.append(('CUNIT1',self.MetaData['Pixel_Size']['Units'],''))
        _hdu_.header.append(('CRVAL1','0.000000',''))
        _hdu_.header.append(('CRPIX1',str(self.MetaData['Sun_Center_Matrix']['X0']),''))
        _hdu_.header.append(('CDELT1',str(self.MetaData['Pixel_Size']['X']),''))
        _hdu_.header.append(('CTYPE2','HPLN-TAN',''))
        _hdu_.header.append(('CUNIT2',self.MetaData['Pixel_Size']['Units'],''))
        _hdu_.header.append(('CRVAL2','0.000000',''))
        _hdu_.header.append(('CRPIX2',str(self.MetaData['Sun_Center_Matrix']['Y0']),''))
        _hdu_.header.append(('CDELT2',str(self.MetaData['Pixel_Size']['Y']),''))

        # Some Solar Parameters
        _hdu_.header.append(('R_SUN',str(self.MetaData['Sun_Coordinates']['R_Sun'].value)[0:7],'arcsecs'))
        _hdu_.header.append(('OBSRAD',str(self.MetaData['limb_off']['R'])[0:7],'arcsecs'))
        _hdu_.header.append(('L0',str(self.MetaData['Sun_Coordinates']['L0'].value)[0:7],'degres'))
        _hdu_.header.append(('B0',str(self.MetaData['Sun_Coordinates']['B0'].value)[0:7],'degres'))
        _hdu_.header.append(('P',str(self.MetaData['Sun_Coordinates']['P'].value)[0:7],'degres'))

        # About the Copyright
        _hdu_.header.append(('comment','COPYRIGHT. Grant of use.',''))
        _hdu_.header.append(('comment','These data are property of Universidade Presbiteriana Mackenzie.'))
        _hdu_.header.append(('comment','The Centro de Radio Astronomia e Astrofisica Mackenzie is reponsible'))
        _hdu_.header.append(('comment','for their distribution. Grant of use permission is given for Academic '))
        _hdu_.header.append(('comment','purposes only. Any questions contact guigue@craam.mackenzie.br'))

        if self.FLAGS['Cal']:
            _hdu_.header.append(('ADC2T',str(self.MetaData['adc2T'])[0:7],'K/ADCu'))
            _hdu_.header.append(('history','Calibrated in Brightness Temperature'))

        if self.FLAGS['Rot_North']:
            _hdu_.header.append(('PARANG',str(self.MetaData['Parallactic_Angle'])[0:7],'degrees'))
            _hdu_.header.append(('history','Rotated, North Up'))


        _hdu_list_ = fits.HDUList([_hdu_])
        try:
            _hdu_list_.writeto(FITSname)
            return True
        except OSError as err:
            print("\n\nWrite Error: {0}\n\n".format(err))
            return False

class MapOff(object):
    
    def __str__(self):
        return 'A Class to determine the off pointing of SST maps.'

    def __init__(self,d,ch=[0,1,2,3],direction='az'):
        self.Version = '20210727T1606BRT'
        self.Channel = {}

        for ich in np.arange(len(ch)):
            chname=str(ch[ich])
            self.Channel.update({chname:Scans(d,ch[ich])})
            self.Channel[chname].Select_Data(direction=direction)
            self.Channel[chname].Extract_Scans(direction=direction)
            self.Channel[chname].Limb_fit()

        return
        
        
    def plotFit(self,ch=0,save=False):

        chname=str(ch)
        if chname in self.Channel.keys():
            angle = np.linspace(0,2*np.pi,500)
            fig  = plt.figure()
            fig.set_size_inches(w=6,h=6)
            Bpos = [0.1,0.1,0.88,0.88]
            ax   = fig.add_subplot(1,1,1, position=Bpos)
            p1   = ax.plot(self.Channel[chname].limb['off'],self.Channel[chname].limb['el'],'.',color='red')
            p2   = ax.plot(self.Channel[chname].limb['X0']+self.Channel[chname].limb['R']*np.cos(angle),
                           self.Channel[chname].limb['Y0']+self.Channel[chname].limb['R']*np.sin(angle),
                           color='black',linestyle='dashed')
            p3   = ax.plot([self.Channel[chname].limb['X0'],self.Channel[chname].limb['X0']],
                           [self.Channel[chname].limb['el'].min(),self.Channel[chname].limb['el'].max()],
                           linestyle='dotted',color='black')
            p3   = ax.plot([self.Channel[chname].limb['off'].min(),self.Channel[chname].limb['off'].max()],
                           [self.Channel[chname].limb['Y0'],self.Channel[chname].limb['Y0']],
                           linestyle='dotted',color='black')
        
            ax.set_aspect('equal')
            ax.text(self.Channel[chname].limb['X0']-600,self.Channel[chname].limb['Y0']+100,
                    'Map Center = (' + str(self.Channel[chname].limb['X0'])[:7]+' , ' + str(self.Channel[chname].limb['Y0'])[:7] + ')  arcsec')
        
            ax.text(self.Channel[chname].limb['X0']-600,self.Channel[chname].limb['Y0']+200,
                    'Solar Radius = ' + str(self.Channel[chname].limb['R'])[:6] + ' +/- ' + str(np.sqrt(np.diag(self.Channel[chname].limb['Cov']))[2])[:4] + ' [arcsec]')

            ax.set_title('Solar Limb : Channel '+chname)
            ax.set_xlabel('Off [arcsec]')
            ax.set_ylabel('El [arcsec]')

            if save:
                fname = 'SST-Limb_Fit_Channel_'+chname+'.pdf'
                plt.savefig(fname, dpi=None, facecolor='w', edgecolor='w',
                            orientation='portrait', papertype=None, format='pdf',
                            transparent=False, bbox_inches=None, pad_inches=0.1,
                            frameon=None, metadata=None)

        else:
            print("\n\n {0:d} is not a valid channel. \n\n".format(ch))

            
        return

    def printFit(self,ch=[0]):

        for ich in np.arange(len(ch)):
            chname = str(ch[ich])
            if chname in self.Channel.keys():
                print("\n\n        Channel {0:s}\n".format(chname))
                print(" off Center = {0:6.2f} , {1:6.2f}  [arcmin]\n".format(self.Channel[chname].limb['X0']/60,self.Channel[chname].limb['Y0']/60))
                print(" Radius = {0:6.2f} +/- {1:5.2f} [arcsec]\n\n".format(self.Channel[chname].limb['R'],np.sqrt(np.diag(self.Channel[chname].limb['Cov']))[2]))

            else:
                print("\n\n {0:s} is not a valid channel. \n\n".format(chname))

            
        return
        
