import numpy as np
import pdb
import scipy.interpolate as intl
from scipy.optimize import curve_fit

import statistics as stat
import datetime as dt

from astropy.io import fits
import astropy.constants as cts

class skydip():

    def Version(self):
        __version = '20230903T1252BRT'
        return __version

    def __str__(self):
        return 'Compute the Optical Depth tau from Intergrated SST Raw Binary Data'

    def __init__(self,d,emin=0,emax=90):

        try:
            d.MetaData['SSTType'] == 'Integration'
        except:
            print('\n  Only SST RBD Integration data accepted\n')
            return

        self.MetaData = d.MetaData
        self.Version  = self.Version()

        self.extractData(d,emin,emax)

        self.Fit = {}
        for chan in np.arange(6):
            chanStr = 'Chan'+str(chan)
            self.Fit.update({chanStr:self.fitTau(chan)})

        self.MeanVal()

        return

    def extractData(self,d,emin,emax):

         x     = (d.Data['opmode'] == 10) & (d.Data['elepos']/1000 >= emin) & (d.Data['elepos']/1000 <= emax)
         ele   = d.Data['elepos'][x] / 1000
         adc   = d.Data['adcval'][x,:]
         azi   = np.mean(d.Data['azipos'][x]) / 1000

         year  = int(d.MetaData['ISODate'][0:4])
         month = int(d.MetaData['ISODate'][5:7])
         day   = int(d.MetaData['ISODate'][8:])
         husec = d.Data['time'][x]

         hours =  husec[0] // 36000000
         minutes = (husec[0] % 36000000) // 600000
         seconds = ((husec[0] % 36000000) % 600000) / 1.0E+04
         seconds_int  = int(seconds)
         seconds_frac = seconds - int(seconds)
         useconds     = int(seconds_frac * 1e6)
         time  = [dt.datetime(year,month,day,hours,minutes,seconds_int,useconds)]

         hours =  husec[-1] // 36000000
         minutes = (husec[-1] % 36000000) // 600000
         seconds = ((husec[-1] % 36000000) % 600000) / 1.0E+04
         seconds_int  = int(seconds)
         seconds_frac = seconds - int(seconds)
         useconds     = int(seconds_frac * 1e6)
         time.append(dt.datetime(year,month,day,hours,minutes,seconds_int,useconds))
         
         self.Data = {'adc': adc,'ele':ele, 'azi':azi, 'time':time}
         return

    def skyModel(self,x,Toff,Ts,tau):
        return Toff + Ts * ( 1 - np.exp(-tau/np.sin(np.radians(x))))

    def fitTau(self,chan=0):
        if (chan < 0) | (chan > 5) :
            print('\n Wrong Channel Number',chan)
            return

        if (chan < 4):
            par0 = [1.0e+04,1.0E+04,0.2]
        else:
            par0 = [1.0e+04,1.0E+04,0.8]

        par, cov = curve_fit(self.skyModel, self.Data['ele'], self.Data['adc'][:,chan], p0 = par0)
        dfit     = self.skyModel(self.Data['ele'],par[0],par[1],par[2])
        perr = (np.sqrt(np.diag(cov)))
        return {'tau':par[2],'stau':perr[2],
                'Toff':par[0],'sToff': perr[0],
                'Ts':par[1],'sTs': perr[1],
                'Model':dfit}

    def MeanVal(self):

        tau212 = np.array([self.Fit['Chan0']['tau'], self.Fit['Chan1']['tau'],
                           self.Fit['Chan2']['tau'], self.Fit['Chan3']['tau'] ])
        stau212 = np.array([self.Fit['Chan0']['stau'], self.Fit['Chan1']['stau'],
                           self.Fit['Chan2']['stau'], self.Fit['Chan3']['stau'] ])
        Toff    = np.array([self.Fit['Chan0']['Toff'], self.Fit['Chan1']['Toff'],
                           self.Fit['Chan2']['Toff'], self.Fit['Chan3']['Toff'],
                           self.Fit['Chan4']['Toff'], self.Fit['Chan5']['Toff']])

        mp212 = 0
        s212 = 0
        for i in np.arange(4):
            mp212 += tau212[i] / stau212[i]**2
            s212 += 1/stau212[i]**2

        mp212 = mp212 / s212
        s212  = np.sqrt(1/s212)

        mp405 = self.Fit['Chan4']['tau']/self.Fit['Chan4']['stau']**2 + \
               self.Fit['Chan5']['tau']/self.Fit['Chan5']['stau']**2

        s405 = 1/self.Fit['Chan4']['stau']**2 + 1/self.Fit['Chan5']['stau']**2

        mp405 = mp405 / s405

        s405 = np.sqrt(1/s405)

        self.Tau = {'t212':mp212, 'st212':s212, 't405': mp405, 'st405':s405}

        return
