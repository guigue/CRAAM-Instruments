import numpy as np
import pdb
import scipy.interpolate as intl
from scipy.optimize import curve_fit

import statistics as stat
import datetime as dt

from astropy.io import fits
import astropy.constants as cts

from CraamTools import contiguo

class TempCal():

    def Version(self):
        __version = '2024-01-22T18:00BRT'
        return __version

    def __str__(self):
        return 'Extract antenna temperature calibration parameters from an SST auxiliary data file.'

    def __init__(self,b):

        try:
            b.MetaData['SSTType'] == 'Auxiliary'
        except:
            print('\n  Only SST RBD auxiliary data accepted\n')
            return

        self.MetaData = b.MetaData
        self.Version  = self.Version()
        self.Data     = self.Process(b)

        return
    

    def Process(self,b):
        
        # Find the cold periods
        xcold, = np.where(b.Data['target']//32 == 1)
        ccold  = contiguo.contiguo(xcold)
        Ncal   = ccold.shape[0]
        Nchan  = 6
        
        # Find the hot periods
        xhot, = np.where(b.Data['target']//32 == 2)
        chot  = contiguo.contiguo(xhot)

        adc2K       = np.zeros((Nchan,Ncal),float)
        Trms        = np.zeros((Nchan,Ncal),float)
        husec       = np.zeros(Ncal,int)
        thot        = np.zeros(Ncal,float)
        tcold       = np.zeros(Ncal,float)
        temperature = np.zeros(Ncal,float)
        pressure    = np.zeros(Ncal,float)
        humidity    = np.zeros(Ncal,float)
        ifrec_temp  = np.zeros(Ncal,float)
        optbox_temp = np.zeros(Ncal,float)
        radome_temp = np.zeros(Ncal,float)
        
        
        for i in np.arange(Ncal):

            tcold[i]       = b.Data['amb_temp'][xcold[ccold[i,0]]:xcold[ccold[i,1]]].mean()
            thot[i]        = b.Data['hot_temp'][xhot[chot[i,0]]:xhot[chot[i,1]]].mean()
            
            husec[i]       = b.Data['time'][xhot[chot[i,0]]:xhot[chot[i,1]]].mean()
            ifrec_temp[i]  = b.Data['if_board'][xhot[chot[i,0]]:xhot[chot[i,1]]].mean()
            optbox_temp[i] = b.Data['opt_temp'][xhot[chot[i,0]]:xhot[chot[i,1]]].mean()
            radome_temp[i] = b.Data['radome_temp'][xhot[chot[i,0]]:xhot[chot[i,1]]].mean()
            temperature[i] = b.Data['temperature'][xhot[chot[i,0]]:xhot[chot[i,1]]].mean()
            pressure[i]    = b.Data['pressure'][xhot[chot[i,0]]:xhot[chot[i,1]]].mean()
            humidity[i]    = b.Data['humidity'][xhot[chot[i,0]]:xhot[chot[i,1]]].mean()
            
            for ch in np.arange(Nchan):

                ADCUcold    = b.Data['adc'][xcold[ccold[i,0]]:xcold[ccold[i,1]],ch].mean()
                ADCUhot     = b.Data['adc'][xhot[chot[i,0]]:xhot[chot[i,1]],ch].mean()
                a           = np.polyfit( [ADCUcold,ADCUhot] , [ tcold[i]+273.15, thot[i]+273.15 ], 1 )
                adc2K[ch,i] = a[0]
                Trms[ch,i]  = - a[1]
                

        return {'adc2K':adc2K,'Trms':Trms,'husec':husec,
                'Thot':thot,'Tcold':tcold,'Temperature':temperature,
                'IFRec_temp':ifrec_temp,'OptBox_temp':optbox_temp,
                'Radome_temp':radome_temp}
                
        

        


