import numpy as np
import datetime as dt

# Global Constants
ColdCod = 1
HotCod = 2
ShiftBits = 5
C2K = 273.15

def Calibrations(b):
    
    if (b.MetaData['SSTType'] != 'Auxiliary'):
        print ('--------------------------------------')
        print (':                                    :')
        print (': Only for Auxiliary (BI) files      :')
        print (':                                    :')
        print ('--------------------------------------')
        return []

    b.CorrectAuxiliary()

    xhot,chot=getCalPos(b,HotCod)
    xcold,ccold=getCalPos(b,ColdCod)

    Ncalibs=np.min([ccold.shape[0],chot.shape[0]])

    Cal={'time'  : np.array(np.zeros(Ncalibs,dtype=dt.datetime)),
         'adc_h' : np.array(np.zeros([Ncalibs,6],dtype=np.float)),
         'adc_c' : np.array(np.zeros([Ncalibs,6],dtype=np.float)),
         'Trec'  : np.array(np.zeros([Ncalibs,6],dtype=np.float)),
         'ADC2K' : np.array(np.zeros([Ncalibs,6],dtype=np.float)),
         'Thot'  : np.array(np.zeros(Ncalibs,dtype=np.float)),
         'Tcold' : np.array(np.zeros(Ncalibs,dtype=np.float)),
         'IF_T'  : np.array(np.zeros(Ncalibs,dtype=np.float)),
         'Opt_T' : np.array(np.zeros(Ncalibs,dtype=np.float)),
         'Rad_T' : np.array(np.zeros(Ncalibs,dtype=np.float))
         }

    y = int(b.MetaData['ISODate'][0:4])
    m = int(b.MetaData['ISODate'][5:7])
    d = int(b.MetaData['ISODate'][8:])

    # Get the Mean time of the Calibrations
    for i in np.arange(Ncalibs):
        ms = np.mean(b.Data['time'][xcold[ccold[i,0]:chot[i,1]]])
        Cal['time'][i] = ms2dt(y,m,d,ms)
        Cal['Thot'][i] = np.mean(b.Data['hot_temp'][xhot[chot[i,0]:chot[i,1]]])+C2K
        Cal['Tcold'][i] = np.mean(b.Data['amb_temp'][xcold[ccold[i,0]:ccold[i,1]]])+C2K
        Cal['IF_T'][i] = np.mean(b.Data['if_board'][xcold[ccold[i,0]:chot[i,1]]])
        Cal['Opt_T'][i] = np.mean(b.Data['opt_temp'][xcold[ccold[i,0]:chot[i,1]]])
        Cal['Rad_T'][i] = np.mean(b.Data['radome_temp'][xcold[ccold[i,0]:chot[i,1]]])
        for ch in np.arange(6):
            Cal['adc_c'][i,ch] = np.mean(b.Data['adc'][xcold[ccold[i,0]:ccold[i,1]],ch])
            Cal['adc_h'][i,ch] = np.mean(b.Data['adc'][xhot[chot[i,0]:chot[i,1]],ch])
            Cal['ADC2K'][i,ch] = (Cal['Thot'][i] - Cal['Tcold'][i]) / (Cal['adc_h'][i,ch]-Cal['adc_c'][i,ch])
            Cal['Trec'][i,ch]  = Cal['adc_h'][i,ch] * Cal['ADC2K'][i,ch] - Cal['Thot'][i]

    return Cal

def getCalPos(b,Cod):
    x=np.where((b.Data['target'] >> ShiftBits) == Cod)
    x=np.asarray(x)
    x=x[0]
    return x,cntgs(x)

def ms2dt(y,m,d,ms):
    
    import datetime as dt

    ms = int(ms)
    hours =  ms // 36000000
    minutes = (ms % 36000000) // 600000
    seconds = ((ms % 36000000) % 600000) / 1.0E+04
    seconds_int  = int(seconds)
    seconds_frac = seconds - int(seconds)
    useconds     = int(seconds_frac * 1.0E+06)

    return dt.datetime(y,m,d,hours,minutes,seconds_int,useconds)

def cntgs(s):

    sDim=s.shape[0]
    if (sDim <= 2):
        return np.array([])

    sDiff = s[1:] - s[0:-1]
    tDisc = np.where(sDiff != 1)
    xDisc = np.asarray(tDisc)
    xDisc = xDisc[0]
    nDisc = xDisc.shape[0]+1
    c     = np.zeros( (nDisc,2), dtype=np.int)
    for i in np.arange(nDisc):
        if (i==0):
            c[0,0] = 0
            c[0,1] = xDisc[0]
        elif (i==nDisc-1):
            c[i,0] = xDisc[i-1]+1
            c[i,1] = s.shape[0]-1
        else:
            c[i,0] = xDisc[i-1]+1
            c[i,1] = xDisc[i]

    return c
            
