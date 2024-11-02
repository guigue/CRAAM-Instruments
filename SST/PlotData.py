#!/usr/bin/python


def getTimeAxis(d):

    ndata = d.Data['time'].shape[0]

    ssttime = np.array(np.empty(ndata),dtype=dt.datetime)
    year  = int(d.MetaData['ISODate'][0:4])
    month = int(d.MetaData['ISODate'][5:7])
    day   = int(d.MetaData['ISODate'][8:])
    for i in np.arange(ndata):
        ms = d.Data['time'][i]
        hours =  ms / 36000000L
        minutes = (ms % 36000000L) / 600000L
        seconds = ((ms % 36000000L) % 600000L) / 1.0E+04
        seconds_int  = int(seconds)
        seconds_frac = seconds - int(seconds)
        useconds     = int(seconds_frac * 1e6)
        ssttime[i] = dt.datetime(year,month,day,hours,minutes,seconds_int,useconds)
                        

    return ssttime
    

if __name__ == "__main__":

    import sys, string, os
    import numpy as np
    import matplotlib.pyplot as plt
    import datetime as dt

    import oRBD

    if len(sys.argv) < 3 :
        print 'Usage: '+sys.argv[0]+' RBDfilename ch'
        sys.exit(1)

    RBDfname=sys.argv[1]
    chn=int(sys.argv[2])
    
    d=oRBD.RBD()
    try:
        d.readRBDinDictionary(RBDfname)
    except Exception, e:
        print "An unexpected exception occurred"
        sys.exit(1)

    if d.MetaData['SSTType'] == 'Auxiliary' :
        fieldname='adc'
    else:
        fieldname='adcval'

    st=getTimeAxis(d)
    plt.plot(st,d.Data[fieldname][:,chn])
    plt.show()

    
    
