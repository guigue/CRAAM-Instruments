import os
import glob
import warnings
import pdb
import datetime as dt
import pickle
from scipy.optimize import curve_fit
from astropy import units as u
from astropy import constants as c
from astropy.time import Time
import datetime as dt

import HATS

############ Global Variables ##########
__version       = "2025-04-11T1530ART"
#######################################

##############################################
#
# HATSTools: a number of tools to monitor HATS data
#
######################################################
#
# HATSday
#
#        Reads data from an entire day.
#
# Parameters:
#        day='YYYY-MM-DD'
#            A string in the ISO format date
# Output:
#        h: A HATS object with the data of the whole day
#        thats : datetime of HATS data
#
######################################################
#
# HATSmonitor
#
#       Plots current data
#
# Parameters:
#       intervalo=nnn
#           The interval to plot in minutes.
#           The interval ends in the current time.
#
#       savefig=False|True
#           If True a PDF file with the plot is created.
#           default is False
#
#####################################################
#
# HATSws
#
#       Retrieves wewather station data for the current day.
#
# Parameters:
#       savefig= False|True
#           If True a PDF file with the plot is created.
#           default is False
#
########################################################

def HATSToolsversion():

    return __version

def HATSday(day=''):
    
    from matplotlib import pyplot as plt
    from matplotlib import dates

    if (day == ''):
        t0 = Time.now()
        day = '{0:04d}-{1:02d}-{2:02d}'.format(t0.datetime.year,
                                               t0.datetime.month,
                                               t0.datetime.day)

    InputPath = os.getenv('HATS_DATA_InputPath')

    lista = glob.glob(InputPath+'/hats-'+day+'T*')
    if lista == []:
        print('\n\n No data for '+ day+'\n\n')
        return

    print('\n\n {0:d} files found \n\n'.format(len(lista)))
    
    lista.sort()
    horas = []
    
    for fn in lista:
        horas.append(fn[-8:-4])
    
    h = HATS.hats(day+' '+horas[0])
    for hh in horas[1:]:
        h=h+HATS.hats(day+' '+hh)

    thats = h.getTimeAxis(h.rbd.Deconv['husec'],h.MetaData)

    fig=plt.figure()
    fig.set_size_inches(((25*u.cm).to(u.imperial.inch)).value,
                        ((15*u.cm).to(u.imperial.inch)).value)
    ax=fig.add_subplot(1,1,1)
    pos=[0.1,0.1,0.85,0.85]
    ax.set_position(pos)
    ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
    ax.plot(thats,h.rbd.Deconv['amplitude'],'-r')
    ax.set_xlabel('UT')
    ax.set_ylabel('mV')
    ax.set_title(day)
                         
    return h,thats
            
def HATSmonitor(intervalo=20,savefig=False):
    
    from matplotlib import pyplot as plt
    from matplotlib import dates

    if (intervalo > 120):
        print('\n\n intervalo should be < 120 minutes\n\n')
        return
    
    interval = intervalo * u.min
    
    t = Time.now()
    arg = '{0:04d}-{1:02d}-{2:02d} {3:02d}00'.format(t.datetime.year,
                                                     t.datetime.month,
                                                     t.datetime.day,
                                                     t.datetime.hour)

    try:
        h = HATS.hats(arg)
    except:
        print("\n\n HATS is not observing\n\n")
        return

    t0 = t - dt.timedelta(seconds=interval.to(u.s).value)
    if (t.datetime.hour > t0.datetime.hour):
        arg0 = '{0:04d}-{1:02d}-{2:02d} {3:02d}00'.format(t0.datetime.year,
                                                          t0.datetime.month,
                                                          t0.datetime.day,
                                                          t0.datetime.hour)
    try:
        h0 = HATS.hats(arg0)
        h = h0 + h
    except:
        pass

    h = h.extract([t0.datetime,t.datetime])
    thats = h.getTimeAxis(h.rbd.Deconv['husec'],h.MetaData)
    
    fig=plt.figure()
    fig.set_size_inches(((25*u.cm).to(u.imperial.inch)).value,
                        ((15*u.cm).to(u.imperial.inch)).value)
    ax=fig.add_subplot(1,1,1)
    pos=[0.1,0.1,0.85,0.85]
    ax.set_position(pos)
    ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
    ax.plot(thats,h.rbd.Deconv['amplitude'],'-r')
    ax.set_xlabel('UT')
    ax.set_ylabel('mV')
    ax.set_title('{0:04d}-{1:02d}-{2:02d} {3:02d} UT'.format(t.datetime.year,
                                                             t.datetime.month,
                                                             t.datetime.day,
                                                             t.datetime.hour))

    if savefig:
        fname = 'Monitor-'+arg+'.pdf'
        plt.savefig(fname, dpi=None, facecolor='w', edgecolor='w',
                  orientation='portrait', papertype=None, format='pdf',
                  transparent=False, bbox_inches=None, pad_inches=0.1,
                  frameon=None, metadata=None)

    return

def HATSws(savefig=False):

    from matplotlib import pyplot as plt
    from matplotlib import dates

    t = Time.now()
    arg = 'hats-{0:04d}-{1:02d}-{2:02d}'.format(t.datetime.year,
                                                  t.datetime.month,
                                                  t.datetime.day)
    ws = HATS.ws(arg)
    
    fig=plt.figure()
    fig.set_size_inches(((20*u.cm).to(u.imperial.inch)).value,
                        ((15*u.cm).to(u.imperial.inch)).value)
    
    ax=fig.add_subplot(3,1,1)
    pos=[0.1,0.6,0.85,0.25]
    ax.set_position(pos)
    ax.set_ylabel('T [°C]')
    ax.set_xticklabels([None]*10)
    ax.plot(ws.data['time'],ws.data['temperature'],'-r')

    ax1=fig.add_subplot(3,1,2)
    pos=[0.1,0.35,0.85,0.25]
    ax1.set_position(pos)
    ax1.plot(ws.data['time'],ws.data['humidity'],'-b')
    ax1.set_ylabel('H [%]')
    ax1.set_xticklabels([None]*10)

    ax2=fig.add_subplot(3,1,3)
    pos=[0.1,0.1,0.85,0.25]
    ax2.set_position(pos)
    ax2.plot(ws.data['time'],ws.data['pressure'],'-k')
    ax2.set_ylabel('P [hPa]')
    ax2.set_xticklabels([None]*10)

    ax2.xaxis.set_major_formatter(dates.DateFormatter('%H'))
 
    ax2.set_xlabel('UT')
    ax.set_title('{0:04d}-{1:02d}-{2:02d}'.format(t.datetime.year,
                                                  t.datetime.month,
                                                  t.datetime.day))

    if savefig:
        fname = 'Weather-'+arg+'.pdf'
        plt.savefig(fname, dpi=None, facecolor='w', edgecolor='w',
                    orientation='portrait', papertype=None, format='pdf',
                    transparent=False, bbox_inches=None, pad_inches=0.1,
                    frameon=None, metadata=None)

    return ws
