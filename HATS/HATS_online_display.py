#
# https://gist.github.com/nikhilkumarsingh/1dcec96a1eb0aeb8975fc13ec5825d43
#

import os
import numpy as np
import matplotlib.pyplot as plt
import pdb
import HATS

_version_ = '2021-12-17T0700BRT'
_record_size_ = 38

class rbd(object):

    def __init__(self,InputPath='/data/HATS/'):
        self.rbd = HATS.rbd(InputPath=InputPath)
        self.MetaData = {'InputPath':InputPath}

    def version(self):
        self.version = __version__
        return self.version

    def loop_read(self,nrecords=4096,fname='hats-2021-11-29T1000.rbd'):

        bytes = os.path.getsize(self.MetaData['InputPath']+'/'+fname)
        nloops = bytes // _record_size_ // nrecords
        print("\n\n N loops = {0:d}\n\n".format(nloops))

        husec = np.empty((0))
        amplitude = np.empty((0))
        time = []
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_ylim(-20,20)
        ax.set_xlim(0,8*nloops)
        for loop in np.arange(nloops):
            self.rbd.from_file(fname,nrecords=nrecords,recoffset=loop*nrecords)
            self.rbd.getFFT(steps=512,window_size=1024)
            #dt_time = list(map(self.rbd.husec2dt,self.rbd.Deconv['husec']))
            #time+=dt_time            
            amplitude=np.concatenate((amplitude,self.rbd.Deconv['amplitude']))
            ax.plot(amplitude,color='black')
            #ax.set_xlim(left=max(0,loop-1000),right=loop+1000)
            
            fig.canvas.draw()


        return
                    
            
            
