#
# https://gist.github.com/nikhilkumarsingh/1dcec96a1eb0aeb8975fc13ec5825d43
#

from multiprocessing import shared_memory as shm
import posix_ipc

import os
import numpy as np
import matplotlib.pyplot as plt
import pdb
import HATS

_version_ = '2021-12-17T0700BRT'
_record_size_ = 38
_ACQUIRE_SHM_NAME_ = 'acquire_shmem'
_ACQUIRE_SEM_NAME_ = 'acquire_semaphore'

class rbd(object):

    def __init__(self,InputPath='/data/HATS/'):
        self.acq_shm = shm.SharedMemory(_ACQUIRE_SHM_NAME_)
        self.acq_sem = posix_ipc.Sempahore(_ACQUIRE_SEM_NAME_)
        self.rbd = HATS.rbd(InputPath=InputPath)
        self.MetaData = {'InputPath':InputPath,
                         'AqcShm':_ACQUIRE_SHM_NAME_,
                         'AcqSem':_ACQUIRE_SEM_NAME_}

    def version(self):
        self.version = __version__
        return self.version

    def loop_read(self):


        return
                    
            
            
