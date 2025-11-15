import glob
import numpy as np
from astropy.io import fits

def read(dirname):

    flist = glob.glob(dirname+'*.fits')
    f=fits.open(flist[0])
    cube = np.ndarray((f[0].data.shape[0],f[0].data.shape[1],len(flist)),dtype=np.float64)
    i=0
    for file in flist:
        f=fits.open(file)
        cube[:,:,i] = f[0].data

    return cube
