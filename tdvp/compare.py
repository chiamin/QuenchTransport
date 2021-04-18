import numpy as np
import pylab as pl

if __name__ == '__main__':
    fname = '/home/chiamin/Dropbox/physics/conductance/jupyterNotebooks/data/N7current_L150_xi25.dat'
    dat = np.loadtxt (fname, dtype=complex, unpack=True)

    pl.plot (dat[0], -dat[1], marker='.')
    pl.plot (dat[0], -dat[2], marker='.')
    pl.plot (dat[0], -dat[3], marker='.')
    pl.plot (dat[0], -dat[4], marker='.')
    pl.plot (dat[0], -dat[5], marker='.')
    pl.plot (dat[0], -dat[6], marker='.')
    pl.plot (dat[0], -dat[7], marker='.')
    pl.show()

