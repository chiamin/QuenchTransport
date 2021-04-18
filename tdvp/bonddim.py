import pylab as pl
from collections import OrderedDict
import sys
from math import pi, acos
import plotsetting as ps
import numpy as np

def exactG (V, tp):
    if V != 0:
        return 0.5 * pi / acos(-0.5*V)
    else:
        return 4*tp**2/(1+tp**2)**2

def get_data (fname):
    with open(fname) as f:
        for line in f:
            if 'dt =' in line:
                dt = float(line.split()[-1])
        print ('dt =',dt)

    ms = OrderedDict()
    terrs = OrderedDict()
    with open(fname) as f:
        for line in f:
            if 'step =' in line:
                step = int(line.split()[-1])
            elif 'Largest link dim' in line:
                m = float(line.split()[-1])
                ms[step] = m
            elif 'Largest truncation error' in line:
                terr = float(line.split()[-1])
                terrs[step] = terr

    tns = ms.keys()
    ts = [dt*tn for tn in tns]
    return ts, ms.values(), terrs.values()

def plot_data (fname, ax1, ax2):
    ts, ms, terrs = get_data (fname)

    ax1.plot (ts, ms, marker='.')
    ax2.plot (ts, terrs, marker='.')

if __name__ == '__main__':
    f1,ax1 = pl.subplots()
    f2,ax2 = pl.subplots()

    files = [i for i in sys.argv[1:] if i[0] != '-']
    for fname in files:
        print (fname)

        plot_data (fname, ax1, ax2)

        if '-pdf' in sys.argv:
            f1.savefig (fname+'_I.pdf')
    ax1.set_ylabel ('$m$')
    ax2.set_ylabel ('truncation error')
    ps.set([ax1,ax2])
    f1.savefig ('m.pdf')
    f2.savefig ('terr.pdf')
    pl.show()
