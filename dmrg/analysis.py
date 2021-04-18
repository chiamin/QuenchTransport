import pylab as pl
from collections import OrderedDict
import sys
from math import pi, acos
import plotsetting as ps
import numpy as np
import fitfun as ff
import matplotlib.colors as colors
import cmasher as cmr

def exactG (V, tp):
    if V != 0:
        return 0.5 * pi / acos(-0.5*V)
    else:
        return 4*tp**2/(1+tp**2)**2

def get_para (fname, key, typ, last=False, n=1):
    with open(fname) as f:
        for line in f:
            if key in line:
                val = list(map(typ,line.split()[-n:]))
                if n == 1: val = val[0]
                if not last:
                    return val
        return val

def get_data (fname):
    L_lead = get_para (fname, 'L_lead', int)
    L_device = get_para (fname, 'L_device', int)
    L = 2*L_lead + L_device
    ns = np.full (L, np.nan)
    Ss = np.full (L, np.nan)
    ms = np.full (L, np.nan)
    with open(fname) as f:
        for line in f:
            line = line.lstrip()
            if line.startswith('n '):
                tmp = line.split()
                i = int(tmp[1])
                n = float(tmp[-1])
                ns[i-1] = n
            elif line.startswith('entang entropy'):
                tmp = line.split()
                i = int(tmp[2])
                S = float(tmp[-1])
                Ss[i-1] = S
            elif 'HS=' in line and 'Bond=' in line:
                tmp = line.split('=')
                bond = int(tmp[-1].split('/')[0])
            elif 'States kept:' in line:
                m = int(line.split('=')[-1])
                ms[bond-1] = m
    return ns, Ss, ms

if __name__ == '__main__':
    files = [i for i in sys.argv[1:] if i[0] != '-']
    for fname in files:
        print (fname)
        ns, Ss, ms = get_data (fname)
        idevL, idevR = get_para (fname, 'device site', int, n=2)

        # density
        f,ax = pl.subplots()
        ax.plot (range(1,len(ns)+1), ns, marker='.')
        ax.set_xlabel ('site')
        ax.set_ylabel ('density')
        ps.set(ax)
        print (np.sum(ns[idevL-1:idevR]))

        # entropy
        f2,ax2 = pl.subplots()
        ax2.plot (range(1,len(ns)+1), Ss, marker='.')
        ax2.set_xlabel ('site')
        ax2.set_ylabel ('entropy')
        ps.set(ax2)

        # entropy
        f3,ax3 = pl.subplots()
        ax3.plot (range(1,len(ns)+1), ms, marker='.')
        ax3.set_xlabel ('site')
        ax3.set_ylabel ('bond dimension')
        ps.set(ax3)
    pl.show()
