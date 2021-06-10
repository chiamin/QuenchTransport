import pylab as pl
from collections import OrderedDict
import sys, glob
from math import pi, acos
import plotsetting as ps
import numpy as np
import fitfun as ff
import matplotlib.colors as colors
import cmasher as cmr
import fitfun as ff

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

def get_hop_t (fname):
    L_lead = get_para (fname, 'L_lead', int)
    L_device = get_para (fname, 'L_device', int)
    L = 2*L_lead + L_device
    ts = np.full (L, np.nan)
    with open(fname) as f:
        for line in f:
            if 'H left lead' in line:
                part = 'L'
            elif 'H right lead' in line:
                part = 'R'
            elif 'H dev' in line:
                part = 'S'
            elif 'Hk, t' in line:
                if part == 'L': offset = 0
                elif part == 'S': offset = L_lead
                elif part == 'R': offset = L_lead + L_device
                tmp = line.split()
                i = int(tmp[-3])
                t = float(tmp[-1])
                ts[i+1+offset] = t
            elif 't_contactL' in line:
                t = float(line.split()[-1])
                ts[L_lead] = t
            elif 't_contactR' in line:
                t = float(line.split()[-1])
                ts[L_lead+L_device] = t
    return ts

def get_data (fname):
    #dmrgdir = get_para (fname, 'input_dir', str)
    #dmrgfile = glob.glob ('../gs/*.out')[0]
    L_lead = get_para (fname, 'L_lead', int)
    L_device = get_para (fname, 'L_device', int)
    L = 2*L_lead + L_device
    Nstep = get_para (fname, 'step =', int, last=True)
    idevL, idevR = get_para (fname, 'device site', int, n=2)
    iC = get_para (fname, 'charge site', int)
    hopt = get_hop_t (fname)

    Is = np.full ((Nstep,L), np.nan)
    Is_spec = np.full ((Nstep,L-1), np.nan)
    ns = np.full ((Nstep,L+1), np.nan)
    Ss = np.full ((Nstep,L+1), np.nan)
    dims = np.full ((Nstep,L), np.nan)
    with open(fname) as f:
        for line in f:
            line = line.lstrip()
            if line.startswith ('step ='):
                tmp = line.split()
                step = int(tmp[-1])
            elif line.startswith('*I'):
                tmp = line.split()
                ilink = int(tmp[-3])

                #I = float(tmp[-1].strip('()').split(',')[0])
                I = float(tmp[-1])
                cc = 2*pi * hopt[ilink]
                Is_spec[step-1,ilink-1] = I * cc
                #Is[step-1,ilink-1] = I * cc
            elif line.startswith('*den '):
                tmp = line.split()
                i = int(tmp[1])
                n = float(tmp[-1])
                ns[step-1,i-1] = n
            elif line.startswith('*entS'):
                tmp = line.split()
                i = int(tmp[1])
                S = float(tmp[-1])
                Ss[step-1,i-1] = S
            elif line.startswith('*m'):
                tmp = line.split()
                i = int(tmp[1])
                m = float(tmp[-1])
                dims[step-1,i-1] = m
    return Is, Is_spec, ns, Ss, dims

def plot_prof (ax, data, dt, label=''):
    Nstep, L = np.shape(data)
    sc = ax.imshow (data, origin='lower', extent=[1, L, dt, Nstep*dt], aspect='auto')
    cb = pl.colorbar (sc)
    ax.set_xlabel ('site')
    ax.set_ylabel ('time')
    cb.ax.set_ylabel (label)

def plot_time_slice (ax, data, n, xs=[], label='', **args):
    Nstep, L = np.shape(data)
    itv = Nstep // n
    if xs == []:
        xs = range(1,L+1)
    for d in data[::itv,:]:
        ax.plot (xs, d, **args)
    ax.plot (xs, data[-1,:], label=label, **args)

def get_basis (fname):
    ens, segs = [],[]
    with open(fname) as f:
        for line in f:
            if 'orbitals, segment, ki, energy' in line:
                for line in f:
                    tmp = line.split()
                    if not tmp[0].isdigit():
                        return np.array(ens), np.array(segs)
                    ens.append (float(tmp[3]))
                    segs.append (tmp[1])

def extrap_current (ts, Il, Ir, plot=False):
    n = 100
    ts = np.reciprocal(ts)[-n:]
    Il = Il[-n:]
    Ir = Ir[-n:]
    if plot:
        f,ax = pl.subplots()
        ax.plot (ts, Il, marker='.', ls='None', label='left')
        ax.plot (ts, Ir, marker='.', ls='None', label='right')
        fitx, fity, stddev, fit = ff.myfit (ts, Il, order=1, ax=ax, refit=True)
        fitx, fity, stddev, fit = ff.myfit (ts, Ir, order=1, ax=ax, refit=True)


if __name__ == '__main__':
    files = [i for i in sys.argv[1:] if i[0] != '-']
    fi,axi = pl.subplots()
    for fname in files:
        print (fname)

        en_basis, segs = get_basis (fname)

        # Get data
        Is, Is_spec, ns, Ss, dims = get_data (fname)
        dt = get_para (fname, 'dt', float)
        m = get_para (fname, 'Largest link dim', int)
        Nstep, L = np.shape(Is)

        # I profile
        '''f,ax = pl.subplots()
        plot_prof (ax, Is_spec, dt, 'current')
        ax.set_title ('$m='+str(m)+'$')
        ps.set(ax)

        f8,ax8 = pl.subplots()
        plot_time_slice (ax8, Is, n=3)
        ps.set(ax8)'''

        # n profile
        '''f2,ax2 = pl.subplots()
        plot_prof (ax2, ns, dt, 'density')
        ax2.set_title ('$m='+str(m)+'$')
        ps.set(ax2)'''

        f,ax = pl.subplots()
        ii = segs == 'L'
        plot_time_slice (ax, ns[:,ii], n=5, marker='.', ls='None', label='L', xs=en_basis[ii])
        ii = segs == 'R'
        plot_time_slice (ax, ns[:,ii], n=5, marker='x', ls='None', label='R', xs=en_basis[ii])
        ii = segs == 'S'
        plot_time_slice (ax, ns[:,ii], n=5, marker='+', ls='None', label='S', xs=en_basis[ii])
        for x in np.where (ii)[0]:
            ax.axvline (en_basis[x], ls='--', c='gray', alpha=0.5)
        ii = segs == 'C'
        plot_time_slice (ax, ns[:,ii], n=5, marker='*', ls='None', label='C', xs=en_basis[ii])
        #ax.plot (range(1,len(en_basis)+1), en_basis)
        ax.set_xlabel ('energy')
        ax.set_ylabel ('occupasion')
        ax.legend()
        ps.set(ax)

        f,ax = pl.subplots()
        sites = np.array(range(1,L+2))
        ii = segs == 'L'
        plot_time_slice (ax, ns[:,ii], n=5, marker='.', ls='None', label='L', xs=sites[ii])
        ii = segs == 'R'
        plot_time_slice (ax, ns[:,ii], n=5, marker='x', ls='None', label='R', xs=sites[ii])
        ii = segs == 'S'
        plot_time_slice (ax, ns[:,ii], n=5, marker='+', ls='None', label='S', xs=sites[ii])
        ii = segs == 'C'
        plot_time_slice (ax, ns[:,ii], n=5, marker='*', ls='None', label='C', xs=sites[ii])
        #ax.plot (range(1,len(en_basis)+1), en_basis)
        ax.set_xlabel ('site')
        ax.set_ylabel ('occupasion')
        ax.legend()
        ps.set(ax)

        # S profile
        '''f5,ax5 = pl.subplots()
        plot_prof (ax5, Ss, dt, 'entropy')
        ax5.set_title ('$m='+str(m)+'$')
        ps.set(ax5)

        f6,ax6 = pl.subplots()
        plot_time_slice (ax6, Ss, n=3)
        ax6.set_xlabel ('site')
        ax6.set_ylabel ('entropy')
        ps.set(ax6)'''

        f,ax = pl.subplots()
        plot_time_slice (ax, dims, n=5)
        ax.set_xlabel ('site')
        ax.set_ylabel ('bond dimension')
        ps.set(ax)

        ts = dt * np.arange(1,Nstep+1)
        idevL, idevR = get_para (fname, 'device site', int, n=2)
        print (idevL, idevR)
        #
        '''N_dev = np.sum (ns[:, idevL-1:idevR], axis=1)
        f3,ax3 = pl.subplots()
        ax3.plot (ts, N_dev, marker='.')
        ps.set(ax3)'''

        # current vs time
        muL = get_para (fname, 'mu_biasL', float)
        muR = get_para (fname, 'mu_biasR', float)
        Vb = muR - muL
        Il = Is_spec[:, idevL-2] / Vb
        Ir = Is_spec[:, idevR-2] / Vb
        axi.plot (ts, Il, label='left')
        axi.plot (ts, Ir, label='right')
        #axi.plot (ts, Is[:, idevL-3] / Vb, marker='.', label='l2')
        #axi.plot (ts, Is[:, idevL] / Vb, marker='.', label='r2')
        #extrap_current (ts, Il, Ir, plot=True)
        axi.set_xlabel ('time')
        axi.set_ylabel ('conductance')
        axi.legend()

        '''dmrgdir = get_para (fname, 'input_dir', str)
        dmrgfile = glob.glob (dmrgdir+'/out')[0]
        tp = get_para (dmrgfile, 't_contactR', float)
        axi.axhline (-exactG(0, tp), ls='--', c='gray')'''

        ps.set(axi)

        if '-pdf' in sys.argv:
            filename = fname.replace('.out','')+'_I.pdf'
            fi.savefig (filename)
        if '-save' in sys.argv:
            with open(fname.replace('.out','')+'_I.txt','w') as f:
                print ('t {link}',file=f)
                for i in range(len(ts)):
                    print (ts[i],*(Is_spec[i,:]/Vb),file=f)

    pl.show()
