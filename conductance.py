import analysis as gt
import plotsetting as ps
import pylab as pl
import sys, glob, os
import numpy as np
from scipy.interpolate import interp1d

def get_current (ts, IL, IR, ax, tbeg, tend):
    Il, errl = gt.get_average_current (ts, IL, ax=ax, tbeg=tbeg, tend=tend, marker='.', label='left contact')
    Ir, errr = gt.get_average_current (ts, IR, ax=ax, tbeg=tbeg, tend=tend, marker='.', label='right contact')
    I = 0.5*(Il+Ir)
    err = 0.5*(errl+errr)
    return I, err

def get_Vb (filename):
    muL = gt.get_para (filename, 'mu_leadL', float)
    muR = gt.get_para (filename, 'mu_leadR', float)
    return muR - muL

def sort_by_first (list1, list2):
    re = zip(*sorted(zip(list1, list2)))
    return map(list,re)

if __name__ == '__main__':
    '''dirrs = [i for i in sys.argv[1:] if i[0] != '-']
    if len(dirrs) == 0:
        print ('Please give the directory')
        exit()
    dirr = dirrs[0]'''
    dirr = '.'
    if not os.path.isdir (dirr):
        print ('Directory not found:',dirr)
        exit()

    files = glob.glob(dirr+'/Vb*.out')
    Vbs = [get_Vb(filename) for filename in files]
    Vbs, files = sort_by_first (Vbs, files)
    #Vbs = sorted ([float(i.split('Vb')[-1]) for i in glob.glob(dirr+'/Vb*')])
    tbeg,tend = 100,200

    Is, errs, Gs, errGs = [],[],[],[]
    for Vb,fname in zip(Vbs,files):
        #fname = 'Vb'+str(Vb)+'/Vb'+str(Vb)+'_m200.out'
        print (fname)
        I_time_space, xss, tss, Iss, nss, ts, IL, IR, Imean, II = gt.get_data (fname)

        if '-fast' in sys.argv:
            I, err = get_current (ts, IL, IR, ax=None, tbeg=tbeg, tend=tend)
        else:
            f1,ax1 = pl.subplots()
            Idata, tmin, tmax, xmin, xmax = gt.to_imag_data (xss, tss, Iss)
            sc = ax1.imshow (Idata, origin='lower', extent=[xmin, xmax, tmin, tmax], aspect='auto')
            #sc = ax1.scatter (xss, tss, c=Iss)
            cb = pl.colorbar (sc)
            ax1.set_xlabel ('site')
            ax1.set_ylabel ('time')
            cb.ax.set_ylabel ('current')
            ps.set (ax1)
            if '-pdf' in sys.argv:
                f1.savefig (fname.replace('.out','_I_tx.pdf'))
            pl.close (f1)

            f,ax = pl.subplots()
            I, err = get_current (ts, IL, IR, ax, tbeg, tend)
            print ('current =', I, err)
            ax.set_xlabel ('time')
            ax.set_ylabel ('current')
            ax.legend()
            ps.set (ax)
            if '-pdf' in sys.argv:
                f.savefig (fname.replace('.out','_I_t.pdf'))
            pl.close (f)
        Is.append (I)
        errs.append (err)

        Gs.append (I/Vb)
        errGs.append (err/Vb)

    # Plot current
    f,ax = pl.subplots()
    ax.errorbar (Vbs, Is, errs, marker='.', label='I')
    ax.set_xlabel ('$V_b$')
    ax.set_ylabel ('conductance')
    # Interpolate
    func = interp1d (Vbs, Is, kind='cubic')
    fitx = np.linspace (min(Vbs), max(Vbs), 100)
    fity = func (fitx)
    ax.plot (fitx, fity, label='I interp')
    # Plot dI/dVb
    slopex = [0.5*(fitx[i]+fitx[i-1]) for i in range(1,len(fitx))]
    slopey = [(fity[i]-fity[i-1])/(fitx[i]-fitx[i-1]) for i in range(1,len(fitx))]
    ax.plot (slopex, slopey, label='dI/dG')
    # Plot I/Vb
    ax.errorbar (Vbs, Gs, errGs, marker='.', label='I/V')
    #
    mu_device = gt.get_para (fname, 'mu_device', float)
    ax.axvline (mu_device, ls='--', c='gray')

    ax.legend()
    ps.set (ax)

    tmp = ''
    if 'ground' in dirr: tmp = '_grnd'
    elif 'symm' in dirr: tmp = '_symm'
    if '-pdf' in sys.argv:
        f.savefig (dirr+'/G_mup'+str(mu_device)+'_tp0.25'+tmp+'.pdf')
    if '-save' in sys.argv:
        with open(dirr+'/G_mup'+str(mu_device)+'_tp0.25'+tmp+'.dat','w') as f:
            print ('Vb I err G errG',file=f)
            for Vb,I,err,G,errG in zip(Vbs,Is,errs,Gs,errGs):
                print (Vb, I, err, G, errG, file=f)
        with open(dirr+'/diffG_mup'+str(mu_device)+'_tp0.25'+tmp+'.dat','w') as f:
            print ('Vb diffG',file=f)
            for Vb,G in zip(slopex,slopey):
                print (Vb,G,file=f)
    pl.show()
