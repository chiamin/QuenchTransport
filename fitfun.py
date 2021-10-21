import numpy as np
import pylab as pl
import random as rand
from math import sqrt
import sys

def get_fitxy (x,fit,order):
  if order == 1: fitx = np.append (x,[0.])
  else:
    dx = max(x) / 100.
    fitx = np.arange (0.,max(x)+dx,dx)

  fity = polyf (fit,fitx)

  return fitx, fity

def myfit (x,y,werr=[],fitpts=0,order=1,refit=False,ax=None,pltargs={},addzero=True,interp='auto',addx=[]):
    if interp == 'auto':
        if order == 1: interp = False
        else: interp = True
    x = x[-fitpts:]
    y = y[-fitpts:]
    werr = werr[-fitpts:]
    ind = np.where (np.isnan(y))
    x = np.delete (x, ind)
    y = np.delete (y, ind)
    if len(werr) != 0:
        werr = np.delete (werr, ind)
    fit = polyfit (x,y,err=werr,order=order)
    fitx = x
    if addzero:
        fitx = sorted(np.append (x, [0.]))
    fitx = sorted(np.append (fitx, addx))
    if interp:
        dx = max(x) / 100.
        fitx = np.arange (0.,max(x)+dx, dx)
    #print x,fitx
    fity = polyf (fit,fitx)
    stddev = standard_deviation_fit (x,y,fit)

    if ax != None:
        ax.plot (fitx,fity,**pltargs)

    if refit: return fitx, fity, stddev, fit
    else: return fitx, fity, stddev

def standard_deviation_fit (x,y,fit):
  sigma = 0.
  for xi,yi in zip(x,y):
    sigma += (yi-polyf(fit,xi))**2
  return (sigma / len(y))**0.5

def polyfit (x,y,order,err=[],plot=False,text=False,c='r',returnerr=0,addpoint=[0],fitpts=0):
  x = x[-fitpts:]
  y = y[-fitpts:]
  def fitf(x): return polyf (fit,x)
  if err == []:
    fit = np.polyfit(x,y,order)
  else:
    for index, item in enumerate(err):
      if item == 0.:
        err[index] = sys.float_info.min
    fit = np.polyfit(x,y,order,w=np.array(err)**(-1))
  if plot:
    pl.plot (x,y,marker='o',ls='None')
    fitx,fity = list(x)+addpoint, map(fitf,list(x)+addpoint)
    pl.plot(fitx,fity,c=c)
    if text: pl.text (0.1,0.8,fitfstr(fit),fontsize=20,transform=pl.gca().transAxes)
    pl.show()
  if not returnerr: return fit#,ferr
  else:
    var = 0
    for xi,yi in zip(x,y):
      var += (yi-fitf(xi))**2
    return fit,sqrt(var)

def polyfiterr (x,y,e,order=1,n=100,plot=False,c='r',text=False,fig=None):
  fit = np.polyfit(x,y,order)
  if fig != None: plt = fig
  else: plt = pl
  if plot:
    def fitf(x): return polyf (fit,x)
    fitx,fity = list(x)+[0], map(fitf,list(x)+[0])
    plt.plot(fitx,fity,c=c)
    if text: plt.text (0.1,0.8,fitfstr(fit),fontsize=20,transform=pl.gca().transAxes)
  ferr = fiterr (x,y,e,order,plot=0,n=n,fig=fig)
  #print fit[-1], ferr
  return fit,ferr

def fiterr (x,y,e,order,n=100,plot=False,fig=None,fitpts=0):
  x = x[-fitpts:]
  y = y[-fitpts:]
  e = e[-fitpts:]
  yex,summ,sum2 = [],0.,0.
  for i in range(n):
    yr = [rand.gauss(yi,ei) for yi,ei in zip(y,e)]
    fit = np.polyfit(x,yr,order)
    yex.append (fit[-1]) # extrapolated y
    summ += fit[-1]
    sum2 += fit[-1]*fit[-1]
    if plot:
      if fig != None: plt = fig
      else: plt = pl
      plt.errorbar(x,y,e,marker='o',ls='None',c='k')
      def fitf(x): return polyf (fit,x)
      fitx,fity = list(x)+[0], map(fitf,list(x)+[0])
      plt.plot(fitx,fity,c='r')
  sigma = sum2/float(n) - (summ/float(n))**2
  if abs(sigma) < 1e-12: return 1e-6
  else: return sqrt(sigma)

def polyf (coef,x):
  x = np.array(x)
  y = coef[-1]
  order = range(len(coef))[::-1]
  for i in range(len(coef)-1):
    y += coef[i]*x**order[i]
  return y

def fitfstr (fit):
  text = "{:.4f}".format(fit[-1])
  text += '+' if fit[-2] > 0 else ''
  text += "{:.4f}".format(fit[-2])+'x'
  order = 2
  for coef in fit[-3::-1]:
    text += '+' if coef > 0 else ''
    text += "{:.4f}".format(coef)+'x^'+str(order)
    order += 1
  return text

