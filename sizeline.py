sradius=False # use solomon radius parameters - if not, use sigma_x
ploterr=False # plot error bars

import pickle
import numpy as np
import matplotlib.pyplot as pl
pl.ion()

data=pickle.load(open("data/sizeline_data.pkl","rb"))

# can't resize existing windows in osx:
# https://github.com/matplotlib/matplotlib/issues/15131
pl.close(0)
pl.figure(0,figsize=(6,8))

pl.subplot(211)

if sradius:
    s=1.92
else:
    s=1.
    
for k,d in data.items():
    x = s*np.array(d['sigx_deconv'])  
    y = np.array(d['sigv_deconv'])
    if len(x)<500:
        fmt='.'
    else:
        fmt=','
    if "e_sigx_deconv" in d and "e_sigv_deconv" in d and ploterr:
        dx = s*np.array(d['e_sigx_deconv'])
        dy = np.array(d['e_sigv_deconv'])
        myplot,=pl.plot(x,y,fmt,label=k)
        pl.errorbar(x,y,xerr=dx,yerr=dy,fmt=',',color=myplot.get_color(),alpha=0.25)
    else:
        myplot,=pl.plot(x,y,fmt,label=k)
    x = s*np.array(d['sigx'])
    y = np.array(d['sigv'])
    pl.plot(x,y,fmt,color=myplot.get_color(),alpha=0.25)

if sradius:
    pl.xlabel(r'R=1.92*$\sigma_x$ [pc]')
else:
    pl.xlabel(r'$\sigma_x$ [pc]')
pl.ylabel(r'$\sigma_v$ [km/s]')
pl.xscale("log")
pl.yscale("log")
pl.legend(loc="best",prop={"size":8})




# grav param plot
pl.subplot(212)

for k,d in data.items():
    if 'MLTE' in d:
        x = d['MLTE']/np.pi/(s*np.array(d['sigx_deconv']))**2
        y = np.array(d['sigv_deconv']) / (s*np.array(d['sigx_deconv']))**0.5

        if len(x)<500:
            fmt='.'
        else:
            fmt=','
#        if "e_sigx_deconv" in d and "e_sigv_deconv" in d and ploterr:

        myplot,=pl.plot(x,y,fmt,label=k)
        x = d['MLTE']/np.pi/(s*np.array(d['sigx_deconv']))**2
        y = np.array(d['sigv_deconv']) / (s*np.array(d['sigx_deconv']))**0.5
        pl.plot(x,y,fmt,color=myplot.get_color(),alpha=0.25)
    else:
        pl.plot(pl.xlim()[1],pl.ylim()[1],',') # to match colors w/upper plot

if sradius:
    pl.xlabel(r'$\Sigma$ = M/R$^2$ [Msun/pc$^2$]')
    pl.ylabel(r'$\sigma_v/\sqrt{R}$ [km/s/sqrt(pc)]')
else:
    pl.xlabel(r'$\Sigma$ = M/$\sigma_x^2$ [Msun/pc$^2$]')
    pl.ylabel(r'$\sigma_v/\sqrt{\sigma_x}$ [km/s/sqrt(pc)]')
pl.xscale("log")
pl.yscale("log")
pl.legend(loc="best",prop={"size":8})

coeff = 0.052 # sigv/sqrt(R) = coeff [km/s/pc0.5] * (Sigma [Msun/pc2])^0.5
xx=np.array(pl.xlim())
pl.plot(xx,coeff* xx**0.5,'k')




pl.subplots_adjust(top=0.95)
