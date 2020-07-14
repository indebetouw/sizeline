from astropy.table import Table
import numpy as np
from astropy import units as u
import pickle

# TODO uncertainties?

sdat={}

heyer=Table.read('heyer09.txt',format='ascii')  # Table 1 from the paper

# beam in pc depends on distance D in kpc
sigbm_pc = 46. *heyer['D']*1000/206265 / 2.355 * u.pc

sdat['Heyer09']={
    'sigx':heyer['R']/1.92 *u.pc, # second spatial moment, pc
    'sigx_deconv':np.sqrt( (heyer['R']/1.92 *u.pc)**2 - sigbm_pc**2 ),
    'sigv':heyer['sigv'] *u.km/u.s, # second velocity moment, km/s
    'sigv_deconv':np.sqrt( heyer['sigv']**2 - (0.21/2.355)**2 ) *u.km/u.s,
    'line':'13CO 1-0', # line used to weight the intensity moments
    'bmaj':46.*u.arcsec, 'bmin':46.*u.arcsec, 'bpa':0.*u.deg,
    'segmentation':'12CO 1-0 SRBY ZI rectangles', # regions used
    'ref':"https://ui.adsabs.harvard.edu/abs/2009ApJ...699.1092H/abstract",
    'LCO':heyer['Lco'] *u.K *u.km/u.s *u.pc**2, # K km/s pc2
    'MLTE':heyer['Mlte'] *u.Msun
    }

sdat['Heyer09 13CO cores']={
    'sigx':heyer['R2']/1.92 *u.pc, # second spatial moment
    'sigx_deconv':np.sqrt( (heyer['R2']/1.92 *u.pc)**2 - sigbm_pc**2 ),
    'sigv':heyer['sigv2'] *u.km/u.s,
    'sigv_deconv':np.sqrt( heyer['sigv']**2 - (0.21/2.355)**2 ) *u.km/u.s,
    'line':'13CO 1-0', # line used to weight the intensity moments
    'bmaj':46.*u.arcsec, 'bmin':46.*u.arcsec, 'bpa':0.*u.deg, # arcsec
    'segmentation':'N(H2) half-max isophot',
    'ref':"https://ui.adsabs.harvard.edu/abs/2009ApJ...699.1092H/abstract",
    'LCO':heyer['Lco2'] *u.K *u.km/u.s *u.pc**2, # K km/s pc2
    'MLTE':heyer['Mlte2'] *u.Msun# Msun
    }




#====================================================================

gc=Table.read('oka_galcen.txt',format='ascii')  # Table 1 from the paper
#i,l,b,v,sigl,sigb,s,sigv,tpk,tmin,lco,mvir,alpha,flag,$
z=np.where(gc['f_Num']=='.')[0] # exclude crowded and confused clouds
sigbm_pc = 34.*8500/206265/2.355 *u.pc
S=np.array(gc['Size'][z])*gc['Size'].unit
sv=np.array(gc['e_VLSR'][z])*gc['e_VLSR'].unit

sdat['Oka 01']={
    'sigx':S, # sqrt(sigl*sigb)
    'sigx_deconv':np.sqrt( S**2 - sigbm_pc**2 ),
    'sigv':sv,
    'sigv_deconv':np.sqrt( sv**2 - (0.21/2.355 *u.km/u.s)**2 ),
    'line':'12CO 1-0',
    'bmaj':34.*u.arcsec, 'bmin':34.*u.arcsec, 'bpa':0.*u.deg, # arcsec 
    'segmentation':'12CO 1-0 5-10K isophot',
    'ref': "https://ui.adsabs.harvard.edu/abs/2001ApJ...562..348O/abstract",
    'LCO': np.array(gc['LumCO'][z])*gc['LumCO'].unit,
    'Mvir':np.array(gc['MassVT'][z])*gc['MassVT'].unit
    }


# oka GC 3-2:
# https://ui.adsabs.harvard.edu/abs/2006JPhCS..54...67O/abstract
# https://ui.adsabs.harvard.edu/abs/2007PASJ...59...15O/abstract
# https://ui.adsabs.harvard.edu/abs/2012ApJS..201...14O/abstract




#====================================================================
# COMPLETE
#shet=Table.read('shetty12.perseus.sizeline.csv',format='ascii')
#pl.plot(shet['col1'].data,shet['col2'],'k.',label="Shetty12 Perseus dendro")


for region in "PerA","OphA":
    root="trim"
    line="13co"
    t=Table.read("complete/"+region+"_"+line+"_physprop.txt",format="ascii")
    if region=="PerA":
        sigbm_pc=0.03/2.355 *u.pc # Per is 250pc
    else:
        sigbm_pc=0.06/2.355 *u.pc # Oph is 131 pc 
        
    sx=np.array(t['rad_pc'])/1.92 *u.pc
    sv=np.array(t['vrms_k']) *u.km/u.s
    
    sdat[region+' 13CO cores']={
        'sigx':sx,
        'e_sigx':t['e_rad_pc']/1.92,
        'sigx_deconv':np.sqrt( sx**2 - sigbm_pc**2 ),
        'e_sigx_deconv':t['e_rad_pc']/1.92,
        'sigv':sv,
        'e_sigv':t['e_vrms_k'],
        'sigv_deconv':np.sqrt( sv**2 - (0.07/2.355 *u.km/u.s)**2 ),
        'e_sigv_deconv':t['e_vrms_k'],
        'line':'13CO 1-0',
        'bmaj':46.*u.arcsec, 'bmin':46.*u.arcsec, 'bpa':0.*u.deg, # arcsec 
        'segmentation':'dendrogram',
        'ref': "https://ui.adsabs.harvard.edu/abs/2001ApJ...562..348O/abstract",
        'MLTE':t['mlte'],
        'e_MLTE':t['e_mlte'],
        'Mvir':t['mvir'],
        'e_Mvir':t['e_mvir']}
 


#------

import pickle
pickle.dump(sdat,open("sizeline_data.pkl","wb"))

