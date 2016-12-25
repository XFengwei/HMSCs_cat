import os
import time
import numpy as np
from astropy import wcs
from astropy import coordinates as coords
from astropy import constants as cons
from astropy import units as u
from astropy.table import Column, Table, join
from astropy.io import ascii, fits
from matplotlib import pyplot as plt

#some constants
mH = cons.m_n
pc = cons.pc
mSun = cons.M_sun

hmscList = ascii.read('../../Tables/hmscList_201612.txt')
disB = ascii.read('../../Tables/BaysianDis.txt')

sedDir = '../../fitsDir/sedFitting/allSixMin/'
hmscList['amin'][hmscList['amin']<20.0] = 20.0
hmscList['amaj'][hmscList['amaj']<21.0] = 21.0

amaj_36 = np.sqrt(hmscList['amaj']**2+36.4**2-19.2**2)/2
amin_36 = np.sqrt(hmscList['amin']**2+36.4**2-19.2**2)/2
#amaj_36 = Column([20/2]*len(hmscList))
#amin_36 = Column([20/2]*len(hmscList))
angSize = np.sqrt(hmscList['amin']*hmscList['amaj'] - 19.5**2)
angSize.name = 'fwhm_eq'
amaj_36.name = 'amaj_36'
amin_36.name = 'amin_36'

hmscList.add_columns([amaj_36, amin_36,angSize])


disB['ExtraInformation'].name = 'Name'
disB['Dist'].name = 'Dist_B'
disB = disB['Name', 'Dist_B', 'sigma_D']

hmscList = join(hmscList, disB, join_type = 'left')

# pixel solid area
pixSr = 11.5**2*u.arcsec**2

massFactor = pixSr.cgs.value*mH*2.8*u.kpc**2/u.cm**2
lumFactor  = 4*np.pi*pixSr.cgs.value*u.kpc**2*u.g/u.s/u.s*u.Hz


# define null lists for some parameters
Tdust     = []
Nh2       = []
mass_1kpc = []
Lum_1kpc  = []
start = time.time()
pfmt = '%6i %18s %6.2f%% finished.'

for isour in range(len(hmscList)):
    sourName = hmscList['Name'][isour]
    sName    = hmscList['sName'][isour]
    ra_0 = hmscList['ra'][isour]
    dec_0 = hmscList['dec'][isour]
    Nh2File    = sedDir+sName+'_Nh2.fits'
    TdustFile  = sedDir+sName+'_Tdust.fits'
    intSEDFile = sedDir+sName+'_intSED.fits'
         
    a1     =  hmscList['amaj_36'][isour]/11.5
    b1     =  hmscList['amin_36'][isour]/11.5
    a2     =  hmscList['amaj_36'][isour]/11.5+0.5
    b2     =  hmscList['amin_36'][isour]/11.5+0.5
    paG    =  hmscList['PA'][isour]*u.deg
    c1 = coords.SkyCoord(ra_0,dec_0,frame = 'fk5',
                         unit = (u.deg,u.deg))
    c2 = coords.SkyCoord(ra_0,dec_0+0.01,frame = 'fk5',
                         unit = (u.deg,u.deg))
    c1G = c1.galactic
    ang = c1G.position_angle(c2.galactic).deg
    pa  = paG+(90-ang)*u.deg
    
    paRad  = pa.cgs.value
    hduNh2    = fits.open(Nh2File)
    hduTdust  = fits.open(TdustFile)
    hduIntSed = fits.open(intSEDFile)
    
    xCenter = np.int(hduNh2[0].header['CRPIX1']-1)
    yCenter = np.int(hduNh2[0].header['CRPIX2']-1)
    
    w = wcs.WCS(hduNh2[0].header)
    centerWorld = np.array([[ra_0, dec_0]])
    centerPix = w.wcs_world2pix(centerWorld, 0)
    x0, y0 = centerPix[0][0], centerPix[0][1]
    y,x = np.indices(hduNh2[0].data.shape)
    inside1 = ((((y-y0)*np.sin(paRad)+(x-x0)*np.cos(paRad))/b1)**2
               +(((y-y0)*np.cos(paRad)-(x-x0)*np.sin(paRad))/a1)**2
               <1)
    inside2 = ((((y-y0)*np.sin(paRad)+(x-x0)*np.cos(paRad))/b2)**2
               +(((y-y0)*np.cos(paRad)-(x-x0)*np.sin(paRad))/a2)**2
               <1)
    
    intInside1 = Column(hduIntSed[0].data[inside1].flatten())
    nh2Inside1 = Column(hduNh2[0].data[inside1].flatten())
    intInside2 = Column(hduIntSed[0].data[inside2].flatten())
    nh2Inside2 = Column(hduNh2[0].data[inside2].flatten())
    #sumInt = np.mean([sum(intInside1[intInside1 == intInside1]),
    #                sum(intInside2[intInside2 == intInside2])])*lumFactor
    #sumNh2 = np.mean([sum(nh2Inside1[nh2Inside1 == nh2Inside1]), 
    #                sum(nh2Inside2[nh2Inside2 == nh2Inside2])])*massFactor
    sumInt = sum(intInside1[intInside1 == intInside1])*lumFactor
    sumNh2 = sum(nh2Inside1[nh2Inside1 == nh2Inside1])*massFactor
    Tdust.append(hduTdust[0].data[yCenter,xCenter])
    Nh2.append(hduNh2[0].data[yCenter,xCenter])
    if len(intInside1[intInside1 == intInside1]) > 0:
        mass_1kpc.append(sumNh2.to(u.M_sun).value)
        Lum_1kpc.append(sumInt.to(u.L_sun).value)
    else:
        mass_1kpc.append(np.nan)
        Lum_1kpc.append(np.nan)
        
        
    print(pfmt %(isour+1, sourName, ((isour+1.0) / len(hmscList)*100)))
    
Nh2 = Column(Nh2, name = 'Nh2')
Tdust = Column(Tdust, name = 'Tdust')
lum = Column(Lum_1kpc, name = 'Lum_1kpc')
mass = Column(mass_1kpc, name = 'mass_1kpc')

L_M_ratio = lum/mass
L_M_ratio.name = 'L_M_ratio'
M_clump = mass*hmscList['Dist_B']**2
L_clump = lum*hmscList['Dist_B']**2

#r_pc_eq = hmscList['fwhm_eq']/3600.0/180*np.pi*hmscList['dis']*1000


r_pc_eq = np.sqrt((hmscList['amaj']*hmscList['amin']-19.2**2))/3600.0/180*np.pi*hmscList['Dist_B']*1000
nH2 = M_clump*mSun.cgs.value/(4/3.0*np.pi*(r_pc_eq*pc.cgs.value)**3)/2.8/mH.cgs

M_clump.name = 'Mclump'
L_clump.name = 'Lclump'
r_pc_eq.name = 'r_pc'
nH2.name = 'n_H2'

hmscList.add_columns([Nh2, Tdust, mass, lum, 
                      M_clump, L_clump, L_M_ratio, r_pc_eq, nH2])

hmscList.write('../../Tables/hmscList_full_20161218.txt', 
                format = 'ascii.ipac')

ending = time.time()
print((ending-start)/60.0, 'minutes have been used.')