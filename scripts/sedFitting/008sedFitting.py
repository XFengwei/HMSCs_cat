import glob
from astropy.io import fits
from astropy import units as u
from higal_sedfitter import smooth, fit, higal_beams
from higal_sedfitter.fit import PixelFitter

# 70, 160, 250, 350 um data downloaded from https://irsa.ipac.caltech.edu/Missions/planck.html
# 870 um data from https://atlasgal.mpifr-bonn.mpg.de/cgi-bin/ATLASGAL_DATASETS.cgi

# renew header of fits file
for Herschel_file in glob.glob('Herschel_data/*.fits'):
    hdu = fits.open(Herschel_file)
    # download data saved in 'Herschel_data' directory and new data saved in 'Herschel_data/Herschel_newhd' directory
    # all files are subfixed by '_newhd.fits'
    fits.writeto('Herschel_data/Herschel_newhd'+Herschel_file[13:-5]+'_newhd.fits',header=hdu[1].header, data=hdu[1].data, overwrite=True)

# add beams info to the header
for fn in glob.glob('Herschel_data/Herschel_newhd/Herschel*_newhd.fits'):
    smooth.add_beam_information_to_higal_header(fn, name_to_um=higal_beams.num_to_um, clobber=True)

# since the longest wavelength corresponds to lowest resolution
# set 500 um as target
target_fn = glob.glob('Herschel_data/Herschel_newhd/Herschel_SPIRE500_newhd.fits')[0]
target_header = fits.getheader(target_fn)

# smooth other images to the resolution of 500 um
# output files are subfixed by 'smooth45.fits' and 'smregrid45.fits' in work directory
smooth.smooth_images_toresolution(45*u.arcsec, skip_existing=False,
                                  globs=['Herschel_data/Herschel_newhd/*070_newhd.fits',
                                         'Herschel_data/Herschel_newhd/*160_newhd.fits',
                                         'Herschel_data/Herschel_newhd/*250_newhd.fits',
                                         'Herschel_data/Herschel_newhd/*350_newhd.fits',
                                         'Herschel_data/Herschel_newhd/*500_newhd.fits'],
                                  regrid=True,
                                  target_header=target_header, clobber=True)

# SED fit
import os
import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import quad
from astropy.io import fits
from astropy.table import Table, Column, vstack
from astropy import units as u
import astropy.coordinates as coords
import astropy.constants as cons
from lmfit import minimize, Parameters, Model

# Frequency interval for integrating the SED
wavLower = 1*u.um
wavUpper = 3000*u.um

freqLower = wavUpper.to(u.Hz, u.spectral()).value
freqUpper = wavLower.to(u.Hz, u.spectral()).value

# Initiate the parameters
iNh2 = 1.0E+22 # in cm^-2
iTdust = 22.0 # in K.
ibeta = 2.0 # 
betaVary = False # If True, Beta will be fitted. If False, beta will be fixed.
betaMin = 1.0
betaMax = 3.0

# 70 micron is left out
wavelengths = np.array([160,250,350,500]) * u.um
frequencies = wavelengths.to(u.Hz, u.spectral()).value.copy()

# Constants
h = cons.h.cgs.value   # Planck constant in CGS unit
k = cons.k_B.cgs.value # Boltzmann constant in CGS unit
c = cons.c.cgs.value # speed of light in CGS unit
mH = cons.m_p.cgs.value # mass of an neutron
muh2 = 2.8 # mean molecular weight adopted from Kauffmann et al. (2008)
rGD = 100.0 # gas-to-dust mass ratio
nu0 = 599.584916E9 # Reference frequency in Hz.
kappa0 = 5.0/1.5 # Dust emissivity at reference frequency

# Define the model
def greybody(nu, nH2=iNh2, Tdust=iTdust, beta=ibeta):

    blackBody = 2*h*nu**3/c**2/(np.exp(h*nu/k/Tdust)-1)
    tau = muh2*mH*nH2*kappa0*(nu/nu0)**beta/rGD
    return blackBody*(1-np.exp(-tau))

gMod = Model(greybody)
gMod.set_param_hint('beta', min = betaMin, max = betaMax, vary = betaVary)
gMod.set_param_hint('Tdust',  max = 80.0)

pars = gMod.make_params()

pfmt = '%6i %18s %6.2f%% finished.'
start = time.time()

# Enter the grandchild dir, choose your own work directory
Herschel_dir = './Herschel_data/Herschel_newhd/'
os.chdir(Herschel_dir)

# Run the model
hdrRaw   = fits.getheader('Herschel_PACS160_newhd_smregrid45.fits') 
hdrNh2   = hdrRaw.copy()
hdrInt  = hdrRaw.copy()
hdrTdust = hdrRaw.copy()
hdrTdust['BUNIT'] = 'K'
hdrInt['BUNIT'] = 'g s^-3'
hdrNh2['BUNIT'] = 'cm^-2'

# get the beam size
# unit is degree
bmaj = fits.getheader('Herschel_SPIRE500_newhd_smregrid45.fits')['BMAJ']
bmin = fits.getheader('Herschel_SPIRE500_newhd_smregrid45.fits')['BMIN']
beam_factor = (bmaj*bmin*60**2)*5.32E-5

# get the 2d data
# unit is MJy/sr
data160 = fits.getdata('Herschel_PACS160_newhd_smregrid45.fits')/beam_factor
data250 = fits.getdata('Herschel_SPIRE250_newhd_smregrid45.fits')/beam_factor
data350 = fits.getdata('Herschel_SPIRE350_newhd_smregrid45.fits')/beam_factor
data500 = fits.getdata('Herschel_SPIRE500_newhd_smregrid45.fits')/beam_factor

# data error constants
errData = np.array([0.076, 5.997, 4.706, 3.733])/beam_factor*u.MJy/u.sr
weis = 1.0/errData

# Initiate the 2d array of three parameters
dataNh2 = data500.copy()
dataTdust = data500.copy()
dataIntegrate = data500.copy()

# For loop
for iCol in range(len(data500[:,1])):
    for iRow in range(len(data500[1,:])):
        data = np.array([data160[iCol,iRow], data250[iCol,iRow], 
                         data350[iCol,iRow], data500[iCol,iRow]])*u.MJy/u.sr
        if ((data160[iCol, iRow] > 0) 
            and (data250[iCol, iRow] > 0) 
            and (data350[iCol, iRow] > 0) 
            and (data500[iCol, iRow] > 0)):
            sedResult = gMod.fit(data.cgs.value, nu = frequencies,weights = weis.cgs.value)
            inTSED = quad(greybody,freqLower, freqUpper, 
                          args = (sedResult.params['nH2'].value, sedResult.params['Tdust'].value, 
                                  sedResult.params['beta'].value))
            dataNh2[iCol,iRow] = sedResult.params['nH2'].value
            dataTdust[iCol, iRow] = sedResult.params['Tdust'].value
            dataIntegrate[iCol,iRow] = inTSED[0]
        elif ((data160[iCol, iRow] <= 0) 
            and (data250[iCol, iRow] > 0.0) 
            and (data350[iCol, iRow] > 0.0) 
            and (data500[iCol, iRow] > 0.0)):
            sedResult = gMod.fit(data.cgs.value[1:], nu = frequencies[1:],weights = weis.cgs.value[1:])
            inTSED = quad(greybody,freqLower, freqUpper, 
                          args = (sedResult.params['nH2'].value, sedResult.params['Tdust'].value, 
                                  sedResult.params['beta'].value))
            dataNh2[iCol,iRow] = sedResult.params['nH2'].value
            dataTdust[iCol, iRow] = sedResult.params['Tdust'].value
            dataIntegrate[iCol,iRow] = inTSED[0]
        elif ((data160[iCol, iRow] > 0) 
              and (data250[iCol, iRow] > 0.0) 
              and (data350[iCol, iRow] > 0.0) 
              and (data500[iCol, iRow] <= 0.0)):
            sedResult = gMod.fit(data.cgs.value[:-1], nu = frequencies[:-1],weights = weis.cgs.value[:-1])
            inTSED = quad(greybody,freqLower, freqUpper, 
                          args = (sedResult.params['nH2'].value, sedResult.params['Tdust'].value, 
                                  sedResult.params['beta'].value))
            dataNh2[iCol,iRow] = sedResult.params['nH2'].value
            dataTdust[iCol, iRow] = sedResult.params['Tdust'].value
            dataIntegrate[iCol,iRow] = inTSED[0]
        else:
            dataNh2[iCol,iRow] = np.nan
            dataTdust[iCol,iRow] = np.nan
            dataIntegrate[iCol,iRow] = np.nan

# write the three parameters into fits files
fits.writeto('Nh2.fits', data = dataNh2, 
    header = hdrNh2, clobber = True)
fits.writeto('Tdust.fits', data = dataTdust, 
    header = hdrTdust, clobber = True)
fits.writeto('intSED.fits', data = dataIntegrate, 
    header = hdrInt, clobber = True)
stop = time.time()
dure = stop - start

# display how much time the code runs
print("Run time = ",dure, "seconds")
# come back to the origin directory
os.chdir('../..')

# visualize the result
import aplpy as ap
# show the distribution of dust temperature
# other two parameters are the same
Tdust_img = ap.FITSFigure('Herschel_data/Herschel_newhd/Tdust.fits')
# the crop the useless boundary
Tdust_img.recenter(x=247.7514630,y=-48.7257175, width=0.14, height=0.14)
# show the colormap
Tdust_img.show_colorscale(stretch='linear',cmap='viridis')
Tdust_img.add_colorbar()
Tdust_img.add_grid()
# give the size scale
# depends on distance of your source
Tdust_img.add_scalebar(0.2*u.pc/(2.92*10**3*u.pc)*180/np.pi*u.deg, corner='top left',color='white')
Tdust_img.scalebar.set_label('0.2 pc')


