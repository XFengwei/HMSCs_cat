
import os
import math
import time
import tarfile
import shutil
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import montage_wrapper as mt
from scipy.integrate import quad
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table, Column, vstack
from astropy import units as u
import astropy.coordinates as coords
import astropy.constants as cons
from lmfit import minimize, Parameters, Model

listFile = '../../Tables/hmscFull_20160405.dat'
sourList = ascii.read(listFile)


sixFITSdir = '../../fitsDir/sedFitting/allSixMin/'

os.chdir(sixFITSdir)

# Frequency interval for integrating the SED
wavLower = 1*u.um
wavUpper = 3000*u.um

freqLower = wavUpper.to(u.Hz, u.spectral()).value
freqUpper = wavLower.to(u.Hz, u.spectral()).value


# Initiate the parameters
iNh2 = 1.0E+22 # in cm^-2
iTdust = 22.0 # in K.
ibeta = 2.0 # 
betaVary = False # If True, Beta will be fitted.
				 # If False, beta will be fixed.
betaMin = 1.0
betaMax = 3.0

wavelengths = np.array([160,250,350, 500, 870]) * u.um
frequencies = wavelengths.to(u.Hz, u.spectral()).value.copy()

errData = np.array([20.0, 20.0, 10.0, 8.0, 8.0])*u.MJy/u.sr
weis = 1.0/errData
# Constants

h = cons.h.cgs.value   # Planck constant in CGS unit
k = cons.k_B.cgs.value # Boltzmann constant in CGS unit
c = cons.c.cgs.value # speed of light in CGS unit
mH = cons.m_n.cgs.value # mass of an neutron
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
gMod.set_param_hint('beta', min = betaMin, max = betaMax, 
	vary = betaVary)
gMod.set_param_hint('Tdust',  max = 80.0)

pars = gMod.make_params()

pfmt = '%6i %18s %6.2f%% finished.'
start = time.time()

for isour in range(400,len(sourList)):
	sName = sourList['sName'][isour]

	if os.path.exists(sName+'_Nh2.fits'):

		hdrRaw   = fits.getheader(sName+'_160_bgRM.fits') 
		hdrNh2   = hdrRaw.copy()
		hdrInt  = hdrRaw.copy()
		hdrTdust = hdrRaw.copy()
		hdrTdust['BUNIT'] = 'K'
		hdrInt['BUNIT'] = 'g s^-3'
		hdrNh2['BUNIT'] = 'cm^-2'
	
		data160 = fits.getdata(sName+'_160_bgRM.fits')
		data250 = fits.getdata(sName+'_250_bgRM.fits')
		data350 = fits.getdata(sName+'_350_bgRM.fits')
		data500 = fits.getdata(sName+'_500_bgRM.fits')
		data870 = fits.getdata(sName+'_870_bgRM.fits')
		dataNh2 = data870.copy()
		dataTdust = data870.copy()
		dataIntegrate = data870.copy()
		
		for iCol in range(len(data870[:,1])):
			for iRow in range(len(data870[1,:])):
				data = np.array([data160[iCol, iRow],
						data250[iCol, iRow],data350[iCol, iRow],data500[iCol, iRow],
						data870[iCol, iRow]])*u.MJy/u.sr
				if ((data160[iCol, iRow] > 60) and
						(data250[iCol, iRow] > 0.0) and
						(data350[iCol, iRow] > 0.0) and
						(data870[iCol, iRow] > 0.0)):
					sedResult = gMod.fit(data.cgs.value, nu = frequencies,
						weights = weis.cgs.value)
					inTSED = quad(greybody,freqLower, freqUpper, 
						args = (sedResult.params['nH2'].value, 
							sedResult.params['Tdust'].value, 
							sedResult.params['beta'].value))
					dataNh2[iCol,iRow] = sedResult.params['nH2'].value
					dataTdust[iCol, iRow] = sedResult.params['Tdust'].value
					dataIntegrate[iCol,iRow] = inTSED[0]
				elif ((data160[iCol, iRow] <= 60) and
						(data250[iCol, iRow] > 0.0) and
						(data350[iCol, iRow] > 0.0) and
						(data870[iCol, iRow] > 0.0)):
					sedResult = gMod.fit(data.cgs.value[1:], nu = frequencies[1:],
						weights = weis.cgs.value[1:])
					inTSED = quad(greybody,freqLower, freqUpper, 
						args = (sedResult.params['nH2'].value, 
							sedResult.params['Tdust'].value, 
							sedResult.params['beta'].value))
					dataNh2[iCol,iRow] = sedResult.params['nH2'].value
					dataTdust[iCol, iRow] = sedResult.params['Tdust'].value
					dataIntegrate[iCol,iRow] = inTSED[0]
				elif ((data160[iCol, iRow] > 60) and
						(data250[iCol, iRow] > 0.0) and
						(data350[iCol, iRow] > 0.0) and
						(data870[iCol, iRow] <= 0.0)):
					sedResult = gMod.fit(data.cgs.value[:-1], nu = frequencies[:-1],
						weights = weis.cgs.value[:-1])
					inTSED = quad(greybody,freqLower, freqUpper, 
						args = (sedResult.params['nH2'].value, 
							sedResult.params['Tdust'].value, 
							sedResult.params['beta'].value))
					dataNh2[iCol,iRow] = sedResult.params['nH2'].value
					dataTdust[iCol, iRow] = sedResult.params['Tdust'].value
					dataIntegrate[iCol,iRow] = inTSED[0]
				else:
					dataNh2[iCol,iRow] = np.nan
					dataTdust[iCol,iRow] = np.nan
					dataIntegrate[iCol,iRow] = np.nan
		
		fits.writeto(sName+'_Nh2.fits', data = dataNh2, 
			header = hdrNh2, clobber = True)
		fits.writeto(sName+'_Tdust.fits', data = dataTdust, 
			header = hdrTdust, clobber = True)
		fits.writeto(sName+'_intSED.fits', data = dataIntegrate, 
			header = hdrInt, clobber = True)
	print(pfmt %(isour+1, sName, ((isour+1.0) / len(sourList)*100)))

stop = time.time()
dure = stop - start

print("Run time = ",dure, "seconds")
