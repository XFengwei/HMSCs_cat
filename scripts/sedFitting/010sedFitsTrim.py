#!/usr/bin python

import os
import math
import time
import tarfile
import shutil
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import montage_wrapper as mt
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table, Column, vstack
from astropy import units as u
import astropy.coordinates as coords

listFile = '../../Tables/hmscList_full_20161218.txt'
sourList = ascii.read(listFile)


sixFITSdir = '../../fitsDir/sedFitting/allSixMin/'


os.chdir(sixFITSdir)

for isour in range(len(sourList)):
	sour_name = sourList['Name'][isour]
	sName = sourList['sName'][isour]

	rawHDU = fits.open(sName+'_Tdust.fits')
	
	data = rawHDU[0].data[7:-7,7:-7]
	data[(data>100) | (data<5)] = np.nan
	hdr = rawHDU[0].header.copy()
	hdr['NAXIS1'] = rawHDU[0].header['NAXIS1'] - 14
	hdr['NAXIS2'] = rawHDU[0].header['NAXIS2'] - 14
	hdr['CRPIX1'] = rawHDU[0].header['CRPIX1'] - 7
	hdr['CRPIX2'] = rawHDU[0].header['CRPIX2'] - 7
	fits.writeto('../sedDir/'+sName+'_Tdust.fits', 
		data = data, header = hdr, clobber = True)

for isour in range(len(sourList)):
	sour_name = sourList['Name'][isour]
	sName = sourList['sName'][isour]

	rawHDU = fits.open(sName+'_Nh2.fits')
	
	data = rawHDU[0].data[7:-7,7:-7]/1.0e22
	data[(data>1e4) | (data<0)] = np.nan
	hdr = rawHDU[0].header.copy()
	hdr['NAXIS1'] = rawHDU[0].header['NAXIS1'] - 14
	hdr['NAXIS2'] = rawHDU[0].header['NAXIS2'] - 14
	hdr['CRPIX1'] = rawHDU[0].header['CRPIX1'] - 7
	hdr['CRPIX2'] = rawHDU[0].header['CRPIX2'] - 7
	fits.writeto('../sedDir/'+sName+'_Nh2.fits', 
		data = data, header = hdr, clobber = True)
