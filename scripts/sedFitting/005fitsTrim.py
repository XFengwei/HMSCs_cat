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

listFile = '../../Tables/hmscFull_20160405.dat'
sourList = ascii.read(listFile)


sixFITSdir = '../../fitsDir/sedFitting/allSixMin/'

os.chdir(sixFITSdir)

pfmt = '%6i %18s %6.2f%% finished.'
start = time.time()

for isour in range(len(sourList)):
	sour_name = sourList['Name'][isour]
	sName = sourList['sName'][isour]

	for bands in ['160', '250', '350', '500', '870']:
		rawHDU = fits.open(sour_name+'_'+bands+'_smregrid36.fits')
		
		data = rawHDU[0].data[7:-7,7:-7]
		hdr = rawHDU[0].header.copy()
		hdr['NAXIS1'] = rawHDU[0].header['NAXIS1'] - 14
		hdr['NAXIS2'] = rawHDU[0].header['NAXIS2'] - 14
		hdr['CRPIX1'] = rawHDU[0].header['CRPIX1'] - 7
		hdr['CRPIX2'] = rawHDU[0].header['CRPIX2'] - 7

		fits.writeto(sName+'_'+bands+'_smregrid36.fits', 
			data = data, header = hdr, clobber = True)

	print(pfmt %(isour+1, sour_name, ((isour+1.0) / len(sourList)*100)))

stop = time.time()
dure = stop - start

print("Run time = ",dure, "seconds")
