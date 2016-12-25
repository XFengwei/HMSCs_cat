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
from astropy.io import ascii,fits
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

	for bands in ['160', '250', '350', '500']:
		rawHDU = fits.open(sName+'_'+bands+'_smregrid36.fits')
		bgHDU  = fits.open((sName+'_'+bands+'_bg.fits'))

		resiData = rawHDU[0].data - bgHDU[0].data

		resiHDR = rawHDU[0].header.copy()

		fits.writeto(sName+'_'+bands+'_bgRM.fits', 
			data = resiData, header = resiHDR, clobber = True)

	shutil.copy(sName+'_870_smregrid36.fits',
			sName+'_870_bgRM.fits')
	
	print(pfmt %(isour+1, sour_name, ((isour+1.0) / len(sourList)*100)))

stop = time.time()
dure = stop - start

print("Run time = ",dure, "seconds")
