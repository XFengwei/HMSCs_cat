from scipy import fftpack
import numpy as np
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
from astropy.table import Table, Column, vstack
from astropy import units as u
import astropy.coordinates as coords
import sys
import glob
from astropy import units as u
from astropy.io import fits
from higal_sedfitter import smooth,fit,higal_beams
from higal_sedfitter.fit import PixelFitter

listFile = '../../Tables/hmscFull_20160405.dat'
sourList = ascii.read(listFile)

sixFITSdir = '../../fitsDir/sedFitting/allSixMin/'
# arguments
#imf = str(sys.argv[1]) # input image file
#stf = str(sys.argv[2]) # output structure file
#bgf = str(sys.argv[3]) # output background file
#frac  = float(sys.argv[4]) # fraction of max power spectrum, to separate hight and low spatial freq
#scale = float(sys.argv[5]) # factor to scale background, usually 1.0
os.chdir(sixFITSdir)
pfmt = '%6i %18s %6.2f%% finished.'
start = time.time()

pixelfitter = PixelFitter(bfixed=True)
for isour in range(len(sourList)):
	sour_name = sourList['Name'][isour]
	higal_field = sour_name
	fmt = {'field':higal_field}

	target_fn = sour_name+'_500.fits'
	target_header = fits.getheader(target_fn)
	#target_header['BMAJ'] = 37.0/3600
	#target_header['BMIN'] = 37.0/3600
	fitsList = [sour_name+'_160.fits', sour_name+'_250.fits',
				sour_name+'_350.fits', sour_name+'_870.fits']
	smooth.smooth_images_toresolution(36.4*u.arcsec, skip_existing=False,
		globs=fitsList,
		target_header=target_header, clobber=True)

	shutil.copy(sour_name+'_500.fits',
			sour_name+'_500_smregrid36.fits')
	bands = ['160', '250', '350', '870']
	for band in bands:
		os.remove(sour_name+'_'+band+'_smooth36.fits')

	print(pfmt %(isour+1, sour_name, ((isour+1.0) / len(sourList)*100)))

stop = time.time()
dure = stop - start

print("Run time = ",dure, "seconds")