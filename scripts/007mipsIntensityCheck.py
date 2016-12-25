#!/usr/bin python
#

from astropy.io import votable as vo
from astropy.io import ascii
from astropy.io import fits
from astropy import wcs
import matplotlib.pyplot as plt
import matplotlib as mpl
import aplpy as ap
import numpy as np
import scipy as sp
import time
import os
from astropy import units as u
import astropy.coordinates as coords
from astropy.table import Table, Column, vstack, join

mpl.rc("font", family="serif", size=19)
mpl.rc("axes", linewidth = 2.0)
mpl.rc("lines", linewidth = 2.0)
mpl.rc("xtick.major", pad = 5, width = 1)
mpl.rc("ytick.major", pad = 5, width = 1)
mpl.rc("xtick.minor", width = 1)
mpl.rc("ytick.minor", width = 1)

listFile = '../Tables/StarlessListBefore24ext.dat'
mipsDir = '../fitsDir/MIPS/'

sourList = ascii.read(listFile, format = 'fixed_width')

pfmt = '%6i %18s %6.2f%% finished.'

start = time.time()

colList = list(sourList.columns)
sourList.remove_columns([colList[0], colList[-1]])
refinedStarlessfList = sourList[sourList['ra'] < 0]
abortedList = sourList[sourList['ra'] < 0]
for isour in range(len(sourList)):
	sour_name = sourList['Name'][isour]
	ra_0 = sourList['ra'][isour]
	dec_0 = sourList['dec'][isour]
	radi = sourList['amaj'][isour] / 3600.0 / 2.0
	
	hdus = fits.open(mipsDir+sour_name+'_mips.fits')
	image = hdus[0].data
	hd = hdus[0].header
	w = wcs.WCS(hd)
	for i in range(image.shape[0]):
		for j in range(image.shape[1]):
			if (image[i][j] != image[i][j]):
				image[i][j] = 100000.0
	centerWorld = np.array([[ra_0, dec_0]])
	centerPix = w.wcs_world2pix(centerWorld, 1)
	radiInPix = radi / hd['CDELT2']
	y,x = np.indices(image.shape)
	r = np.hypot(x - centerPix[0][0], y - centerPix[0][1])
	ind = np.argsort(r.flat)
	r_sorted = r.flat[ind]
	i_sorted = image.flat[ind]
	subInd = np.where(r_sorted < radiInPix)
	subData = i_sorted[subInd[0]]
	lenVeryBright = len(subData[np.where(subData > 1000.0)[0]])
	if lenVeryBright < 10:
		refinedStarlessfList.add_row(sourList[isour])
	else:
		abortedList.add_row(sourList[isour])

	print  pfmt %(isour+1, sour_name, ((isour+1.0) / len(sourList)*100))

print len(refinedStarlessfList)
refinedStarlessfList.write("../Tables/StarlessAfter24ext.dat", 
	format = "ascii.fixed_width")
abortedList.write("../Tables/abortedList.dat", 
	format = "ascii.fixed_width")

	
