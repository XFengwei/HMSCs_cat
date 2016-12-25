#!/usr/bin python
#

from astropy.io import votable as vo
from astropy.io import ascii
from astropy.io import fits
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

mpl.rc("font", family="Complex", size=17)
mpl.rc("axes", linewidth = 2.0)
mpl.rc("lines", linewidth = 2.0)
mpl.rc("xtick.major", pad = 5, width = 1)
mpl.rc("ytick.major", pad = 5, width = 1)
mpl.rc("xtick.minor", width = 1)
mpl.rc("ytick.minor", width = 1)


listFile = '../Tables/hmscList_full_20161218.txt'
outDir = '../epsFigs/pixSEDfigs/'
sourList = ascii.read(listFile)

atlasgalDir = "../fitsDir/ATLASGAL/"
sedDir = '../fitsDir/sedFitting/sedDir/'

pfmt = '%6i %18s %6.2f%% finished.'

start = time.time()

for isour in range(len(sourList)):
	sour_name = sourList['Name'][isour]
	sName = sourList['sName'][isour]
	amaj = sourList['amaj'][isour]*u.arcsec
	amin = sourList['amin'][isour]*u.arcsec
	paG   =  sourList['PA'][isour]

	ra_0 = sourList['ra'][isour]
	dec_0 = sourList['dec'][isour]
	c1 = coords.SkyCoord(ra_0,dec_0,frame = 'fk5',
		unit = (u.deg,u.deg))
	c2 = coords.SkyCoord(ra_0,dec_0+0.01,frame = 'fk5',
		unit = (u.deg,u.deg))
	c1G = c1.galactic
	ang = c1G.position_angle(c2.galactic).deg
	pa  = paG+90-ang

	fig = plt.figure(figsize = (18,18))

	hdu = fits.open(sedDir+sName+'_Nh2.fits')
	data = hdu[0].data.flatten()
	data = Column(data)
	data = data[data == data]
	data = data[(data > 0) & (data < 1000)]
	n,bins = np.histogram(data, 1000, normed = 1)
	cumu = np.cumsum(n)
	factor = 1.0/cumu.max()
	cumu = factor*cumu
	cdf = Table([bins[:-1],cumu])
	cdf = cdf[cdf['col1']>0.02]
	cdf = cdf[cdf['col1']<0.995]
	maxData = np.max(cdf['col0'])
	minData = np.min(cdf['col0'])
	#logcdf = np.log10(cdf["col0"])
	ticks = [np.int(10*(minData+0.1*(maxData-minData)))/10.0,
		np.int(10*(minData+0.3*(maxData-minData)))/10.0,
		np.int(10*(minData+0.5*(maxData-minData)))/10.0,
		np.int(10*(minData+0.7*(maxData-minData)))/10.0,
		np.int(10*(minData+0.9*(maxData-minData)))/10.0]

	fb3 = ap.FITSFigure(sedDir+sName+'_Nh2.fits',
		figure = fig, subplot = [0.506, 0.1, 0.2, 0.2], 
		aspect="auto")
	fb3.ticks.set_color("black")
	fb3.ticks.set_length(10)
	fb3.show_colorscale(vmin = minData, vmax = maxData,
		cmap = 'viridis')
	conFactor = np.pi*(np.deg2rad(5.061775000E-03))**2/4/np.log(2)
	levs = list(np.array([0.3,0.4,0.5,0.7,0.9,1.3,1.8,2.5,4,5.5,7])*1e-6/conFactor)
	fb3.show_contour(atlasgalDir+sour_name+'_870.fits',
		levels = levs, colors = "red",zorder = 1)
	fb3.show_markers(ra_0,dec_0,marker = "+", s = 300, 
		linewidth = 3, c="white", zorder = 100)
	#fb3.show_ellipses(ra_0,dec_0,amin.to(u.deg).value, amaj.to(u.deg).value,
	#	angle = pa, linewidth = 2, facecolors = "None", edgecolors = 'white', zorder = 100)
	fb3.recenter(ra_0, dec_0, radius = 0.025)
	#fb3.axis_labels.hide()
	#fb3.tick_labels.hide()
	fb3.tick_labels.set_xformat("hhmmss")
	fb3.tick_labels.set_yformat("ddmm")
	fb3.add_colorbar(box = [0.506, 0.303, 0.2,0.01], 
		box_orientation = "horizontal",location = 'top', 
		axis_label_text = 'N$_{H_2}$ ($10^{22}$ cm$^{-2}$)',
		ticks = ticks)

	fb3.add_label(0.07, 0.9, "(a)", relative = True, 
		horizontalalignment='left', color = "white",
		zorder = 1010)
	
	fb3.add_label(-0.32, 1.165,  sour_name, relative = True,
		horizontalalignment='left', color = "black")
	

	#hdu = fits.open(sedDir+sName+'_Tdust.fits')
	#data = hdu[0].data.flatten()
	#data = Column(data)
	#data = data[data == data]
	#data = data[data > 0]
	#n,bins = np.histogram(data, 1000, normed = 1)
	#cumu = np.cumsum(n)
	#factor = 1.0/cumu.max()
	#cumu = factor*cumu
	#cdf = Table([bins[:-1],cumu])
	#cdf = cdf[cdf['col1']>0.2]
	#cdf = cdf[cdf['col1']<0.995]
	#maxData = np.max(cdf['col0'])
	#minData = np.min(cdf['col0'])
	##logcdf = np.log10(cdf["col0"])
	#ticks = [np.int(10*(minData+0.1*(maxData-minData)))/10.0,
	#	np.int(10*(minData+0.3*(maxData-minData)))/10.0,
	#	np.int(10*(minData+0.5*(maxData-minData)))/10.0,
	#	np.int(10*(minData+0.7*(maxData-minData)))/10.0,
	#	np.int(10*(minData+0.9*(maxData-minData)))/10.0]

	fb4 = ap.FITSFigure(sedDir+sName+'_Tdust.fits',
		figure = fig, subplot = [0.709, 0.1, 0.2, 0.2], 
		aspect="auto")
	fb4.ticks.set_color("black")
	fb4.ticks.set_length(10)
	fb4.axis_labels.hide()
	fb4.tick_labels.hide()
	fb4.show_colorscale(vmin = 13, vmax = 27 , 
		cmap = 'inferno')
	fb4.set_nan_color("white")
	fb4.show_markers(ra_0,dec_0,marker = "+", s = 300, 
		linewidth = 3, c="white")
	#fb4.show_circles(ra_0,dec_0,radius = radi, 
	#	linewidth = 2, color = "white")
	fb4.recenter(ra_0, dec_0, radius = 0.025)
	fb4.tick_labels.set_xformat("hhmmss")
	fb4.tick_labels.set_yformat("ddmm")
	fb4.add_colorbar(box = [0.709, 0.303, 0.2,0.01], 
		box_orientation = "horizontal",
		location = 'top', axis_label_text = 'T$_{dust}$ (K)',
		ticks = [14,16,18,20,22,24,26])
	#fb3.add_colorbar(location='top', 
	#	box_orientation = "horizontal", ticks = ticks,
	#	axis_label_text = "Myr sr$^{-1}$")
	fb4.add_label(0.07, 0.9, "(b)", relative = True, 
		horizontalalignment='left', color = "white")
	#fb4.add_beam( facecolor = 'cyan',
	#	corner = 'bottom right',
	#	linewidth = 0.5)
	#fb3.add_label(0.95, 0.9,  "70 $\mu$m", relative = True,
	#	horizontalalignment='right', color = "white")


	plt.savefig(outDir+sour_name+'_sedImshow.eps', bbox_inches='tight',
		papertype='a2')
	plt.savefig(outDir+sour_name+'_sedImshow.pdf', bbox_inches='tight',
		papertype='a2')
	plt.savefig(outDir+sour_name+'_sedImshow.png', bbox_inches='tight', dpi = 150,
		papertype='a2')  
	#plt.savefig(outDir+sour_name+'_mulPlot.png', bbox_inches='tight', dpi = 300,
	#	papertype='a2')
	fig.clf()

	print(pfmt %(isour+1, sour_name, ((isour+1.0) / len(sourList)*100)))

ending = time.time()
print((ending-start)/60.0, 'minutes have been used.')