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

mpl.rc("font", family="Complex", size=18)
mpl.rc("axes", linewidth = 2.0)
mpl.rc("lines", linewidth = 2.0)
mpl.rc("xtick.major", pad = 5, width = 1)
mpl.rc("ytick.major", pad = 5, width = 1)
mpl.rc("xtick.minor", width = 1)
mpl.rc("ytick.minor", width = 1)

listFile = '../Tables/hmscList_201612.txt'
outDir = '../epsFigs/multiPlot/'
sourList = ascii.read(listFile)

atlasgalDir = "../fitsDir/ATLASGAL/"
higalDir = '../fitsDir/HiGAL/'
iracDir =  '../fitsDir/IRAC/'
mipsDir =  '../fitsDir/MIPS/'
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
	coor = coords.SkyCoord(ra_0, dec_0, frame = 'icrs',
	 	unit=(u.deg, u.deg))
	ra_0 = coor.fk5.ra.deg
	dec_0 = coor.fk5.dec.deg
	c1 = coords.SkyCoord(ra_0,dec_0,frame = 'fk5',
		unit = (u.deg,u.deg))
	c2 = coords.SkyCoord(ra_0,dec_0+0.01,frame = 'fk5',
		unit = (u.deg,u.deg))
	c1G = c1.galactic
	ang = c1G.position_angle(c2.galactic).deg
	pa  = paG+90-ang
	#radi = sourList['amaj'][isour] / 3600.0 / 2

	ap.make_rgb_cube([iracDir+sour_name+'/'+sour_name+'_I4.fits',
		iracDir+sour_name+'/'+sour_name+'_I2.fits',
		iracDir+sour_name+'/'+sour_name+'_I1.fits'],
		iracDir+sour_name+'/'+sour_name+'Spitzer_cube.fits', 
		system='EQUJ', north='Ture')

	ap.make_rgb_image(iracDir+sour_name+'/'+sour_name+'Spitzer_cube.fits',
		iracDir+sour_name+'/'+sour_name+'_Spitzer_rgb.png', 
		pmin_r=0.5, pmax_r=99.6, pmin_g=0.5, pmax_g=99.6, 
		pmin_b=0.5, pmax_b=99.6, stretch_r='log', 
		stretch_g='log', stretch_b='log', embed_avm_tags=False)

	fig = plt.figure(figsize = (18,18))

	fb1 = ap.FITSFigure(iracDir+sour_name+'/'+sour_name+'_I4.fits',
		figure = fig, subplot = [0.1, 0.1, 0.2, 0.2], 
		aspect="auto")
	fb1.ticks.set_length(10)
	fb1.show_rgb(iracDir+sour_name+'/'+sour_name+'_Spitzer_rgb.png')
	#fb1.show_markers(ra_0,dec_0,marker = "+", s = 300, 
	#	linewidth = 3, c="white")
	#fb1.show_circles(ra_0,dec_0,radius = radi, 
	#	linewidth = 2, color = "blue")
	fb1.show_ellipses(ra_0,dec_0,amin.to(u.deg).value, amaj.to(u.deg).value,
		angle = pa, linewidth = 2, facecolors = "None", edgecolors = 'white', zorder = 100)
	fb1.recenter(ra_0, dec_0, radius = 0.025)
	fb1.tick_labels.set_xformat("hhmmss")
	fb1.tick_labels.set_yformat("ddmm")

	fb1.add_label(0.93, 0.92, 'red: 8.0 $\mu$m', relative=True, 
		horizontalalignment='right', color="white")
	fb1.add_label(0.93, 0.85, 'green: 4.5 $\mu$m', relative=True, 
		horizontalalignment='right', color="white")
	fb1.add_label(0.93, 0.78, 'blue: 3.6 $\mu$m', relative=True, 
		horizontalalignment='right', color="white")
	fb1.add_label(0.07, 0.9, "(a)", relative = True,
		horizontalalignment='left', color = "white")
	fb1.add_label(0.0, 1.05,  sour_name, relative = True,
		horizontalalignment='left', color = "black")

	hdu = fits.open(mipsDir+sour_name+'_mips.fits')
	data = hdu[0].data.flatten()
	data = Column(data)
	data = data[data == data]
	data = data[data > 0]
	n,bins = np.histogram(data, 1000, normed = 1)
	cumu = np.cumsum(n)
	factor = 1.0/cumu.max()
	cumu = factor*cumu
	cdf = Table([bins[:-1],cumu])
	cdf = cdf[cdf['col1']>0.2]
	cdf = cdf[cdf['col1']<0.995]
	maxData = np.max(cdf['col0'])
	minData = np.min(cdf['col0'])
	logcdf = np.log10(cdf["col0"])
	ticks = [np.int(10**(logcdf.min()+0.2*(logcdf.max()-logcdf.min()))),
		np.int(10**(logcdf.min()+0.5*(logcdf.max()-logcdf.min()))),
		np.int(10**(logcdf.min()+0.8*(logcdf.max()-logcdf.min())))]
	fb2 = ap.FITSFigure(mipsDir+sour_name+'_mips.fits',
		figure = fig, subplot = [0.303, 0.1, 0.2, 0.2], 
		aspect="auto")
	#fb2.ticks.set_color("black")
	fb2.ticks.set_length(10)
	fb2.axis_labels.hide()
	fb2.tick_labels.hide()
	fb2.tick_labels.set_xformat("hhmmss")
	fb2.tick_labels.set_yformat("ddmm")
	fb2.show_colorscale(vmin = minData, vmax = maxData,
		cmap = 'viridis', stretch = 'log')
	fb2.set_nan_color("khaki")
	#fb2.show_markers(ra_0,dec_0,marker = "+", s = 300, 
	#	linewidth = 3, c="white")
	fb2.show_ellipses(ra_0,dec_0,amin.to(u.deg).value, amaj.to(u.deg).value,
		angle = pa, linewidth = 2, facecolors = "None", edgecolors = 'white', zorder = 100)
	fb2.recenter(ra_0, dec_0, radius = 0.025)
	fb2.tick_labels.set_xformat("hhmmss")
	fb2.tick_labels.set_yformat("ddmm")
	fb2.add_colorbar(box = [0.303, 0.303, 0.2,0.01], 
		box_orientation = "horizontal",
		location = 'top', ticks = ticks,
		axis_label_text = "Myr sr$^{-1}$")
	fb2.add_label(0.07, 0.9, "(b)", relative = True,
		horizontalalignment='left', color = "white")

	fb2.add_label(0.93, 0.90, '24 $\mu$m', relative=True, 
		horizontalalignment='right', color="white")


	hdu = fits.open(higalDir+sour_name+'/'+sour_name+'_070.fits')
	data = hdu[0].data.flatten()
	data = Column(data)
	data = data[data == data]
	data = data[data > 0]
	n,bins = np.histogram(data, 1000, normed = 1)
	cumu = np.cumsum(n)
	factor = 1.0/cumu.max()
	cumu = factor*cumu
	cdf = Table([bins[:-1],cumu])
	cdf = cdf[cdf['col1']>0.2]
	cdf = cdf[cdf['col1']<0.995]
	maxData = np.max(cdf['col0'])
	minData = np.min(cdf['col0'])
	logcdf = np.log10(cdf["col0"])
	ticks = [np.int(10**(logcdf.min()+0.2*(logcdf.max()-logcdf.min()))),
		np.int(10**(logcdf.min()+0.5*(logcdf.max()-logcdf.min()))),
		np.int(10**(logcdf.min()+0.8*(logcdf.max()-logcdf.min())))]

	fb3 = ap.FITSFigure(higalDir+sour_name+'/'+sour_name+'_070.fits',
		figure = fig, subplot = [0.506, 0.1, 0.2, 0.2], 
		aspect="auto")
	fb3.ticks.set_length(10)
	fb3.show_colorscale(vmin = minData, vmax = maxData,
		cmap = 'viridis', stretch = 'log')
	#conFactor = np.pi*(np.deg2rad(5.061775000E-03))**2/4/np.log(2)
	#levs = list(np.array([0.3,0.4,0.5,0.7,0.9,1.3,1.8,2.5,4,5.5,7])*1e-6/conFactor)
	#fb3.show_contour(atlasgalDir+sour_name+'_870.fits',
	#	levels = levs, colors = "red",zorder = 1)
	#fb3.show_markers(ra_0,dec_0,marker = "+", s = 300, 
	#	linewidth = 3, c="white")
	fb3.show_ellipses(ra_0,dec_0,amin.to(u.deg).value, amaj.to(u.deg).value,
		angle = pa, linewidth = 2, facecolors = "None", edgecolors = 'white', zorder = 100)
	fb3.recenter(ra_0, dec_0, radius = 0.025)
	fb3.axis_labels.hide()
	fb3.tick_labels.hide()
	fb3.tick_labels.set_xformat("hhmmss")
	fb3.tick_labels.set_yformat("ddmm")
	fb3.add_colorbar(box = [0.506, 0.303, 0.2,0.01], 
		box_orientation = "horizontal",location = 'top', 
		axis_label_text = 'Myr sr$^{-1}$',
		ticks = ticks)

	fb3.add_label(0.07, 0.9, "(c)", relative = True,
		horizontalalignment='left', color = "white",
		zorder = 1010)
	
	fb3.add_label(0.93, 0.90, '70 $\mu$m', relative=True, 
		horizontalalignment='right', color="white")

	hdu = fits.open(higalDir+sour_name+'/'+sour_name+'_250.fits')
	data = hdu[0].data.flatten()
	data = Column(data)
	data = data[data == data]
	data = data[data > 0]
	n,bins = np.histogram(data, 1000, normed = 1)
	cumu = np.cumsum(n)
	factor = 1.0/cumu.max()
	cumu = factor*cumu
	cdf = Table([bins[:-1],cumu])
	cdf = cdf[cdf['col1']>0.2]
	cdf = cdf[cdf['col1']<0.995]
	maxData = np.max(cdf['col0'])
	minData = np.min(cdf['col0'])
	logcdf = np.log10(cdf["col0"])
	ticks = [np.int(10**(logcdf.min()+0.2*(logcdf.max()-logcdf.min()))),
		np.int(10**(logcdf.min()+0.5*(logcdf.max()-logcdf.min()))),
		np.int(10**(logcdf.min()+0.8*(logcdf.max()-logcdf.min())))]

	fb4 = ap.FITSFigure(higalDir+sour_name+'/'+sour_name+'_250.fits',
		figure = fig, subplot = [0.709, 0.1, 0.2, 0.2], 
		aspect="auto")
	fb4.ticks.set_length(10)
	fb4.show_colorscale(vmin = minData, vmax = maxData,
		cmap = 'viridis', stretch = 'log')
	conFactor = np.pi*(np.deg2rad(5.061775000E-03))**2/4/np.log(2)
	levs = list(np.array([0.3,0.4,0.5,0.7,0.9,1.3,1.8,2.5,4,5.5,7,9,11.5,14])*1e-6/conFactor)
	fb4.show_contour(atlasgalDir+sour_name+'_870.fits',
		levels = levs, colors = "red",zorder = 1)
	fb4.show_markers(ra_0,dec_0,marker = "+", s = 300, 
		linewidth = 3, c="white")
	fb4.show_ellipses(ra_0,dec_0,amin.to(u.deg).value, amaj.to(u.deg).value,
		angle = pa, linewidth = 2, facecolors = "None", edgecolors = 'white', zorder = 100)
	fb4.recenter(ra_0, dec_0, radius = 0.025)
	fb4.axis_labels.hide()
	fb4.tick_labels.hide()
	fb4.tick_labels.set_xformat("hhmmss")
	fb4.tick_labels.set_yformat("ddmm")
	fb4.add_colorbar(box = [0.709, 0.303, 0.2,0.01], 
		box_orientation = "horizontal",location = 'top', 
		axis_label_text = 'Myr sr$^{-1}$',
		ticks = ticks)
	fb4.add_label(0.07, 0.9, "(d)", relative = True,
		horizontalalignment='left', color = "white")
	fb4.add_label(0.93, 0.90, '250 $\mu$m', relative=True, 
		horizontalalignment='right', color="white")
	
	
	plt.savefig(outDir+sour_name+'_mulPlot.eps', bbox_inches='tight',
		papertype='a2')
	plt.savefig(outDir+sour_name+'_mulPlot.pdf', bbox_inches='tight',
		papertype='a2')
	plt.savefig(outDir+sour_name+'_mulPlot.png', bbox_inches='tight', dpi = 150,
		papertype='a2')
	#plt.savefig(outDir+sour_name+'_mulPlot.png', bbox_inches='tight', dpi = 300,
	#	papertype='a2')
	fig.clf()

	print(pfmt %(isour+1, sour_name, ((isour+1.0) / len(sourList)*100)))

    
ending = time.time()
print((ending-start)/60.0, 'minutes have been used.')