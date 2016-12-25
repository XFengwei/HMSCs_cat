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


mpl.rc("font", family="Complex", size=15)
mpl.rc("axes", linewidth = 1.0)
mpl.rc("lines", linewidth = 1.0)
mpl.rc("xtick.major", pad = 5, width = 1)
mpl.rc("ytick.major", pad = 5, width = 1)
mpl.rc("xtick.minor", width = 1)
mpl.rc("ytick.minor", width = 1)


listFile = '../Tables/hmscList_201612.txt'
sourList = ascii.read(listFile)

ra = sourList['ra'].data
dec = sourList['dec'].data

coor = coords.SkyCoord(ra, dec, frame = 'icrs',
        unit=(u.degree, u.degree))

l = Column(coor.galactic.l.deg, name = 'l')
b = Column(coor.galactic.b.deg, name = 'b')

PlanckCO = '../fitsDir/PlanckCO_Type3_car.fits'
#PlanckCO = '/Users/yuan/work/Python_Demo/PlanckCO/PlanckCO_Type3_car.fits'

fig = plt.figure(figsize = (14,14))


height = 2.3529411764705883
fb1 = ap.FITSFigure(PlanckCO,
		figure = fig, subplot = [0.13, 0.85, 0.85, 0.1], 
		aspect="auto")
fb1.recenter(50,0,width = 20, height = height)
fb1.show_colorscale(cmap = 'viridis', vmin = 10, vmax = 1000, stretch = 'log')
fb1.tick_labels.set_xformat("ddd")
fb1.tick_labels.set_yformat("dd.dd")
fb1.ticks.set_yspacing(0.5)

fb1.ticks.set_xspacing(5)
fb1.set_axis_labels('l','b')
fb1.hide_xaxis_label()
fb1.show_markers(l,b)
fb1.show_regions('../Tables/filaments_boxed_region.reg')
fb1.add_colorbar(box = [0.5, 0.955, 0.48,0.01], 
		location='top', ticks=[20,50,100,200,500,1000],
		box_orientation = "horizontal", 
		axis_label_text = "K km s$^{-1}$")

fb2 = ap.FITSFigure(PlanckCO,
	figure = fig, subplot = [0.13, 0.7, 0.85, 0.1], 
	aspect="auto")
fb2.recenter(30,0,width = 20, height = height)
fb2.show_colorscale(cmap = 'viridis', vmin = 10, vmax = 1000, stretch = 'log')
fb2.tick_labels.set_xformat("ddd")
fb2.tick_labels.set_yformat("dd.dd")
fb2.ticks.set_yspacing(0.5)
fb2.ticks.set_xspacing(5)
fb2.set_axis_labels('l','b')
fb2.hide_xaxis_label()
fb2.show_markers(l,b)
fb2.show_regions('../Tables/filaments_boxed_region.reg')

fb3 = ap.FITSFigure(PlanckCO,
	figure = fig, subplot = [0.13, 0.55, 0.85, 0.1], 
	aspect="auto")
fb3.recenter(10,0,width = 20, height = height)
fb3.show_colorscale(cmap = 'viridis', vmin = 10, vmax = 1000, stretch = 'log')
fb3.tick_labels.set_xformat("ddd")
fb3.tick_labels.set_yformat("dd.dd")
fb3.ticks.set_yspacing(0.5)
fb3.ticks.set_xspacing(5)
fb3.set_axis_labels('l','b')
fb3.hide_xaxis_label()
fb3.show_markers(l,b)
fb3.show_regions('../Tables/filaments_boxed_region.reg')


#fb4 = ap.FITSFigure(PlanckCO,
#	figure = fig, subplot = [0.13, 0.4, 0.85, 0.1], 
#	aspect="auto")
#fb4.recenter(-10,0,width = 20, height = height)
#fb4.show_colorscale(cmap = 'viridis', vmin = 10, vmax = 1000, stretch = 'log')
#fb4.tick_labels.set_xformat("ddd")
#fb4.tick_labels.set_yformat("dd.dd")
#fb4.ticks.set_yspacing(0.5)
#fb4.ticks.set_xspacing(5)
#fb4.set_axis_labels('l','b')
#fb4.hide_xaxis_label()
#fb4.show_markers(l,b)
#fb4.show_regions('../Tables/filaments_boxed_region.reg')

fb4 = ap.FITSFigure(PlanckCO,
	figure = fig, subplot = [0.13, 0.4, 0.85, 0.1], 
	aspect="auto")
fb4.recenter(-10,0,width = 20, height = height)
fb4.show_colorscale(cmap = 'viridis', vmin = 10, vmax = 1000, stretch = 'log')
fb4.tick_labels.set_xformat("ddd")
fb4.tick_labels.set_yformat("dd.dd")
fb4.ticks.set_yspacing(0.5)
fb4.ticks.set_xspacing(5)
fb4.set_axis_labels('l','b')
fb4.hide_xaxis_label()
fb4.show_markers(l,b)
fb4.show_regions('../Tables/filaments_boxed_region.reg')

fb5 = ap.FITSFigure(PlanckCO,
	figure = fig, subplot = [0.13, 0.25, 0.85, 0.1], 
	aspect="auto")
fb5.recenter(-30,0,width = 20, height = height)
fb5.show_colorscale(cmap = 'viridis', vmin = 10, vmax = 1000, stretch = 'log')
fb5.tick_labels.set_xformat("ddd")
fb5.tick_labels.set_yformat("dd.dd")
fb5.ticks.set_yspacing(0.5)
fb5.ticks.set_xspacing(5)
fb5.set_axis_labels('l','b')
fb5.hide_xaxis_label()
fb5.show_markers(l,b)
fb5.show_regions('../Tables/filaments_boxed_region.reg')

fb6 = ap.FITSFigure(PlanckCO,
	figure = fig, subplot = [0.13, 0.1, 0.85, 0.1], 
	aspect="auto")
fb6.recenter(-50,0,width = 20, height = height)
fb6.show_colorscale(cmap = 'viridis', vmin = 10, vmax = 1000, stretch = 'log')
fb6.tick_labels.set_xformat("ddd")
fb6.tick_labels.set_yformat("dd.dd")
fb6.ticks.set_yspacing(0.5)
fb6.ticks.set_xspacing(5)
fb6.set_axis_labels('l','b')
fb6.show_markers(l,b)
fb6.show_regions('../Tables/filaments_boxed_region.reg')

plt.savefig('../epsFigs/starLessDistribution.eps', bbox_inches='tight',
		papertype='a2')
plt.savefig('../epsFigs/starLessDistribution.pdf', bbox_inches='tight',
		papertype='a2')