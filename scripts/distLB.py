import os
import math
import time
import tarfile
import shutil
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy import constants as cons
from astropy.table import Table, Column, vstack, join
import astropy.coordinates as coords
from astropy import units as u
from scipy import stats
import matplotlib.mlab as mlab
from matplotlib.ticker import MultipleLocator

mpl.rc("font", family="Complex", size=15)
mpl.rc("axes", linewidth = 1.0)
mpl.rc("lines", linewidth = 1.0)
mpl.rc("xtick.major", pad = 8, size = 12, width = 1)
mpl.rc("ytick.major", pad = 8, size = 12, width = 1)
mpl.rc("xtick.minor", size = 6, width = 1)
mpl.rc("ytick.minor", size = 6, width = 1)

sourListFile = '../Tables/hmscList_201612.txt'
sourList = ascii.read(sourListFile)

coord = coords.SkyCoord(sourList['ra'].data, 
	sourList['dec'].data, frame='fk5', 
	unit = (u.deg, u.deg))

l = Column(coord.galactic.l.deg, name = 'l')
b = Column(coord.galactic.b.deg, name = 'b')

for i in range(len(l)):
	if l[i]>180:
		l[i] -= 360 


# Distributions of excitation temperatures and column densities
fig1 = plt.figure(1, figsize = (10,4.5))
fig1.subplots_adjust(left = 0.07, right = 0.91, wspace = 0.38,
		bottom = 0.07, top = 0.91)
ax1 = fig1.add_subplot(111)
#ax1.set_yticks([0.05, 0.1, 0.15, 0.2])
ax1.set_ylim(0,105)
ax1.set_yticks([ 20, 40,60,80,100])
ax1.set_xticks([ -60, -40, -20, 0, 20, 40,60])
ax1.set_xticklabels([r'300${\degree}$', r'320${\degree}$', 
	'340${\degree}$', '0${\degree}$', '20${\degree}$', 
	'40${\degree}$', '60${\degree}$'])
xmin = -65
xmax = 65
ax1.set_xlim(xmax,xmin)
ax1.hist(l.data, 
	60, edgecolor='black', facecolor = 'yellow',
	alpha = 1, range = (-60,60),  linewidth = 1)

ax1.minorticks_on()
ax1.xaxis.set_minor_locator(MultipleLocator(5))
ax1.set_xlabel(r"Galactic Longitude (deg)")
ax1.set_ylabel('Number')

compCoords = [49.4857, 30.7593, 23.0, 15.041, 11.3, -7.8, -23.    , -27., -33., -55.]
compNames = ['W51', 'W43', 'G23', 'M17', 'W31/W33', 'NGC6357/6334', 'G337', 'G333', 'G327',  'G305']

for icom in range(len(compCoords)):
	comL = compCoords[icom]
	comN = compNames[icom]
	ax1.plot([comL, comL],[45, 60],'--r')
	ax1.text(comL, 62, comN, rotation=90,
		verticalalignment='bottom', horizontalalignment='center')

ax1.text(0, 85, 'GC', rotation=90,
	verticalalignment='bottom', horizontalalignment='center')

#ymin,ymax = plt.ylim()
#ax1.text(xmin+(xmax-xmin)*0.07, 
#		ymax-(ymax-ymin)*0.12, '(a)')

fig1.savefig('../epsFigs/distL.eps' ,dpi = 300, 
	bbox_inches='tight', papertype='a2')
fig1.savefig('../epsFigs/distL.pdf' ,dpi = 300, 
	bbox_inches='tight', papertype='a2')
fig1.clear()



# Distributions of excitation temperatures and column densities
fig2 = plt.figure(2, figsize = (6,4.2))
fig2.subplots_adjust(left = 0.07, right = 0.91, wspace = 0.38,
		bottom = 0.07, top = 0.91)
ax1 = fig2.add_subplot(111)
#ax1.set_yticks([0.05, 0.1, 0.15, 0.2])
ax1.set_ylim(0,85)
ax1.set_yticks([ 20, 40,60,80])
ax1.set_xticks([ -1,-0.5,0,0.5,1])
ax1.set_xticklabels(['-1.0$^\circ$', '-0.5$^\circ$', 
	'0.0$^\circ$', '0.5$^\circ$', '1.0$^\circ$'])
xmin = -1.4
xmax = 1.4
ax1.set_xlim(xmin,xmax)
ax1.hist(b.data, 
	20, edgecolor='black', facecolor = 'yellow',
	alpha = 1, range = (-1,1),  linewidth = 2)

ax1.minorticks_on()
#ax1.xaxis.set_minor_locator(MultipleLocator(2))
ax1.set_xlabel(r"Galactic Latitude (deg)")
ax1.set_ylabel('Number')
#ymin,ymax = plt.ylim()
#ax1.text(xmin+(xmax-xmin)*0.07, 
#		ymax-(ymax-ymin)*0.12, '(a)')

fig2.savefig('../epsFigs/distB.eps' ,dpi = 300, 
	bbox_inches='tight', papertype='a2')


