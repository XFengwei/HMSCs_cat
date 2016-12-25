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
from mpl_toolkits.axes_grid1 import make_axes_locatable

mpl.rc("font", family="Complex", size=16)
mpl.rc("axes", linewidth = 1)
mpl.rc("lines", linewidth = 1)
mpl.rc("xtick.major", pad = 8, size = 8, width = 1)
mpl.rc("ytick.major", pad = 8, size = 8, width = 1)
mpl.rc("xtick.minor", size = 4, width = 1)
mpl.rc("ytick.minor", size = 4, width = 1)

#read the data 
hmscList = ascii.read('../Tables/hmscList_full_20161218.txt')
hiiList  = ascii.read('../Tables/hiiList_20161218.dat')
mmbList  = ascii.read('/Users/yuan/Desktop/MaserUTas/MMB/Tables/mmbGP_with_full_para_20161218.txt')
hiiList = hiiList[hiiList['Sp']>0.5]
mmbList = mmbList[mmbList['Sp_870']>0.5]
mmbList = mmbList[mmbList['L_M_ratio']<10000]
hiiList = hiiList[hiiList['L_M_ratio']<10000]
mmbList = mmbList[mmbList['Lclump'].mask == False]
hmscList = hmscList[hmscList['Lclump'].mask == False]
hmscList['L_M_ratio'].mask[hmscList['Mclump'].mask == True] = True

#%matplotlib inline
fig1 = plt.figure(1, figsize = (6.0,4.8))
fig1.subplots_adjust(left = 0.15, right = 0.8,
		bottom = 0.15, top = 0.8)
fig1.clf()
ax1 = fig1.add_subplot(111)
ax1.minorticks_on()

mlReDir = '../Tables/MLrelation/'
for idot in ['01','02','03','04']:
	data = ascii.read(mlReDir+'curve'+idot+'.dat')
	ax1.plot(data['Mass'], data['Lum'], 'k-', zorder = 0)

for idot in ['01','02','03','04','05']:
	data = ascii.read(mlReDir+'hori'+idot+'.dat')
	ax1.plot(data['Mass'], data['Lum'], marker = 'd',
		linestyle=':', color = 'black',
		markersize=2.5, zorder = 0, 
		markerfacecolor = 'None')

for idot in ['01','02','03','04','05']:
	data = ascii.read(mlReDir+'vert'+idot+'.dat')
	ax1.plot(data['Mass'], data['Lum'], marker = 'o',
		linestyle='-', color = 'black',
		markersize=2.5, zorder = 0)

a1,b1 = 1.27, 1.16
a2,b2 = 1.13, 0.95
x = np.logspace(np.log10(2),6,1000)
y1 = 10**(a1*np.log10(x))*10**b1
y2 = 10**(a2*np.log10(x))*10**b2
ax1.plot(x,y1,'k-')
ax1.plot(x,y2,'k--')

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylim(1, 1e7)
ax1.set_xlim(1,1e5)


ax1.scatter(hiiList['Mclump'],hiiList['Lclump'], marker = '.', s = 50, 
            color = 'blue')
ax1.scatter(mmbList['Mclump'],mmbList['Lclump'], marker = '+', s = 12, 
            color = 'lime')
ax1.scatter(hmscList['Mclump'],hmscList['Lclump'], marker = 'x', s = 12, 
            color = 'red')
ax1.set_ylabel('$L_{clump}$ ($L_\odot$)')
ax1.set_xlabel('$M_{clump}$ ($M_\odot$)')


fig1.savefig('../epsFigs/relationML.eps', papertype='a2',
				bbox_inches='tight')
fig1.savefig('../epsFigs/relationML.pdf', papertype='a2',
				bbox_inches='tight')