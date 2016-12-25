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

mpl.rc("font", family="Complex", size=16)
mpl.rc("axes", linewidth =  1 )
mpl.rc("lines", linewidth = 1 )
mpl.rc("xtick.major", pad = 8, size = 8, width = 1)
mpl.rc("ytick.major", pad = 8, size = 8, width = 1)
mpl.rc("xtick.minor", size = 4, width = 1 )
mpl.rc("ytick.minor", size = 4, width = 1 )
hmscList = ascii.read('../Tables/hmscList_full_20161218.txt')
hiiList  = ascii.read('../Tables/hiiList_20161218.dat')
mmbList  = ascii.read('/Users/yuan/Desktop/MaserUTas'+
    '/MMB/Tables/mmbGP_with_full_para_20161218.txt')
hiiList = hiiList[hiiList['Sp']>0.5]
mmbList = mmbList[mmbList['Sp_870']>0.5]
mmbList = mmbList[mmbList['L_M_ratio']<10000]
hiiList = hiiList[hiiList['L_M_ratio']<10000]
hmscList['L_M_ratio'].mask[hmscList['Mclump'].mask == True] = True

xmin = 0.04
xmax = 5
ymin = 20 
ymax = 5e4
LumDense = hmscList['Lclump']/(np.pi*hmscList['r_pc']**2)
fig8 = plt.figure(2, figsize = (4.5,9.5))
fig8.subplots_adjust(left = 0.07, right = 0.91, hspace=0.02,
                     bottom = 0.09, top = 0.96)

ax1 = fig8.add_subplot(311)
ax1.scatter(hmscList['r_pc'], LumDense,marker = 'x', s = 22, 
            color = 'red')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(xmin,xmax)
ax1.set_ylim([ymin,ymax])
ax1.minorticks_on()
ax1.set_xticklabels([])
#ax1.set_yticks([0,10,20,30,40,50])

#ax1.set_ylabel(r'$L_{clump}/4\pi r_{eq}^{2}$ ($L_\odot$ pc$^{-2}$)')
ax1.text(10**(np.log10(xmax)-(np.log10(xmax)-np.log10(xmin))*0.08), 
         10**(np.log10(ymax)-(np.log10(ymax)-np.log10(ymin))*0.12), 
         'starless', horizontalalignment = 'right')

ymin = 150 
ymax = 2e6
LumDense = mmbList['Lclump']/(np.pi*mmbList['r_pc']**2)
ax2 = fig8.add_subplot(312)
ax2.scatter(mmbList['r_pc'], LumDense, marker = '+', s = 40, 
            color = 'lime')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim(xmin,xmax)
ax2.set_ylim([ymin,ymax])
ax2.minorticks_on()
ax2.set_xticklabels([])
#ax1.set_yticks([0,10,20,30,40,50])

ax2.set_ylabel(r'Luminosity Surface Density ($L_\odot$ pc$^{-2}$)')
ax2.text(10**(np.log10(xmax)-(np.log10(xmax)-np.log10(xmin))*0.08), 
         10**(np.log10(ymax)-(np.log10(ymax)-np.log10(ymin))*0.12), 
         'CH$_3$OH', horizontalalignment = 'right')

ymin = 150 
ymax = 2e6
LumDense = hiiList['Lclump']/(np.pi*hiiList['r_pc']**2)
ax3 = fig8.add_subplot(313)
ax3.scatter(hiiList['r_pc'], LumDense,marker = '.', s = 80, 
            color = 'blue')
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xlim(xmin,xmax)
ax3.set_ylim([ymin,ymax])
ax3.minorticks_on()
#ax1.set_yticks([0,10,20,30,40,50])
ax3.set_xlabel(r"r$_{eq}$ (pc)")
#ax3.set_ylabel(r'Luminosity Surface Density ($L_\odot$ pc$^{-2}$)')
ax3.text(10**(np.log10(xmax)-(np.log10(xmax)-np.log10(xmin))*0.08), 
         10**(np.log10(ymax)-(np.log10(ymax)-np.log10(ymin))*0.12), 
         'HII', horizontalalignment = 'right')

fig8.savefig('../epsFigs/LumSurfaceDense.eps' ,dpi = 300, 
             bbox_inches='tight', papertype='a2')
fig8.savefig('../epsFigs/LumSurfaceDense.pdf' ,dpi = 300, 
             bbox_inches='tight', papertype='a2')


