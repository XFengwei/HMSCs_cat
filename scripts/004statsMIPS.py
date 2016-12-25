#!/usr/bin python
#
# This script is used to identity and classify YSOs based on 
# the scheme of Koenig et al. (2012), ApJ, 744, 130
# 
# Originally written by Jinghua Yuan.
# 
# V0: Nov. 14th 2014
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
from astroquery.irsa import Irsa
from astroquery import irsa



# Extinction laws

sourList = ascii.read('../Tables/CsengeriList.dat',
	format = 'fixed_width')
result_dir = '../Tables/'
# Rieke & Lebofsky (1985, ApJ, 288, 618)
AJ = 0.282 # in Av
AH = 0.175
AK = 0.112

# Flaherty et al. (2007, ApJ, 663, 1069)
AW1 = 0.632*AK
AW2 = 0.53*AK
outfile = result_dir+'mipsStats.dat'
ysoStats = open(outfile, 'w')
fmt = '%18s %6i \n'
pfmt = '%18s %6i %6.2f%% finished.' 

start = time.time()
Irsa.ROW_LIMIT = 100000000
Irsa.TIMEOUT = 10000
for isou in range(len(sourList)):
	source = sourList['Name'][isou]
	ra = sourList['ra'][isou]
	dec = sourList['dec'][isou]
	radius = sourList['amaj'][isou]/2.0

	mipsTable = Irsa.query_region(coords.SkyCoord(ra,dec,
		unit = (u.degree,u.degree), frame = "icrs"), 
			catalog = 'mipsgalc',
			spatial = "Cone",
			radius = radius * u.arcsec)
	
	print pfmt %(source, len(mipsTable), ((isou +1.0 ) / len(sourList))*100)
	ysoStats.write(fmt %(source, len(mipsTable)))

ending = time.time()
print (ending-start)/60.0, 'minutes have been used.'
# Phase 1: 
