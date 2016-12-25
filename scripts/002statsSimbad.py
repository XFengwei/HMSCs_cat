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
from astroquery.simbad import Simbad


# Extinction laws
listFile = '../Tables/CsengeriList.dat'
sourCoor = ascii.read(listFile, format = 'fixed_width')
result_dir = '../Tables/'

outfile = result_dir+'simbadStats.dat'
simStats = open(outfile, 'w')
fmt = '%18s %6i \n'
pfmt = '%6i %18s %6i %6.2f%% finished.' 

Simbad.TIMEOUT = 1000000
customSimbad = Simbad()
customSimbad.add_votable_fields('ra(d)','dec(d)', 'otype')
customSimbad.remove_votable_fields('coordinates')
customSimbad.TIMEOUT = 100000
start = time.time()
for isou in range(len(sourCoor)):
	ra = sourCoor['ra'][isou]
	dec = sourCoor['dec'][isou]
	c = coords.SkyCoord(ra, dec, unit=('deg','deg'), frame='icrs')
	source = sourCoor['Name'][isou]
	radius = sourCoor['amaj'][isou]/2.0
	result = customSimbad.query_region(c, radius = radius * u.arcsec)
	lenSF = 0
	if result:
		for isim in range(len(result)):
			if (result['OTYPE'][isim] == 'Radio(cm)' or
				 result['OTYPE'][isim] == 'Radio' or
				 result['OTYPE'][isim] == 'HII' or
				 result['OTYPE'][isim] == 'IR' or
				 result['OTYPE'][isim] == 'IR>30um' or
				 result['OTYPE'][isim] == 'IR<30um' or
				 result['OTYPE'][isim] == 'Maser' or
				 result['OTYPE'][isim] == 'HH' or
				 result['OTYPE'][isim] == 'Outflow' or
				 result['OTYPE'][isim] == 'outflow?' or
				 result['OTYPE'][isim] == 'Candidate_YSO' or 
				 result['OTYPE'][isim] == 'Candidate_pMS*' or
				 result['OTYPE'][isim] == 'Candidate_TTau*' or
				 result['OTYPE'][isim] == 'Candidate_Ae*' or
				 result['OTYPE'][isim] == 'YSO' or result['OTYPE'][isim] == 'Ae*' or
				 result['OTYPE'][isim] == 'pMS*' or result['OTYPE'][isim] == 'TTau*'):
				 lenSF += 1 

	print(pfmt %(isou, source, lenSF, ((isou+1.0) / len(sourCoor)*100)))
	simStats.write(fmt %(source, lenSF))	

ending = time.time()
print((ending-start)/60.0, 'minutes have been used.')
# Phase 1: 
