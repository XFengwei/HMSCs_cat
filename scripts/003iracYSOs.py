#!/usr/bin python
#
# This script is used to identity and classify YSOs based on 
# the scheme of Gutermuth et al. (2009), ApJS, 184, 18
# 
# Originally written by Jinghua Yuan.
# 
# V0: Dec. 30th 2013
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



sourList = ascii.read('../Tables/hmscList_full_20160628.txt',
	format = 'fixed_width')
result_dir = '../Tables/'

outfile = result_dir+'iracysosStats.dat'
ysoStats = open(outfile, 'w')
fmt = '%18s %6i \n'
pfmt = '%6i %18s %6i %6.2f%% finished.' 

start = time.time()
Irsa.ROW_LIMIT = 100000000
Irsa.TIMEOUT = 10000
for isou in range(len(sourList)):
	source = sourList['Name'][isou]
	ra = sourList['ra'][isou]
	dec = sourList['dec'][isou]
	radius = sourList['amaj'][isou]/2.0
	coor = coords.SkyCoord(ra,dec, unit = (u.degree,u.degree),
	 	frame = "icrs")
	if (coor.galactic.l.deg > 10.8 and coor.galactic.l.deg < 349.2): 
		raw_psc = Irsa.query_region(coor, catalog = 'glimpse_s07',
				spatial = "Cone", radius = radius * u.arcsec)
	else :
		raw_psc = Irsa.query_region(coor, catalog = 'glimpse2_v2cat',
				spatial = "Cone", radius = radius * u.arcsec)

	raw_psc['mag3_6'].mask[raw_psc['d3_6m']>=0.2]=True
	raw_psc['mag4_5'].mask[raw_psc['d4_5m']>=0.2]=True
	raw_psc['color36_45'].mask[raw_psc['d3_6m']>=0.2]=True
	raw_psc['color36_45'].mask[raw_psc['d4_5m']>=0.2]=True
	raw_psc['color36_58'].mask[raw_psc['d3_6m']>=0.2]=True
	raw_psc['color45_58'].mask[raw_psc['d4_5m']>=0.2]=True
	raw_psc['color36_80'].mask[raw_psc['d3_6m']>=0.2]=True
	raw_psc['color45_80'].mask[raw_psc['d4_5m']>=0.2]=True
	raw_psc['mag5_8'].mask[raw_psc['d5_8m']>=0.2]=True
	raw_psc['mag8_0'].mask[raw_psc['d8_0m']>=0.2]=True
	raw_psc['color36_58'].mask[raw_psc['d5_8m']>=0.2]=True
	raw_psc['color36_80'].mask[raw_psc['d8_0m']>=0.2]=True
	raw_psc['color58_80'].mask[raw_psc['d5_8m']>=0.2]=True
	raw_psc['color58_80'].mask[raw_psc['d8_0m']>=0.2]=True
	raw_psc['d3_6m'].mask[raw_psc['d3_6m']>=0.2]=True
	raw_psc['d4_5m'].mask[raw_psc['d4_5m']>=0.2]=True
	
	sig_I1_I2 = np.sqrt(raw_psc['d3_6m']**2+raw_psc['d4_5m']**2)
	
	sig_I1_I2.name = 'sig_I1_I2'
	sig_I2_I3 = np.sqrt(raw_psc['d4_5m']**2+ raw_psc['d5_8m']**2)
	
	sig_I2_I3.name = 'sig_I2_I3'
	sig_I2_I4 = np.sqrt(raw_psc['d4_5m']**2+raw_psc['d8_0m']**2)
	
	sig_I2_I4.name = 'sig_I2_I4'
	sig_I1_I3 = np.sqrt(raw_psc['d3_6m']**2+raw_psc['d5_8m']**2)
	
	sig_I1_I3.name = 'sig_I1_I3'

	sig_I3_I4 =  np.sqrt(raw_psc['d5_8m']**2+raw_psc['d8_0m']**2)
	
	sig_I3_I4.name = 'sig_I3_I4'
	
	yso_type = raw_psc['designation'][raw_psc['ra']>0]
	yso_type.name = 'yso_type'
	yso_type[0:]= 'field'
	# Add new clumns to the table
	raw_psc.add_columns([sig_I1_I2, sig_I2_I3, sig_I2_I4, 
			sig_I1_I3, sig_I3_I4, yso_type])
	
	irac_classI = raw_psc[raw_psc['ra']<0]
	irac_classII = raw_psc[raw_psc['ra']<0]
	
	for iraw in range(len(raw_psc)):
		if ((not raw_psc['sig_I3_I4'].mask[iraw]) and
			raw_psc['color36_45'][iraw]>0.7 and
			raw_psc['color45_58'][iraw]>0.7):
			raw_psc['yso_type'][iraw] = 'I'
			irac_classI.add_row(raw_psc[iraw])
	
		elif ((not raw_psc['sig_I3_I4'].mask[iraw]) and
			raw_psc['color45_80'][iraw]-raw_psc['sig_I2_I4'][iraw]>0.5 and
			raw_psc['color36_58'][iraw]-raw_psc['sig_I1_I3'][iraw]>0.35 and
			raw_psc['color36_58'][iraw]+
			raw_psc['sig_I1_I3'][iraw]<=0.14/0.04*
			(raw_psc['color45_80'][iraw]
			-raw_psc['sig_I2_I4'][iraw]-0.5)+0.5 and
			raw_psc['color36_45'][iraw]-raw_psc['sig_I1_I2'][iraw]>0.15):
			raw_psc['yso_type'][iraw] = 'II'
			irac_classII.add_row(raw_psc[iraw])
	
	all_ysos = vstack([irac_classI, irac_classII])
	print(pfmt %(isou+1, source, len(all_ysos), ((isou+1.0) / len(sourList)*100)))
	ysoStats.write(fmt %(source, len(all_ysos)))
	
ending = time.time()

print((ending-start)/60.0, 'minutes have been used.')

# Phase 1: 
