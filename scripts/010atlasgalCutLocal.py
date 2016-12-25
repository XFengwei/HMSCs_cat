#!/usr/bin python

import os
import math
import time
import tarfile
import shutil
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import montage_wrapper as mt
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table, Column, vstack
from astropy import units as u
import astropy.coordinates as coords

listFile = '../Tables/StarlessAfter24ext.dat'
outDir = '../fitsDir/ATLASGAL/'
sourList = ascii.read(listFile, format = 'fixed_width')

data_dir = '/Volumes/Sailor_1T/ATLAS_GAL/all_data/'
hdr_dir = '/Volumes/Sailor_1T/ATLAS_GAL/hd_dir/'
list_dir = '/Volumes/Sailor_1T/ATLAS_GAL/list_dir/'
#sour_name = raw_input('Please Specify the source name:\n')
#ra_0 = float(raw_input('Please input the RA (in degree):\n'))
#dec_0 = float(raw_input('Please input the Dec (in degree):\n'))
#im_size = float(raw_input('Please input the box size (in degree):\n'))
pfmt = '%6i %18s %6.2f%% finished.'
os.chdir(outDir)

start = time.time()
for isour in range(len(sourList)):
	sour_name = sourList['Name'][isour]
	ra_0 = sourList['ra'][isour]
	dec_0 = sourList['dec'][isour]
	coor = coords.SkyCoord(ra_0, dec_0, frame = 'icrs',
		unit=(u.degree, u.degree))
	l_0 = coor.galactic.l.deg
	b_0 = coor.galactic.b.deg
	im_size = 0.1  # in degree
	pix_num = int(im_size*3600/6.0+1)
	center_pix = pix_num/2

	l_min = l_0-(im_size/2)
	l_max = l_0+(im_size/2)
	b_min = b_0-(im_size/2)
	b_max = b_0-(im_size/2)
	if not os.path.exists(sour_name+'_870.fits'):
		# Modify header files 
		hdr_raw = hdr_dir+'870.hdr'
		hdr_in = open(hdr_raw, 'r')
		hdr_out = open('870.hdr', 'w')
		for line in hdr_in:
			words = line.split()
			if words[0] == 'NAXIS1':
				line = 'NAXIS1  =                 '+str(pix_num)+'\n'
			elif words[0] == 'NAXIS2':
				line = 'NAXIS2  =                 '+str(pix_num)+'\n'
			elif words[0] == 'CRVAL1':
				line = 'CRVAL1  =        '+str(ra_0)+'\n'
			elif words[0] == 'CRVAL2':
				line = 'CRVAL2  =        '+str(dec_0)+'\n'
			elif words[0] == 'CRPIX1':
				line = 'CRPIX1  =            '+str(center_pix)+'\n'
			elif words[0] == 'CRPIX2':
				line = 'CRPIX2  =            '+str(center_pix)+'\n'
			hdr_out.write(line)
		hdr_in.close()
		hdr_out.close()
		
		os.mkdir('870')
		os.mkdir('870/raw')
		
		coor = ascii.read(list_dir+'list.txt')
		for icoor in range(len(coor)):
			if (l_0+3.2 > coor['l'][icoor] and l_0-3.2 < coor['l'][icoor]):
				shutil.copy(data_dir+coor['Name'][icoor], 
					'870/raw/'+coor['Name'][icoor])
		
		
		for bands in ['870']:
			os.chdir(bands)
			os.mkdir('projected')
			os.mkdir('diffdir')
			os.mkdir('corrdir')
			mt.mImgtbl('raw','rimages.tbl')
			mt.mProjExec('rimages.tbl', '../'+bands+'.hdr', 'projected', 'stats.tbl',
				raw_dir='raw')
			mt.mImgtbl('projected', 'pimages.tbl')
			len_dir = len(os.listdir('projected'))
			if len_dir < 3 :
				mt.mAdd('pimages.tbl', '../'+bands+'.hdr', 
					'../'+sour_name+'_'+bands+'.fits', 
					img_dir='projected')
			else:
				mt.mOverlaps('pimages.tbl', 'diffs.tbl')
				mt.mDiffExec('diffs.tbl',  '../'+bands+'.hdr', 'diffdir', 
					proj_dir = 'projected')
				mt.mFitExec('diffs.tbl', 'fits.tbl', 'diffdir')
				if ((len(os.listdir("diffdir")) > 1 and 
					os.path.getsize("diffdir/"+os.listdir("diffdir")[1]) < 10000) or 
					(len(os.listdir("diffdir")) < 1)) :
					listPro = os.listdir('projected')
					listPro = np.array(listPro)
					fileSize = np.array(range(len(listPro)))
					for ifile in range(len(listPro)):
						fileSize[ifile] = os.path.getsize('projected/'+listPro[ifile])
					proTable = Table([listPro, fileSize], names = ["Name", "size"])
					proTable.sort(["size", "Name"])
					shutil.copy('projected/'+proTable['Name'][-2],
						'../'+sour_name+'_'+bands+'.fits')
					shutil.copy('projected/'+proTable['Name'][-1],
						'../'+sour_name+'_'+bands+'_area.fits')
				else:
					mt.mBgModel('pimages.tbl', 'fits.tbl', 'corrections.tbl')
					mt.mBgExec('pimages.tbl', 'corrections.tbl', 'corrdir',
						proj_dir = 'projected')
					mt.mAdd('pimages.tbl', '../'+bands+'.hdr', 
						'../'+sour_name+'_'+bands+'.fits', img_dir = 'corrdir')
			os.chdir('..')
			
		for bands in ['870']:
			shutil.rmtree(bands)
			os.remove(bands+'.hdr')
			os.remove(sour_name+'_'+bands+'_area.fits')

		#convert Jy/beam to MJy/sr
		conFactor = np.pi*(np.deg2rad(5.061775000E-03))**2/4/np.log(2)
		hdu = fits.open(sour_name+'_'+bands+'.fits')
		hdu[0].data = hdu[0].data*1e-6/conFactor
		hdu[0].header['BUNIT'] = 'MJy/sr'
		hdu.writeto(sour_name+'_'+bands+'.fits', clobber = True)

	print  pfmt %(isour+1, sour_name, ((isour+1.0) / len(sourList)*100))

stop = time.time()
dure = stop - start

print "Run time = ",dure, "seconds"

