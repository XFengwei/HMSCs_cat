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
from astropy.table import Table, Column, vstack
from astropy import units as u
import astropy.coordinates as coords
from astroquery.vizier import Vizier

rawList = ascii.read("../Tables/CsengeriList.dat",
	format = 'fixed_width')
ysoList = ascii.read("../Tables/iracysosStats.dat")
mipsList  = ascii.read("../Tables/mipsStats.dat")
simList = ascii.read("../Tables/simbadStats.dat")
rawList.sort("Name")
ysoList.sort("Name")
mipsList.sort("Name")
simList.sort("Name")
sfIndi = ysoList['nYSOs'] + mipsList['nMIPS'] + simList['nSIMBAD']
sfIndi.name = "nSF"
rawList.add_columns([ysoList['nYSOs'], mipsList['nMIPS'],
	simList['nSIMBAD'], sfIndi])
sub1 = rawList[rawList['nSIMBAD']==0]
sub2 = sub1[sub1['nYSOs']==0]

starlessList = rawList[rawList['nSF']==0]
sfList = rawList[rawList['nSF'] > 0]
starlessList.write("../Tables/StarlessListBefore24ext.dat", 
	format = "ascii.fixed_width")

sfList.write("../Tables/sfList.dat", 
	format = "ascii.fixed_width")

print(len(sub1), len(sub2), len(starlessList))