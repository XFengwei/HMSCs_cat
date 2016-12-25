#!/usr/bin python

import os
import math
import time
import tarfile
import shutil
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
#import montage_wrapper as mt
from astropy.io import ascii
from astropy.table import Table, Column, vstack
from astropy import units as u
import astropy.coordinates as coords
from astroquery.vizier import Vizier

Vizier.ROW_LIMIT = -1
catalog_list = Vizier.find_catalogs('Csengeri 2014')

catalogs = Vizier.get_catalogs(catalog_list.keys())

atlasgalList = catalogs['J/A+A/565/A75/table1']
print len(atlasgalList)
atlasgalList.rename_column("RAJ2000", "ra")
atlasgalList.rename_column("DEJ2000", "dec")

coor = coords.SkyCoord(atlasgalList['ra'], atlasgalList['dec'], 
	frame = 'icrs', unit = {u.deg,u.deg})
l = coor.galactic.l.deg
b = coor.galactic.b.deg
l = Column(l, name = "l")
l[l>180] = l[l>180]-360
b = Column(b, name = "b")
atlasgalList.add_columns([l,b])
atlasgalList = atlasgalList[atlasgalList['l']>-60.0]
atlasgalList = atlasgalList[atlasgalList['b']>-1.0]
atlasgalList = atlasgalList[atlasgalList['b']<1.0]
print len(atlasgalList)
atlasgalList = atlasgalList[atlasgalList['Sp']>0.5]
atlasgalList.remove_columns(["_RAJ2000","_DEJ2000", 'l', 'b'])

atlasgalList.write("../Tables/CsengeriList.dat", 
	format = "ascii.fixed_width")
print(len(atlasgalList))