import os
import glob
import time
from astropy.io import fits
from astropy import units as u
from astropy.table import Table, Column, vstack
import astropy.coordinates as coords
import astropy.constants as cons
from scipy.integrate import quad
import numpy as np
from higal_sedfitter import smooth, fit, higal_beams
import matplotlib.pyplot as plt
import matplotlib as mpl
from lmfit import minimize, Parameters, Model
import FITS_tools.header_tools as FITS_header_tools
import aplpy as ap
import shutil

def make_new_dir(dirname):
    if os.path.exists(dirname):
        shutil.rmtree(dirname)
        os.mkdir(dirname)
    else:
        os.mkdir(dirname)

def greybody(nu, nH2=iNh2, Tdust=iTdust, beta=ibeta):
    blackBody = 2*h*nu**3/c**2/(np.exp(h*nu/k/Tdust)-1)
    tau = muh2*mH*nH2*kappa0*(nu/nu0)**beta/rGD
    return blackBody*(1-np.exp(-tau))

# test the sedFitting result of single point
def point_test(x, y):
    # x: pixel num at x coordinate, y: pixel num at y coordinate
    os.chdir('/Users/xufengwei/ATOMS/I16272/ApexHerschel/converted.fits/')
    Tdust_point = fits.getdata('fitresult.fits/Tdust.fits')[x][y]
    Nh2_point = fits.getdata('fitresult.fits/Nh2.fits')[x][y]
    cont070_point = fits.getdata('Herschel_PACS070_converted_smregrid45.fits')[x][y]
    cont160_point = fits.getdata('Herschel_PACS160_converted_smregrid45.fits')[x][y]
    cont250_point = fits.getdata('Herschel_PACS250_converted_smregrid45.fits')[x][y]
    cont350_point = fits.getdata('Herschel_PACS350_converted_smregrid45.fits')[x][y]
    cont500_point = fits.getdata('Herschel_PACS500_converted_smregrid45.fits')[x][y]
    cont870_point = fits.getdata('Apex_ATLASGAL870_converted_smregrid45.fits')[x][y]
    datapoint = [cont070_point,cont160_point,cont250_point,cont350_point,cont500_point,cont870_point]
    wavelengths = np.array([70,160,250,350,500,870])*u.um
    frequencies = wavelengths.to(u.Hz, u.spectral()).value.copy()
    plt.scatter(frequencies, datapoint)
    
    wavLower = 60*u.um
    wavUpper = 1000*u.um
    freqLower = wavUpper.to(u.Hz, u.spectral()).value
    freqUpper = wavLower.to(u.Hz, u.spectral()).value

    freq_range = np.linspace(freqLower,freqUpper,200)
    fitflux = 10**17*greybody(freq_range,Nh2_point,Tdust_point,beta=2)
#     fitflux2 = 10**17*greybody(freq_range,Nh2_point,23.1,beta=2)
    print(Nh2_point)
    print(Tdust_point)
    plt.plot(freq_range, fitflux,color='b')
#     plt.plot(freq_range, fitflux2,'r')
    plt.xscale('log')
    plt.yscale('log')
    #plt.ylim(500,2000)
    
def plot_result(parameter, distance, sourcename, center_x, center_y, crop_size=0.16, recenter=True, add_grid=True,
                add_beam=True, beam_size=True, add_scalebar=True):
    os.chdir('/Users/xufengwei/ATOMS/I16272/ApexHerschel/converted.fits/fitresult.fits/')
    parameter_dir = {'Tdust': 'Tdust.fits', 'Nh2':'Nh2.fits', 'intSED':'intSED.fits'}
    img = ap.FITSFigure(parameter_dir[parameter])
    img.show_colorscale(stretch='log',cmap='hot')
    img.add_colorbar()
    img.colorbar.set_width(0.3)
    if recenter:
        img.recenter(x=center_x, y=center_y, width=crop_size, height=crop_size)
    if add_grid:
        img.add_grid()
    if add_beam:
        img.add_beam(color='black',hatch='',linewidth=3, alpha=0.1)
        img.beam.set_edgecolor('white')
        img.add_label(0.1,0.15,'beam size',relative='axes',color='k',fontsize=15)
    if add_scalebar:
        img.add_scalebar(0.2*u.pc/(distance*10**3*u.pc)*180/np.pi*u.deg, corner='top left',color='white')
        img.scalebar.set_label('0.2 pc')
    img.set_title(parameter+' distribution of '+sourcename)
    img.savefig('Tdust.png',dpi=600)
#     img.close()    

# record the time
start = time.time()

# renew header of fits file
# no need for Apex header
os.chdir('/Users/xufengwei/ATOMS/I16272/ApexHerschel/')
make_new_dir('newhd.fits')
for Herschel_file in glob.glob('Herschel*.fits'):
    hdu = fits.open(Herschel_file)
    fits.writeto('newhd.fits/'+Herschel_file[:-5]+'_newhd.fits', 
                 header=hdu[1].header, data=hdu[1].data, overwrite=True)

# convert unit to MJy/beam
# get the beam size of six bands and calculate the beam factor of each of them
bmaj = np.array(
        [FITS_header_tools.header_to_platescale(fits.open('Herschel_PACS070.fits')[1].header,use_units=True).value,
        FITS_header_tools.header_to_platescale(fits.open('Herschel_PACS160.fits')[1].header,use_units=True).value,
        fits.getheader('Herschel_SPIRE250.fits')['BMAJ'],
        fits.getheader('Herschel_SPIRE350.fits')['BMAJ'],
        fits.getheader('Herschel_SPIRE500.fits')['BMAJ'],
        fits.getheader('Apex_ATLASGAL870.fits')['BMAJ']])
bmin = np.array(
        [FITS_header_tools.header_to_platescale(fits.open('Herschel_PACS070.fits')[1].header,use_units=True).value,
        FITS_header_tools.header_to_platescale(fits.open('Herschel_PACS160.fits')[1].header,use_units=True).value,
        fits.getheader('Herschel_SPIRE250.fits')['BMIN'],
        fits.getheader('Herschel_SPIRE350.fits')['BMIN'],
        fits.getheader('Herschel_SPIRE500.fits')['BMIN'],
        fits.getheader('Apex_ATLASGAL870.fits')['BMIN']])

beam_factor = bmaj*bmin*3600**2*np.array([2.35E-5,2.35E-5,2.66E-5,2.66E-5,2.66E-5,2.66E-5])
beam_factor_dir = {'070':beam_factor[0],'160':beam_factor[1],'250':beam_factor[2],
                   '350':beam_factor[3],'500':beam_factor[4],'870':beam_factor[5]}

# data unit convertion
make_new_dir('converted.fits')
for newhd_file in glob.glob('newhd.fits/*_newhd.fits'):
    band = newhd_file[-14:-11]
    hdu = fits.open(newhd_file)
    fits.writeto('converted.fits/'+newhd_file[11:-11]+'_converted.fits', 
                 header=hdu[0].header, data=hdu[0].data/beam_factor_dir[band], overwrite=True)
for Apex_file in glob.glob('Apex_ATLASGAL870.fits'):
    band = Apex_file[-8:-5]
    hdu = fits.open(Apex_file)
    fits.writeto('converted.fits/'+Apex_file[0:-5]+'_converted.fits', 
                 header=hdu[0].header, data=hdu[0].data/beam_factor_dir[band], overwrite=True)
    
# smooth
os.chdir('/Users/xufengwei/ATOMS/I16272/ApexHerschel/converted.fits/')
for fn in glob.glob('*converted.fits'):
    if fn != 'Apex_ATLASGAL870_converted.fits':
        smooth.add_beam_information_to_higal_header(fn, name_to_um=higal_beams.num_to_um, clobber=True)

# since the longest wavelength corresponds to lowest resolution
# Apex resoultion is about 19 arcsec
# All bands are smoothed to 500 um
target_fn = glob.glob('*500_converted.fits')[0]
target_header = fits.getheader(target_fn)
smooth.smooth_images_toresolution(45*u.arcsec, skip_existing=False,
                                  globs=['*_converted.fits'],
                                  regrid=True,
                                  target_header=target_header, clobber=True)

# start to fit
# Frequency interval
wavLower = 1*u.um
wavUpper = 3000*u.um

freqLower = wavUpper.to(u.Hz, u.spectral()).value
freqUpper = wavLower.to(u.Hz, u.spectral()).value

# Initiate the parameters
iNh2 = 1.0E+22 # in cm^-2
iTdust = 22.0 # in K.
ibeta = 2.0 # 
betaVary = False # If True, Beta will be fitted. If False, beta will be fixed.
betaMin = 1.0
betaMax = 3.0

wavelengths = np.array([70,160,250,350,500,870])*u.um
frequencies = wavelengths.to(u.Hz, u.spectral()).value.copy()

# Constants
h = cons.h.cgs.value   # Planck constant in CGS unit
k = cons.k_B.cgs.value # Boltzmann constant in CGS unit
c = cons.c.cgs.value # speed of light in CGS unit
mH = cons.m_p.cgs.value # mass of an neutron
muh2 = 2.8 # mean molecular weight adopted from Kauffmann et al. (2008)
rGD = 100.0 # gas-to-dust mass ratio
nu0 = 599.584916E9 # Reference frequency in Hz.
kappa0 = 5.0/1.5 # Dust emissivity at reference frequency

# Define the model
gMod = Model(greybody)
gMod.set_param_hint('beta', min = betaMin, max = betaMax, vary = betaVary)
gMod.set_param_hint('Tdust',  max = 80.0)
pars = gMod.make_params()

# initiate the model
# copy the header format and define the unit
hdrRaw = fits.getheader('Herschel_PACS160_converted_smregrid45.fits') 
hdrNh2 = hdrRaw.copy()
hdrInt = hdrRaw.copy()
hdrTdust = hdrRaw.copy()
hdrNh2['BUNIT'] = 'cm^-2'
hdrInt['BUNIT'] = 'g s^-3'
hdrTdust['BUNIT'] = 'K'

# get the 2d data
# unit is MJy/sr
data070 = fits.getdata('Herschel_PACS070_converted_smregrid45.fits')
data160 = fits.getdata('Herschel_PACS160_converted_smregrid45.fits')
data250 = fits.getdata('Herschel_SPIRE250_converted_smregrid45.fits')
data350 = fits.getdata('Herschel_SPIRE350_converted_smregrid45.fits')
data500 = fits.getdata('Herschel_SPIRE500_converted_smregrid45.fits')
data870 = fits.getdata('Apex_ATLASGAL870_converted_smregrid45.fits')

# data error, unit: MJy/beam
# manual copy from ds9
# choose a specific region first
# rms: root mean square, flatfield: mean value or baseline, stddev: standard deviation
rms = np.array([0.100218/beam_factor[0],
                0.318015/beam_factor[1], 
                6.21757/beam_factor[2], 
                5.05017/beam_factor[3], 
                3.97634/beam_factor[4],
                0.535837/beam_factor[5]])*u.MJy/u.sr
flatfield = np.array([0.0992824/beam_factor[0],
                      0.315178/beam_factor[1], 
                      6.09269/beam_factor[2], 
                      4.93347/beam_factor[3], 
                      3.87736/beam_factor[4],
                      0.532882/beam_factor[5]])*u.MJy/u.sr
stddev = np.array([0.0136641/beam_factor[0],
                   0.0423781/beam_factor[1],
                   1.23993/beam_factor[2], 
                   1.07941/beam_factor[3], 
                   0.881712/beam_factor[4],
                   0.0561999/beam_factor[5]])*u.MJy/u.sr
solid_signal = flatfield + stddev*3

weis = 1.0/stddev

# Initiate the 2d array of three parameters
dataNh2 = data500.copy()
dataTdust = data500.copy()
dataIntegrate = data500.copy()

# Pixel-by-pixel fitting
for iCol in range(len(data500[:,1])):
    for iRow in range(len(data500[1,:])):
        data = np.array([data070[iCol,iRow], data160[iCol,iRow], data250[iCol,iRow], 
                         data350[iCol,iRow], data500[iCol,iRow], data870[iCol,iRow]])*u.MJy/u.sr
        if ((data070[iCol,iRow] > solid_signal[0].value)
            and (data160[iCol,iRow] > solid_signal[1].value)
            and (data250[iCol,iRow] > 0) 
            and (data350[iCol,iRow] > 0) 
            and (data500[iCol,iRow] > 0)
            and (data870[iCol,iRow] > 0)):
            sedResult = gMod.fit(data.cgs.value, nu = frequencies, weights = weis.cgs.value)
            inTSED = quad(greybody, freqLower, freqUpper, 
                          args = (sedResult.params['nH2'].value, sedResult.params['Tdust'].value, 
                                  sedResult.params['beta'].value))
            dataNh2[iCol,iRow] = sedResult.params['nH2'].value
            dataTdust[iCol,iRow] = sedResult.params['Tdust'].value
            dataIntegrate[iCol,iRow] = inTSED[0]
        else:
            dataNh2[iCol,iRow] = np.nan
            dataTdust[iCol,iRow] = np.nan
            dataIntegrate[iCol,iRow] = np.nan

# write the three parameters into fits files
make_new_dir('fitresult.fits')
fits.writeto('fitresult.fits/Nh2.fits', data = dataNh2, 
    header = hdrNh2, clobber = True)
fits.writeto('fitresult.fits/Tdust.fits', data = dataTdust, 
    header = hdrTdust, clobber = True)
fits.writeto('fitresult.fits/intSED.fits', data = dataIntegrate, 
    header = hdrInt, clobber = True)

# plot the result
plot_result('Tdust', 2.92, 'I16272', center_x=247.7514630, center_y=-48.7257175)

# record the time
stop = time.time()
dure = stop - start

# display how much time the code runs
print("Run time = ",dure, "seconds")
