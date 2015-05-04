#########################################################
#SCRIPT TO EXECUTE THE KMOS PIPELINE AUTOMATICALLY 
#BUILDING THE .SOF FILES EN ROUTE. MUST BE EXECUTED FROM 
#WITHIN THE $KMOS_DYN_CALIBRATIONS DIRECTORY 
########################################################




#import the relevant modules
import os, sys, numpy as np, random, matplotlib.mlab as mlab, matplotlib.pyplot as plt, math, matplotlib.cm as cm
import pyraf
import numpy.polynomial.polynomial as poly
import lmfit
import scipy
import math
from lmfit.models import GaussianModel, ExponentialModel, LorentzianModel, VoigtModel, PolynomialModel
from scipy import stats
from scipy import optimize
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
from pylab import *
from matplotlib.colors import LogNorm

#add the class file to the PYTHONPATH
sys.path.append('/Users/owenturner/Documents/PhD/KMOS/Analysis_Pipeline/Python_code/Class')

#import the class 
from pipelineClass import pipelineOps
from cubeClass import cubeOps

#Create an instance of the class 
pipe_methods = pipelineOps()


#Writing robust pipeline that will generate the .sof files at each stage. 
#Put each of the different wavebands in a separate folder and execute this recipe 
#still need to figure out how to navigate the directory structure from the raw files 
#Probably need to have the environment variables set each time. 
#Check to see if log.txt file exists 

if os.path.isfile('log.txt'):
	os.system('rm log.txt')

#Check to see if the directory variables have been set in the first place 
try:
	#If any are False - need to exit. Variables aren't set 
	if os.environ['KMOS_SCIENCE'] and os.environ['KMOS_DYN_CALIBRATIONS'] and os.environ['KMOS_RAW']:
		print 'Directory Variables are set'
except KeyError:
	print 'The Directory variables are not set'
	raise

#Check to see if the user has changed the $KMOS_SCIENCE, $KMOS_RAW and $KMOS_DYN_CALIBRATIONS environment variables 
if raw_input('Are you happy with the values of the directory variables?\nKMOS_SCIENCE: %s\nKMOS_RAW: %s\nKMOS_DYN_CALIBRATIONS: %s\n(Y/N): ' % \
											(os.environ['KMOS_SCIENCE'], os.environ['KMOS_RAW'],os.environ['KMOS_DYN_CALIBRATIONS'])) == 'Y':
	print 'Continuing Calibrations and Reduction'
else:
	raise ValueError('Did Not Change directory variables... Quitting')

#First find the grating ID of one of the object,sky files in the $KMOS_RAW directory
#This will be used later on in the script to construct the correct .sof files  
#Write out the types of the raw files to the log 
os.system('dfits $KMOS_RAW/*.fits | fitsort dpr.type ins.grat1.id >> log.txt')

#reload the log file and only keep those with OBJECT,SKY type 
Table = np.loadtxt('log.txt', dtype='str')
tupe = zip(Table[:,0], Table[:,1], Table[:,2])
zip_names = []
obj_names = []
#Loop around searching for the keyword OBJECT,SKY 
for entry in tupe:
    for i in range(len(entry)):
        if entry[i].find('OBJECT,SKY') != -1:
            zip_names.append(entry)
            #Create the object_names.txt file
            #obj_names_tup = (entry[0], entry[3])
            #obj_names.append(obj_names_tup)

grating_ID = zip_names[0][2]
os.system('rm log.txt')

#Write out to the object.txt file
#if os.path.isfile('object_names.txt'):
#	os.system('rm object_names.txt')
#with open('object_names.txt', 'a') as f:
#	for entry in obj_names:
#		f.write('%s\t%s\n' % (entry[0], entry[1]))


#Happy that directory variables are there and that they are pointing to the right place 
#Now change directory to the Calibrations folder and construct the dark sof 
##########################################
#STEP 1: DARK FRAMES AND BAD PIXEL MAP
#
#Created sof name:	dark.sof
##########################################

#Write out the types of the raw files to the log 
os.system('dfits $KMOS_RAW/*.fits | fitsort dpr.type >> log.txt')

#reload the log file and only keep those with DARK type 
Table = np.loadtxt('log.txt', dtype='str')
tupe = zip(Table[:,0], Table[:,1])
zip_names = []

#Loop around searching for the keyword DARK
for entry in tupe:
    for i in range(len(entry)):
        if entry[i].find('DARK') != -1:
            zip_names.append(entry)

#Check for existence of dark.sof 
if os.path.isfile('dark.sof'):
	os.system('rm dark.sof')

with open('dark.sof', 'a') as f:
	for entry in zip_names:
		f.write('%s\t%s\n' % (entry[0], entry[1]))
os.system('rm log.txt')

print 'Computing dark frames'

#os.system('esorex kmos_dark dark.sof')

#Generates the two files master_dark.fits and badpixel_dark.fits that are used 
#By the flatfield recipe

#We want to grow the badpixel mask, to account for the 
#pixels surrounding the bad ones which are also saturated 

print 'Extending bad pixel mask'

pipe_methods.badPixelextend(badpmap='badpixel_dark.fits')

###################################################
#STEP 2: FLATFIELDING
#
#Required .sof filename: flatfield.sof
#in the sof require 3 flat,off dpr.type file s
#and 18 flat,lamp files which correspond to 
#3 files for each of the 6 different rotator angles
###################################################

#Play the same game with generating the .sof file
#Write out the types of the raw files to the log 
os.system('dfits $KMOS_RAW/*.fits | fitsort dpr.type >> log.txt')

#reload the log file and only keep those with DARK type 
Table = np.loadtxt('log.txt', dtype='str')
tupe = zip(Table[:,0], Table[:,1])
zip_names = []

#Loop around searching for the keyword DARK
for entry in tupe:
    for i in range(len(entry)):
        if entry[i].find('FLAT,OFF') != -1:
			new_entry = (entry[0], str(entry[i].replace('FLAT,OFF', 'FLAT_OFF')))
			zip_names.append(new_entry)

        elif entry[i].find('FLAT,LAMP') != -1:
			new_entry = (entry[0], str(entry[i].replace('FLAT,LAMP', 'FLAT_ON')))
			zip_names.append(new_entry)

#Check for existence of flat.sof 
if os.path.isfile('flat.sof'):
	os.system('rm flat.sof')

with open('flat.sof', 'a') as f:
	for entry in zip_names:
		f.write('%s\t%s\n' % (entry[0], entry[1]))
os.system('rm log.txt')
#Also add the badpixel_dark_added.fits file to the flat.sof
with open('flat.sof', 'a') as f:
	f.write('badpixel_dark_Added.fits\tBADPIXEL_DARK')


print 'Computing flat fields for the different potato angles'

#os.system('esorex kmos_flat flat.sof')

####################################################
#STEP 3: ARCS
#
#Created .sof name: wave.sof
#Require the products from both the flat and dark recipes 
#badpixel_dark_Added.fits, badpixel_flat_HHH.fits, 
#flat_edge_HHH.fits, ycal_HHH.fits, xcal_HHH.fits, 
#master_flat_HHH.fits
#AND 6 ARC_ON, 1 ARC_OFF and 3 static calibration 
#files listed in the document 
##################################################

#Play the same game with generating the .sof file
#Write out the types of the raw files to the log 
os.system('dfits $KMOS_RAW/*.fits | fitsort dpr.type >> log.txt')

#reload the log file and only keep those with DARK type 
Table = np.loadtxt('log.txt', dtype='str')
tupe = zip(Table[:,0], Table[:,1])
zip_names = []

#Loop around searching for the keyword DARK
for entry in tupe:
    for i in range(len(entry)):
        if entry[i].find('WAVE,OFF') != -1:
			new_entry = (entry[0], str(entry[i].replace('WAVE,OFF', 'ARC_OFF')))
			zip_names.append(new_entry)

        elif entry[i].find('WAVE,LAMP') != -1:
			new_entry = (entry[0], str(entry[i].replace('WAVE,LAMP', 'ARC_ON')))
			zip_names.append(new_entry)

#Check for existence of flat.sof 
if os.path.isfile('wave.sof'):
	os.system('rm wave.sof')

with open('wave.sof', 'a') as f:
	for entry in zip_names:
		f.write('%s\t%s\n' % (entry[0], entry[1]))
os.system('rm log.txt')

#Also add the calibration products to the wave.sof
calib_dir = os.environ['KMOS_CALIB']
#This is the first time the grating dependence comes into play 
if grating_ID == 'H': 
	neon_name = calib_dir + '/kmos_ar_ne_list_h.fits'
	bpixel_flat_name = 'badpixel_flat_HHH.fits'
	flat_edge_name = 'flat_edge_HHH.fits'
	master_flat_name = 'master_flat_HHH.fits'
	xcal_name = 'xcal_HHH.fits'
	ycal_name = 'ycal_HHH.fits'
	lcal_name = 'lcal_HHH.fits'
elif grating_ID == 'YJ': 
	neon_name = calib_dir + '/kmos_ar_ne_list_yj.fits'
	bpixel_flat_name = 'badpixel_flat_YJYJYJ.fits'
	flat_edge_name = 'flat_edge_YJYJYJ.fits'
	master_flat_name = 'master_flat_YJYJYJ.fits'
	xcal_name = 'xcal_YJYJYJ.fits'
	ycal_name = 'ycal_YJYJYJ.fits'
	lcal_name = 'lcal_YJYJYJ.fits'
elif grating_ID == 'IZ': 
	neon_name = calib_dir + '/kmos_ar_ne_list_iz.fits'
	bpixel_flat_name = 'badpixel_flat_IZIZIZ.fits'
	flat_edge_name = 'flat_edge_IZIZIZ.fits'
	master_flat_name = 'master_flat_IZIZIZ.fits'
	xcal_name = 'xcal_IZIZIZ.fits'
	ycal_name = 'ycal_IZIZIZ.fits'
	lcal_name = 'lcal_IZIZIZ.fits'
elif grating_ID == 'HK': 
	neon_name = calib_dir + '/kmos_ar_ne_list_hk.fits'
	bpixel_flat_name = 'badpixel_flat_HKHKHK.fits'
	flat_edge_name = 'flat_edge_HKHKHK.fits'
	master_flat_name = 'master_flat_HKHKHK.fits'
	xcal_name = 'xcal_HKHKHK.fits'
	ycal_name = 'ycal_HKHKHK.fits'
	lcal_name = 'lcal_HKHKHK.fits'
elif grating_ID == 'K': 
	neon_name = calib_dir + '/kmos_ar_ne_list_k.fits'
	bpixel_flat_name = 'badpixel_flat_KKK.fits'
	flat_edge_name = 'flat_edge_KKK.fits'
	master_flat_name = 'master_flat_KKK.fits'
	xcal_name = 'xcal_KKK.fits'
	ycal_name = 'ycal_KKK.fits'
	lcal_name = 'lcal_KKK.fits'
else:
	raise ValueError('Did not recognise grating ID')


waveband_name = calib_dir + '/kmos_wave_band.fits'
reftable_name = calib_dir + '/kmos_wave_ref_table.fits'	
with open('wave.sof', 'a') as f:
	f.write('%s\tBADPIXEL_FLAT\n' % bpixel_flat_name)
	f.write('%s\tFLAT_EDGE\n' % flat_edge_name)
	f.write('%s\tMASTER_FLAT\n' % master_flat_name)
	f.write('%s\tXCAL\n' % xcal_name)
	f.write('%s\tYCAL\n' % ycal_name)
	f.write('%s\tARC_LIST\n' % neon_name)
	f.write('%s\tWAVE_BAND\n' % waveband_name)
	f.write('%s\tREF_LINES\n' % reftable_name)
	


print 'Computing wavelength calibration'

#os.system('esorex kmos_wave_cal wave.sof')

#####################################################
#STEP 4: ILLUMINATION CORRECTION (Optional)
#
#Required .sof filename: illum_cor.sof
#####################################################

#Leaving out the illumination correction at the moment 
#Not an important part of the pipeline
#print 'Attempting to create illumination correction'

#os.system('esorex kmos_illumination %s' % illum_cor_sof)

#######################################################
#STEP 5: STANDARD STARS
#
#Created .sof file star.sof
#######################################################

#print 'Calibrating Standard Stars'

#os.system('esorex kmo_std_star %s' % std_star_sof)

#At this stage we've reached the end of the calibration. Ready for correcting  
#the readout column noise and performing the sub-pixel alignment

#For the most basic run of the pipeline, don't need to process 
#standard stars. Will skip this for now.

#########################################################
#STEP 6: READOUT COLUMN NOISE
#
#Required - list of object and skyfiles piped from dfits 
#called object_names.txt
#########################################################

print 'Correcting for readout column bias'



pipe_methods.applyCorrection('object_names.txt', 'badpixel_dark_Added.fits', lcal_name)

########################################################
#STEP 7: SHIFT AND ALIGN
#
#Required - list of corrected files from above 
#called corrected_object_names.txt 
#Be sure before going this far - this step will take ages
#There are some specifiable things in this function. 
#Later on I'll look at adding some options and variables 
#to this part of the code using an argument parser
########################################################

#print 'Shifting and alligning all images - this will take some time'

#pipe_methods.applyShiftAllExtensionsMin(fileList = corrected_object_names, badpmap='badpixel_dark_Added.fits', \
#  	 vertSegments=1, horSegments=1, interp_type='spline3')

###########################################################

#######################################################################
#SCIENCE REDUCTION: Use everything we've made and create the data cubes
#
#Required: reduction sof sci_reduc.sof
#This file will contain all the corrected and shifted observations 
#produced by the above two manual recipes. Make sure they have the correct 
#names in the .sof file.
#Also will move directory into science folder 
#Must have a folder called Science_Output in the directory above 
#the callibration directory
########################################################################

#print 'Starting Data Cube Reconstruction'

#os.system('up')
#os.system('cd Science_Output')

#os.system('esorex kmo_sci_red %s' % sci_reduc_sof)

#Should all work given that the correct .sof files are given and that these contain the right things 
#Need to think about automatically generating these and the object names later. 