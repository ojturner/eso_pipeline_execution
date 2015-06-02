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
import fnmatch
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
sys.path.append('/disk1/turner/PhD/KMOS/Analysis_Pipeline/Python_code/Class')

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

#Get the values of the directory variables 
calib_dir = os.environ['KMOS_CALIB']
raw_dir = os.environ['KMOS_RAW']
science_dir = os.environ['KMOS_SCIENCE']
dyn_dir = os.environ['KMOS_DYN_CALIBRATIONS']

#Check for the existence of any Subtracted files in the raw directory and delete if they exist 
for file in os.listdir(os.environ['KMOS_RAW']):
    if fnmatch.fnmatch(file, '*Subtracted*.fits'):
        os.system('rm $KMOS_RAW/*Subtracted*.fits')
#Should clear up the problems with creating object_names.txt files with subtracted files present
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
        if entry[i].find('OBJECT,SKY') != -1 and entry[i].find('OBJECT,SKY,STD,FLUX') == -1 :
			zip_names.append(entry)
            #Append only the names of these files to the obj_names vector
			obj_names.append(entry[0])
#create a new temporary directory for the object,sky files 
if os.path.isdir('temp/'):
	os.system('rm -rf temp/')
os.system('mkdir temp/')
#make sure that only the raw images make it into the object_names.txt file
corr_list = []
#loop around the obj_names and check for the word 'Corrected'
for item in obj_names:
	if item.find('Corrected') != -1:
		corr_list.append(item)
#Now remove from the obj_names array any of the flagged files 
for item in corr_list:
	obj_names.remove(item)
#obj_names should now contain only the raw object,sky files 
#copy the object,sky files into this temporary directory
for item in obj_names:
	os.system('cp %s temp/' % item)
#Now simply pipe the output from dfits to the object_names.txt file
if os.path.isfile('temp.txt'):
	os.system('rm temp.txt')

os.system('dfits temp/*.fits | fitsort ocs.arm1.type >> temp.txt')

#Stupidly this also leaves the directory in the output
#Need to reload the temporary file, change the occurence of to the $KMOS_RAW directory
object_names_Table = np.loadtxt('temp.txt', dtype='str')
new_list = []
sky_list = []
for i in range(len(object_names_Table)):
	temp_new_list = []
	for j in range(len(object_names_Table[i])):
		newString = object_names_Table[i][j].replace('temp/', raw_dir + '/')
		temp_new_list.append(newString)
		if object_names_Table[i][j] == 'S':
			sky_list.append(object_names_Table[i][0])
	new_list.append(temp_new_list)

for i in range(len(sky_list)):
	sky_list[i] = sky_list[i][5:]
print '[INFO]: This is the list of sky_frames: %s' % sky_list
#Finally ready for writing to the object_names.txt file 
if os.path.isfile('object_names.txt'):
	os.system('rm object_names.txt')
with open('object_names.txt', 'a') as f:
	for entry in new_list:
		f.write('%s\t%s\n' % (entry[0], entry[1]))
#clean up by removing the temporary directory 
os.system('rm -rf temp/')
os.system('rm temp.txt')

grating_ID = zip_names[0][2]
os.system('rm log.txt')

#Apply initial subtraction, placing files in a sub-directory within the raw folder 
sub_dir_name = raw_dir + '/Subtracted'
if os.path.isdir(sub_dir_name):
	os.system('rm -rf %s' % sub_dir_name)
os.system('mkdir %s' % sub_dir_name) 
#Apply the subtraction method 
pipe_methods.applySubtraction('object_names.txt')
#Move the created files into the subtracted directory 
os.system('mv $KMOS_RAW/*Subtracted.fits %s' % sub_dir_name)

#After defining the grating ID, set the values of the created files 
#conditionally, depending on which grating is being used

if grating_ID == 'H': 
	neon_name = calib_dir + '/kmos_ar_ne_list_h.fits'
	oh_name = calib_dir + '/kmos_oh_spec_h.fits'
	bpixel_flat_name = 'badpixel_flat_HHH.fits'
	flat_edge_name = 'flat_edge_HHH.fits'
	master_flat_name = 'master_flat_HHH.fits'
	xcal_name = 'xcal_HHH.fits'
	ycal_name = 'ycal_HHH.fits'
	lcal_name = 'lcal_HHH.fits'
	telluric_name = 'telluric_HHH.fits'
elif grating_ID == 'YJ': 
	neon_name = calib_dir + '/kmos_ar_ne_list_yj.fits'
	oh_name = calib_dir + '/kmos_oh_spec_yj.fits'
	bpixel_flat_name = 'badpixel_flat_YJYJYJ.fits'
	flat_edge_name = 'flat_edge_YJYJYJ.fits'
	master_flat_name = 'master_flat_YJYJYJ.fits'
	xcal_name = 'xcal_YJYJYJ.fits'
	ycal_name = 'ycal_YJYJYJ.fits'
	lcal_name = 'lcal_YJYJYJ.fits'
	telluric_name = 'telluric_YJYJYJ.fits'
elif grating_ID == 'IZ': 
	neon_name = calib_dir + '/kmos_ar_ne_list_iz.fits'
	oh_name = calib_dir + '/kmos_oh_spec_iz.fits'
	bpixel_flat_name = 'badpixel_flat_IZIZIZ.fits'
	flat_edge_name = 'flat_edge_IZIZIZ.fits'
	master_flat_name = 'master_flat_IZIZIZ.fits'
	xcal_name = 'xcal_IZIZIZ.fits'
	ycal_name = 'ycal_IZIZIZ.fits'
	lcal_name = 'lcal_IZIZIZ.fits'
	telluric_name = 'telluric_IZIZIZ.fits'
elif grating_ID == 'HK': 
	neon_name = calib_dir + '/kmos_ar_ne_list_hk.fits'
	oh_name = calib_dir + '/kmos_oh_spec_hk.fits'
	bpixel_flat_name = 'badpixel_flat_HKHKHK.fits'
	flat_edge_name = 'flat_edge_HKHKHK.fits'
	master_flat_name = 'master_flat_HKHKHK.fits'
	xcal_name = 'xcal_HKHKHK.fits'
	ycal_name = 'ycal_HKHKHK.fits'
	lcal_name = 'lcal_HKHKHK.fits'
	telluric_name = 'telluric_HKHKHK.fits'
elif grating_ID == 'K': 
	neon_name = calib_dir + '/kmos_ar_ne_list_k.fits'
	oh_name = calib_dir + '/kmos_oh_spec_k.fits'
	bpixel_flat_name = 'badpixel_flat_KKK.fits'
	flat_edge_name = 'flat_edge_KKK.fits'
	master_flat_name = 'master_flat_KKK.fits'
	xcal_name = 'xcal_KKK.fits'
	ycal_name = 'ycal_KKK.fits'
	lcal_name = 'lcal_KKK.fits'
	telluric_name = 'telluric_KKK.fits'
else:
	raise ValueError('Did not recognise grating ID')


waveband_name = calib_dir + '/kmos_wave_band.fits'
reftable_name = calib_dir + '/kmos_wave_ref_table.fits'
spec_lookup_name = calib_dir + '/kmos_spec_type.fits'

#Happy that directory variables are there and that they are pointing to the right place 
#Now change directory to the Calibrations folder and construct the dark sof 
#########################################
#STEP 1: DARK FRAMES AND BAD PIXEL MAP
#Created sof name:	dark.sof
##########################################
##Write out the types of the raw files to the log 
#os.system('dfits $KMOS_RAW/*.fits | fitsort dpr.type >> log.txt')#
##reload the log file and only keep those with DARK type 
#Table = np.loadtxt('log.txt', dtype='str')
#tupe = zip(Table[:,0], Table[:,1])
#zip_names = []#
##Loop around searching for the keyword DARK
#for entry in tupe:
#    for i in range(len(entry)):
#        if entry[i].find('DARK') != -1:
#            zip_names.append(entry)#
##Check for existence of dark.sof 
#if os.path.isfile('dark.sof'):
#	os.system('rm dark.sof')#
#with open('dark.sof', 'a') as f:
#	for entry in zip_names:
#		f.write('%s\t%s\n' % (entry[0], entry[1]))
#os.system('rm log.txt')#
#print '[INFO]: Computing dark frames'#
#os.system('esorex kmos_dark dark.sof')#
##Generates the two files master_dark.fits and badpixel_dark.fits that are used 
##By the flatfield recipe#
##We want to grow the badpixel mask, to account for the 
##pixels surrounding the bad ones which are also saturated #
#print '[INFO]: Extending bad pixel mask'#
#pipe_methods.badPixelextend(badpmap='badpixel_dark.fits')#
####################################################
##STEP 2: FLATFIELDING
##
##Required .sof filename: flatfield.sof
##in the sof require 3 flat,off dpr.type file s
##and 18 flat,lamp files which correspond to 
##3 files for each of the 6 different rotator angles
#####################################################
##Play the same game with generating the .sof file
##Write out the types of the raw files to the log 
#os.system('dfits $KMOS_RAW/*.fits | fitsort dpr.type >> log.txt')#
##reload the log file and only keep those with DARK type 
#Table = np.loadtxt('log.txt', dtype='str')
#tupe = zip(Table[:,0], Table[:,1])
#zip_names = []#
##Loop around searching for the keyword DARK
#for entry in tupe:
#    for i in range(len(entry)):
#        if entry[i].find('FLAT,OFF') != -1:
#			new_entry = (entry[0], str(entry[i].replace('FLAT,OFF', 'FLAT_OFF')))
#			zip_names.append(new_entry)#
#        elif entry[i].find('FLAT,LAMP') != -1:
#			new_entry = (entry[0], str(entry[i].replace('FLAT,LAMP', 'FLAT_ON')))
#			zip_names.append(new_entry)#
##Check for existence of flat.sof 
#if os.path.isfile('flat.sof'):
#	os.system('rm flat.sof')#
#with open('flat.sof', 'a') as f:
#	for entry in zip_names:
#		f.write('%s\t%s\n' % (entry[0], entry[1]))
#os.system('rm log.txt')
##Also add the badpixel_dark_added.fits file to the flat.sof
#with open('flat.sof', 'a') as f:
#	f.write('badpixel_dark_Added.fits\tBADPIXEL_DARK')##
#print '[INFO]: Computing flat fields for the different potato angles'#
#os.system('esorex kmos_flat flat.sof')#
#####################################################
##STEP 3: ARCS
##
##Created .sof name: wave.sof
##Require the products from both the flat and dark recipes 
##badpixel_dark_Added.fits, badpixel_flat_HHH.fits, 
##flat_edge_HHH.fits, ycal_HHH.fits, xcal_HHH.fits, 
##master_flat_HHH.fits
##AND 6 ARC_ON, 1 ARC_OFF and 3 static calibration 
##files listed in the document 
####################################################
##Play the same game with generating the .sof file
##Write out the types of the raw files to the log 
#os.system('dfits $KMOS_RAW/*.fits | fitsort dpr.type >> log.txt')#
##reload the log file and only keep those with DARK type 
#Table = np.loadtxt('log.txt', dtype='str')
#tupe = zip(Table[:,0], Table[:,1])
#zip_names = []#
##Loop around searching for the keyword DARK
#for entry in tupe:
#    for i in range(len(entry)):
#        if entry[i].find('WAVE,OFF') != -1:
#			new_entry = (entry[0], str(entry[i].replace('WAVE,OFF', 'ARC_OFF')))
#			zip_names.append(new_entry)#
#        elif entry[i].find('WAVE,LAMP') != -1:
#			new_entry = (entry[0], str(entry[i].replace('WAVE,LAMP', 'ARC_ON')))
#			zip_names.append(new_entry)#
##Check for existence of flat.sof 
#if os.path.isfile('wave.sof'):
#	os.system('rm wave.sof')#
#with open('wave.sof', 'a') as f:
#	for entry in zip_names:
#		f.write('%s\t%s\n' % (entry[0], entry[1]))
#os.system('rm log.txt')#
##Also add the calibration products to the wave.sof
##This is the first time the grating dependence comes into play 
#with open('wave.sof', 'a') as f:
#	f.write('%s\tBADPIXEL_FLAT\n' % bpixel_flat_name)
#	f.write('%s\tFLAT_EDGE\n' % flat_edge_name)
#	f.write('%s\tMASTER_FLAT\n' % master_flat_name)
#	f.write('%s\tXCAL\n' % xcal_name)
#	f.write('%s\tYCAL\n' % ycal_name)
#	f.write('%s\tARC_LIST\n' % neon_name)
#	f.write('%s\tWAVE_BAND\n' % waveband_name)
#	f.write('%s\tREF_LINES\n' % reftable_name)
#	##
#print '[INFO]: Computing wavelength calibration'#
#os.system('esorex kmos_wave_cal wave.sof')#
######################################################
##STEP 4: ILLUMINATION CORRECTION (Optional)
##
##Required .sof filename: illum_cor.sof
#######################################################
##Leaving out the illumination correction at the moment 
##Not an important part of the pipeline
##print 'Attempting to create illumination correction'#
##os.system('esorex kmos_illumination %s' % illum_cor_sof)#
########################################################
##STEP 5: STANDARD STARS
##
##Created .sof file star.sof
#########################################################
#print '[INFO]: Calibrating Standard Stars'#
##Check the raw frames directory for the standard star files 
#os.system('dfits $KMOS_RAW/*.fits | fitsort dpr.type >> log.txt')#
##reload the log file and only keep those with OBJECT,SKY,STD,FLUX type 
#Table = np.loadtxt('log.txt', dtype='str')
#tupe = zip(Table[:,0], Table[:,1])
#zip_names = []#
##Loop around searching for the keyword OBJECT,SKY,STD,FLUX
#for entry in tupe:
#    for i in range(len(entry)):
#        if entry[i].find('OBJECT,SKY,STD,FLUX') != -1:
#			new_entry = (entry[0], str(entry[i].replace('OBJECT,SKY,STD,FLUX', 'STD')))
#			zip_names.append(new_entry)##
##Check for existence of star.sof 
#if os.path.isfile('star.sof'):
#	os.system('rm star.sof')#
#with open('star.sof', 'a') as f:
#	for entry in zip_names:
#		f.write('%s\t%s\n' % (entry[0], entry[1]))
#os.system('rm log.txt')##
#with open('star.sof', 'a') as f:
#	f.write('%s\tMASTER_FLAT\n' % master_flat_name)
#	f.write('%s\tXCAL\n' % xcal_name)
#	f.write('%s\tYCAL\n' % ycal_name)
#	f.write('%s\tLCAL\n' % lcal_name)
#	f.write('%s\tWAVE_BAND\n' % waveband_name)
#	f.write('%s\tSPEC_TYPE_LOOKUP\n' % spec_lookup_name)
#	##
#print '[INFO]: Computing Standard Star calibration'
#os.system('esorex kmos_std_star star.sof')#
##At this stage we've reached the end of the calibration. Ready for correcting  
##the readout column noise and performing the sub-pixel alignment
##########################################################
##STEP 6: READOUT COLUMN NOISE
##
##Required - list of object and skyfiles piped from dfits 
##called object_names.txt
###########################################################
#print '[INFO]: Correcting for readout column bias'
#if os.path.isfile('$KMOS_RAW/*Corrected.fits'):
#	os.system('rm $KMOS_RAW/*Corrected.fits')
#pipe_methods.applyCorrection('object_names.txt', 'badpixel_dark_Added.fits', lcal_name)
##Create corrected directory within the $KMOS_RAW directort
#Corrected_dir_name = raw_dir + '/Corrected'
#if os.path.isdir(Corrected_dir_name):
#	os.system('rm -rf %s' % Corrected_dir_name)
#os.system('mkdir %s' % Corrected_dir_name)
##Move all of the newly created Corrected files into this directory 
#os.system('mv $KMOS_RAW/*Corrected.fits %s' % Corrected_dir_name)
##Also copy in the sky frames
#for entry in sky_list:
#	#Doing this for ease of corrected_object file creation
#	os.system('cp $KMOS_RAW/%s %s' % (entry, Corrected_dir_name))
##Create the corrected object names file 
##Now simple to create the corrected_object_names.txt file
#if os.path.isfile('corrected_object_names.txt'):
#	os.system('rm corrected_object_names.txt')
#os.system('dfits $KMOS_RAW/Corrected/*.fits | fitsort ocs.arm1.type >> corrected_object_names.txt')
##NOTE THIS WORKS WELL SO LONG AS THE DIRECTORY STRUCTURE ISN'T TOO LONG
##Now also apply the subtraction to populate the /Corrected directory with these files
#pipe_methods.applySubtraction('corrected_object_names.txt')
########################################################
##STEP 7: SHIFT AND ALIGN#
##Required - list of corrected files from above 
##called corrected_object_names.txt 
##Be sure before going this far - this step will take ages
##There are some specifiable things in this function. 
##Later on I'll look at adding some options and variables 
##to this part of the code using an argument parser
##########################################################
#print '[INFO]: Shifting and alligning all images - this will take some time'
##remove files which could cause issues should they exist
#if os.path.isfile('temp_shift.fits'):
#	os.system('rm temp*')
#if os.path.isfile('$KMOS_RAW/*Shifted.fits'):
#	os.system('rm $KMOS_RAW/*Shifted.fits')
##Apply the shift 
#pipe_methods.applyShiftAllExtensionsMin(fileList = 'corrected_object_names.txt', badpmap='badpixel_dark_Added.fits', \
#  	 vertSegments=1, horSegments=1, interp_type='spline3')
##This has created the shifted files in the /Corrected directory
##Follow the exact same steps as above to create the shifted directory 
##And move the create files to here, then apply the subtraction 
#Shifted_dir_name = raw_dir + '/Shifted'
#if os.path.isdir(Shifted_dir_name):
#	os.system('rm -rf %s' % Shifted_dir_name)
#os.system('mkdir %s' % Shifted_dir_name)
##Move all of the newly created Shifted files into this directory 
#os.system('mv %s/*Shifted.fits %s' % (Corrected_dir_name, Shifted_dir_name))#
##Also copy in the sky frames
#for entry in sky_list:
#	#Doing this for ease of Shifted_object file creation
#	os.system('cp $KMOS_RAW/%s %s' % (entry, Shifted_dir_name))
##Create the Shifted object names file 
##Now simple to create the Shifted_object_names.txt file
#if os.path.isfile('shifted_object_names.txt'):
#	os.system('rm shifted_object_names.txt')#
#os.system('dfits $KMOS_RAW/Shifted/*.fits | fitsort ocs.arm1.type >> shifted_object_names.txt')
##NOTE THIS WORKS WELL SO LONG AS THE DIRECTORY STRUCTURE ISN'T TOO LONG
##Now also apply the subtraction to populate the /shifted directory with these files
#pipe_methods.applySubtraction('shifted_object_names.txt')#
##Make a plot of the shift results 
#pipe_methods.shiftPlot('Shift_Coords.txt')
##########################################################
###All lists of files now made and objects have been shifted. 
###Need to think of a way to get the output from frameCheck into 
###the science directory #
######################################################################
###SCIENCE REDUCTION: Use everything we've made and create the data cubes
###Required: reduction sof sci_reduc.sof
###This file will contain all the corrected and shifted observations 
###produced by the above two manual recipes. Make sure they have the correct 
###names in the .sof file.
###Also will move directory into science folder 
###Must have a folder called Science_Output in the directory above 
###the callibration directory
########################################################################
###Wipe out the science directory before beginning, this is meant as a complete reset 
#os.system('rm -rf $KMOS_SCIENCE/*')
#print '[INFO]: Creating %s Grating sci_reduc.sof file' % grating_ID
#if os.path.isfile('sci_reduc.sof'):
#	os.system('rm sci_reduc.sof')#
#with open('sci_reduc.sof', 'a') as f:
#	f.write('%s\tXCAL\n' % xcal_name)
#	f.write('%s\tYCAL\n' % ycal_name)
#	f.write('%s\tLCAL\n' % lcal_name)
#	f.write('%s\tMASTER_FLAT\n' % master_flat_name)
#	f.write('%s\tTELLURIC\n' % telluric_name)
#	f.write('%s\tWAVE_BAND\n' % waveband_name)
#	f.write('%s\tOH_SPEC\n' % oh_name)
##Grating sof now created - need to create the sky datacubes in the science directory 
##Have access to the names of the sky files via sky_list 
#sky_frame_first = raw_dir + '/' + sky_list[0]
#sky_frame_last = raw_dir + '/' + sky_list[-1]
##Create a sky-reconstruct sof file 
#if os.path.isfile('sky_reconstruct.sof'):
#	os.system('rm sky_reconstruct.sof')
#with open('sky_reconstruct.sof', 'a') as f:
#	f.write('%s\tSCIENCE\n' % sky_frame_first)
#	f.write('%s\tSCIENCE\n' % sky_frame_last)
#	f.write('%s\tXCAL\n' % xcal_name)
#	f.write('%s\tYCAL\n' % ycal_name)
#	f.write('%s\tLCAL\n' % lcal_name)
#	f.write('%s\tMASTER_FLAT\n' % master_flat_name)
#	f.write('%s\tWAVE_BAND\n' % waveband_name)
#	f.write('%s\tOH_SPEC\n' % oh_name)
#print '[INFO]: Reconstructing Both sky cubes'
##Execute the esorex reconstruct recipe to create this raw sky multi cube 
#os.system('esorex --output-dir=$KMOS_SCIENCE kmos_sci_red --no_subtract sky_reconstruct.sof')
##The outputs from this are two file called sci_reconstructed_***.fits in the science directory
##os.system('mv $KMOS_SCIENCE/cube_science.fits $KMOS_SCIENCE/cube_science_1.fits')
##remove the sky_reconstruct file 
#os.system('rm sky_reconstruct.sof')
##Now have two reconstructed files in the science directory which we combine to give the individual IFU cubes
##first set the names of the two files 
#sky_first_combine_name = science_dir + '/sci_reconstructed_' + sky_list[0] 
#sky_last_combine_name = science_dir + '/sci_reconstructed_' + sky_list[-1]
##check for existence of combine.sof
#if os.path.isfile('sky_combine.sof'):
#	os.system('rm sky_combine.sof')
#with open('sky_combine.sof', 'a') as f:
#	f.write('%s\n' % sky_first_combine_name)
#	f.write('%s' % sky_last_combine_name)
#print '[INFO]: Combining sky cubes'
##Execute the KMOS recipe for combining the multi-cubes - save the output to the science directory
#os.system('esorex --output-dir=$KMOS_SCIENCE kmo_combine --method="header" --edge_nan=TRUE sky_combine.sof')
##Everything Now in place for the multiExtractSpec method - which checks the performance of the individual frames 
##This should deposit everything in the science directory
pipe_methods.multiExtractSpec(sci_dir=science_dir, frameNames='shifted_object_names.txt', tracked_name='C_STARS_7656')

