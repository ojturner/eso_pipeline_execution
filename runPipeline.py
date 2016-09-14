#########################################################
# SCRIPT TO EXECUTE THE KMOS PIPELINE AUTOMATICALLY
# BUILDING THE .SOF FILES EN ROUTE. MUST BE EXECUTED FROM
# WITHIN THE $KMOS_DYN_CALIBRATIONS DIRECTORY
########################################################

# import the relevant modules
import os
import sys
import numpy as np
import fnmatch


# add the class file to the PYTHONPATH
sys.path.append('/home/oturner/disk1/turner/' +
                'PhD/KMOS/Analysis_Pipeline/Python_code/Class')

from pipelineClass import pipelineOps


def run_pipeline(tracked_list,
                 quick=False,
                 shift=False,
                 telluric_method='standard',
                 illum_cor='None',
                 pix_scale=0.2):

    """
    Def:
    run the reduction pipeline, including all of the major steps.
    also includes additional processing, bad pixel mask growth,
    column bias correction and optional pre-science reduction shifting
    of the two dimensional data, controlled with the shift keyword.

    Note this function requires the following directory structure
    and associated environment variables pointing to them:

    $KMOS_RAW
              - contains all raw files. This includes the raw OSO object
                and sky frames as well as all the raw calibration frames
                including darks, flats, wave cal and standard star

    $KMOS_DYN_CALIBRATIONS
              - empty directory to house the calibration and .sof files
                produced as the pipeline works

    $KMOS_SCIENCE
              - empty directory to house the science ouptut products
                from running the pipeline

    Input:
            tracked_star - name of one of the stars tracked throughout the 
                           night by one of the ifus. This is used to look 
                           at the seeing profile and create the optimal 
                           integrated spectrum extraction 

            quick - set to true to go straight through to the
                    multi_extract_spec method at the end of the file

            shift - set to true to include the 2D data shift prior to
                    feeding into the science reduction pipeline. Currently
                    there doesn't seem to be a significant advantage
                    to this over the Davies algorithm

            telluric_method - 'flat' - use only unit conversion, don't correct
                                       the shape
                              'poly' - Use polynomial fit to the full telluric
                                        spectrum to correct the shape a bit
                              'standard' - Use full telluric spectrum for corr

            illum_cor - If illumination correction files are available
                        specify either flat_sky or flat_sky_flat as the
                        method. The Final science products are corrected
                        for the illumination gradients. 

    Output:
            reduced data cubes for all of the raw object files
    """
    # Create an instance of the class
    pipe_methods = pipelineOps()

    # Writing robust pipeline that will generate the .sof files at each stage.
    # Put each of the different wavebands in a separate folder and
    # execute this recipe
    # still need to figure out how to navigate the directory
    # structure from the raw files
    # Probably need to have the environment variables set each time.
    # Check to see if log.txt file exists

    if os.path.isfile('log.txt'):

        os.system('rm log.txt')

    # Check to see if the directory variables have been set in the first place
    try:

        # If any are False - need to exit. Variables aren't set
        if os.environ['KMOS_SCIENCE'] \
                and os.environ['KMOS_DYN_CALIBRATIONS']\
                and os.environ['KMOS_RAW']:

            print '[INFO]: Directory Variables are set'

    except KeyError:

        print '[INFO]: The Directory variables are not set'

        raise

    # Check to see if the user has changed the $KMOS_SCIENCE,
    # $KMOS_RAW and $KMOS_DYN_CALIBRATIONS environment variables
    if raw_input('Are you happy with the values of the directory variables?' +
                 '\nKMOS_SCIENCE: %s' % os.environ['KMOS_SCIENCE'] +
                 '\nKMOS_RAW: %s' % os.environ['KMOS_RAW'] +
                 '\nKMOS_DYN_CALIBRATIONS: %s\n(Y/N): '
                 % os.environ['KMOS_DYN_CALIBRATIONS']) == 'Y':

        print '[INFO]: Continuing Calibrations and Reduction'

    else:

        raise ValueError('Did Not Change directory variables... Quitting')

    # Get the values of the directory variables
    calib_dir = os.environ['KMOS_CALIB']

    raw_dir = os.environ['KMOS_RAW']

    science_dir = os.environ['KMOS_SCIENCE']

    dyn_dir = os.environ['KMOS_DYN_CALIBRATIONS']

    if not quick:

        # Check for the existence of any Subtracted
        # files in the raw directory and delete if they exist
        for file in os.listdir(os.environ['KMOS_RAW']):

            if fnmatch.fnmatch(file, '*Subtracted*.fits'):

                os.system('rm $KMOS_RAW/*Subtracted*.fits')

        # First find the grating ID of one of the object,
        # sky files in the $KMOS_RAW directory
        # This will be used later on in the
        # script to construct the correct .sof files
        # Write out the types of the raw files to the log
        os.system('dfits $KMOS_RAW/*.fits |'
                  + ' fitsort dpr.type ins.grat1.id >> log.txt')

        # reload the log file and only keep those with OBJECT,SKY type
        Table = np.loadtxt('log.txt', dtype='str')

        tupe = zip(Table[:, 0], Table[:, 1], Table[:, 2])

        zip_names = []

        obj_names = []

        # Loop around searching for the keyword OBJECT,SKY

        for entry in tupe:

            for i in range(len(entry)):

                if entry[i].find('OBJECT,SKY') != -1 \
                        and entry[i].find('OBJECT,SKY,STD,FLUX') == -1:

                    zip_names.append(entry)

                    # Append only the names of these
                    # files to the obj_names vector
                    obj_names.append(entry[0])

        # create a new temporary directory for the object,sky files
        if os.path.isdir('temp/'):

            os.system('rm -rf temp/')

        os.system('mkdir temp/')

        # make sure that only the raw images make
        # it into the object_names.txt file
        corr_list = []

        # loop around the obj_names and check for the word 'Corrected'
        for item in obj_names:

            if item.find('Corrected') != -1:

                corr_list.append(item)

        # Now remove from the obj_names array any of the flagged files
        for item in corr_list:

            obj_names.remove(item)

        # obj_names should now contain only the raw object,sky files
        # copy the object,sky files into this temporary directory
        for item in obj_names:

            os.system('cp %s temp/' % item)

        # Now simply pipe the output from dfits to the object_names.txt file
        if os.path.isfile('temp.txt'):

            os.system('rm temp.txt')

        os.system('dfits temp/*.fits | fitsort ocs.arm1.type >> temp.txt')

        # Stupidly this also leaves the directory in the output
        # Need to reload the temporary file,
        # change the occurence of to the $KMOS_RAW directory
        object_names_Table = np.loadtxt('temp.txt', dtype='str')

        new_list = []

        sky_list = []

        for i in range(len(object_names_Table)):

            temp_new_list = []

            for j in range(len(object_names_Table[i])):

                newString = object_names_Table[i][j].replace('temp/',
                                                             raw_dir + '/')

                temp_new_list.append(newString)

                if object_names_Table[i][j] == 'S':

                    sky_list.append(object_names_Table[i][0])

            new_list.append(temp_new_list)

        for i in range(len(sky_list)):

            sky_list[i] = sky_list[i][5:]

        print '[INFO]: This is the list of sky_frames: %s' % sky_list

        # Finally ready for writing to the object_names.txt file
        if os.path.isfile('object_names.txt'):

            os.system('rm object_names.txt')

        with open('object_names.txt', 'a') as f:

            for entry in new_list:

                f.write('%s\t%s\n' % (entry[0], entry[1]))
        # clean up by removing the temporary directory
        os.system('rm -rf temp/')

        os.system('rm temp.txt')

        grating_ID = zip_names[0][2]

        os.system('rm log.txt')

        # place files in a sub-directory within the raw folder
        sub_dir_name = raw_dir + '/Subtracted'

        if os.path.isdir(sub_dir_name):

            os.system('rm -rf %s' % sub_dir_name)

        os.system('mkdir %s' % sub_dir_name)

        # Apply the subtraction method
        pipe_methods.applySubtraction('object_names.txt')

        # Move the created files into the subtracted directory
        os.system('mv $KMOS_RAW/*Subtracted.fits %s' % sub_dir_name)

        # After defining the grating ID, set the values of the created files
        # conditionally, depending on which grating is being used

        if grating_ID == 'H':
            neon_name = calib_dir + '/kmos_ar_ne_list_h.fits'
            oh_name = calib_dir + '/kmos_oh_spec_h.fits'
            bpixel_flat_name = 'BADPIXEL_FLAT_HHH.fits'
            flat_edge_name = 'FLAT_EDGE_HHH.fits'
            master_flat_name = 'MASTER_FLAT_HHH.fits'
            xcal_name = 'XCAL_HHH.fits'
            ycal_name = 'YCAL_HHH.fits'
            lcal_name = 'LCAL_HHH.fits'
            telluric_name = 'TELLURIC_HHH.fits'
            atmos_name = calib_dir + '/kmos_atmos_h.fits'

        elif grating_ID == 'YJ':
            neon_name = calib_dir + '/kmos_ar_ne_list_yj.fits'
            oh_name = calib_dir + '/kmos_oh_spec_yj.fits'
            bpixel_flat_name = 'BADPIXEL_FLAT_YJYJYJ.fits'
            flat_edge_name = 'FLAT_EDGE_YJYJYJ.fits'
            master_flat_name = 'MASTER_FLAT_YJYJYJ.fits'
            xcal_name = 'XCAL_YJYJYJ.fits'
            ycal_name = 'YCAL_YJYJYJ.fits'
            lcal_name = 'LCAL_YJYJYJ.fits'
            telluric_name = 'TELLURIC_YJYJYJ.fits'
            atmos_name = calib_dir + '/kmos_atmos_yj.fits'

        elif grating_ID == 'IZ':
            neon_name = calib_dir + '/kmos_ar_ne_list_iz.fits'
            oh_name = calib_dir + '/kmos_oh_spec_iz.fits'
            bpixel_flat_name = 'BADPIXEL_FLAT_IZIZIZ.fits'
            flat_edge_name = 'FLAT_EDGE_IZIZIZ.fits'
            master_flat_name = 'MASTER_FLAT_IZIZIZ.fits'
            xcal_name = 'XCAL_IZIZIZ.fits'
            ycal_name = 'YCAL_IZIZIZ.fits'
            lcal_name = 'LCAL_IZIZIZ.fits'
            telluric_name = 'TELLURIC_IZIZIZ.fits'
            atmos_name = calib_dir + '/kmos_atmos_iz.fits'

        elif grating_ID == 'HK':
            neon_name = calib_dir + '/kmos_ar_ne_list_hk.fits'
            oh_name = calib_dir + '/kmos_oh_spec_hk.fits'
            bpixel_flat_name = 'BADPIXEL_FLAT_HKHKHK.fits'
            flat_edge_name = 'FLAT_EDGE_HKHKHK.fits'
            master_flat_name = 'MASTER_FLAT_HKHKHK.fits'
            xcal_name = 'XCAL_HKHKHK.fits'
            ycal_name = 'YCAL_HKHKHK.fits'
            lcal_name = 'LCAL_HKHKHK.fits'
            telluric_name = 'TELLURIC_HKHKHK.fits'
            atmos_name = calib_dir + '/kmos_atmos_hk.fits'

        elif grating_ID == 'K':
            neon_name = calib_dir + '/kmos_ar_ne_list_k.fits'
            oh_name = calib_dir + '/kmos_oh_spec_k.fits'
            bpixel_flat_name = 'BADPIXEL_FLAT_KKK.fits'
            flat_edge_name = 'FLAT_EDGE_KKK.fits'
            master_flat_name = 'MASTER_FLAT_KKK.fits'
            xcal_name = 'XCAL_KKK.fits'
            ycal_name = 'YCAL_KKK.fits'
            lcal_name = 'LCAL_KKK.fits'
            telluric_name = 'TELLURIC_KKK.fits'
            atmos_name = calib_dir + '/kmos_atmos_k.fits'

        else:

            raise ValueError('Did not recognise grating ID')

        waveband_name = calib_dir + '/kmos_wave_band.fits'

        reftable_name = calib_dir + '/kmos_wave_ref_table.fits'

        spec_lookup_name = calib_dir + '/kmos_spec_type.fits'

        master_dark_name = 'MASTER_DARK.fits'

        # Happy that directory variables are
        # there and that they are pointing to the right place
        # Now change directory to the Calibrations
        # folder and construct the dark sof

        # STEP 1: DARK FRAMES AND BAD PIXEL MAP
        # Created sof name:  dark.sof
        # Write out the types of the raw files to the log
        os.system('dfits $KMOS_RAW/*.fits | fitsort dpr.type >> log.txt')

        # reload the log file and only keep those with DARK type
        Table = np.loadtxt('log.txt', dtype='str')

        tupe = zip(Table[:, 0], Table[:, 1])

        zip_names = []

        # Loop around searching for the keyword DARK
        for entry in tupe:

            for i in range(len(entry)):

                if entry[i].find('DARK') != -1:

                    zip_names.append(entry)

        # Check for existence of dark.sof
        if os.path.isfile('dark.sof'):

            os.system('rm dark.sof')

        with open('dark.sof', 'a') as f:

            for entry in zip_names:

                f.write('%s\t%s\n' % (entry[0], entry[1]))

        os.system('rm log.txt')

        print '[INFO]: Computing dark frames'

        os.system('esorex kmos_dark dark.sof')

        # Generates the two files master_dark.fits
        # and BADPIXEL_DARK.fits that are used
        # By the flatfield recipe#
        # We want to grow the badpixel mask, to account for the
        # pixels surrounding the bad ones which are also saturated
        print '[INFO]: Extending bad pixel mask'

        pipe_methods.badPixelextend(badpmap='BADPIXEL_DARK.fits')

        ###################################################
        # STEP 2: FLATFIELDING
        # Required .sof filename: flatfield.sof
        # in the sof require 3 flat,off dpr.type file s
        # and 18 flat,lamp files which correspond to
        # 3 files for each of the 6 different rotator angles
        #####################################################

        # Play the same game with generating the .sof file
        # Write out the types of the raw files to the log

        os.system('dfits $KMOS_RAW/*.fits | fitsort dpr.type >> log.txt')

        # reload the log file and only keep those with DARK type

        Table = np.loadtxt('log.txt', dtype='str')

        tupe = zip(Table[:, 0], Table[:, 1])

        zip_names = []

        # Loop around searching for the keyword DARK

        for entry in tupe:

            for i in range(len(entry)):

                if entry[i].find('FLAT,OFF') != -1:

                    new_entry = (entry[0],
                                 str(entry[i].replace('FLAT,OFF', 'FLAT_OFF')))

                    zip_names.append(new_entry)

                elif entry[i].find('FLAT,LAMP') != -1:

                    new_entry = (entry[0],
                                 str(entry[i].replace('FLAT,LAMP', 'FLAT_ON')))

                    zip_names.append(new_entry)

        # Check for existence of flat.sof
        if os.path.isfile('flat.sof'):

            os.system('rm flat.sof')

        with open('flat.sof', 'a') as f:

            for entry in zip_names:

                f.write('%s\t%s\n' % (entry[0], entry[1]))

        os.system('rm log.txt')

        # Also add the BADPIXEL_DARK_added.fits file to the flat.sof

        with open('flat.sof', 'a') as f:

            f.write('BADPIXEL_DARK_Added.fits\tBADPIXEL_DARK')

        print '[INFO]: Computing flat fields for the different potato angles'

        os.system('esorex kmos_flat flat.sof')

        ####################################################
        # STEP 3: ARCS
        #
        # Created .sof name: wave.sof
        # Require the products from both the flat and dark recipes
        # BADPIXEL_DARK_Added.fits, badpixel_flat_HHH.fits,
        # flat_edge_HHH.fits, ycal_HHH.fits, xcal_HHH.fits,
        # master_flat_HHH.fits
        # AND 6 ARC_ON, 1 ARC_OFF and 3 static calibration
        # files listed in the document
        #####################################################

        # Play the same game with generating the .sof file
        # Write out the types of the raw files to the log

        os.system('dfits $KMOS_RAW/*.fits | fitsort dpr.type >> log.txt')

        # reload the log file and only keep those with WAVEOFF type

        Table = np.loadtxt('log.txt', dtype='str')

        tupe = zip(Table[:, 0], Table[:, 1])

        zip_names = []

        # Loop around searching for the keyword DARK

        for entry in tupe:

            for i in range(len(entry)):

                if entry[i].find('WAVE,OFF') != -1:

                    new_entry = (entry[0],
                                 str(entry[i].replace('WAVE,OFF', 'ARC_OFF')))

                    zip_names.append(new_entry)

                elif entry[i].find('WAVE,LAMP') != -1:

                    new_entry = (entry[0],
                                 str(entry[i].replace('WAVE,LAMP', 'ARC_ON')))

                    zip_names.append(new_entry)

        # Check for existence of flat.sof

        if os.path.isfile('wave.sof'):

            os.system('rm wave.sof')

        with open('wave.sof', 'a') as f:

            for entry in zip_names:

                f.write('%s\t%s\n' % (entry[0], entry[1]))

        os.system('rm log.txt')

        # Also add the calibration products to the wave.sof
        # This is the first time the grating dependence comes into play

        with open('wave.sof', 'a') as f:
            f.write('%s\tBADPIXEL_FLAT\n' % bpixel_flat_name)
            f.write('%s\tFLAT_EDGE\n' % flat_edge_name)
            f.write('%s\tMASTER_FLAT\n' % master_flat_name)
            f.write('%s\tXCAL\n' % xcal_name)
            f.write('%s\tYCAL\n' % ycal_name)
            f.write('%s\tARC_LIST\n' % neon_name)
            f.write('%s\tWAVE_BAND\n' % waveband_name)
            f.write('%s\tREF_LINES\n' % reftable_name)

        print '[INFO]: Computing wavelength calibration'

        os.system('esorex kmos_wave_cal wave.sof')

        #####################################################
        # STEP 4: ILLUMINATION CORRECTION (Optional)
        #
        # Required .sof filename: illum_cor.sof
        ######################################################

        # it may be that sky flats were not taken

        if illum_cor == 'None':

            print 'No illumination correction frames available'

        # do the same thing with preparing the sof file

        elif illum_cor == 'flat_sky':

            print '[INFO]: Processing FLAT_SKY illumination correction'

            os.system('dfits $KMOS_RAW/*.fits | fitsort dpr.type >> log.txt')

            # reload the log file and only keep those with FLAT,SKY type

            Table = np.loadtxt('log.txt', dtype='str')

            tupe = zip(Table[:, 0], Table[:, 1])

            zip_names = []

            # Loop around searching for the keyword FLAT,SKY

            for entry in tupe:

                for i in range(len(entry)):

                    if entry[i].find('FLAT,SKY') != -1:

                        new_entry = (entry[0],
                                     str(entry[i].replace('FLAT,SKY', 'FLAT_SKY')))

                        zip_names.append(new_entry)

            # Check for existence of flat.sof

            if os.path.isfile('illum_corr.sof'):

                os.system('rm illum_corr.sof')

            with open('illum_corr.sof', 'a') as f:

                for entry in zip_names:

                    f.write('%s\t%s\n' % (entry[0], entry[1]))

            os.system('rm log.txt')

            # Also add the calibration products to the wave.sof
            # This is the first time the grating dependence comes into play

            with open('illum_corr.sof', 'a') as f:
                f.write('%s\tMASTER_DARK\n' % master_dark_name)
                f.write('%s\tFLAT_EDGE\n' % flat_edge_name)
                f.write('%s\tMASTER_FLAT\n' % master_flat_name)
                f.write('%s\tXCAL\n' % xcal_name)
                f.write('%s\tYCAL\n' % ycal_name)
                f.write('%s\tLCAL\n' % lcal_name)
                f.write('%s\tWAVE_BAND\n' % waveband_name)

            print '[INFO]: Computing illumination correction'

            os.system('esorex kmos_illumination --pix_scale=%s illum_corr.sof' % pix_scale)

        elif illum_cor == 'flat_sky_flat':

            print '[INFO]: Processing FLAT_SKY_FLAT illumination correction'

            os.system('dfits $KMOS_RAW/*.fits | fitsort dpr.type >> log.txt')

            # reload the log file and only keep those with FLAT,SKY type

            Table = np.loadtxt('log.txt', dtype='str')

            tupe = zip(Table[:, 0], Table[:, 1])

            zip_names = []

            # Loop around searching for the keyword FLAT,SKY

            for entry in tupe:

                for i in range(len(entry)):

                    if entry[i].find('FLAT,SKY,FLAT') != -1:

                        new_entry = (entry[0],
                                     str(entry[i].replace('FLAT,SKY,FLAT', 'FLAT_SKY_FLAT')))

                        zip_names.append(new_entry)

            # Check for existence of flat.sof

            if os.path.isfile('illum_corr.sof'):

                os.system('rm illum_corr.sof')

            with open('illum_corr.sof', 'a') as f:

                for entry in zip_names:

                    f.write('%s\t%s\n' % (entry[0], entry[1]))

            os.system('rm log.txt')

            # Also add the calibration products to the wave.sof
            # This is the first time the grating dependence comes into play

            with open('illum_corr.sof', 'a') as f:
                f.write('%s\tXCAL\n' % xcal_name)
                f.write('%s\tYCAL\n' % ycal_name)
                f.write('%s\tLCAL\n' % lcal_name)
                f.write('%s\tWAVE_BAND\n' % waveband_name)

            print '[INFO]: Computing illumination correction'

            os.system('esorex kmos_illumination_flat --pix_scale=%s illum_corr.sof' % pix_scale)

        else:

            raise TypeError('Specified incorrect illumination corr format')

        ########################################################
        # STEP 5: STANDARD STARS
        #
        # Created .sof file star.sof
        ##########################################################

        print '[INFO]: Calibrating Standard Stars'

        # Check the raw frames directory for the standard star files

        os.system('dfits $KMOS_RAW/*.fits | fitsort dpr.type >> log.txt')

        # reload the log file and only keep those with OBJECT,SKY,STD,FLUX type

        Table = np.loadtxt('log.txt', dtype='str')

        tupe = zip(Table[:, 0], Table[:, 1])

        zip_names = []

        # Loop around searching for the keyword OBJECT,SKY,STD,FLUX

        for entry in tupe:

            for i in range(len(entry)):

                if entry[i].find('OBJECT,SKY,STD,FLUX') != -1:

                    new_entry = (entry[0],
                                 str(entry[i].replace('OBJECT,SKY,STD,FLUX',
                                                      'STD')))
                    zip_names.append(new_entry)

        # Check for existence of star.sof

        if os.path.isfile('star.sof'):

            os.system('rm star.sof')

        with open('star.sof', 'a') as f:

            for entry in zip_names:

                f.write('%s\t%s\n' % (entry[0], entry[1]))

        os.system('rm log.txt')

        with open('star.sof', 'a') as f:
            f.write('%s\tMASTER_FLAT\n' % master_flat_name)
            f.write('%s\tXCAL\n' % xcal_name)
            f.write('%s\tYCAL\n' % ycal_name)
            f.write('%s\tLCAL\n' % lcal_name)
            f.write('%s\tWAVE_BAND\n' % waveband_name)
            f.write('%s\tSPEC_TYPE_LOOKUP\n' % spec_lookup_name)
            f.write('%s\tATMOS_MODEL\n' % atmos_name)
        # Now there are also solar spectra
        # which can be used to aid the correction.
        # This is useful in the K, H and HK bands.
        # So add these conditionally now
        # depending on the grating that is being passed through the pipeline

        if grating_ID == 'H':

            solar_name = calib_dir + '/kmos_solar_h_2400.fits'

            with open('star.sof', 'a') as f:

                f.write('%s\tSOLAR_SPEC\n' % solar_name)

        elif grating_ID == 'K':

            solar_name = calib_dir + '/kmos_solar_k_1700.fits'

            with open('star.sof', 'a') as f:

                f.write('%s\tSOLAR_SPEC\n' % solar_name)

        elif grating_ID == 'HK':

            solar_name = calib_dir + '/kmos_solar_hk_1100.fits'

            with open('star.sof', 'a') as f:

                f.write('%s\tSOLAR_SPEC\n' % solar_name)

        # star.sof now ready compute the telluric correction

        print '[INFO]: Computing Standard Star calibration'

        os.system('esorex kmos_std_star star.sof')

        # now correct the telluric spectrum with the polynomial fit
        # or with the flat telluric spectrum

        if telluric_method == 'flat':

            print '[INFO:] Using flat telluric spectrum'

            pipe_methods.telluric_correct_flat(grating_ID, dyn_dir)

        elif telluric_method == 'poly':

            print '[INFO:] Using polynomial fit for telluric'

            pipe_methods.telluric_correct(grating_ID, dyn_dir)

        elif telluric_method == 'standard':

            print '[INFO:] Using full spectrum for telluric'

        else:

            raise TypeError('Invalid telluric method specified')

        # At this stage we've reached the
        # end of the calibration. Ready for correcting
        # the readout column noise and performing the sub-pixel alignment

        #########################################################
        # STEP 6: READOUT COLUMN NOISE
        #
        # Required - list of object and skyfiles piped from dfits
        # called object_names.txt
        #############################################################

        print '[INFO]: Correcting for readout column bias'

        if os.path.isfile('$KMOS_RAW/*Corrected.fits'):

            os.system('rm $KMOS_RAW/*Corrected.fits')

        pipe_methods.applyCorrection('object_names.txt',
                                     'BADPIXEL_DARK_Added.fits',
                                     lcal_name)

        # Create corrected directory within the $KMOS_RAW directort
        Corrected_dir_name = raw_dir + '/Corrected'

        if os.path.isdir(Corrected_dir_name):

            os.system('rm -rf %s' % Corrected_dir_name)

        os.system('mkdir %s' % Corrected_dir_name)

        # Move all of the newly created Corrected files into this directory

        os.system('mv $KMOS_RAW/*Corrected.fits %s' % Corrected_dir_name)

        # Also copy in the sky frames

        for entry in sky_list:

            # Doing this for ease of corrected_object file creation

            os.system('cp $KMOS_RAW/%s %s' % (entry, Corrected_dir_name))

        # Create the corrected object names file
        # Now simple to create the corrected_object_names.txt file

        if os.path.isfile('corrected_object_names.txt'):

            os.system('rm corrected_object_names.txt')

        # this is a bit hacky - but have to do it because the file names
        # are becoming too long and dfits cuts them off. Again, this requires
        # that we are in the Calibrations directory to run this script.
        # otherwise it will not work.

        os.system('dfits ../raw_frames/Corrected/*.fits |'
                  + ' fitsort ocs.arm1.type >> corrected_object_names.txt')

        # now open the corrected_object_names.txt filename and fill in
        # with the full path to the files

        f = open('corrected_object_names.txt', 'r')

        filedata = f.read()

        f.close()

        newdata = filedata.replace('../raw_frames','%s' % raw_dir)

        f = open('corrected_object_names.txt', 'w')

        f.write(newdata)

        f.close()

        # everything should be well

        # Now also apply the subtraction to populate
        # the /Corrected directory with these files

        pipe_methods.applySubtraction('corrected_object_names.txt')

        #####################################################
        # STEP 7: SHIFT AND ALIGN
        # Required - list of corrected files from above
        # called corrected_object_names.txt
        # Be sure before going this far - this step will take ages
        # There are some specifiable things in this function.
        # Later on I'll look at adding some options and variables
        # to this part of the code using an argument parser
        ########################################################

        # conditional execution of this section
        # dependent upon the shift keyword

        if shift:

            print '[INFO]: Shifting and aligning' \
                  + ' all images - this will take some time'

            # remove files which could cause issues should they exist

            if os.path.isfile('temp_shift.fits'):

                os.system('rm temp*')

            if os.path.isfile('$KMOS_RAW/*Shifted.fits'):

                os.system('rm $KMOS_RAW/*Shifted.fits')

            # FUNCTION TO DO ADDITIONAL BAD PIXEL/COSMIC RAY REMOVAL
            # Apply the shift#

            pipe_methods.applyShiftAllExtensionsMin(fileList='corrected_object'
                                                             + '_names.txt',
                                                    badpmap='BADPIXEL_DARK'
                                                            + '_Added.fits',
                                                    vertSegments=1,
                                                    horSegments=1,
                                                    interp_type='spline3')

            # This has created the shifted files in the /Corrected directory
            # Follow the exact same steps as
            # above to create the shifted directory
            # And move the create files to here, then apply the subtraction

            Shifted_dir_name = raw_dir + '/Shifted'

            if os.path.isdir(Shifted_dir_name):

                os.system('rm -rf %s' % Shifted_dir_name)

            os.system('mkdir %s' % Shifted_dir_name)

            # Move all of the newly created Shifted files into this directory
            os.system('mv %s/*Shifted.fits %s' % (Corrected_dir_name,
                                                  Shifted_dir_name))
            # Also copy in the sky frames

            for entry in sky_list:

                # Doing this for ease of Shifted_object file creation

                os.system('cp $KMOS_RAW/%s %s' % (entry, Shifted_dir_name))

            # Create the Shifted object names file
            # Now simple to create the Shifted_object_names.txt file

            if os.path.isfile('shifted_object_names.txt'):

                os.system('rm shifted_object_names.txt')

            os.system('dfits ../raw_frames/Shifted/*.fits |'
                      + ' fitsort ocs.arm1.type >> shifted_object_names.txt')

            # now open the corrected_object_names.txt filename and fill in
            # with the full path to the files

            f = open('shifted_object_names.txt', 'r')

            filedata = f.read()

            f.close()

            newdata = filedata.replace('../raw_frames','%s' % raw_dir)

            f = open('shifted_object_names.txt', 'w')

            f.write(newdata)

            f.close()

            # Now also apply the subtraction to
            # populate the /shifted directory with these files

            pipe_methods.applySubtraction('shifted_object_names.txt')

            # Make a plot of the shift results

            pipe_methods.shiftPlot('Shift_Coords.txt')

        ########################################################
        # All lists of files now made and objects have been shifted
        ####################################################################

        # SCIENCE REDUCTION: Use everything we've
        # made and create the data cubes
        # Required: reduction sof sci_reduc.sof
        # This file will contain all the corrected and shifted observations
        # produced by the above two manual recipes.
        # Make sure they have the correct
        # names in the .sof file.
        # Also will move directory into science folder
        # Must have a folder called Science_Output in the directory above
        # the callibration directory
        ######################################################################

        print '[INFO]: Creating %s Grating sci_reduc.sof file' % grating_ID

        if os.path.isfile('sci_reduc.sof'):

            os.system('rm sci_reduc.sof')

        with open('sci_reduc.sof', 'a') as f:

            f.write('%s\tXCAL\n' % xcal_name)
            f.write('%s\tYCAL\n' % ycal_name)
            f.write('%s\tLCAL\n' % lcal_name)
            f.write('%s\tMASTER_FLAT\n' % master_flat_name)
            f.write('%s\tTELLURIC\n' % telluric_name)
            f.write('%s\tWAVE_BAND\n' % waveband_name)
            f.write('%s\tOH_SPEC\n' % oh_name)

        # Grating sof now created - need to
        # create the sky datacubes in the science directory
        # Have access to the names of the sky files via sky_list

        sky_frame_first = raw_dir + '/' + sky_list[0]

        sky_frame_last = raw_dir + '/' + sky_list[-1]

        # Create a sky-reconstruct sof file
        if os.path.isfile('sky_reconstruct.sof'):

            os.system('rm sky_reconstruct.sof')

        with open('sky_reconstruct.sof', 'a') as f:

            f.write('%s\tSCIENCE\n' % sky_frame_first)
            f.write('%s\tSCIENCE\n' % sky_frame_last)
            f.write('%s\tXCAL\n' % xcal_name)
            f.write('%s\tYCAL\n' % ycal_name)
            f.write('%s\tLCAL\n' % lcal_name)
            f.write('%s\tMASTER_FLAT\n' % master_flat_name)
            f.write('%s\tWAVE_BAND\n' % waveband_name)
            f.write('%s\tOH_SPEC\n' % oh_name)

        print '[INFO]: Reconstructing Both sky cubes'

        # Execute the esorex reconstruct recipe
        # to create this raw sky multi cube
        os.system('esorex --output-dir=$KMOS_SCIENCE'
                  + ' kmos_sci_red --no_subtract=TRUE --oscan=FALSE'
                  + ' sky_reconstruct.sof')

        # The outputs from this are two file called
        # sci_reconstructed_***.fits in the science directory

        # remove the sky_reconstruct file
        os.system('rm sky_reconstruct.sof')

        # Now have two reconstructed files in the science
        # directory which we combine to give the individual IFU cubes
        # first set the names of the two files
        sky_first_combine_name = science_dir + \
            '/SCI_RECONSTRUCTED_' + sky_list[0]

        sky_last_combine_name = science_dir + \
            '/SCI_RECONSTRUCTED_' + sky_list[-1]

        # check for existence of combine.sof
        if os.path.isfile('sky_combine.sof'):

            os.system('rm sky_combine.sof')

        with open('sky_combine.sof', 'a') as f:

            f.write('%s\tSCI_RECONSTRUCTED\n' % sky_first_combine_name)
            f.write('%s\tSCI_RECONSTRUCTED' % sky_last_combine_name)

        print '[INFO]: Combining sky cubes'

        # Execute the KMOS recipe for combining the
        # multi-cubes - save the output to the science directory
        os.system('esorex --output-dir=$KMOS_SCIENCE kmos_combine'
                  + ' --method="header" --edge_nan=TRUE sky_combine.sof')

    # Everything Now in place for the multiExtractSpec
    # method - which checks the performance of the individual frames
    # This should deposit everything in the science directory

    # if shift has been set to true use the shifted_object_names.txt file
    if shift:

        pipe_methods.multiExtractSpec(sci_dir=science_dir,
                                      frameNames='shifted_object_names.txt',
                                      tracked_list=tracked_list,
                                      pix_scale=pix_scale)

    else:

        pipe_methods.multiExtractSpec(sci_dir=science_dir,
                                      frameNames='corrected_object_names.txt',
                                      tracked_list=tracked_list,
                                      pix_scale=pix_scale)

    # names of the stars for different pointings
    # C_STARS_7656 - GP1
    # C_STARS_4833
    # C_STARS_4266
    # C_STARS_26263 - GP2
    #
    # SSA_P1
    # R60634 - Det_1
    # R39837 - Det_2
    # R47740 - Det_3
    #
    # SSA_P2
    # R33643 - Det_1
    # R26948 - Det_2
    # R33270 - Det_3
    #
    # UDS - 20772

ssa_p1_list = ['R60634', 'R39837', 'R47740']
ssa_p2_list = ['R33643', 'R26948', 'R33270']
goods_p1_h_list = ['C_STARS_7656', 'C_STARS_7656', 'C_STARS_7656']
goods_p2_h_list = ['C_STARS_23991', 'C_STARS_23991', 'C_STARS_23991']
goods_p1_list = ['C_STARS_7656', 'C_STARS_7656', 'C_STARS_7656']
goods_p1_list_v2 = ['C_STARS_7656', 'C_STARS_4833', 'C_STARS_4266']
goods_p2_list = ['C_STARS_23991', 'C_STARS_23991', 'C_STARS_23991']
goods_p2_list_v2 = ['C_STARS_24357', 'C_STARS_23991', 'C_STARS_26263']
est_list = ['C_STARS_7656', 'C_STARS_12280', 'C_STARS_7656']
mason_list = ['Star_MCID_1564', 'Star_MCID_1564', 'Star_MCID_1564']
lee_list_one = ['n55_19', 'n55_19', 'n55_19']
lee_list_two = ['n55_19', 'n55_19', 'n55_19']
clash_list = ['CLASH0000415', 'CLASH0000415','CLASH0000415']


# now run the pipeline
run_pipeline(tracked_list=mason_list,
             quick=True,
             shift=False,
             telluric_method='standard',
             illum_cor='flat_sky',
             pix_scale=0.2)