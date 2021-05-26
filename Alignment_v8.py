# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 13:23:21 2018

@author: mawgui

This is the current Production Program
Source Files are those written to /python/filesnew/<FITS file object name> by the Reduction_v1.py
Source File object name must be manually provided
All of exposures through a given filter are aligned and stacked together.  If this is not desired, exposures must be manually sorted
Output is written to /python/filesnew/<FITS file object name>
"""
###

# Not optimized, runs multiple variantions of combination.  This is a fully functioning and usable tool!! - 022020!!

### Load the Astropy and other functions and programs - not all used
objid = 'ZOW2021aab'
date = '20210504'

def main(date, objid):
    import numpy as np
    
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    #   %matplotlib inline
    
    from astropy.stats import sigma_clip
    from astropy.stats import sigma_clipped_stats
    # from numpy import mean
    
    import glob
    import os
    from astropy.io import fits
    from scipy.ndimage import interpolation as interp
    
    from skimage.feature.register_translation import (register_translation, _upsampled_dft)
    
    ## This turns off warnings: not a great way to code
    ## But when we show the images, sometimes we're taking the logarithm of zero and it doesn't like that
    ## Which would matter if we were doing math, but we're just taking a look at images, so we can ignore it. 
    import warnings
    warnings.filterwarnings('ignore')
    import time
    import shutil
    
    ## Create a function to plot an image cube
    ## for this code a "cube" means a stack of image data arrays
    def plot_grid(datacube,imagenames):
        no_A = len(datacube) ## number of individual images in the cube
        xplots = int(np.around(np.sqrt(no_A))) ## number of image grid columns
        yplots = xplots + 1 ## number of image grid rows--sometimes there are one fewer, but that's okay
    #     print no_A, xplots, yplots ## this line is for troubleshooting
        gs = gridspec.GridSpec(yplots, xplots) ## define the image grid
        plt.figure(figsize=(15,15)) ## set the figure size
        for i in range(no_A): 
            ## plots each individual image within the image grid: 
            B = datacube[i]
            plt.subplot(gs[i])
            plt.imshow(np.log10(B), origin='lower', cmap='gray')
            plt.title(imagenames[i])
    
    #TargetObject = input("Please provide the Name of your target as described in your File Directory name.  You have 5 seconds.  ")
    TargetObject = objid
    date = date
    ## Here we define the path to our Data
    path = '/home/david/ZowadaTransients/Lights/' + TargetObject + "/" + date
    
    processedpath = '/home/david/ZowadaTransients/Processed'
    
    
    #TargetObject = '2MASS J18014991+1100284'
    #TargetObject = 'GD 246'
    
    # TargetObject = 'NGC 2682 FBC 1616' # For Test
    
    #time.sleep(5)
    
    targetdir = processedpath +'/'+ str(TargetObject)+'/'+date
    if not os.path.isdir(processedpath +'/'+ str(TargetObject)):
        os.mkdir(processedpath +'/'+ str(TargetObject))
    if not os.path.isdir(targetdir):
        os.mkdir(targetdir)
    
    
    print(targetdir)
    
    #TOGGLESTACK=input("What stacking method do you wish to perform during alignment?  Answer: 'Summed', 'Averaged', 'Mean' or 'Median'?  You have 3 seconds.  ")
    #TOGGLESTACK = 'Summed'
    #TOGGLESTACK='Averaged'
    
    #print("No worries: this module is not yet functional")
    print(" ")
    print(" ")
    print(" ")
    
    print("You can speed up processing by electing to NOT see selected statistics and graphics.")
    #TOGGLESIGMA=input("Answer 'Yes' or 'No' to whether you wish to see these statistics and data. You have 3 seconds.  ")
    TOGGLESIGMA = 'No'
    time.sleep(3)
    
    #TOGGLESIGMA='Yes'
       
    #$$$ This section identifies the Target Object. Updated code needs to be tested!
    
    # This section creates the directory of "writeto" files for the Target Directory
    '''
    dirlist_rproc=glob.glob(targetdir + '/' + '*-r *.fit')+glob.glob(targetdir + '/' + '*-r *.fts')+glob.glob(targetdir + '/' + '*-r-*.fit')+glob.glob(targetdir + '/' + '*-r-*.fts')
    dirlist1_rproc=dirlist_rproc+glob.glob(targetdir + '/' + '*-r *.fits')+glob.glob(targetdir + '/' + '*-r-*.fits')
    # print(dirlist1_rproc)
    
    dirlist_gproc=glob.glob(targetdir + '/' + '*-g *.fit')+glob.glob(targetdir + '/' + '*-g *.fts')+glob.glob(targetdir + '/' + '*-g-*.fit')+glob.glob(targetdir + '/' + '*-g-*.fts')
    dirlist1_gproc=dirlist_gproc+glob.glob(targetdir + '/' + '*-g *.fits')+glob.glob(targetdir + '/' + '*-g-*.fits')
    # print(dirlist1_gproc)
    
    dirlist_uproc=glob.glob(targetdir + '/' + '*-u *.fit')+glob.glob(targetdir + '/' + '*-u *.fts')+glob.glob(targetdir + '/' + '*-u-*.fit')+glob.glob(targetdir + '/' + '*-u-*.fts')
    dirlist1_uproc=dirlist_uproc+glob.glob(targetdir + '/' + '*-u *.fits')+glob.glob(targetdir + '/' + '*-u-*.fits')
    # print(dirlist1_uproc)
    
    dirlist_iproc=glob.glob(targetdir + '/' + '*-i *.fit')+glob.glob(targetdir + '/' + '*-i *.fts')+glob.glob(targetdir + '/' + '*-i-*.fit')+glob.glob(targetdir + '/' + '*-i-*.fts')
    dirlist1_iproc=dirlist_iproc+glob.glob(targetdir + '/' + '*-i *.fits')+glob.glob(targetdir + '/' + '*-i-*.fits')
    # print(dirlist1_iproc)
    
    dirlist_zproc=glob.glob(targetdir + '/' + '*-z *.fit')+glob.glob(targetdir + '/' + '*-z *.fts')+glob.glob(targetdir + '/' + '*-z-*.fit')+glob.glob(targetdir + '/' + '*-z-*.fts')
    dirlist1_zproc=dirlist_zproc+glob.glob(targetdir + '/' + '*-z *.fits')+glob.glob(targetdir + '/' + '*-z-*.fits')
    # print(dirlist1_zproc)
    
    dirlist_redproc=glob.glob(targetdir + '/' + '*-red *.fit')+glob.glob(targetdir + '/' + '*-red *.fts')+glob.glob(targetdir + '/' + '*-red-*.fit')+glob.glob(targetdir + '/' + '*-red-*.fts')
    dirlist1_redproc=dirlist_redproc+glob.glob(targetdir + '/' + '*-red *.fits')+glob.glob(targetdir + '/' + '*-red-*.fits')
    # print(dirlist1_redproc)
    
    dirlist_blueproc=glob.glob(targetdir + '/' + '*-blue *.fit')+glob.glob(targetdir + '/' + '*-blue *.fts')+glob.glob(targetdir + '/' + '*-blue-*.fit')+glob.glob(targetdir + '/' + '*-blue-*.fts')
    dirlist1_blueproc=dirlist_blueproc+glob.glob(targetdir + '/' + '*-blue *.fits')+glob.glob(targetdir + '/' + '*-blue-*.fits')
    # print(dirlist1_blueproc)
    
    dirlist_greenproc=glob.glob(targetdir + '/' + '*-green *.fit')+glob.glob(targetdir + '/' + '*-green *.fts')+glob.glob(targetdir + '/' + '*-green-*.fit')+glob.glob(targetdir + '/' + '*-green-*.fts')
    dirlist1_greenproc=dirlist_greenproc+glob.glob(targetdir + '/' + '*-green *.fits')+glob.glob(targetdir + '/' + '*-green-*.fits')
    # print(dirlist1_greenproc)
    
    dirlist_Haproc=glob.glob(targetdir + '/' + '*-Ha *.fit')+glob.glob(targetdir + '/' + '*-Ha *.fts')+glob.glob(targetdir + '/' + '*-Ha-*.fit')+glob.glob(targetdir + '/' + '*-Ha-*.fts')
    dirlist1_Haproc=dirlist_Haproc+glob.glob(targetdir + '/' + '*-Ha *.fits')+glob.glob(targetdir + '/' + '*-Ha-*.fits')
    # print(dirlist1_Haproc)
    
    dirlist_lumproc=glob.glob(targetdir + '/' + '*-lum *.fit')+glob.glob(targetdir + '/' + '*-lum *.fts')+glob.glob(targetdir + '/' + '*-lum-*.fit')+glob.glob(targetdir + '/' + '*-lum-*.fts')
    dirlist1_lumproc=dirlist_lumproc+glob.glob(targetdir + '/' + '*-lum *.fits')+glob.glob(targetdir + '/' + '*-lum-*.fits')
    # print(dirlist1_lumproc)
    '''
    dirlist1_rproc=glob.glob(targetdir + '/' + '*-r.fit')+glob.glob(targetdir + '/' + '*-r.fts')+glob.glob(targetdir + '/' + '*-r.fits')
    dirlist1_gproc=glob.glob(targetdir + '/' + '*-g.fit')+glob.glob(targetdir + '/' + '*-g.fts')+glob.glob(targetdir + '/' + '*-g.fits')
    dirlist1_uproc=glob.glob(targetdir + '/' + '*-u.fit')+glob.glob(targetdir + '/' + '*-u.fts')+glob.glob(targetdir + '/' + '*-u.fits')
    dirlist1_iproc=glob.glob(targetdir + '/' + '*-i.fit')+glob.glob(targetdir + '/' + '*-i.fts')+glob.glob(targetdir + '/' + '*-i.fits')
    dirlist1_zproc=glob.glob(targetdir + '/' + '*-z.fit')+glob.glob(targetdir + '/' + '*-z.fts')+glob.glob(targetdir + '/' + '*-z.fits')
    dirlist1_redproc=glob.glob(targetdir + '/' + '*-red.fit')+glob.glob(targetdir + '/' + '*-red.fts')+glob.glob(targetdir + '/' + '*-red.fits')
    dirlist1_blueproc=glob.glob(targetdir + '/' + '*-blue.fit')+glob.glob(targetdir + '/' + '*-blue.fts')+glob.glob(targetdir + '/' + '*-blue.fits')
    dirlist1_greenproc=glob.glob(targetdir + '/' + '*-green.fit')+glob.glob(targetdir + '/' + '*-green.fts')+glob.glob(targetdir + '/' + '*-green.fits')
    dirlist1_Haproc=glob.glob(targetdir + '/' + '*-Ha.fit')+glob.glob(targetdir + '/' + '*-Ha.fts')+glob.glob(targetdir + '/' + '*-Ha.fits')
    dirlist1_lumproc=glob.glob(targetdir + '/' + '*-lum.fit')+glob.glob(targetdir + '/' + '*-lum.fts')+glob.glob(targetdir + '/' + '*-lum.fits')
    
    dirlist1_proc=dirlist1_blueproc+dirlist1_redproc+dirlist1_greenproc+dirlist1_rproc+dirlist1_gproc \
    +dirlist1_iproc+dirlist1_uproc+dirlist1_zproc+dirlist1_lumproc+dirlist1_Haproc
    print(dirlist1_proc)
    
    newraw_image_data = {}
    for image_name in dirlist1_proc: newraw_image_data[image_name] = fits.getdata(image_name)
    
     
    for image in dirlist1_proc: print (image, newraw_image_data[image].shape)
    
    # This section is to build information to describe the stacking performed
    numexpr = len(dirlist1_rproc)
    numexpg = len(dirlist1_gproc)
    numexpu = len(dirlist1_uproc)
    numexpi = len(dirlist1_iproc)
    numexpz = len(dirlist1_zproc)
    numexpred = len(dirlist1_redproc)
    numexpblue = len(dirlist1_blueproc)
    numexpgreen = len(dirlist1_greenproc)
    numexpHa = len(dirlist1_Haproc)
    numexplum = len(dirlist1_lumproc)
    
    
    ##  This is to Display Frames Before Alignment and Average or Summation
    try:
        rimagecube = np.stack([newraw_image_data[science_frame] for science_frame in dirlist1_rproc],axis=0)
        plot_grid(rimagecube,dirlist1_rproc)
    except:
        pass
    
    try: 
        gimagecube = np.stack([newraw_image_data[science_frame] for science_frame in dirlist1_gproc],axis=0)
        # plot_grid(gimagecube,dirlist1_gproc)
    except:
        pass
    
    try: 
        uimagecube = np.stack([newraw_image_data[science_frame] for science_frame in dirlist1_uproc],axis=0)
        # plot_grid(uimagecube,dirlist1_uproc)
    except:
        pass
    
    try: 
        iimagecube = np.stack([newraw_image_data[science_frame] for science_frame in dirlist1_iproc],axis=0)
       # plot_grid(iimagecube,dirlist1_iproc)
    except:
        pass
    
    try: 
        zimagecube = np.stack([newraw_image_data[science_frame] for science_frame in dirlist1_zproc],axis=0)
        # plot_grid(zimagecube,dirlist1_zproc)
    except:
        pass
    
    try: 
        redimagecube = np.stack([newraw_image_data[science_frame] for science_frame in dirlist1_redproc],axis=0)
        # plot_grid(redimagecube,dirlist1_redproc)
    except:
        pass
    
    try: 
        greenimagecube = np.stack([newraw_image_data[science_frame] for science_frame in dirlist1_greenproc],axis=0)
        # plot_grid(greenimagecube,dirlist1_greenproc)
    except:
        pass
    
    try: 
        blueimagecube = np.stack([newraw_image_data[science_frame] for science_frame in dirlist1_blueproc],axis=0)
        # plot_grid(blueimagecube,dirlist1_blueproc)
    except:
        pass
    
    try: 
        Haimagecube = np.stack([newraw_image_data[science_frame] for science_frame in dirlist1_Haproc],axis=0)
        # plot_grid(Haimagecube,dirlist1_Haproc)
    except:
        pass
    
    try: 
        lumimagecube = np.stack([newraw_image_data[science_frame] for science_frame in dirlist1_lumproc],axis=0)
        # plot_grid(lumimagecube,dirlist1_lumproc)
    except:
        pass
    
    try: 
        allimagecube = np.stack([newraw_image_data[science_frame] for science_frame in dirlist1_proc],axis=0)
        if(TOGGLESIGMA == 'Yes'):
            plot_grid(allimagecube,dirlist1_proc)
    except:
        pass
    
    plt.show()
    
    
    
    ##  This is to Display What Stacked Frames would look like before Alignment and Average or Summation
    try: 
        rimage_stacked = np.average(rimagecube, axis=0)
        if(TOGGLESIGMA == 'Yes'):
            plt.figure(figsize=(15,15));
            plt.title('Image: R Unaligned and Stacked');
            plt.imshow(np.log10(rimage_stacked), origin='lower', cmap='gray', vmin=1.0, vmax=5.0)
    except:
        pass
    
    try:
        gimage_stacked = np.average(gimagecube, axis=0)
    #    if(TOGGLESIGMA == 'Yes'):
            #plt.figure(figsize=(15,15));
            #plt.title('Image: G Unaligned and Stacked');
            #plt.imshow(np.log10(gimage_stacked), origin='lower', cmap='gray', vmin=1.0, vmax=5.0)
            #plt.show()
            #mean, median, std = sigma_clipped_stats(gimage_stacked)
            #print(mean, median, std)
    except:
        pass
    
    try: 
        uimage_stacked = np.average(uimagecube, axis=0)
        if(TOGGLESIGMA == 'Yes'):
            plt.figure(figsize=(15,15));
            plt.title('Image: U Unaligned and Stacked');
            plt.imshow(np.log10(uimage_stacked), origin='lower', cmap='gray', vmin=1.0, vmax=5.0)
    except:
        pass
    
    try: 
        iimage_stacked = np.average(iimagecube, axis=0)
    #    if(TOGGLESIGMA == 'Yes'):
            #plt.figure(figsize=(15,15));
            #plt.title('Image: I Unaligned and Stacked');
            #plt.imshow(np.log10(iimage_stacked), origin='lower', cmap='gray', vmin=1.0, vmax=5.0)
    except:
        pass
    
    try:
        zimage_stacked = np.average(zimagecube, axis=0)
    #    if(TOGGLESIGMA == 'Yes'):
            #plt.figure(figsize=(15,15));
            #plt.title('Image: Z Unaligned and Stacked');
            #plt.imshow(np.log10(zimage_stacked), origin='lower', cmap='gray', vmin=1.0, vmax=5.0)
    except:
        pass
    
    try:
        redimage_stacked = np.average(redimagecube, axis=0)
    #    if(TOGGLESIGMA == 'Yes'):
            #plt.figure(figsize=(15,15));
            #plt.title('Image: Red Unaligned and Stacked');
            #plt.imshow(np.log10(redimage_stacked), origin='lower', cmap='gray', vmin=1.0, vmax=5.0)
    except:
        pass
        
    try:
        greenimage_stacked = np.average(greenimagecube, axis=0)
    #    if(TOGGLESIGMA == 'Yes'):
            #plt.figure(figsize=(15,15));
            #plt.title('Image: Green Unaligned and Stacked');
            #plt.imshow(np.log10(greenimage_stacked), origin='lower', cmap='gray', vmin=1.0, vmax=5.0)
    except:
        pass
    
    try:
        blueimage_stacked = np.average(blueimagecube, axis=0)
    #    if(TOGGLESIGMA == 'Yes'):
            #plt.figure(figsize=(15,15));
            #plt.title('Image: Blue Unaligned and Stacked');
            #plt.imshow(np.log10(blueimage_stacked), origin='lower', cmap='gray', vmin=1.0, vmax=5.0)
    except:
        pass
    
    try:
        Haimage_stacked = np.average(Haimagecube, axis=0)
    #    if(TOGGLESIGMA == 'Yes'):
            #plt.figure(figsize=(15,15));
            #plt.title('Image: Ha Unaligned and Stacked');
            #plt.imshow(np.log10(Haimage_stacked), origin='lower', cmap='gray', vmin=1.0, vmax=5.0)
            # plt.show()
    except:
        pass
    
    try:
        lumimage_stacked = np.average(lumimagecube, axis=0)
    #    if(TOGGLESIGMA == 'Yes'):
            #plt.figure(figsize=(15,15));
            #plt.title('Image: Lum Unaligned and Stacked');
            #plt.imshow(np.log10(lumimage_stacked), origin='lower', cmap='gray', vmin=1.0, vmax=5.0)
            # plt.show()
    except:
        pass
    
    try:
        allimage_stacked = np.average(allimagecube, axis=0)
        if(TOGGLESIGMA == 'Yes'):
            plt.figure(figsize=(15,15));
            plt.title('Image: All Unaligned and Stacked');
            plt.imshow(np.log10(allimage_stacked), origin='lower', cmap='gray', vmin=1.0, vmax=5.0)
    except:
        pass 
    
    plt.show()
    
    ##  THIS IS THE ALIGNMENT MODULE per se
    
    ## choose an image to be the baseline of zero shift
    # These need to be unique to each filter because the WCS coordinates are set uniquely in the Headers
    try:
        zerofiler=dirlist1_rproc[0]
        zero_shift_imager = newraw_image_data[zerofiler]
        print(zerofiler, zero_shift_imager)
    except:
        pass
    
    try:
        zerofileg=dirlist1_gproc[0]
        zero_shift_imageg = newraw_image_data[zerofileg]
        print(zerofileg, zero_shift_imageg)
    except:
        pass
    
    try:    
        zerofilei=dirlist1_iproc[0]
        zero_shift_imagei = newraw_image_data[zerofilei]
        print(zerofilei, zero_shift_imagei)
    except:
        pass
    
    try:
        zerofilez=dirlist1_zproc[0]
        zero_shift_imagez = newraw_image_data[zerofilez]
        print(zerofilez, zero_shift_imagez)
    except:
        pass
    
    try:
        zerofileu=dirlist1_uproc[0]
        zero_shift_imageu = newraw_image_data[zerofileu]
        print(zerofileu, zero_shift_imageu)
    except:
        pass
    
    try:
        zerofilered=dirlist1_redproc[0]
        zero_shift_imagered = newraw_image_data[zerofilered]
        print(zerofilered, zero_shift_imagered)
    except:
        pass
    
    try:
        zerofileblue=dirlist1_blueproc[0]
        zero_shift_imageblue = newraw_image_data[zerofileblue]
        print(zerofileblue, zero_shift_imageblue)
    except:
        pass
    
    try:
        zerofilegreen=dirlist1_greenproc[0]
        zero_shift_imagegreen = newraw_image_data[zerofilegreen]
        print(zerofilegreen, zero_shift_imagegreen)
    except:
        pass
    
    try:
        zerofileHa=dirlist1_Haproc[0]
        zero_shift_imageHa = newraw_image_data[zerofileHa]
        print(zerofileHa, zero_shift_imageHa)
    except:
        pass
    
    try:
        zerofilelum=dirlist1_lumproc[0]
        zero_shift_imagelum = newraw_image_data[zerofilelum]
        print(zerofilelum, zero_shift_imagelum)
    except:
        pass
    
    
    ## Imitated format is available in Reduction_Alignment Module Only.py
    
    ## Create dictionaries to hold shifts for other images: 
    imshifts = {} # dictionary to hold the x and y shift pairs for each image
    rimshifts = {} # dictionary to hold the x and y shift pairs for each image
    gimshifts = {} # dictionary to hold the x and y shift pairs for each image
    uimshifts = {} # dictionary to hold the x and y shift pairs for each image
    iimshifts = {} # dictionary to hold the x and y shift pairs for each image
    zimshifts = {} # dictionary to hold the x and y shift pairs for each image
    redimshifts = {} # dictionary to hold the x and y shift pairs for each image
    greenimshifts = {} # dictionary to hold the x and y shift pairs for each image
    blueimshifts = {} # dictionary to hold the x and y shift pairs for each image
    Haimshifts = {} # dictionary to hold the x and y shift pairs for each image
    lumimshifts = {} # dictionary to hold the x and y shift pairs for each image
    allimshifts = {} # dictionary to hold the x and y shift pairs for each image
    
    ## new list for shifted image names: 
    shifted_rimages_out  = ['shifted_' + im for im in dirlist1_rproc]
    shifted_gimages_out  = ['shifted_' + im for im in dirlist1_gproc]
    shifted_uimages_out  = ['shifted_' + im for im in dirlist1_uproc]
    shifted_iimages_out  = ['shifted_' + im for im in dirlist1_iproc]
    shifted_zimages_out  = ['shifted_' + im for im in dirlist1_zproc]
    shifted_redimages_out  = ['shifted_' + im for im in dirlist1_redproc]
    shifted_greenimages_out  = ['shifted_' + im for im in dirlist1_greenproc]
    shifted_blueimages_out  = ['shifted_' + im for im in dirlist1_blueproc]
    shifted_Haimages_out  = ['shifted_' + im for im in dirlist1_Haproc]
    shifted_lumimages_out  = ['shifted_' + im for im in dirlist1_lumproc]
    shifted_allimages_out  = ['shifted_' + im for im in dirlist1_proc]
    
    ###&&&& Through Here Has Been Proofread  -- Maybe rimshifts --> imshifts everywhere, just customize second section
    
    ## new dictionaries and shifts for shifted image data: 
    '''for image in newraw_image_data:
        ## register_translation is a function that calculates shifts by comparing 2-D arrays
        result, error, diffphase = register_translation(
            newraw_image_data[zerofile], 
            newraw_image_data[image], 1000)
        imshifts[image] = result
    
    for image in dirlist1_rproc:
        ## register_translation is a function that calculates shifts by comparing 2-D arrays
        result, error, diffphase = register_translation(
            newraw_image_data[zerofiler], 
            newraw_image_data[image], 1000)
        rimshifts[image] = result
    '''
    for i in range(len(dirlist1_rproc)):
        ## register_translation is a function that calculates shifts by comparing 2-D arrays
        result, error, diffphase = register_translation(
            newraw_image_data[zerofiler], 
            newraw_image_data[dirlist1_rproc[i]], 1000)
        rimshifts[dirlist1_rproc[i]] = result
    
    shifted_r_images = {}
    for i in range(len(dirlist1_rproc)):
        print(dirlist1_rproc[i])
        print(shifted_rimages_out[i])
        ## interp.shift is the function doing the heavy lifting here,
        ## it's reinterpolating each array into the new, shifted one
        shifted_r_images[shifted_rimages_out[i]] = interp.shift(
            newraw_image_data[dirlist1_rproc[i]], 
            rimshifts[dirlist1_rproc[i]])
        print(shifted_r_images[shifted_rimages_out[i]])
        print(newraw_image_data[dirlist1_rproc[i]])
    
    for i in range(len(dirlist1_gproc)):
        ## register_translation is a function that calculates shifts by comparing 2-D arrays
        result, error, diffphase = register_translation(
            newraw_image_data[zerofileg], 
            newraw_image_data[dirlist1_gproc[i]], 1000)
        gimshifts[dirlist1_gproc[i]] = result
    
    shifted_g_images = {}
    for i in range(len(dirlist1_gproc)):
        print(dirlist1_gproc[i])
        print(shifted_gimages_out[i])
        ## interp.shift is the function doing the heavy lifting here,
        ## it's reinterpolating each array into the new, shifted one
        shifted_g_images[shifted_gimages_out[i]] = interp.shift(
            newraw_image_data[dirlist1_gproc[i]], 
            gimshifts[dirlist1_gproc[i]])
        print(shifted_g_images[shifted_gimages_out[i]])
        print(newraw_image_data[dirlist1_gproc[i]])
    
    for i in range(len(dirlist1_iproc)):
        ## register_translation is a function that calculates shifts by comparing 2-D arrays
        result, error, diffphase = register_translation(
            newraw_image_data[zerofilei], 
            newraw_image_data[dirlist1_iproc[i]], 1000)
        iimshifts[dirlist1_iproc[i]] = result
    
    shifted_i_images = {}
    for i in range(len(dirlist1_iproc)):
        print(dirlist1_iproc[i])
        print(shifted_iimages_out[i])
        ## interp.shift is the function doing the heavy lifting here,
        ## it's reinterpolating each array into the new, shifted one
        shifted_i_images[shifted_iimages_out[i]] = interp.shift(
            newraw_image_data[dirlist1_iproc[i]], 
            iimshifts[dirlist1_iproc[i]])
        print(shifted_i_images[shifted_iimages_out[i]])
        print(newraw_image_data[dirlist1_iproc[i]])
    
    for i in range(len(dirlist1_zproc)):
        ## register_translation is a function that calculates shifts by comparing 2-D arrays
        result, error, diffphase = register_translation(
            newraw_image_data[zerofilez], 
            newraw_image_data[dirlist1_zproc[i]], 1000)
        zimshifts[dirlist1_zproc[i]] = result
    
    shifted_z_images = {}
    for i in range(len(dirlist1_zproc)):
        print(dirlist1_zproc[i])
        print(shifted_zimages_out[i])
        ## interp.shift is the function doing the heavy lifting here,
        ## it's reinterpolating each array into the new, shifted one
        shifted_z_images[shifted_zimages_out[i]] = interp.shift(
            newraw_image_data[dirlist1_zproc[i]], 
            zimshifts[dirlist1_zproc[i]])
        print(shifted_z_images[shifted_zimages_out[i]])
        print(newraw_image_data[dirlist1_zproc[i]])
    
    for i in range(len(dirlist1_uproc)):
        ## register_translation is a function that calculates shifts by comparing 2-D arrays
        result, error, diffphase = register_translation(
            newraw_image_data[zerofileu], 
            newraw_image_data[dirlist1_uproc[i]], 1000)
        uimshifts[dirlist1_uproc[i]] = result
    
    shifted_u_images = {}
    for i in range(len(dirlist1_uproc)):
        print(dirlist1_uproc[i])
        print(shifted_uimages_out[i])
        ## interp.shift is the function doing the heavy lifting here,
        ## it's reinterpolating each array into the new, shifted one
        shifted_u_images[shifted_uimages_out[i]] = interp.shift(
            newraw_image_data[dirlist1_uproc[i]], 
            uimshifts[dirlist1_uproc[i]])
        print(shifted_u_images[shifted_uimages_out[i]])
        print(newraw_image_data[dirlist1_uproc[i]])
    
    for i in range(len(dirlist1_blueproc)):
        ## register_translation is a function that calculates shifts by comparing 2-D arrays
        result, error, diffphase = register_translation(
            newraw_image_data[zerofileblue], 
            newraw_image_data[dirlist1_blueproc[i]], 1000)
        blueimshifts[dirlist1_blueproc[i]] = result
    
    shifted_blue_images = {}
    for i in range(len(dirlist1_blueproc)):
        print(dirlist1_blueproc[i])
        print(shifted_blueimages_out[i])
        ## interp.shift is the function doing the heavy lifting here,
        ## it's reinterpolating each array into the new, shifted one
        shifted_blue_images[shifted_blueimages_out[i]] = interp.shift(
            newraw_image_data[dirlist1_blueproc[i]], 
            blueimshifts[dirlist1_blueproc[i]])
        print(shifted_blue_images[shifted_blueimages_out[i]])
        print(newraw_image_data[dirlist1_blueproc[i]])
        
    for i in range(len(dirlist1_greenproc)):
        ## register_translation is a function that calculates shifts by comparing 2-D arrays
        result, error, diffphase = register_translation(
            newraw_image_data[zerofilegreen], 
            newraw_image_data[dirlist1_greenproc[i]], 1000)
        greenimshifts[dirlist1_greenproc[i]] = result
    
    shifted_green_images = {}
    for i in range(len(dirlist1_greenproc)):
        print(dirlist1_greenproc[i])
        print(shifted_greenimages_out[i])
        ## interp.shift is the function doing the heavy lifting here,
        ## it's reinterpolating each array into the new, shifted one
        shifted_green_images[shifted_greenimages_out[i]] = interp.shift(
            newraw_image_data[dirlist1_greenproc[i]], 
            greenimshifts[dirlist1_greenproc[i]])
        print(shifted_green_images[shifted_greenimages_out[i]])
        print(newraw_image_data[dirlist1_greenproc[i]])
    
    for i in range(len(dirlist1_redproc)):
        ## register_translation is a function that calculates shifts by comparing 2-D arrays
        result, error, diffphase = register_translation(
            newraw_image_data[zerofilered], 
            newraw_image_data[dirlist1_redproc[i]], 1000)
        redimshifts[dirlist1_redproc[i]] = result
    
    shifted_red_images = {}
    for i in range(len(dirlist1_redproc)):
        print(dirlist1_redproc[i])
        print(shifted_redimages_out[i])
        ## interp.shift is the function doing the heavy lifting here,
        ## it's reinterpolating each array into the new, shifted one
        shifted_red_images[shifted_redimages_out[i]] = interp.shift(
            newraw_image_data[dirlist1_redproc[i]], 
            redimshifts[dirlist1_redproc[i]])
        print(shifted_red_images[shifted_redimages_out[i]])
        print(newraw_image_data[dirlist1_redproc[i]])
    
    for i in range(len(dirlist1_Haproc)):
        ## register_translation is a function that calculates shifts by comparing 2-D arrays
        result, error, diffphase = register_translation(
            newraw_image_data[zerofileHa], 
            newraw_image_data[dirlist1_Haproc[i]], 1000)
        Haimshifts[dirlist1_Haproc[i]] = result
    
    shifted_Ha_images = {}
    for i in range(len(dirlist1_Haproc)):
        print(dirlist1_Haproc[i])
        print(shifted_Haimages_out[i])
        ## interp.shift is the function doing the heavy lifting here,
        ## it's reinterpolating each array into the new, shifted one
        shifted_Ha_images[shifted_Haimages_out[i]] = interp.shift(
            newraw_image_data[dirlist1_Haproc[i]], 
            Haimshifts[dirlist1_Haproc[i]])
        print(shifted_Ha_images[shifted_Haimages_out[i]])
        print(newraw_image_data[dirlist1_Haproc[i]])
    
    for i in range(len(dirlist1_lumproc)):
        ## register_translation is a function that calculates shifts by comparing 2-D arrays
        result, error, diffphase = register_translation(
            newraw_image_data[zerofilelum], 
            newraw_image_data[dirlist1_lumproc[i]], 1000)
        lumimshifts[dirlist1_lumproc[i]] = result
    
    shifted_lum_images = {}
    for i in range(len(dirlist1_lumproc)):
        print(dirlist1_lumproc[i])
        print(shifted_lumimages_out[i])
        ## interp.shift is the function doing the heavy lifting here,
        ## it's reinterpolating each array into the new, shifted one
        shifted_lum_images[shifted_lumimages_out[i]] = interp.shift(
            newraw_image_data[dirlist1_lumproc[i]], 
            lumimshifts[dirlist1_lumproc[i]])
        print(shifted_lum_images[shifted_lumimages_out[i]])
        print(newraw_image_data[dirlist1_lumproc[i]])
    
    for i in range(len(dirlist1_blueproc)):
        ## register_translation is a function that calculates shifts by comparing 2-D arrays
        result, error, diffphase = register_translation(
            newraw_image_data[zerofileblue], 
            newraw_image_data[dirlist1_blueproc[i]], 1000)
        blueimshifts[dirlist1_blueproc[i]] = result
    
    '''
    shifted_all_images = {}
    for i in range(len(dirlist1_proc)):
        print(dirlist1_proc[i])
        print(shifted_allimages_out[i])
        ## interp.shift is the function doing the heavy lifting here,
        ## it's reinterpolating each array into the new, shifted one
        shifted_all_images[shifted_allimages_out[i]] = interp.shift(
            newraw_image_data[dirlist1_proc[i]], 
            imshifts[dirlist1_proc[i]])
    #    print(shifted_all_images[shifted_allimages_out[i]])
    #    print(newraw_image_data[dirlist1_proc[i]])
    '''
    #In [29]:
    ## Create aligned arrays: 
    try:
        blueimagencube  = np.stack(shifted_blue_images.values(),axis=0)
    except:
        pass
    
    try:
        redimagencube  = np.stack(shifted_red_images.values(),axis=0)
    except:
        pass
    
    try:
        greenimagencube  = np.stack(shifted_green_images.values(),axis=0)
    except:
        pass
    
    try:
        rimagencube  = np.stack(shifted_r_images.values(),axis=0)
    except:
        pass
    
    try:
        gimagencube  = np.stack(shifted_g_images.values(),axis=0)
    except:
        pass
    
    try:
        uimagencube  = np.stack(shifted_u_images.values(),axis=0)
    except:
        pass
    
    try:
        iimagencube  = np.stack(shifted_i_images.values(),axis=0)
    except:
        pass
    
    try:
        zimagencube  = np.stack(shifted_z_images.values(),axis=0)
    except:
        pass
    
    try:
        lumimagencube  = np.stack(shifted_lum_images.values(),axis=0)
    except:
        pass
    
    try:
        Haimagencube  = np.stack(shifted_Ha_images.values(),axis=0)
    except:
        pass
    '''
    try:
        allimagencube  = np.stack(shifted_all_images.values(),axis=0)
    except:
        pass
    '''
    ## average combined final image: 
    try:
        blueimagen_stacked = np.average(blueimagencube, axis=0)
    except:
        pass
    
    try:
        redimagen_stacked = np.average(redimagencube, axis=0)
    except:
        pass
    
    try:
        greenimagen_stacked = np.average(greenimagencube, axis=0)
    except:
        pass
    
    try:
        rimagen_stacked = np.average(rimagencube, axis=0)
    except:
        pass
    
    try:
        gimagen_stacked = np.average(gimagencube, axis=0)
    except:
        pass
    
    try:
        uimagen_stacked = np.average(uimagencube, axis=0)
    except:
        pass
    
    try:
        iimagen_stacked = np.average(iimagencube, axis=0)
    except:
        pass
    
    try:
        zimagen_stacked = np.average(zimagencube, axis=0)
    except:
        pass
    
    try:
        lumimagen_stacked = np.average(lumimagencube, axis=0)
    except:
        pass
    
    try:
        Haimagen_stacked = np.average(Haimagencube, axis=0)
    except:
        pass
    '''
    try:
        allimagen_stacked = np.average(allimagencube, axis=0)
    except:
        pass
    '''
    # imagen_cube.shape  
    
    #mean, median, std = sigma_clipped_stats(gimagen_stacked)
    #print(mean, median, std)
    
    ## Sum combined final image: 
    try:
        blueimagen_summed = np.sum(blueimagencube, axis=0)
    except:
        pass
    
    try:
        redimagen_summed = np.sum(redimagencube, axis=0)
    except:
        pass
    
    try:
        greenimagen_summed = np.sum(greenimagencube, axis=0)
    except:
        pass
    
    try:
        rimagen_summed = np.sum(rimagencube, axis=0)
    except:
        pass
    
    try:
        gimagen_summed = np.sum(gimagencube, axis=0)
    except:
        pass
    
    try:
        uimagen_summed = np.sum(uimagencube, axis=0)
    except:
        pass
    
    try:
        iimagen_summed = np.sum(iimagencube, axis=0)
    except:
        pass
    
    try:
        zimagen_summed = np.sum(zimagencube, axis=0)
    except:
        pass
    
    try:
        lumimagen_summed = np.sum(lumimagencube, axis=0)
    except:
        pass
    
    try:
        Haimagen_summed = np.sum(Haimagencube, axis=0)
    except:
        pass
    '''
    try:
        allimagen_summed = np.sum(allimagencube, axis=0)
    except:
        pass
    '''
    #mean, median, std = sigma_clipped_stats(rimagen_summed)
    #print(mean, median, std)
    
    
    #mean, median, std = sigma_clipped_stats(gimagen_summed)
    #print(mean, median, std)
    
    ## mean combined final image, ignoring nan: 
    try:
        blueimagen_nanmean = np.nanmean(blueimagencube, axis=0)
    except:
        pass
    
    try:
        redimagen_nanmean = np.nanmean(redimagencube, axis=0)
    except:
        pass
    
    try:
        greenimagen_nanmean = np.nanmean(greenimagencube, axis=0)
    except:
        pass
    
    try:
        rimagen_nanmean = np.nanmean(rimagencube, axis=0)
    except:
        pass
    
    try:
        gimagen_nanmean = np.nanmean(gimagencube, axis=0)
    except:
        pass
    
    try:
        uimagen_nanmean = np.nanmean(uimagencube, axis=0)
    except:
        pass
    
    try:
        iimagen_nanmean = np.nanmean(iimagencube, axis=0)
    except:
        pass
    
    try:
        zimagen_nanmean = np.nanmean(zimagencube, axis=0)
    except:
        pass
    
    try:
        lumimagen_nanmean = np.nanmean(lumimagencube, axis=0)
    except:
        pass
    
    try:
        Haimagen_nanmean = np.nanmean(Haimagencube, axis=0)
    except:
        pass
    
    '''
    try:
        allimagen_nanmean = np.nanmean(allimagencube, axis=0)
    except:
        pass
    '''
    ## median combined final image, ignoring nan: 
    try:
        blueimagen_nanmedian = np.nanmedian(blueimagencube, axis=0)
    except:
        pass
    
    try:
        redimagen_nanmedian = np.nanmedian(redimagencube, axis=0)
    except:
        pass
    
    try:
        greenimagen_nanmedian = np.nanmedian(greenimagencube, axis=0)
    except:
        pass
    
    try:
        rimagen_nanmedian = np.nanmedian(rimagencube, axis=0)
    except:
        pass
    
    try:
        gimagen_nanmedian = np.nanmedian(gimagencube, axis=0)
    except:
        pass
    
    try:
        uimagen_nanmedian = np.nanmedian(uimagencube, axis=0)
    except:
        pass
    
    try:
        iimagen_nanmedian = np.nanmedian(iimagencube, axis=0)
    except:
        pass
    
    try:
        zimagen_nanmedian = np.nanmedian(zimagencube, axis=0)
    except:
        pass
    
    try:
        lumimagen_nanmedian = np.nanmedian(lumimagencube, axis=0)
    except:
        pass
    
    try:
        Haimagen_nanmedian = np.nanmedian(Haimagencube, axis=0)
    except:
        pass
    
    '''
    try:
        allimagen_nanmedian = np.nanmedian(allimagencube, axis=0)
    except:
        pass
    '''
    
    #In [30]:
    ## show the final summed arrays as images: 
    graymin = 2
    graymax = 5
    
    try:
        if(TOGGLESIGMA == 'Yes'):
            plt.figure(figsize=(15,15));
            plt.title('Image: R Aligned and Stacked');
            plt.imshow(np.log10(rimagen_summed), origin='lower', cmap='gray', vmin=graymin, vmax=graymax)
    except: 
        pass
    
    
    try:
        if(TOGGLESIGMA == 'Yes'):
            plt.figure(figsize=(15,15));
            plt.title('Image: U Aligned and Stacked');
            plt.imshow(np.log10(uimagen_summed), origin='lower', cmap='gray', vmin=graymin, vmax=graymax)
    except:
        pass
    
    
    
    try:
        if(TOGGLESIGMA == 'Yes'):
            plt.figure(figsize=(15,15));
            plt.title('Image: G Aligned and Stacked');
            plt.imshow(np.log10(gimagen_summed), origin='lower', cmap='gray', vmin=graymin, vmax=graymax)
    except:
        pass
    
    try:
        if(TOGGLESIGMA == 'Yes'):
            plt.figure(figsize=(15,15));
            plt.title('Image: I Aligned and Stacked');
            plt.imshow(np.log10(iimagen_summed), origin='lower', cmap='gray', vmin=graymin, vmax=graymax)
    except:
        pass
    
    try:
        if(TOGGLESIGMA == 'Yes'):
            plt.figure(figsize=(15,15));
            plt.title('Image: Z Aligned and Stacked');
            plt.imshow(np.log10(zimagen_summed), origin='lower', cmap='gray', vmin=graymin, vmax=graymax)
    except:
        pass
    
    try:
        if(TOGGLESIGMA == 'Yes'):
            plt.figure(figsize=(15,15));
            plt.title('Image: Red Aligned and Stacked');
            plt.imshow(np.log10(redimagen_summed), origin='lower', cmap='gray', vmin=graymin, vmax=graymax)
    except:
        pass
    
    try:
        if(TOGGLESIGMA == 'Yes'):
            plt.figure(figsize=(15,15));
            plt.title('Image: Green Aligned and Stacked');
            plt.imshow(np.log10(greenimagen_summed), origin='lower', cmap='gray', vmin=graymin, vmax=graymax)
    except:
        pass
    
    try:
        if(TOGGLESIGMA == 'Yes'):
            plt.figure(figsize=(15,15));
            plt.title('Image: Blue Aligned and Stacked');
            plt.imshow(np.log10(blueimagen_summed), origin='lower', cmap='gray', vmin=graymin, vmax=graymax)
    except:
        pass
    
    try:
        if(TOGGLESIGMA == 'Yes'):
            plt.figure(figsize=(15,15));
            plt.title('Image: Ha Aligned and Stacked');
            plt.imshow(np.log10(Haimagen_summed), origin='lower', cmap='gray', vmin=graymin, vmax=graymax)
    except:
        pass
    
    try:
        if(TOGGLESIGMA == 'Yes'):
            plt.figure(figsize=(15,15));
            plt.title('Image: Lum Aligned and Stacked');
            plt.imshow(np.log10(lumimagen_summed), origin='lower', cmap='gray', vmin=graymin, vmax=graymax)
    except:
        pass
    
    '''
    try:
        if(TOGGLESIGMA == 'Yes'):
            plt.figure(figsize=(15,15));
            plt.title('Image: All Aligned and Stacked');
            plt.imshow(np.log10(allimagen_summed), origin='lower', cmap='gray', vmin=graymin, vmax=graymax)
    except:
        pass
    '''
    
    plt.show()
    
    
    ### This section writes the Averaged and Summed Data to New Files
    
    # This section creates the Target Directories for Ave and Sum and Nanmean and Nanmedian
    avetargetdir = targetdir + '/' + 'averaged'
    print(avetargetdir)
    
    if not os.path.isdir(avetargetdir):
        os.mkdir(avetargetdir)
    
    
    sumtargetdir = targetdir + '/' + 'summed'
    print(sumtargetdir)
    
    if not os.path.isdir(sumtargetdir):
        os.mkdir(sumtargetdir)
    
    
    nanmeantargetdir = targetdir + '/' + 'nanmean'
    print(nanmeantargetdir)
    
    if not os.path.isdir(nanmeantargetdir):
        os.mkdir(nanmeantargetdir)
    
        
    nanmediantargetdir = targetdir + '/' + 'nanmedian'
    print(nanmediantargetdir)
    
    if not os.path.isdir(nanmediantargetdir):
        os.mkdir(nanmediantargetdir)
    
    # Here we form the base files for overwriting
    
    # These are for the Averaged files
    try:
        source = dirlist1_rproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = avetargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_uproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = avetargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_gproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = avetargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_iproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = avetargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_zproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = avetargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_redproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = avetargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_blueproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = avetargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_greenproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = avetargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_Haproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = avetargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_lumproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = avetargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    # These are for the Summed files
    try:
        source = dirlist1_rproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = sumtargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_uproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = sumtargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_gproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = sumtargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_iproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = sumtargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_zproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = sumtargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_redproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = sumtargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_blueproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = sumtargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_greenproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = sumtargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_Haproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = sumtargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_lumproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = sumtargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    # These are for the Nanmean files
    try:
        source = dirlist1_rproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = nanmeantargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_uproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = nanmeantargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_gproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = nanmeantargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_iproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = nanmeantargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_zproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = nanmeantargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_redproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = nanmeantargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_blueproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = nanmeantargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_greenproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = nanmeantargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_Haproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = nanmeantargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_lumproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = nanmeantargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    # These are for the Nanmedian files
    try:
        source = dirlist1_rproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = nanmediantargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_uproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = nanmediantargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_gproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = nanmediantargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_iproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = nanmediantargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_zproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = nanmediantargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_redproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = nanmediantargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_blueproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = nanmediantargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_greenproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = nanmediantargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_Haproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = nanmediantargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    try:
        source = dirlist1_lumproc[0]
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = nanmediantargetdir + '/' + name1
        shutil.copy2(source, destination)
    except:
        pass
    
    # Here we create the FITS file for the Average Frames
    
    path=avetargetdir
    
    # negative values
    #    baselined_image[baselined_image<0] = 0  # To remove negative valuees
    #    darkdeb_gimdata_out['d_' + i] = baselined_image
    #    threshold_image = darkdeb_gimdata_out['d_' + i] 
    #    threshold_image[threshold_image<min] = 0.00001 # To filter small valuees
     #   darkdeb_gimdata_out['d_' + i] = threshold_image
    
    i = 0
    
    try:
        source = dirlist1_rproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = rimagen_stacked
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and average combined image')
    #    procHDU.header['STACK_COUNT'] = 25
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpr
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_uproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = uimagen_stacked
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and average combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpu
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_gproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = gimagen_stacked
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and average combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpg
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_iproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = iimagen_stacked
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and average combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpi
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_zproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = zimagen_stacked
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and average combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpz
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_redproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = redimagen_stacked
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and average combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpred
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_blueproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = blueimagen_stacked
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and average combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpblue
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_greenproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = greenimagen_stacked
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and average combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpgreen
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_Haproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = Haimagen_stacked
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and average combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpHa
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_lumproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = lumimagen_stacked
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and average combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexplum
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    # Here we create the FITS file for the Summed Frames
    
    path=sumtargetdir
    
    try:
        source = dirlist1_rproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = rimagen_summed
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and summed image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpr
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_uproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = uimagen_summed
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and summed image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpu
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_gproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = gimagen_summed
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and summed image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpg
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_iproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = iimagen_summed
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and summed image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpi
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_zproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = zimagen_summed
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and summed image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpz
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_redproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = redimagen_summed
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and summed image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpred
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_blueproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = blueimagen_summed
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and summed image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpblue
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_greenproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = greenimagen_summed
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and summed image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpgreen
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_Haproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = Haimagen_summed
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and summed image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpHa
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_lumproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = lumimagen_summed
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and summed image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexplum
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    ## This is for NanMean combination
    
    path=nanmeantargetdir
    
    try:
        source = dirlist1_rproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = rimagen_nanmean
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and nanmean combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpr
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_uproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = uimagen_nanmean
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and nanmean combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpu
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_gproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = gimagen_nanmean
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and nanmean combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpg
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_iproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = iimagen_nanmean
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and nanmean combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpi
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_zproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = zimagen_nanmean
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and nanmean combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpz
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_redproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = redimagen_nanmean
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and nanmean combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpred
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_blueproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = blueimagen_nanmean
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and nanmean combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpblue
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_greenproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = greenimagen_nanmean
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and nanmean combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpgreen
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_Haproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = Haimagen_nanmean
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and nanmean combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpHa
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_lumproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = lumimagen_nanmean
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and nanmean combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexplum
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    ## Finally, let's look at NanMedian combination
    
    path=nanmediantargetdir
    
    try:
        source = dirlist1_rproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = rimagen_nanmedian
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and nanmedian combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpr
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_uproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = uimagen_nanmedian
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and nanmedian combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpu
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_gproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = gimagen_nanmedian
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and nanmedian combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpg
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_iproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = iimagen_nanmedian
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and nanmedian combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpi
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_zproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = zimagen_nanmedian
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and nanmedian combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpz
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_redproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = redimagen_nanmedian
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and nanmedian combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpred
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_blueproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = blueimagen_nanmedian
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and nanmedian combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpblue
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_greenproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = greenimagen_nanmedian
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and nanmedian combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpgreen
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_Haproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = Haimagen_nanmedian
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and nanmedian combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexpHa
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass
    
    try:
        source = dirlist1_lumproc[i]
        print(newraw_image_data[source])
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(newraw_image_data[source])
            print(mean, median, std)
        name = str(source)
        name1 = name.replace(targetdir,'')
        destination = path + '/' + name1
        print(destination)
        os.chdir(path)
        HDUList=fits.open(destination)
        primaryHeader=HDUList[0].header
        HDUList.close()
        procData = lumimagen_nanmedian
        procData = np.uint16(procData)
        print(procData)
        procHDU = fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Aligned and nanmedian combined image')
        procHDU.header.append('Stack_Count', end=True)
        procHDU.header['Stack_Count'] = numexplum
        procHDU.writeto(destination,overwrite=True)
        os.chdir(path)
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(procData)
            print(mean, median, std)
    except:
        pass    
    
    
    #  Let's see what has happened!
        ## R Filter
    try:
        newraw_image_data = {}
    #    
        for image_name in dirlist1_rproc: 
            newraw_image_data[image_name] = fits.getdata(image_name)
            postproc = newraw_image_data[image_name]
            if(TOGGLESIGMA == 'Yes'):
                mean, median, std = sigma_clipped_stats(postproc)
                print(mean, median, std)
            print(postproc)
    #    
        if(TOGGLESIGMA == 'Yes'):
             mean, median, std = sigma_clipped_stats(rimagen_stacked)
             print(mean, median, std)
             print(rimagen_stacked)
    #    
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(rimagen_summed)
            print(mean, median, std)
            print(rimagen_summed)
    #    
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(rimagen_nanmean)
            print(mean, median, std)
            print(rimagen_nanmean)
    except:
        pass
    
    try:
        mean, median, std = sigma_clipped_stats(rimagen_nanmedian)
        print(mean, median, std)
        print(rimagen_nanmedian)
    except:
        pass
    
    
        ## U Filter
    try:
        newraw_image_data = {}
    #
        for image_name in dirlist1_uproc: 
            newraw_image_data[image_name] = fits.getdata(image_name)
            postproc = newraw_image_data[image_name]
            if(TOGGLESIGMA == 'Yes'):
                mean, median, std = sigma_clipped_stats(postproc)
                print(mean, median, std)
                print(postproc)
    #    
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(uimagen_stacked)
            print(mean, median, std)
            print(uimagen_stacked)
    #    
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(uimagen_summed)
            print(mean, median, std)
            print(uimagen_summed)
    #    
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(uimagen_nanmean)
            print(mean, median, std)
            print(uimagen_nanmean)
    #    
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(uimagen_nanmedian)
            print(mean, median, std)
            print(uimagen_nanmedian)
    except:
        pass
    
    print("Image alignment program has concluded.  Your images have been successfully processed.")
    
    ##  END  ###
if __name__ == '__main__':
    main(date, objid)