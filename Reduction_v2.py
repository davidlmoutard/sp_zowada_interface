# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 13:23:21 2018

@author: mawgui

This is the current Production Program
Source Files must be contained in /python/filesnew
Output is written to /python/filesnew/<FITS file object name>
"""
###

### Load the Astropy and other functions and programs - not all used

def main(date,objid):

    import numpy as np

    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    #   %matplotlib inline
    

    from astropy.stats import sigma_clipped_stats
    # from numpy import mean
    
    import glob
    import os
    from astropy.io import fits
    from scipy.ndimage import interpolation as interp
    
    from skimage.feature.register_translation import (register_translation, _upsampled_dft)
    import gc
    ## This turns off warnings: not a great way to code
    ## But when we show the images, sometimes we're taking the logarithm of zero and it doesn't like that
    ## Which would matter if we were doing math, but we're just taking a look at images, so we can ignore it. 
    import warnings
    warnings.filterwarnings('ignore')
    import time
    gc.enable()
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
    
    print("You can significantly speed up processing by electing to NOT see selected statistics and graphics.")
    #TOGGLESIGMA=input("Answer 'Yes' or 'No' to whether you wish to see these statistics and data. You have 3 seconds.  ")
    TOGGLESIGMA = 'No'
    time.sleep(3)
    
    targetname = objid
    date = date
    ## Here we define the path to our Data
    path = '/home/david/ZowadaTransients/Lights/' + targetname + "/" + date
    calibrationpath = '/home/david/ZowadaTransients/Calibration'
    flatspath = '/home/david/ZowadaTransients/Flats'
    os.chdir(path)
    # print(path)
    
    ## First, we create Directories of Science Images by Filter 
    '''        
    dirlist_blueimages=glob.glob(path + '/' + '*-blue *.fit')+glob.glob(path + '/' + '*-blue *.fts')+glob.glob(path + '/' + '*-blue-*.fit')+glob.glob(path + '/' + '*-blue-*.fts')
    dirlist1_blueimages=dirlist_blueimages+glob.glob(path + '/' + '*-blue *.fits')+glob.glob(path + '/' + '*-blue-*.fits')
    # print(dirlist1_blueimages)
    dirlist_redimages=glob.glob(path + '/' + '*-red *.fit')+glob.glob(path + '/' + '*-red *.fts')+glob.glob(path + '/' + '*-red-*.fit')+glob.glob(path + '/' + '*-red-*.fts')
    dirlist1_redimages=dirlist_redimages+glob.glob(path + '/' + '*-red *.fits')+glob.glob(path + '/' + '*-red-*.fits')
    # print(dirlist1_redimages)
    dirlist_greenimages=glob.glob(path + '/' + '*-green *.fit')+glob.glob(path + '/' + '*-green *.fts')+glob.glob(path + '/' + '*-green-*.fit')+glob.glob(path + '/' + '*-green-*.fts')
    dirlist1_greenimages=dirlist_greenimages+glob.glob(path + '/' + '*-green *.fits')+glob.glob(path + '/' + '*-green-*.fits')
    # print(dirlist1_greenimages)
    dirlist_rimages=glob.glob(path + '/' + '*-r *.fit')+glob.glob(path + '/' + '*-r *.fts')+glob.glob(path + '/' + '*-r-*.fit')+glob.glob(path + '/' + '*-r-*.fts')
    dirlist1_rimages=dirlist_rimages+glob.glob(path + '/' + '*-r *.fits')+glob.glob(path + '/' + '*-r-*.fits')
    # print(dirlist1_rimages)
    dirlist_gimages=glob.glob(path + '/' + '*-g *.fit')+glob.glob(path + '/' + '*-g *.fts')+glob.glob(path + '/' + '*-g-*.fit')+glob.glob(path + '/' + '*-g-*.fts')
    dirlist1_gimages=dirlist_gimages+glob.glob(path + '/' + '*-g *.fits')+glob.glob(path + '/' + '*-g-*.fits')
    # print(dirlist1_gimages)
    dirlist_iimages=glob.glob(path + '/' + '*-i *.fit')+glob.glob(path + '/' + '*-i *.fts')+glob.glob(path + '/' + '*-i-*.fit')+glob.glob(path + '/' + '*-i-*.fts')
    dirlist1_iimages=dirlist_iimages+glob.glob(path + '/' + '*-i *.fits')+glob.glob(path + '/' + '*-i-*.fits')
    # print(dirlist1_iimages)
    dirlist_uimages=glob.glob(path + '/' + '*-u*.fit')+glob.glob(path + '/' + '*-u*.fts')+glob.glob(path + '/' + '*-u-*.fit')+glob.glob(path + '/' + '*-u-*.fts')
    dirlist1_uimages=dirlist_uimages+glob.glob(path + '/' + '*-u*.fits')+glob.glob(path + '/' + '*-u-*.fits')
    print(dirlist1_uimages)
    dirlist_zimages=glob.glob(path + '/' + '*-z *.fit')+glob.glob(path + '/' + '*-z *.fts')+glob.glob(path + '/' + '*-z-*.fit')+glob.glob(path + '/' + '*-z-*.fts')
    dirlist1_zimages=dirlist_zimages+glob.glob(path + '/' + '*-z *.fits')+glob.glob(path + '/' + '*-z-*.fits')
    # print(dirlist1_zimages)
    dirlist_lumimages=glob.glob(path + '/' + '*-lum *.fit')+glob.glob(path + '/' + '*-lum *.fts')+glob.glob(path + '/' + '*-lum-*.fit')+glob.glob(path + '/' + '*-lum-*.fts')
    dirlist1_lumimages=dirlist_lumimages+glob.glob(path + '/' + '*-lum *.fits')+glob.glob(path + '/' + '*-lum-*.fits')
    # print(dirlist1_lumimages)
    dirlist_Haimages=glob.glob(path + '/' + '*-Ha *.fit')+glob.glob(path + '/' + '*-Ha *.fts')+glob.glob(path + '/' + '*-Ha-*.fit')+glob.glob(path + '/' + '*-Ha-*.fts')
    dirlist1_Haimages=dirlist_Haimages+glob.glob(path + '/' + '*-Ha *.fits')+glob.glob(path + '/' + '*-Ha-*.fits')
    # print(dirlist1_Haimages)
    '''
    dirlist1_blueimages=glob.glob(path + '/' + '*-blue.fit')+glob.glob(path + '/' + '*-blue.fts')+glob.glob(path + '/' + '*-blue.fits')
    dirlist1_redimages=glob.glob(path + '/' + '*-red.fit')+glob.glob(path + '/' + '*-red.fts')+glob.glob(path + '/' + '*-red.fits')
    dirlist1_greenimages=glob.glob(path + '/' + '*-green.fit')+glob.glob(path + '/' + '*-green.fts')+glob.glob(path + '/' + '*-green.fits')
    dirlist1_rimages=glob.glob(path + '/' + '*-r.fit')+glob.glob(path + '/' + '*-r.fts')+glob.glob(path + '/' + '*-r.fits')
    dirlist1_uimages=glob.glob(path + '/' + '*-u.fit')+glob.glob(path + '/' + '*-u.fts')+glob.glob(path + '/' + '*-u.fits')
    dirlist1_gimages=glob.glob(path + '/' + '*-g.fit')+glob.glob(path + '/' + '*-g.fts')+glob.glob(path + '/' + '*-g.fits')
    dirlist1_iimages=glob.glob(path + '/' + '*-i.fit')+glob.glob(path + '/' + '*-i.fts')+glob.glob(path + '/' + '*-i.fits')
    dirlist1_zimages=glob.glob(path + '/' + '*-z.fit')+glob.glob(path + '/' + '*-z.fts')+glob.glob(path + '/' + '*-z.fits')
    dirlist1_Haimages=glob.glob(path + '/' + '*-Ha.fit')+glob.glob(path + '/' + '*-Ha.fts')+glob.glob(path + '/' + '*-Ha.fits')
    dirlist1_lumimages=glob.glob(path + '/' + '*-lum.fit')+glob.glob(path + '/' + '*-lum.fts')+glob.glob(path + '/' + '*-lum.fits')
    dirlist1_images=dirlist1_blueimages+dirlist1_redimages+dirlist1_greenimages+dirlist1_rimages+dirlist1_gimages \
    +dirlist1_iimages+dirlist1_uimages+dirlist1_zimages+dirlist1_lumimages+dirlist1_Haimages
    print(dirlist1_images)
    
    ## Then we create Directories for Bias and Dark frames
    dirlist_bias=glob.glob(calibrationpath + '/' + '*Bia*.fits')+glob.glob(calibrationpath + '/' + '*Bia*.fts')
    dirlist1_bias=dirlist_bias+glob.glob(calibrationpath + '/' + '*Bia*.fit')
    # print(dirlist1_bias)
    
    dirlist_dark=glob.glob(calibrationpath + '/' + '*Dar*.fits')+glob.glob(calibrationpath + '/' + '*Dar*.fts')
    dirlist1_dark=dirlist_dark+glob.glob(calibrationpath + '/' + '*Dar*.fit')
    # print(dirlist1_dark)
    
    ## Then we create Directories for Flat frames by Filter
    dirlist_blueflat=glob.glob(flatspath+'/'+'*Blue*.fits')+glob.glob(flatspath+'/'+'*Blue*.fts')
    dirlist1_blueflat=dirlist_blueflat+glob.glob(flatspath+'/'+'*Blue*.fit')
    # print(dirlist1_blueflat)
    dirlist_redflat=glob.glob(flatspath+'/'+'*Red*.fits')+glob.glob(flatspath+'/'+'*Red*.fts')
    dirlist1_redflat=dirlist_redflat+glob.glob(flatspath+'/'+'*Red*.fit')
    # print(dirlist1_redflat)
    dirlist_greenflat=glob.glob(flatspath+'/'+'*Green*.fits')+glob.glob(flatspath+'/'+'*Green*.fts')
    dirlist1_greenflat=dirlist_greenflat+glob.glob(flatspath+'/'+'*Green*.fit')
    # print(dirlist1_greenflat)
    dirlist_lumflat=glob.glob(flatspath+'/'+'*Luminance*.fits')+glob.glob(flatspath+'/'+'*Luminance*.fts')
    dirlist1_lumflat=dirlist_lumflat+glob.glob(flatspath+'/'+'*Luminance*.fit')
    # print(dirlist1_lumflat)
    dirlist_Haflat=glob.glob(flatspath+'/'+'*Ha-*.fits')+glob.glob(flatspath+'/'+'*Ha-*.fts')
    dirlist1_Haflat=dirlist_Haflat+glob.glob(flatspath+'/'+'*Ha-*.fit')
    # print(dirlist1_Haflat)
    dirlist_gflat=glob.glob(flatspath+'/'+'*-g-*.fits')+glob.glob(flatspath+'/'+'*-g-*.fts')
    dirlist1_gflat=dirlist_gflat+glob.glob(flatspath+'/'+'*-g-*.fit')
    # print(dirlist1_gflat)
    dirlist_iflat=glob.glob(flatspath+'/'+'*-i-*.fits')+glob.glob(flatspath+'/'+'*-i-*.fts')
    dirlist1_iflat=dirlist_iflat+glob.glob(flatspath+'/'+'*-i-*.fit')
    # print(dirlist1_iflat)
    dirlist_rflat=glob.glob(flatspath+'/'+'*-r-*.fits')+glob.glob(flatspath+'/'+'*-r-*.fts')
    dirlist1_rflat=dirlist_rflat+glob.glob(flatspath+'/'+'*-r-*.fit')
    # print(dirlist1_rflat)
    dirlist_uflat=glob.glob(flatspath+'/'+'*-u-*.fits')+glob.glob(flatspath+'/'+'*-u-*.fts')
    dirlist1_uflat=dirlist_uflat+glob.glob(flatspath+'/'+'*-u-*.fit')
    # print(dirlist1_uflat)
    dirlist_zflat=glob.glob(flatspath+'/'+'*-z-*.fits')+glob.glob(flatspath+'/'+'*-z-*.fts')
    dirlist1_zflat=dirlist_zflat+glob.glob(flatspath+'/'+'*-z-*.fit')
    # print(dirlist1_zflat)
    dirlist1_flat=dirlist1_blueflat+dirlist1_redflat+dirlist1_greenflat+dirlist1_rflat+dirlist1_gflat \
    +dirlist1_iflat+dirlist1_uflat+dirlist1_zflat+dirlist1_lumflat+dirlist1_Haflat
    # print(dirlist1_flat)
    
    # This is our Directory of ALL Frames
    all_images_list=dirlist1_images+dirlist1_bias+dirlist1_dark+dirlist1_blueflat+dirlist1_redflat \
    +dirlist1_greenflat+dirlist1_rflat+dirlist1_zflat+dirlist1_uflat+dirlist1_iflat+dirlist1_gflat \
    +dirlist1_lumflat+dirlist1_Haflat
    # print(all_images_list)
    
    ## We check to be sure that all exposures are the same size.  This also checks to ensure that we are working okay so far.
    raw_image_data = {}
    for image_name in all_images_list: raw_image_data[image_name] = fits.getdata(image_name)
    
    for image in all_images_list: print (image, raw_image_data[image].shape)
    
    ## Here, we look at size, shape and appearance of the Bias frames
    biascube = np.stack([raw_image_data[bias_frame] for bias_frame in dirlist1_bias],axis=0)
    # print(biascube.shape) ## this line is for troubleshooting
    
    # plot_grid(biascube,dirlist1_bias)
    # plt.show()
    
    ## Next, look at size, shape and appearance of the Dark frames
    darkcube = np.stack([raw_image_data[dark_frame] for dark_frame in dirlist1_dark],axis=0)
    # print(darkcube.shape) ## this line is for troubleshooting
    
    # plot_grid(darkcube,dirlist1_dark)
    # plt.show()
    
    ## Finally, look at Flat frames 
    ## Create Arrays of Flats 
    
    try:
        blueflatcube = np.stack([raw_image_data[vflat_frame] for vflat_frame in dirlist1_blueflat],axis=0)
    except:
        pass
    
    try:
        redflatcube = np.stack([raw_image_data[vflat_frame] for vflat_frame in dirlist1_redflat],axis=0)
    except:
        pass
    
    try:
        rflatcube = np.stack([raw_image_data[vflat_frame] for vflat_frame in dirlist1_rflat],axis=0)
    except:
        pass
    
    try:
        gflatcube = np.stack([raw_image_data[vflat_frame] for vflat_frame in dirlist1_gflat],axis=0)
    except:
        pass
    
    try:
        uflatcube = np.stack([raw_image_data[vflat_frame] for vflat_frame in dirlist1_uflat],axis=0)
    except:
        pass
    
    try:
        iflatcube = np.stack([raw_image_data[vflat_frame] for vflat_frame in dirlist1_iflat],axis=0)
    except:
        pass
    
    try:
        zflatcube = np.stack([raw_image_data[vflat_frame] for vflat_frame in dirlist1_zflat],axis=0)
    except:
        pass
    
    try:
        lumflatcube = np.stack([raw_image_data[vflat_frame] for vflat_frame in dirlist1_lumflat],axis=0)
    except:
        pass
    
    try:
        Haflatcube = np.stack([raw_image_data[vflat_frame] for vflat_frame in dirlist1_Haflat],axis=0)
    except:
        pass
    
    try:
        greenflatcube = np.stack([raw_image_data[vflat_frame] for vflat_frame in dirlist1_greenflat],axis=0)
    except:
        pass    
    
    ## show the Arrays of Flats
    # plot_grid(blueflatcube,dirlist1_blueflat)
    # plot_grid(redflatcube,dirlist1_redflat)
    # plot_grid(greenflatcube,dirlist1_greenflat)
    # plot_grid(rflatcube,dirlist1_rflat)
    # plot_grid(gflatcube,dirlist1_gflat)
    # plot_grid(uflatcube,dirlist1_uflat)
    # plot_grid(iflatcube,dirlist1_iflat)
    # plot_grid(zflatcube,dirlist1_zflat)
    # plot_grid(lumflatcube,dirlist1_lumflat)
    # plot_grid(Haflatcube,dirlist1_Haflat)
    
    
    # plt.show()
    
    
    ## Now we can show the Science Images(s).  We are not printing because we do not need to.
    ## create an array of raw Science Images
    try:
        blueimagecube = np.stack([raw_image_data[science_frame] for science_frame in dirlist1_blueimages],axis=0)
        # plot_grid(blueimagecube,dirlist1_blueimages)
    except:
        pass
    
    try:
        redimagecube = np.stack([raw_image_data[science_frame] for science_frame in dirlist1_redimages],axis=0)
        # plot_grid(redimagecube,dirlist1_redimages)
    except:
        pass
    
    try:
        greenimagecube = np.stack([raw_image_data[science_frame] for science_frame in dirlist1_greenimages],axis=0)
        # plot_grid(greenimagecube,dirlist1_greenimages)
    except:
        pass
    
    try:
        uimagecube = np.stack([raw_image_data[science_frame] for science_frame in dirlist1_uimages],axis=0)
        plot_grid(uimagecube,dirlist1_uimages)
    except:
        pass
    
    try:
        iimagecube = np.stack([raw_image_data[science_frame] for science_frame in dirlist1_iimages],axis=0)
        # plot_grid(iimagecube,dirlist1_iimages)
    except:
        pass
    
    try:
        gimagecube = np.stack([raw_image_data[science_frame] for science_frame in dirlist1_gimages],axis=0)
        # plot_grid(gimagecube,dirlist1_gimages)
    except:
        pass
    
    try:
        rimagecube = np.stack([raw_image_data[science_frame] for science_frame in dirlist1_rimages],axis=0)
        # plot_grid(rimagecube,dirlist1_rimages)
    except:
        pass
    
    try:
        zimagecube = np.stack([raw_image_data[science_frame] for science_frame in dirlist1_zimages],axis=0)
        # plot_grid(zimagecube,dirlist1_zimages)
    except:
        pass
    
    try:
        lumimagecube = np.stack([raw_image_data[science_frame] for science_frame in dirlist1_lumimages],axis=0)
        # plot_grid(lumimagecube,dirlist1_lumimages)
    except:
        pass
    
    try:
        Haimagecube = np.stack([raw_image_data[science_frame] for science_frame in dirlist1_Haimages],axis=0)
        # plot_grid(Haimagecube,dirlist1_Haimages)
    except:
        pass
    
    plt.show()
    
    ## Create Master Dark and Master Bias Frames.  
    master_bias_med = np.median(biascube, axis=0) ## to median combine them instead
    master_bias = master_bias_med
    
    
    ## Debias the Flats
    ## Create Debias Directories for Flats
    debias_bluelist_in = dirlist1_blueflat
    debias_bluelist_out = ['d_' + im for im in debias_bluelist_in]
    debias_redlist_in = dirlist1_redflat
    debias_redlist_out = ['d_' + im for im in debias_redlist_in]
    debias_greenlist_in = dirlist1_greenflat
    debias_greenlist_out = ['d_' + im for im in debias_greenlist_in]
    debias_rlist_in = dirlist1_rflat
    debias_rlist_out = ['d_' + im for im in debias_rlist_in]
    debias_glist_in = dirlist1_gflat
    debias_glist_out = ['d_' + im for im in debias_glist_in]
    debias_ulist_in = dirlist1_uflat
    debias_ulist_out = ['d_' + im for im in debias_ulist_in]
    debias_ilist_in = dirlist1_iflat
    debias_ilist_out = ['d_' + im for im in debias_ilist_in]
    debias_zlist_in = dirlist1_zflat
    debias_zlist_out = ['d_' + im for im in debias_zlist_in]
    debias_lumlist_in = dirlist1_lumflat
    debias_lumlist_out = ['d_' + im for im in debias_lumlist_in]
    debias_Halist_in = dirlist1_Haflat
    debias_Halist_out = ['d_' + im for im in debias_Halist_in]
    
    ##  Following lines generate directory of exposure times to go with original Frames.  This is also a good WIP check so I keep the printout
    Expo={}
    for i in all_images_list: 
        fitsdata = {}
        fitsdata[i]=fits.getdata(i)
        hdulist = fits.open(i)
        data = hdulist[0].data.astype(float)
        header = hdulist[0].header
        Exposure = hdulist[0].header['EXPTIME']
        Target = hdulist[0].header['OBJECT'] 
        Expo[i]=Exposure
        print(i)
        print(Exposure)
        print(Target)
    
    # for i in all_images_list: 
    #     print(i)
    #     print(Expo[i])
    
    ##  Following lines generate directory of exposure times for use with Debias Frames
    Expo1 = {}
    debias_all_images_list_in = all_images_list  
    debias_all_images_list_out = ['d_' + im for im in debias_all_images_list_in]
    # print(debias_all_images_list_out)
    
    for i in debias_all_images_list_in:
        Expo1['d_' + i]=Expo[i]
    
    # for i in debias_all_images_list_out:    
    #     print(i)
    #     print(Expo1[i])
    
    ## Create empty Dictionaries for Debiased Flat Frames
    debias_bluedata_out = {}    
    debias_reddata_out = {}
    debias_greendata_out = {}
    debias_rdata_out = {}
    debias_gdata_out = {}
    debias_udata_out = {}
    debias_idata_out = {}
    debias_zdata_out = {}
    debias_lumdata_out = {}
    debias_Hadata_out = {}
    
    ## Subtract Bias Master
    for i in range(len(debias_bluelist_in)):  
        debias_bluedata_out[debias_bluelist_out[i]] = raw_image_data[debias_bluelist_in[i]] - master_bias
    
    for i in range(len(debias_redlist_in)):  
        debias_reddata_out[debias_redlist_out[i]] = raw_image_data[debias_redlist_in[i]] - master_bias
    
    for i in range(len(debias_greenlist_in)):  
        debias_greendata_out[debias_greenlist_out[i]] = raw_image_data[debias_greenlist_in[i]] - master_bias
    
    for i in range(len(debias_rlist_in)):  
        debias_rdata_out[debias_rlist_out[i]] = raw_image_data[debias_rlist_in[i]] - master_bias
    #    print(raw_image_data[debias_rlist_in[i]])
    #    print(debias_rlist_in[i])
    #    print(master_bias)
    #    print(debias_rdata_out[debias_rlist_out[i]])
    #    print(debias_rlist_out[i])
    
    for i in range(len(debias_glist_in)):  
        debias_gdata_out[debias_glist_out[i]] = raw_image_data[debias_glist_in[i]] - master_bias
    #    print(raw_image_data[debias_glist_in[i]])
    #    print(debias_glist_in[i])
    #    print(master_bias)
    #    print(debias_gdata_out[debias_glist_out[i]])
    #    print(debias_glist_out[i])
    
    for i in range(len(debias_ulist_in)):  
        debias_udata_out[debias_ulist_out[i]] = raw_image_data[debias_ulist_in[i]] - master_bias
        print(raw_image_data[debias_ulist_in[i]])
        print(debias_ulist_in[i])
        print(master_bias)
        print(debias_udata_out[debias_ulist_out[i]])
        print(debias_ulist_out[i])
    
    for i in range(len(debias_ilist_in)):  
        debias_idata_out[debias_ilist_out[i]] = raw_image_data[debias_ilist_in[i]] - master_bias
    #    print(raw_image_data[debias_ilist_in[i]])
    #    print(debias_ilist_in[i])
    #    print(master_bias)
    #    print(debias_idata_out[debias_ilist_out[i]])
    #    print(debias_ilist_out[i])
    
    for i in range(len(debias_zlist_in)):  
        debias_zdata_out[debias_zlist_out[i]] = raw_image_data[debias_zlist_in[i]] - master_bias
    
    for i in range(len(debias_lumlist_in)):  
        debias_lumdata_out[debias_lumlist_out[i]] = raw_image_data[debias_lumlist_in[i]] - master_bias
    
    for i in range(len(debias_Halist_in)):  
        debias_Hadata_out[debias_Halist_out[i]] = raw_image_data[debias_Halist_in[i]] - master_bias
    #     print(raw_image_data[debias_Halist_in[i]])
    #     print(debias_Halist_in[i])
    #     print(master_bias)
    #     print(debias_Hadata_out[debias_Halist_out[i]])
    #     print(debias_Halist_out[i])
    
    ## Create Arrays of debiased Flats 
    try:
        bluedebiascube = np.stack([debias_bluedata_out[vdebias_frame] for vdebias_frame in debias_bluelist_out],axis=0)
    except:
        pass
    
    try:
        reddebiascube = np.stack([debias_reddata_out[vdebias_frame] for vdebias_frame in debias_redlist_out],axis=0)
    except:
        pass
    
    try:
        rdebiascube = np.stack([debias_rdata_out[vdebias_frame] for vdebias_frame in debias_rlist_out],axis=0)
    except:
        pass
    
    try:
        gdebiascube = np.stack([debias_gdata_out[vdebias_frame] for vdebias_frame in debias_glist_out],axis=0)
    except:
        pass
    
    try:
        udebiascube = np.stack([debias_udata_out[vdebias_frame] for vdebias_frame in debias_ulist_out],axis=0)
    except:
        pass
    
    try:
        idebiascube = np.stack([debias_idata_out[vdebias_frame] for vdebias_frame in debias_ilist_out],axis=0)
    except:
        pass
    
    try:
        zdebiascube = np.stack([debias_zdata_out[vdebias_frame] for vdebias_frame in debias_zlist_out],axis=0)
    except:
        pass
    
    try:
        lumdebiascube = np.stack([debias_lumdata_out[vdebias_frame] for vdebias_frame in debias_lumlist_out],axis=0)
    except:
        pass
    
    try:
        Hadebiascube = np.stack([debias_Hadata_out[vdebias_frame] for vdebias_frame in debias_Halist_out],axis=0)
    except:
        pass
    
    try:
        greendebiascube = np.stack([debias_greendata_out[vdebias_frame] for vdebias_frame in debias_greenlist_out],axis=0)
    except:
        pass
        
    ## Show the Arrays of Debiased Flats
    # plot_grid(bluedebiascube,debias_bluelist_out)
    # plot_grid(reddebiascube,debias_redlist_out)
    # plot_grid(greendebiascube,debias_greenlist_out)
    # plot_grid(rdebiascube,debias_rlist_out)
    # plot_grid(gdebiascube,debias_glist_out)
    # plot_grid(udebiascube,debias_ulist_out)
    # plot_grid(idebiascube,debias_ilist_out)
    # plot_grid(zdebiascube,debias_zlist_out)
    # plot_grid(lumdebiascube,debias_lumlist_out)
    # plot_grid(Hadebiascube,debias_Halist_out)
    
    # plt.show()
    
    ## Debias the Science Frames
    
    ## Create Debias Directories for Science Frames
    debias_blueimages_in = dirlist1_blueimages
    debias_blueimages_out = ['d_' + im for im in debias_blueimages_in]
    debias_redimages_in = dirlist1_redimages
    debias_redimages_out = ['d_' + im for im in debias_redimages_in]
    debias_greenimages_in = dirlist1_greenimages
    debias_greenimages_out = ['d_' + im for im in debias_greenimages_in]
    debias_rimages_in = dirlist1_rimages
    debias_rimages_out = ['d_' + im for im in debias_rimages_in]
    debias_gimages_in = dirlist1_gimages
    debias_gimages_out = ['d_' + im for im in debias_gimages_in]
    debias_uimages_in = dirlist1_uimages
    debias_uimages_out = ['d_' + im for im in debias_uimages_in]
    debias_iimages_in = dirlist1_iimages
    debias_iimages_out = ['d_' + im for im in debias_iimages_in]
    debias_zimages_in = dirlist1_zimages
    debias_zimages_out = ['d_' + im for im in debias_zimages_in]
    debias_lumimages_in = dirlist1_lumimages
    debias_lumimages_out = ['d_' + im for im in debias_lumimages_in]
    debias_Haimages_in = dirlist1_Haimages
    debias_Haimages_out = ['d_' + im for im in debias_Haimages_in]
    
    ## Create empty Dictionaries
    debias_blueimdata_out = {}
    debias_redimdata_out = {}
    debias_greenimdata_out = {}
    debias_rimdata_out = {}
    debias_gimdata_out = {}
    debias_uimdata_out = {}
    debias_iimdata_out = {}
    debias_zimdata_out = {}
    debias_lumimdata_out = {}
    debias_Haimdata_out = {}
    
    ## Subtract Bias Master
    print("## Subtract Bias Master")
    
    
    for i in range(len(debias_blueimages_in)):  
        debias_blueimdata_out[debias_blueimages_out[i]] = raw_image_data[debias_blueimages_in[i]] - master_bias
    
    for i in range(len(debias_redimages_in)):  
        debias_redimdata_out[debias_redimages_out[i]] = raw_image_data[debias_redimages_in[i]] - master_bias
    
    for i in range(len(debias_greenimages_in)):  
        debias_greenimdata_out[debias_greenimages_out[i]] = raw_image_data[debias_greenimages_in[i]] - master_bias
    
    for i in range(len(debias_rimages_in)):  
        debias_rimdata_out[debias_rimages_out[i]] = raw_image_data[debias_rimages_in[i]] - master_bias
    #    print(debias_rimages_in[i], debias_rimages_out[i])
    #    print(raw_image_data[debias_rimages_in[i]])
    #    print(debias_rimdata_out[debias_rimages_out[i]])
    #    print(master_bias)
    
    for i in range(len(debias_gimages_in)):  
        debias_gimdata_out[debias_gimages_out[i]] = raw_image_data[debias_gimages_in[i]] - master_bias
    #    print(debias_gimages_in[i], debias_gimages_out[i])
    #    print(raw_image_data[debias_gimages_in[i]])
    #    print(debias_gimdata_out[debias_gimages_out[i]])
    #    print(master_bias)
    
    #  Need to investigate negative values for master bias subtracted u image
    for i in range(len(debias_uimages_in)):  
        debias_uimdata_out[debias_uimages_out[i]] = raw_image_data[debias_uimages_in[i]] - master_bias
    #    print(debias_uimages_in[i], debias_uimages_out[i])
    #    print(raw_image_data[debias_uimages_in[i]])
    #    print(debias_uimdata_out[debias_uimages_out[i]])
    #    print(master_bias)
    
    for i in range(len(debias_iimages_in)):  
        debias_iimdata_out[debias_iimages_out[i]] = raw_image_data[debias_iimages_in[i]] - master_bias
    #     print(debias_iimages_in[i], debias_iimages_out[i])
    #     print(raw_image_data[debias_iimages_in[i]])
    #     print(debias_iimdata_out[debias_iimages_out[i]])
    #     print(master_bias)
    
    for i in range(len(debias_zimages_in)):  
        debias_zimdata_out[debias_zimages_out[i]] = raw_image_data[debias_zimages_in[i]] - master_bias
    
    for i in range(len(debias_lumimages_in)):  
        debias_lumimdata_out[debias_lumimages_out[i]] = raw_image_data[debias_lumimages_in[i]] - master_bias
    
    for i in range(len(debias_Haimages_in)):  
        debias_Haimdata_out[debias_Haimages_out[i]] = raw_image_data[debias_Haimages_in[i]] - master_bias
    #     print(debias_Haimages_in[i], debias_Haimages_out[i])
    #     print(raw_image_data[debias_Haimages_in[i]])
    #     print(debias_Haimdata_out[debias_Haimages_out[i]])
    #     print(master_bias)
    
    ## Create Arrays of Debiased Science Images
    try:
        blueimdebiascube = np.stack([debias_blueimdata_out[vdebias_frame] for vdebias_frame in debias_blueimages_out],axis=0)
    except:
        pass
    
    try:
        redimdebiascube = np.stack([debias_redimdata_out[vdebias_frame] for vdebias_frame in debias_redimages_out],axis=0)
    except:
        pass
    
    try:
        rimdebiascube = np.stack([debias_rimdata_out[vdebias_frame] for vdebias_frame in debias_rimages_out],axis=0)
    except:
        pass
    
    try:
        gimdebiascube = np.stack([debias_gimdata_out[vdebias_frame] for vdebias_frame in debias_gimages_out],axis=0)
    except:
        pass
    
    try:
        uimdebiascube = np.stack([debias_uimdata_out[vdebias_frame] for vdebias_frame in debias_uimages_out],axis=0)
    except:
        pass
    
    try:
        iimdebiascube = np.stack([debias_iimdata_out[vdebias_frame] for vdebias_frame in debias_iimages_out],axis=0)
    except:
        pass
    
    try:
        zimdebiascube = np.stack([debias_zimdata_out[vdebias_frame] for vdebias_frame in debias_zimages_out],axis=0)
    except:
        pass
    
    try:
        lumimdebiascube = np.stack([debias_lumimdata_out[vdebias_frame] for vdebias_frame in debias_lumimages_out],axis=0)
    except:
        pass
    
    try:
        Haimdebiascube = np.stack([debias_Haimdata_out[vdebias_frame] for vdebias_frame in debias_Haimages_out],axis=0)
    except:
        pass
    
    try:
        greenimdebiascube = np.stack([debias_greenimdata_out[vdebias_frame] for vdebias_frame in debias_greenimages_out],axis=0)
    except:
        pass
    
    ## Show the Arrays of Debiased Science Frames
    # plot_grid(blueimdebiascube,debias_blueimages_out)
    # plot_grid(redimdebiascube,debias_redimages_out)
    # plot_grid(greenimdebiascube,debias_greenimages_out)
    # plot_grid(rimdebiascube,debias_rimages_out)
    # plot_grid(gimdebiascube,debias_gimages_out)
    # plot_grid(uimdebiascube,debias_uimages_out)
    # plot_grid(iimdebiascube,debias_iimages_out)
    # plot_grid(zimdebiascube,debias_zimages_out)
    # plot_grid(lumimdebiascube,debias_lumimages_out)
    # plot_grid(Haimdebiascube,debias_Haimages_out)
        
    # plt.show()
    
    
    ## Create Dark Subtracted Flats
    # Create Dark Master for for Calculations
    Darkframe=dirlist1_dark[0]
    print(Darkframe) # This is helpful to look at to see just how light the darks are!
    Expo_dark=Expo[Darkframe]
    master_dark_med = np.median(darkcube, axis=0) ## to median combine them instead
    # filtered_data = sigma_clip(darkcube, sigma=3, cenfunc=mean, copy=False)
    # master_dark_sigma = np.median(filtered_data, axis=0) ## to median combine them instead
    
    master_dark = master_dark_med
    # print(master_dark_med) 
    ## Following steps may be useful if I revisit the dark frame calculation
    # print(filtered_data)  
    # delta = master_dark_sigma - master_dark_med
    # print(delta)
    
    # Create Exposure Times for Debiased Frames for Dark Adjustment
    Expo2 = {}
    debias_all_images_list_in2 = debias_all_images_list_out  
    debias_all_images_list_out2 = ['d_' + im for im in debias_all_images_list_in2]
    # print(debias_all_images_list_out2)
    
    for i in debias_all_images_list_in2:
        Expo2['d_' + i]=Expo1[i]
    
    # for i in debias_all_images_list_out2:    
    #    print(i)
    #    print(Expo2[i])
    
    ## Create Directories for Dark Subtracted Flats
    dark_bluelist_in = debias_bluelist_out
    dark_bluelist_out = ['d_' + im for im in dark_bluelist_in]
    dark_redlist_in = debias_redlist_out
    dark_redlist_out = ['d_' + im for im in dark_redlist_in]
    dark_greenlist_in = debias_greenlist_out
    dark_greenlist_out = ['d_' + im for im in dark_greenlist_in]
    dark_rlist_in = debias_rlist_out
    dark_rlist_out = ['d_' + im for im in dark_rlist_in]
    dark_glist_in = debias_glist_out
    dark_glist_out = ['d_' + im for im in dark_glist_in]
    dark_ulist_in = debias_ulist_out
    dark_ulist_out = ['d_' + im for im in dark_ulist_in]
    dark_ilist_in = debias_ilist_out
    dark_ilist_out = ['d_' + im for im in dark_ilist_in]
    dark_zlist_in = debias_zlist_out
    dark_zlist_out = ['d_' + im for im in dark_zlist_in]
    dark_lumlist_in = debias_lumlist_out
    dark_lumlist_out = ['d_' + im for im in dark_lumlist_in]
    dark_Halist_in = debias_Halist_out
    dark_Halist_out = ['d_' + im for im in dark_Halist_in]
    
    
    ## Subtract Scaled Dark Master from Flats
    ## Create empty Dictionaries for Dark Debiased Flats
    darkdeb_bluedata_out = {}
    darkdeb_reddata_out = {}
    darkdeb_greendata_out = {}
    darkdeb_rdata_out = {}
    darkdeb_gdata_out = {}
    darkdeb_udata_out = {}
    darkdeb_idata_out = {}
    darkdeb_zdata_out = {}
    darkdeb_lumdata_out = {}
    darkdeb_Hadata_out = {}
    
    ## Subtract Exposure Adjusted Dark Master from Debiased Flat Frames
    min=0
    
    for i in dark_bluelist_in:
        new_master_dark = {}
        new_master_dark = Expo1[i]/Expo_dark * (master_dark - master_bias)
        baselined_flat = debias_bluedata_out[i] - new_master_dark
    #    baselined_flat[baselined_flat<0] = 0  # To remove negative valuees
        darkdeb_bluedata_out['d_' + i] = baselined_flat
        threshold_flat = darkdeb_bluedata_out['d_' + i] 
        threshold_flat[threshold_flat<min] = 0.00001 # To filter small valuees
        darkdeb_bluedata_out['d_' + i] = threshold_flat
    #     print(darkdeb_bluedata_out['d_' + i])    
    #     print(debias_bluedata_out[i])    
    #     print(master_dark)
    #     print(master_bias)
    #     print(new_master_dark)
    #     print(Expo1[i],Expo_dark,Expo1[i]/Expo_dark)
    
    for i in dark_redlist_in:
        new_master_dark = {}
        new_master_dark = Expo1[i]/Expo_dark * (master_dark - master_bias)
        baselined_flat = debias_reddata_out[i] - new_master_dark # To remove negative values
    #    baselined_flat[baselined_flat<0] = 0  # To remove negative valuees
        darkdeb_reddata_out['d_' + i] = baselined_flat
        threshold_flat = darkdeb_reddata_out['d_' + i] 
        threshold_flat[threshold_flat<min] = 0.00001 # To filter small valuees
        darkdeb_reddata_out['d_' + i] = threshold_flat
    #    print(darkdeb_reddata_out['d_' + i])    
    #    print(debias_reddata_out[i])    
    #    print(master_dark)
    #    print(new_master_dark)
    #    print(Expo1[i],Expo_dark,Expo1[i]/Expo_dark)
    
    for i in dark_greenlist_in:
        new_master_dark = {}
        new_master_dark = Expo1[i]/Expo_dark * (master_dark - master_bias)
        baselined_flat = debias_greendata_out[i] - new_master_dark # To remove negative values
    #    baselined_flat[baselined_flat<0] = 0  # To remove negative valuees
        darkdeb_greendata_out['d_' + i] = baselined_flat
        threshold_flat = darkdeb_greendata_out['d_' + i] 
        threshold_flat[threshold_flat<min] = 0.00001 # To filter small valuees
        darkdeb_greendata_out['d_' + i] = threshold_flat
    #    print(darkdeb_greendata_out['d_' + i])    
    #    print(debias_greendata_out[i])    
    #    print(master_dark)
    #    print(new_master_dark)
    #    print(Expo1[i],Expo_dark,Expo1[i]/Expo_dark)
    
    for i in dark_rlist_in:
        new_master_dark = {}
        new_master_dark = Expo1[i]/Expo_dark * (master_dark - master_bias)
        baselined_flat = debias_rdata_out[i] - new_master_dark # To remove negative values
    #    baselined_flat[baselined_flat<0] = 0  # To remove negative valuees
        darkdeb_rdata_out['d_' + i] = baselined_flat
        threshold_flat = darkdeb_rdata_out['d_' + i] 
        threshold_flat[threshold_flat<min] = 0.00001 # To filter small valuees
        darkdeb_rdata_out['d_' + i] = threshold_flat
    #    print(darkdeb_rdata_out['d_' + i])    
    #    print(debias_rdata_out[i])    
    #    print(master_dark)
    #    print(new_master_dark)
    #    print(Expo1[i],Expo_dark,Expo1[i]/Expo_dark)
    
    for i in dark_glist_in:
        new_master_dark = {}
        new_master_dark = Expo1[i]/Expo_dark * (master_dark - master_bias)
        baselined_flat = debias_gdata_out[i] - new_master_dark # To remove negative values
    #    baselined_flat[baselined_flat<0] = 0  # To remove negative valuees
        darkdeb_gdata_out['d_' + i] = baselined_flat
        threshold_flat = darkdeb_gdata_out['d_' + i] 
        threshold_flat[threshold_flat<min] = 0.00001 # To filter small valuees
        darkdeb_gdata_out['d_' + i] = threshold_flat
    #    print(darkdeb_gdata_out['d_' + i])    
    #    print(debias_gdata_out[i])    
    #    print(master_dark)
    #    print(new_master_dark)
    #    print(Expo1[i],Expo_dark,Expo1[i]/Expo_dark)
    
    #$$$ Note slightly different organization -- Check to ensure these are both okay!
    for i in dark_ulist_in:
        new_master_dark = {}
        new_master_dark = Expo1[i]/Expo_dark * (master_dark - master_bias)
        baselined_flat = debias_udata_out[i] - new_master_dark # To remove negative values
    #    baselined_flat[baselined_flat<0] = 0  # To remove negative valuees
        darkdeb_udata_out['d_' + i] = baselined_flat
    #     threshold_flat = darkdeb_udata_out['d_' + i] 
    #     threshold_flat[threshold_flat<min] = 0.00001 # To filter small valuees
    #    darkdeb_udata_out['d_' + i] = threshold_flat
    #    print(baselined_flat)
    #    print('Baselined Flat')
    #    print(darkdeb_udata_out['d_' + i])    
    #    print(['d_' + i])    
    #    print(debias_udata_out[i])    
    #    print([i])    
    #    print(master_dark)
    #    print('Master Dark')
    #    print(master_bias)
    #    print('Master Bias')    
    #    print(new_master_dark)
    #    print('New Master Dark')
    #    print(Expo1[i],Expo_dark,Expo1[i]/Expo_dark)
    
    for i in dark_ilist_in:
        new_master_dark = {}
        new_master_dark = Expo1[i]/Expo_dark * (master_dark - master_bias)
        baselined_flat = debias_idata_out[i] - new_master_dark # To remove negative values
    #    baselined_flat[baselined_flat<0] = 0  # To remove negative valuees
        darkdeb_idata_out['d_' + i] = baselined_flat
        threshold_flat = darkdeb_idata_out['d_' + i] 
        threshold_flat[threshold_flat<min] = 0.00001 # To filter small valuees
        darkdeb_idata_out['d_' + i] = threshold_flat
    #     print(darkdeb_idata_out['d_' + i])    
    #     print(debias_idata_out[i])    
    #     print(master_dark)
    #     print(new_master_dark)
    #     print(Expo1[i],Expo_dark,Expo1[i]/Expo_dark)
    
    for i in dark_zlist_in:
        new_master_dark = {}
        new_master_dark = Expo1[i]/Expo_dark * (master_dark - master_bias)
    #    baselined_flat = debias_zdata_out[i] - new_master_dark # To remove negative values
        baselined_flat[baselined_flat<0] = 0  # To remove negative valuees
        darkdeb_zdata_out['d_' + i] = baselined_flat
        threshold_flat = darkdeb_zdata_out['d_' + i] 
        threshold_flat[threshold_flat<min] = 0.00001 # To filter small valuees
        darkdeb_zdata_out['d_' + i] = threshold_flat
    #     print(darkdeb_zdata_out['d_' + i])    
    #     print(debias_zdata_out[i])    
    #     print(master_dark)
    #     print(new_master_dark)
    #     print(Expo1[i],Expo_dark,Expo1[i]/Expo_dark)
    
    for i in dark_lumlist_in:
        new_master_dark = {}
        new_master_dark = Expo1[i]/Expo_dark * (master_dark - master_bias)
        baselined_flat = debias_lumdata_out[i] - new_master_dark # To remove negative values
    #     baselined_flat[baselined_flat<0] = 0  # To remove negative valuees
        darkdeb_lumdata_out['d_' + i] = baselined_flat
        threshold_flat = darkdeb_lumdata_out['d_' + i] 
        threshold_flat[threshold_flat<min] = 0.00001 # To filter small valuees
        darkdeb_lumdata_out['d_' + i] = threshold_flat
    #     print(darkdeb_lumdata_out['d_' + i])    
    #     print(debias_lumdata_out[i])    
    #     print(master_dark)
    #     print(new_master_dark)
    #     print(Expo1[i],Expo_dark,Expo1[i]/Expo_dark)
    
    for i in dark_Halist_in:
        min = 0
        new_master_dark = {}
        new_master_dark = Expo1[i]/Expo_dark * (master_dark - master_bias)
        baselined_flat = debias_Hadata_out[i] - new_master_dark # To remove negative values
    #     baselined_flat[baselined_flat<0] = 0  # To remove negative valuees
        darkdeb_Hadata_out['d_' + i] = baselined_flat
        threshold_flat = darkdeb_Hadata_out['d_' + i] 
    #     threshold_flat[threshold_flat<min] = 0.00001 # To filter small valuees
        darkdeb_Hadata_out['d_' + i] = threshold_flat
    #     print(darkdeb_Hadata_out['d_' + i])    
    #     print(i)
    #     print(debias_Hadata_out[i])    
    #     print(master_dark)
    #     print(master_bias)
    #     print(new_master_dark)
    #     print(Expo1[i],Expo_dark,Expo1[i]/Expo_dark)
    
    ## Create Arrays of Dark Debiased Flats 
    try:
        bluedarkcube = np.stack([darkdeb_bluedata_out[vdark_frame] for vdark_frame in dark_bluelist_out],axis=0)
    except:
        pass
    
    try:
        reddarkcube = np.stack([darkdeb_reddata_out[vdark_frame] for vdark_frame in dark_redlist_out],axis=0)
    except:
        pass
    
    try:
        rdarkcube = np.stack([darkdeb_rdata_out[vdark_frame] for vdark_frame in dark_rlist_out],axis=0)
    except:
        pass
    
    try:
        gdarkcube = np.stack([darkdeb_gdata_out[vdark_frame] for vdark_frame in dark_glist_out],axis=0)
    except:
        pass
    
    try:
        udarkcube = np.stack([darkdeb_udata_out[vdark_frame] for vdark_frame in dark_ulist_out],axis=0)
    except:
        pass
    
    try:
        idarkcube = np.stack([darkdeb_idata_out[vdark_frame] for vdark_frame in dark_ilist_out],axis=0)
    except:
        pass
    
    try:
        zdarkcube = np.stack([darkdeb_zdata_out[vdark_frame] for vdark_frame in dark_zlist_out],axis=0)
    except:
        pass
    
    try:
        lumdarkcube = np.stack([darkdeb_lumdata_out[vdark_frame] for vdark_frame in dark_lumlist_out],axis=0)
    except:
        pass
    
    try:
        Hadarkcube = np.stack([darkdeb_Hadata_out[vdark_frame] for vdark_frame in dark_Halist_out],axis=0)
    except:
        pass
    
    try:
        greendarkcube = np.stack([darkdeb_greendata_out[vdark_frame] for vdark_frame in dark_greenlist_out],axis=0)
    except:
        pass
    
    ## Show the Arrays of Dark Debiased Flats
    # plot_grid(bluedarkcube,dark_bluelist_out)
    # plot_grid(reddarkcube,dark_redlist_out)
    # plot_grid(greendarkcube,dark_greenlist_out)
    # plot_grid(rdarkcube,dark_rlist_out)
    # plot_grid(gdarkcube,dark_glist_out)
    #plot_grid(udarkcube,dark_ulist_out)
    # plot_grid(idarkcube,dark_ilist_out)
    # plot_grid(zdarkcube,dark_zlist_out)
    # plot_grid(lumdarkcube,dark_lumlist_out)
    # plot_grid(Hadarkcube,dark_Halist_out)
    
    # plt.show()
    
    
    ## Create empty Dictionaries for Dark Debiased Science Frames   
    darkdeb_blueimdata_out = {}
    darkdeb_redimdata_out = {}
    darkdeb_greenimdata_out = {}
    darkdeb_rimdata_out = {}
    darkdeb_gimdata_out = {}
    darkdeb_uimdata_out = {}
    darkdeb_iimdata_out = {}
    darkdeb_zimdata_out = {}
    darkdeb_lumimdata_out = {}
    darkdeb_Haimdata_out = {}
    
    #In [15]:
    ## Create Dark Debiased Directories for Science Frames
    ## This is probably inefficient but it will do for now...
    print("## Create Dark Debiased Directories for Science Frames")
    
    dark_blueimages_in = debias_blueimages_out
    dark_blueimages_out = ['d_' + im for im in dark_blueimages_in]
    dark_redimages_in = debias_redimages_out
    dark_redimages_out = ['d_' + im for im in dark_redimages_in]
    dark_greenimages_in = debias_greenimages_out
    dark_greenimages_out = ['d_' + im for im in dark_greenimages_in]
    dark_rimages_in = debias_rimages_out
    dark_rimages_out = ['d_' + im for im in dark_rimages_in]
    dark_gimages_in = debias_gimages_out
    dark_gimages_out = ['d_' + im for im in dark_gimages_in]
    dark_uimages_in = debias_uimages_out
    dark_uimages_out = ['d_' + im for im in dark_uimages_in]
    dark_iimages_in = debias_iimages_out
    dark_iimages_out = ['d_' + im for im in dark_iimages_in]
    dark_zimages_in = debias_zimages_out
    dark_zimages_out = ['d_' + im for im in dark_zimages_in]
    dark_lumimages_in = debias_lumimages_out
    dark_lumimages_out = ['d_' + im for im in dark_lumimages_in]
    dark_Haimages_in = debias_Haimages_out
    dark_Haimages_out = ['d_' + im for im in dark_Haimages_in]
    
    
    ## Subtract Exposure Adjusted Dark Master from Debiased Science Frames
    min=0
    
    for i in dark_blueimages_in:
        new_master_dark = {}
        new_master_dark = Expo1[i]/Expo_dark * (master_dark - master_bias)
        darkdeb_blueimdata_out['d_' + i] = debias_blueimdata_out[i] - new_master_dark 
        baselined_image = debias_blueimdata_out[i] - new_master_dark # To remove negative values
    #    baselined_image[baselined_image<0] = 0  # To remove negative valuees
        darkdeb_blueimdata_out['d_' + i] = baselined_image
        threshold_image = darkdeb_blueimdata_out['d_' + i] 
        threshold_image[threshold_image<min] = 0.00001 # To filter small valuees
        darkdeb_blueimdata_out['d_' + i] = threshold_image
    #    print(Expo1[i],Expo_dark,Expo1[i]/Expo_dark)
    #    print(new_master_dark)
    #    print('new_master_dark')
    #    print(master_dark)    
    #    print('master dark')
    #    print(master_bias)    
    #    print('master bias')
    #    print(darkdeb_blueimdata_out['d_' + i])    
    #    print(i)
    #    print(debias_blueimdata_out[i])    
    
    for i in dark_redimages_in:
        new_master_dark = {}
        new_master_dark = Expo1[i]/Expo_dark * (master_dark - master_bias)
        darkdeb_redimdata_out['d_' + i] = debias_redimdata_out[i] - new_master_dark 
        baselined_image = debias_redimdata_out[i] - new_master_dark # To remove negative values
    #    baselined_image[baselined_image<0] = 0  # To remove negative valuees
        darkdeb_redimdata_out['d_' + i] = baselined_image
        threshold_image = darkdeb_redimdata_out['d_' + i] 
        threshold_image[threshold_image<min] = 0.00001 # To filter small valuees
        darkdeb_redimdata_out['d_' + i] = threshold_image
    #    print(Expo1[i],Expo_dark,Expo1[i]/Expo_dark)
    #    print(new_master_dark)
    #    print('new_master_dark')
    #    print(master_dark)    
    #    print('master dark')
    #    print(darkdeb_redimdata_out['d_' + i])    
    #    print(i)
    #    print(debias_redimdata_out[i])    
    
    for i in dark_greenimages_in:
        new_master_dark = {}
        new_master_dark = Expo1[i]/Expo_dark * (master_dark - master_bias)
        darkdeb_greenimdata_out['d_' + i] = debias_greenimdata_out[i] - new_master_dark 
        baselined_image = debias_greenimdata_out[i] - new_master_dark # To remove negative values
    #    baselined_image[baselined_image<0] = 0  # To remove negative valuees
        darkdeb_greenimdata_out['d_' + i] = baselined_image
        threshold_image = darkdeb_greenimdata_out['d_' + i] 
        threshold_image[threshold_image<min] = 0.00001 # To filter small valuees
        darkdeb_greenimdata_out['d_' + i] = threshold_image
    #    print(Expo1[i],Expo_dark,Expo1[i]/Expo_dark)
    #    print(new_master_dark)
    #    print('new_master_dark')
    #    print(master_dark)    
    #    print('master dark')
    #    print(darkdeb_greenimdata_out['d_' + i])    
    #    print(i)
    #    print(debias_greenimdata_out[i])    
    
    for i in dark_rimages_in:
        new_master_dark = {}
        new_master_dark = Expo1[i]/Expo_dark * (master_dark - master_bias)
        darkdeb_rimdata_out['d_' + i] = debias_rimdata_out[i] - new_master_dark 
        baselined_image = debias_rimdata_out[i] - new_master_dark # To remove negative values
    #    baselined_image[baselined_image<0] = 0  # To remove negative valuees
        darkdeb_rimdata_out['d_' + i] = baselined_image
        threshold_image = darkdeb_rimdata_out['d_' + i] 
        threshold_image[threshold_image<min] = 0.00001 # To filter small valuees
        darkdeb_rimdata_out['d_' + i] = threshold_image
        print(Expo1[i],Expo_dark,Expo1[i]/Expo_dark)
        print(new_master_dark)
        print('new_master_dark')
        print(master_dark)    
        print('master dark')
        print(darkdeb_rimdata_out['d_' + i])    
        print(i)
        print(debias_rimdata_out[i])    
    
    for i in dark_gimages_in:
        new_master_dark = {}
        new_master_dark = Expo1[i]/Expo_dark * (master_dark - master_bias)
        darkdeb_gimdata_out['d_' + i] = debias_gimdata_out[i] - new_master_dark 
        baselined_image = debias_gimdata_out[i] - new_master_dark # To remove negative values
    #    baselined_image[baselined_image<0] = 0  # To remove negative valuees
        darkdeb_gimdata_out['d_' + i] = baselined_image
        threshold_image = darkdeb_gimdata_out['d_' + i] 
        threshold_image[threshold_image<min] = 0.00001 # To filter small valuees
        darkdeb_gimdata_out['d_' + i] = threshold_image
    #    print(Expo1[i],Expo_dark,Expo1[i]/Expo_dark)
    #    print(new_master_dark)
    #    print('new_master_dark')
    #    print(master_dark)    
    #    print('master dark')
    #    print(darkdeb_gimdata_out['d_' + i])    
    #    print(i)
    #    print(debias_gimdata_out[i])    
    
    #$$$ Note that this one is different from the other!
    for i in dark_uimages_in:
        new_master_dark = {}
        new_master_dark = Expo1[i]/Expo_dark * (master_dark - master_bias)
        darkdeb_uimdata_out['d_' + i] = debias_uimdata_out[i] - new_master_dark 
    #    baselined_image = debias_uimdata_out[i] - new_master_dark # To remove negative values
    #    baselined_image[baselined_image<0] = 0.1 # To remove negative valuees
    #    darkdeb_uimdata_out['d_' + i] = baselined_image
        threshold_image = darkdeb_uimdata_out['d_' + i] 
        threshold_image[threshold_image<min] = 0.01 # To filter small valuees
        darkdeb_uimdata_out['d_' + i] = threshold_image
        print(Expo1[i],Expo_dark,Expo1[i]/Expo_dark)
        print(new_master_dark)
        print('new_master_dark')
        print(master_dark)    
        print('master dark')
        print(darkdeb_uimdata_out['d_' + i])    
        print(['d_' + i])
        print(debias_uimdata_out[i])    
        print(i)
        print('Reminder',Expo1[i],Expo_dark,Expo1[i]/Expo_dark)
    
    for i in dark_iimages_in:
        new_master_dark = {}
        new_master_dark = Expo1[i]/Expo_dark * (master_dark - master_bias)
        darkdeb_iimdata_out['d_' + i] = debias_iimdata_out[i] - new_master_dark 
        baselined_image = debias_iimdata_out[i] - new_master_dark # To remove negative values
    #    baselined_image[baselined_image<0] = 0  # To remove negative valuees
        darkdeb_iimdata_out['d_' + i] = baselined_image
        threshold_image = darkdeb_iimdata_out['d_' + i] 
        threshold_image[threshold_image<min] = 0.00001 # To filter small valuees
        darkdeb_iimdata_out['d_' + i] = threshold_image
    #    print(Expo1[i],Expo_dark,Expo1[i]/Expo_dark)
    #    print(new_master_dark)
    #    print('new_master_dark')
    #    print(master_dark)    
    #    print('master dark')
    #    print(darkdeb_iimdata_out['d_' + i])    
    #    print(i)
    #    print(debias_iimdata_out[i])    
    
    for i in dark_zimages_in:
        new_master_dark = {}
        new_master_dark = Expo1[i]/Expo_dark * (master_dark - master_bias)
        darkdeb_zimdata_out['d_' + i] = debias_zimdata_out[i] - new_master_dark 
        baselined_image = debias_zimdata_out[i] - new_master_dark # To remove negative values
    #    baselined_image[baselined_image<0] = 0  # To remove negative valuees
        darkdeb_zimdata_out['d_' + i] = baselined_image
        threshold_image = darkdeb_zimdata_out['d_' + i] 
        threshold_image[threshold_image<min] = 0.00001 # To filter small valuees
        darkdeb_zimdata_out['d_' + i] = threshold_image
    #    print(Expo1[i],Expo_dark,Expo1[i]/Expo_dark)
    #    print(new_master_dark)
    #    print('new_master_dark')
    #    print(master_dark)    
    #    print('master dark')
    #    print(darkdeb_zimdata_out['d_' + i])    
    #    print(i)
    #    print(debias_zimdata_out[i])    
    
    for i in dark_lumimages_in:
        new_master_dark = {}
        new_master_dark = Expo1[i]/Expo_dark * (master_dark - master_bias)
        darkdeb_lumimdata_out['d_' + i] = debias_lumimdata_out[i] - new_master_dark 
        baselined_image = debias_lumimdata_out[i] - new_master_dark # To remove negative values
    #    baselined_image[baselined_image<0] = 0  # To remove negative valuees
        darkdeb_lumimdata_out['d_' + i] = baselined_image
        threshold_image = darkdeb_lumimdata_out['d_' + i] 
        threshold_image[threshold_image<min] = 0.00001 # To filter small valuees
        darkdeb_lumimdata_out['d_' + i] = threshold_image
    #    print(Expo1[i],Expo_dark,Expo1[i]/Expo_dark)
    #    print(new_master_dark)
    #    print('new_master_dark')
    #    print(master_dark)    
    #    print('master dark')
    #    print(darkdeb_lumimdata_out['d_' + i])    
    #    print(i)
    #    print(debias_lumimdata_out[i])    
    
    for i in dark_Haimages_in:
        new_master_dark = {}
        new_master_dark = Expo1[i]/Expo_dark * (master_dark - master_bias)
        baselined_image = debias_Haimdata_out[i] - new_master_dark # To filter negative values
        baselined_image[baselined_image < 0] = 0.0001  # To filter negative values
        darkdeb_Haimdata_out['d_' + i] = baselined_image
        threshold_image = darkdeb_Haimdata_out['d_' + i] 
        threshold_image[threshold_image<min] = 0.00001 # To filter small valuees
        darkdeb_Haimdata_out['d_' + i] = threshold_image
    #     print(Expo1[i],Expo_dark,Expo1[i]/Expo_dark)
    #     print(new_master_dark)
    #     print('new_master_dark')
    #     print(master_dark)    
    #     print('master dark')
    #     print(master_bias)    
    #     print('master bias')
    #     print(darkdeb_Haimdata_out['d_' + i])    
    #     print(i)
    #     print(debias_Haimdata_out[i])    
    
    ## Create Arrays of Dark Debiased Science Images
    try:
        blueimdarkcube = np.stack([darkdeb_blueimdata_out[vdark_frame] for vdark_frame in dark_blueimages_out],axis=0)
    except:
        pass
    
    try:
        redimdarkcube = np.stack([darkdeb_redimdata_out[vdark_frame] for vdark_frame in dark_redimages_out],axis=0)
    except:
        pass
    
    try:
        rimdarkcube = np.stack([darkdeb_rimdata_out[vdark_frame] for vdark_frame in dark_rimages_out],axis=0)
    except:
        pass
    
    try:
        gimdarkcube = np.stack([darkdeb_gimdata_out[vdark_frame] for vdark_frame in dark_gimages_out],axis=0)
    except:
        pass
    
    try:
        uimdarkcube = np.stack([darkdeb_uimdata_out[vdark_frame] for vdark_frame in dark_uimages_out],axis=0)
    except:
        pass
    
    try:
        iimdarkcube = np.stack([darkdeb_iimdata_out[vdark_frame] for vdark_frame in dark_iimages_out],axis=0)
    except:
        pass
    
    try:
        zimdarkcube = np.stack([darkdeb_zimdata_out[vdark_frame] for vdark_frame in dark_zimages_out],axis=0)
    except:
        pass
    
    try:
        lumimdarkcube = np.stack([darkdeb_lumimdata_out[vdark_frame] for vdark_frame in dark_lumimages_out],axis=0)
    except:
        pass
    
    try:
        Haimdarkcube = np.stack([darkdeb_Haimdata_out[vdark_frame] for vdark_frame in dark_Haimages_out],axis=0)
    except:
        pass
    
    try:
        greenimdarkcube = np.stack([darkdeb_greenimdata_out[vdark_frame] for vdark_frame in dark_greenimages_out],axis=0)
    except:
        pass 
    
    ## Show the Arrays of Dark Debiased Science Images
    # plot_grid(blueimdarkcube,dark_blueimages_out)
    # plot_grid(redimdarkcube,dark_redimages_out)
    # plot_grid(greenimdarkcube,dark_greenimages_out)
    # plot_grid(rimdarkcube,dark_rimages_out)
    # plot_grid(gimdarkcube,dark_gimages_out)
    # plot_grid(uimdarkcube,dark_uimages_out)
    # plot_grid(iimdarkcube,dark_iimages_out)
    # plot_grid(zimdarkcube,dark_zimages_out)
    # plot_grid(lumimdarkcube,dark_lumimages_out)
    # plot_grid(Haimdarkcube,dark_Haimages_out)
    
    # plt.show()
    
    ## Create Reduced Science Images - Divide Dark Debiased Images by Dark Debiased Flats
    ## Create Normalized Master Flats
    #  Create Master Flats by Filter Based upon Median (Why Median?)
    try:
        master_flat_blue = np.median(bluedarkcube, axis=0) ## to median combine them 
    except:
        pass
    
    try:
        master_flat_red = np.median(reddarkcube, axis=0) ## to median combine them 
    except:
        pass
    
    try:
        master_flat_green = np.median(greendarkcube, axis=0) ## to median combine them 
    except:
        pass
    
    try:
        master_flat_r = np.median(rdarkcube, axis=0) ## to median combine them 
    except:
        pass
    
    try:
        master_flat_g = np.median(gdarkcube, axis=0) ## to median combine them 
    except:
        pass
    
    try:
        master_flat_u = np.median(udarkcube, axis=0) ## to median combine them 
    except:
        pass
    
    try:
        master_flat_i = np.median(idarkcube, axis=0) ## to median combine them 
    except:
        pass
    
    try:
        master_flat_z = np.median(zdarkcube, axis=0) ## to median combine them 
    except:
        pass
    
    try:
        master_flat_lum = np.median(lumdarkcube, axis=0) ## to median combine them 
    except:
        pass
    
    try:
        master_flat_Ha = np.median(Hadarkcube, axis=0) ## to median combine them 
    except:
        pass
    
    # Create Master Normalized Flats By Filter 
    normalized_master_flat_red = master_flat_red/np.median(master_flat_red)
    normalized_master_flat_blue = master_flat_blue/np.median(master_flat_blue)
    normalized_master_flat_green = master_flat_green/np.median(master_flat_green)
    normalized_master_flat_r = master_flat_r/np.median(master_flat_r)
    # print(master_flat_r)
    # print(np.median(master_flat_r))
    # print(normalized_master_flat_r)
    # plt.title('Master Flat r')
    # plt.figure(figsize=(15,15)) 
    # plt.imshow(master_flat_u, origin='lower', cmap='gray');
    # plt.title('Master Flat r')
    # plt.show()
    # print(normalized_master_flat_r)
    # plt.title('Normalized Master Flat r')
    # plt.figure(figsize=(15,15)) 
    # plt.imshow(normalized_master_flat_r, origin='lower', cmap='gray');
    # plt.title('Normalized Master Flat r')
    # plt.show()
    normalized_master_flat_g = master_flat_g/np.median(master_flat_g)
    normalized_master_flat_i = master_flat_i/np.median(master_flat_i)
    normalized_master_flat_u = master_flat_u/np.median(master_flat_u)
    # print(master_flat_u)
    # print(np.median(master_flat_u))
    # print(normalized_master_flat_u)
    normalized_master_flat_z = master_flat_z/np.median(master_flat_z)
    normalized_master_flat_lum = master_flat_lum/np.median(master_flat_lum)
    normalized_master_flat_Ha = master_flat_Ha/np.median(master_flat_Ha)
    # print(master_flat_u)
    # print(np.median(master_flat_u))
    # plt.title('Master Flat u')
    # plt.figure(figsize=(15,15)) 
    # plt.imshow(master_flat_u, origin='lower', cmap='gray');
    # plt.title('Master Flat u')
    # plt.show()
    # print(normalized_master_flat_u)
    # plt.title('Normalized Master Flat u')
    # plt.figure(figsize=(15,15)) 
    # plt.imshow(normalized_master_flat_u, origin='lower', cmap='gray');
    # plt.title('Normalized Master Flat u')
    # plt.show()
    
    ## Perform Final Stage of Image Reduction
    # Create Empty Directories for Reduced Images
    reduced_blueimdata_out = {}
    reduced_redimdata_out = {}
    reduced_greenimdata_out = {}
    reduced_rimdata_out = {}
    reduced_gimdata_out = {}
    reduced_uimdata_out = {}
    reduced_iimdata_out = {}
    reduced_zimdata_out = {}
    reduced_lumimdata_out = {}
    reduced_Haimdata_out = {}
    
    # Create List of Files (Check)
    flat_blueimages_in = dark_blueimages_out
    reduced_blueimages_out = ['f_' + im for im in flat_blueimages_in]
    flat_redimages_in = dark_redimages_out
    reduced_redimages_out = ['f_' + im for im in flat_redimages_in]
    flat_greenimages_in = dark_greenimages_out
    reduced_greenimages_out = ['f_' + im for im in flat_greenimages_in]
    flat_rimages_in = dark_rimages_out
    reduced_rimages_out = ['f_' + im for im in flat_rimages_in]
    flat_gimages_in = dark_gimages_out
    reduced_gimages_out = ['f_' + im for im in flat_gimages_in]
    flat_iimages_in = dark_iimages_out
    reduced_iimages_out = ['f_' + im for im in flat_iimages_in]
    # print(reduced_iimages_out)
    flat_uimages_in = dark_uimages_out
    reduced_uimages_out = ['f_' + im for im in flat_uimages_in]
    # print(reduced_uimages_out)
    flat_zimages_in = dark_zimages_out
    reduced_zimages_out = ['f_' + im for im in flat_zimages_in]
    flat_lumimages_in = dark_lumimages_out
    reduced_lumimages_out = ['f_' + im for im in flat_lumimages_in]
    flat_Haimages_in = dark_Haimages_out
    reduced_Haimages_out = ['f_' + im for im in flat_Haimages_in]
    
    # Create Reduced Images # CHECK FORMULA!
    for i in dark_blueimages_out:
        reduced_blueimdata_out['f_' + i] = darkdeb_blueimdata_out[i]/normalized_master_flat_blue 
    
    for i in dark_redimages_out:
        reduced_redimdata_out['f_' + i] = darkdeb_redimdata_out[i]/normalized_master_flat_red 
    
    for i in dark_greenimages_out:
        reduced_greenimdata_out['f_' + i] = darkdeb_greenimdata_out[i]/normalized_master_flat_green 
    
    for i in dark_rimages_out:
        reduced_rimdata_out['f_' + i] = darkdeb_rimdata_out[i]/normalized_master_flat_r 
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(reduced_rimdata_out['f_' + i])
            print(i)
            print(mean, median, std)
    #    print(darkdeb_rimdata_out[i])    
    #    print('Before')
    #    print(reduced_rimdata_out['f_' + i])    
    #    print('After')
    #    print(master_flat_r)
    #    print('master flat r')  
    #    print(normalized_master_flat_r)
    #    print('normalized master flat r')  
    
    ## Looks at the original image
    for i in dirlist1_rimages:
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(raw_image_data[i])
            print(i)
            print(mean, median, std)
    
    for i in dark_gimages_out:
        reduced_gimdata_out['f_' + i] = darkdeb_gimdata_out[i]/normalized_master_flat_g 
        # if(TOGGLESIGMA == 'Yes'):
    #    mean, median, std = sigma_clipped_stats(reduced_gimdata_out['f_' + i])
    #    print(mean, median, std)
    #    print(darkdeb_gimdata_out[i])    
    #    print('Before')
    #    print(reduced_gimdata_out['f_' + i])    
    #    print('After')
    #    print(master_flat_g)
    #    print('master flat g')  
    #    print(normalized_master_flat_g)
    #    print('normalized master flat g')  
    
    for i in dark_iimages_out:
        reduced_iimdata_out['f_' + i] = darkdeb_iimdata_out[i]/normalized_master_flat_i 
        # if(TOGGLESIGMA == 'Yes'):
    #    mean, median, std = sigma_clipped_stats(reduced_iimdata_out['f_' + i])
    #    print(mean, median, std)
    
    for i in dark_uimages_out:
        reduced_uimdata_out['f_' + i] = darkdeb_uimdata_out[i]/normalized_master_flat_u 
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(reduced_uimdata_out['f_' + i])
            print(i)
            print(mean, median, std)
        print(darkdeb_uimdata_out[i])    
        print('Before')
        print(reduced_uimdata_out['f_' + i])    
        print('After')
        print(master_flat_u)
        print('master flat u')  
        print(normalized_master_flat_u)
        print('normalized master flat u')  
    
    ## Looks at the original image
    for i in dirlist1_uimages:
        if(TOGGLESIGMA == 'Yes'):
            mean, median, std = sigma_clipped_stats(raw_image_data[i])
            print(i)
            print(mean, median, std)
    
    for i in dark_zimages_out:  
        reduced_zimdata_out['f_' + i] = darkdeb_zimdata_out[i]/normalized_master_flat_z 
        #if(TOGGLESIGMA == 'Yes'):
    #    mean, median, std = sigma_clipped_stats(reduced_zimdata_out['f_' + i])
    #    print(mean, median, std)
    
    for i in dark_lumimages_out:
        reduced_lumimdata_out['f_' + i] = darkdeb_lumimdata_out[i]/normalized_master_flat_lum 
    
    for i in dark_Haimages_out:
        reduced_Haimdata_out['f_' + i] = darkdeb_Haimdata_out[i]/normalized_master_flat_Ha 
    #    print(darkdeb_Haimdata_out[i])    
    #    print('Before')
    #    print(reduced_Haimdata_out['f_' + i])    
    #    print('After')
    #    print(master_flat_Ha)
    #    print('master flat Ha')  
    #    print(normalized_master_flat_Ha)
    #    print('normalized master flat Ha    `  
        
    ## Create Arrays of Reduced Science Images
    try:
        blueredimdarkcube = np.stack([reduced_blueimdata_out[reduced_frame] for reduced_frame in reduced_blueimages_out],axis=0)
    except:
        pass
    
    try:
        redredimdarkcube = np.stack([reduced_redimdata_out[reduced_frame] for reduced_frame in reduced_redimages_out],axis=0)
    except:
        pass
    
    try:
        rredimdarkcube = np.stack([reduced_rimdata_out[reduced_frame] for reduced_frame in reduced_rimages_out],axis=0)
    except:
        pass
    
    try:
        gredimdarkcube = np.stack([reduced_gimdata_out[reduced_frame] for reduced_frame in reduced_gimages_out],axis=0)
    except:
        pass
    
    try:
        uredimdarkcube = np.stack([reduced_uimdata_out[reduced_frame] for reduced_frame in reduced_uimages_out],axis=0)
    except:
        pass
    
    try:
        iredimdarkcube = np.stack([reduced_iimdata_out[reduced_frame] for reduced_frame in reduced_iimages_out],axis=0)
    except:
        pass
    
    try:
        zredimdarkcube = np.stack([reduced_zimdata_out[reduced_frame] for reduced_frame in reduced_zimages_out],axis=0)
    except:
        pass
    
    try:
        lumredimdarkcube = np.stack([reduced_lumimdata_out[reduced_frame] for reduced_frame in reduced_lumimages_out],axis=0)
    except:
        pass
    
    try:
        Haredimdarkcube = np.stack([reduced_Haimdata_out[reduced_frame] for reduced_frame in reduced_Haimages_out],axis=0)
    except:
        pass
    
    try:    
        greenredimdarkcube = np.stack([reduced_greenimdata_out[reduced_frame] for reduced_frame in reduced_greenimages_out],axis=0)
    except:
        pass
        
    ## Show the Arrays of Reduced Science Images
    # plot_grid(blueredimdarkcube,reduced_blueimages_out)
    # plot_grid(redredimdarkcube,reduced_redimages_out)
    # plot_grid(greenredimdarkcube,reduced_greenimages_out)
    # plot_grid(rredimdarkcube,reduced_rimages_out)
    # plot_grid(gredimdarkcube,reduced_gimages_out)
    # plot_grid(uredimdarkcube,reduced_uimages_out)
    # plot_grid(iredimdarkcube,reduced_iimages_out)
    # plot_grid(zredimdarkcube,reduced_zimages_out)
    # plot_grid(lumredimdarkcube,reduced_lumimages_out)
    # plot_grid(Haredimdarkcube,reduced_Haimages_out)
    
    # plt.show()
    
    ## This section writes data to files
    
    #  In case path definitions are required for troubleshooting
    # path = '/python/filesnew'
    # calibrationpath = '/python/CalibrationNew'
    # flatspath = '/python/FlatsNew'
    
    # This section creates the Processed Directory
    print("# This section creates the Processed Directory")
    
    
    processedpath =  '/home/david/ZowadaTransients/Processed'
    
    processeddir = processedpath
    print(processeddir)
                   
    if not os.path.isdir(processeddir):
        os.mkdir(processeddir)
        
    # Again, in case path changes are required
    # os.chdir(path)
    # print(path)
           
    #$$$ This section identifies the Target Object. Updated code needs to be tested!
    Target={}
    
    for i in dirlist1_images: 
        fitsdata = {}
        fitsdata[i]=fits.getdata(i)
        hdulist = fits.open(i)
        data = hdulist[0].data.astype(float)
        header = hdulist[0].header
    #    Exposure = hdulist[0].header['EXPTIME']
        TargetObject = hdulist[0].header['OBJECT'] 
    #    Expo[i]=Exposure
    #    Target[i]=TargetObject
        print(i)
    #    print(Exposure)
        print(TargetObject)
    
    
    numrfiles=len(dirlist1_rimages)
    # print(numrfiles)
    numgfiles=len(dirlist1_gimages)
    # print(numgfiles)
    numufiles=len(dirlist1_uimages)
    # print(numufiles)
    numifiles=len(dirlist1_iimages)
    # print(numifiles)
    numzfiles=len(dirlist1_zimages)
    # print(numzfiles)
    numredfiles=len(dirlist1_redimages)
    # print(numredfiles)
    numgreenfiles=len(dirlist1_greenimages)
    # print(numgreenfiles)
    numbluefiles=len(dirlist1_blueimages)
    # print(numbluefiles)
    numHafiles=len(dirlist1_Haimages)
    # print(numHafiles)
    numlumfiles=len(dirlist1_lumimages)
    # print(numlumfiles)
    
    # This section creates the Target Directory
    targetdir = processedpath +'/'+ str(TargetObject) +'/'+date
    print(targetdir)
    
                   
    if not os.path.isdir(processedpath +'/'+ str(TargetObject)):
        os.mkdir(processedpath +'/'+ str(TargetObject))
    if not os.path.isdir(targetdir):
        os.mkdir(targetdir)
        
    
    import shutil
    
    # This section creates the "writeto" files
    # This may or may not be necessary
    for i in range(numrfiles):
        source = dirlist1_rimages[i]
        source
        name = str(source)
        name1 = name.replace(path,'')
        name1
        destination = targetdir + '/' + name1
        shutil.copy2(source, destination)
    
    
    for i in range(numufiles):
        source = dirlist1_uimages[i]
        source
        name = str(source)
        name1 = name.replace(path,'')
        name1
        destination = targetdir + '/' + name1
        shutil.copy2(source, destination)
    
        
    for i in range(numifiles):
        source = dirlist1_iimages[i]
        source
        name = str(source)
        name1 = name.replace(path,'')
        name1
        destination = targetdir + '/' + name1
        shutil.copy2(source, destination)
    
    
    for i in range(numzfiles):
        source = dirlist1_zimages[i]
        source
        name = str(source)
        name1 = name.replace(path,'')
        name1
        destination = targetdir + '/' + name1
        shutil.copy2(source, destination)
    
    
    for i in range(numgfiles):
        source = dirlist1_gimages[i]
        source
        name = str(source)
        name1 = name.replace(path,'')
        name1
        destination = targetdir + '/' + name1
        shutil.copy2(source, destination)
        
    
    for i in range(numredfiles):
        source = dirlist1_redimages[i]
        source
        name = str(source)
        name1 = name.replace(path,'')
        name1
        destination = targetdir + '/' + name1
        shutil.copy2(source, destination)
        
    
    for i in range(numbluefiles):
        source = dirlist1_blueimages[i]
        source
        name = str(source)
        name1 = name.replace(path,'')
        name1
        destination = targetdir + '/' + name1
        shutil.copy2(source, destination)
    
        
    for i in range(numgreenfiles):
        source = dirlist1_greenimages[i]
        source
        name = str(source)
        name1
        name1 = name.replace(path,'')
        destination = targetdir + '/' + name1
        shutil.copy2(source, destination)
        
    
    for i in range(numHafiles):
        source = dirlist1_Haimages[i]
        source
        name = str(source)
        name1 = name.replace(path,'')
        name1
        destination = targetdir + '/' + name1
        shutil.copy2(source, destination)
       
    
    for i in range(numlumfiles):
        source = dirlist1_lumimages[i]
        source
        name = str(source)
        name1 = name.replace(path,'')
        name1
        destination = targetdir + '/' + name1
        shutil.copyfile(source, destination)
        
        
    # This section creates the directory of "writeto" files for the Target Directory
    '''
    dirlist_rproc=glob.glob(targetdir + '/' + '*-r *.fit')+glob.glob(targetdir + '/' + '*-r *.fts')+glob.glob(targetdir + '/' + '*-r-*.fit')+glob.glob(targetdir + '/' + '*-r-*.fts')
    dirlist1_rproc=dirlist_rproc+glob.glob(targetdir + '/' + '*-r *.fits')+glob.glob(targetdir + '/' + '*-r-*.fits')
    # print(dirlist1_rproc)
    
    dirlist_gproc=glob.glob(targetdir + '/' + '*-g *.fit')+glob.glob(targetdir + '/' + '*-g *.fts')+glob.glob(targetdir + '/' + '*-g-*.fit')+glob.glob(targetdir + '/' + '*-g-*.fts')
    dirlist1_gproc=dirlist_gproc+glob.glob(targetdir + '/' + '*-g *.fits')+glob.glob(targetdir + '/' + '*-g-*.fits')
    # print(dirlist1_gproc)
    
    dirlist_uproc=glob.glob(targetdir + '/' + '*-u*.fit')+glob.glob(targetdir + '/' + '*-u*.fts')+glob.glob(targetdir + '/' + '*-u-*.fit')+glob.glob(targetdir + '/' + '*-u-*.fts')
    dirlist1_uproc=dirlist_uproc+glob.glob(targetdir + '/' + '*-u*.fits')+glob.glob(targetdir + '/' + '*-u-*.fits')
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
    
    #$$$ These are inconsistent because of limited testing past this point.  Some of the ordering and file names may be off because of testing runs
    # At last, we write to the Target Directory
    # We print out data (at least for now) because of the 32768 error
    
    for i in range(numgfiles):
        os.chdir(path)
        HDUList=fits.open(dirlist1_gimages[i])
        primaryHeader=HDUList[0].header
        HDUList.close()
        file=reduced_gimages_out[i]
        print(file)
        print(path)
    #   The 32768 fudge is to remove an unknown error adding 32768 to all values
    #    procData=reduced_gimdata_out[file] - reduced_gimdata_out[file]
    #    procData=procData.astype(int)
        procData = reduced_gimdata_out[file]# - 32768
        procData = np.uint16(procData)
        print(i)
        if(TOGGLESIGMA == 'Yes'):
           mean, median, std = sigma_clipped_stats(procData)
           print(mean, median, std)
        print(procData)
        os.chdir(targetdir)
        procHDU=fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Reduced and Processed Image')
        procHDU.writeto(dirlist1_gproc[i],overwrite=True)
        os.chdir(path)
    
    newraw_image_data = {}
    for image_name in dirlist1_gproc: 
        newraw_image_data[image_name] = fits.getdata(image_name)
        postproc = newraw_image_data[image_name]
        if(TOGGLESIGMA == 'Yes'):
           mean, median, std = sigma_clipped_stats(postproc)
           print(mean, median, std)
        print(i)
        print(postproc)
    
    for i in range(numrfiles):
        os.chdir(path)
        HDUList=fits.open(dirlist1_rimages[i])
        primaryHeader=HDUList[0].header
        HDUList.close()
        file=reduced_rimages_out[i]
        print(file)
        print(path)
    #   The 32768 fudge is to remove an unknown error adding 32768 to all values
        print(reduced_rimdata_out[file])
        procData=reduced_rimdata_out[file]# - 32768
        procData = np.uint16(procData)
        if(TOGGLESIGMA == 'Yes'):
           mean, median, std = sigma_clipped_stats(procData)
           print(mean, median, std)
        print(procData)
        os.chdir(targetdir)
        procHDU=fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Reduced and Processed Image')
        procHDU.writeto(dirlist1_rproc[i],overwrite=True)
        os.chdir(path)
    
    
    for i in range(numufiles):
        os.chdir(path)
        HDUList=fits.open(dirlist1_uimages[i])
        primaryHeader=HDUList[0].header
        HDUList.close()
        file=reduced_uimages_out[i]
        print(file)
        print(path)
    #   The 32768 fudge is to remove an unknown error adding 32768 to all values
        print(reduced_uimdata_out[file])
        procData=reduced_uimdata_out[file]# - 32768
        procData = np.uint16(procData)
        if(TOGGLESIGMA == 'Yes'):
           mean, median, std = sigma_clipped_stats(procData)
           print(mean, median, std)
        print(procData)
        os.chdir(targetdir)
        procHDU=fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Reduced and Processed Image')
        procHDU.writeto(dirlist1_uproc[i],overwrite=True)
        os.chdir(path)
    
    ## We check to be sure that the exposures are the correct size and that before and after data ARE THE SAME!!  These are test print-outs because of the 32768 Error
    #newraw_image_data = {}
    #for image_name in dirlist1_uproc: newraw_image_data[image_name] = fits.getdata(image_name)
    
    #for image in dirlist1_uproc: print (image, newraw_image_data[image].shape)
    
    #for i in range(len(dirlist1_uproc)):  
    #    print(newraw_image_data[dirlist1_uproc[i]])
    
    #for i in range(numufiles):
    #    file=reduced_uimages_out[i]
    #    print(file)
    #    print(reduced_uimdata_out[file])
    
    for i in range(numifiles):
        os.chdir(path)
        HDUList=fits.open(dirlist1_iimages[i])
        primaryHeader=HDUList[0].header
        HDUList.close()
        file=reduced_iimages_out[i]
        print(file)
        print(path)
    #   The 32768 fudge is to remove an unknown error adding 32768 to all values
        procData=reduced_iimdata_out[file]# -32768
        procData = np.uint16(procData)
        print(procData)
        if(TOGGLESIGMA == 'Yes'):
           mean, median, std = sigma_clipped_stats(procData)
           print(mean, median, std)
        os.chdir(targetdir)
        procHDU=fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Reduced and Processed Image')
        procHDU.writeto(dirlist1_iproc[i],overwrite=True)
        os.chdir(path)
    
    for i in range(numzfiles):
        os.chdir(path)
        HDUList=fits.open(dirlist1_zimages[i])
        primaryHeader=HDUList[0].header
        HDUList.close()
        file=reduced_zimages_out[i]
        print(file)
        print(path)
    #   The 32768 fudge is to remove an unknown error adding 32768 to all values
        procData=reduced_zimdata_out[file]#-32768
        procData = np.uint16(procData)
        print(procData)
        if(TOGGLESIGMA == 'Yes'):
           mean, median, std = sigma_clipped_stats(procData)
           print(mean, median, std)
        os.chdir(targetdir)
        procHDU=fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Reduced and Processed Image')
        procHDU.writeto(dirlist1_zproc[i],overwrite=True)
        os.chdir(path)
    
    for i in range(numredfiles):
        os.chdir(path)
        HDUList=fits.open(dirlist1_redimages[i])
        primaryHeader=HDUList[0].header
        HDUList.close()
        file=reduced_redimages_out[i]
        print(file)
        print(path)
    #   The 32768 fudge is to remove an unknown error adding 32768 to all values
        procData=reduced_redimdata_out[file]#-32768
        procData = np.uint16(procData)
        print(procData)
        if(TOGGLESIGMA == 'Yes'):
           mean, median, std = sigma_clipped_stats(procData)
           print(mean, median, std)
        os.chdir(targetdir)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Reduced and Processed Image')
        procHDU.writeto(dirlist1_redproc[i],overwrite=True)
        os.chdir(path)
    
    for i in range(numbluefiles):
        os.chdir(path)
        HDUList=fits.open(dirlist1_blueimages[i])
        primaryHeader=HDUList[0].header
        HDUList.close()
        file=reduced_blueimages_out[i]
        print(file)
        print(path)
    #   The 32768 fudge is to remove an unknown error adding 32768 to all values
        procData=reduced_blueimdata_out[file]#-32768
        procData = np.uint16(procData)
        print(procData)
        if(TOGGLESIGMA == 'Yes'):
           mean, median, std = sigma_clipped_stats(procData)
           print(mean, median, std)
        os.chdir(targetdir)
        procHDU=fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Reduced and Processed Image')
        procHDU.writeto(dirlist1_blueproc[i],overwrite=True)
        os.chdir(path)
    
    for i in range(numgreenfiles):
        os.chdir(path)
        HDUList=fits.open(dirlist1_greenimages[i])
        primaryHeader=HDUList[0].header
        HDUList.close()
        file=reduced_greenimages_out[i]
        print(file)
        print(path)
    #   The 32768 fudge is to remove an unknown error adding 32768 to all values
        procData=reduced_greenimdata_out[file]#-32768
        procData = np.uint16(procData)
        print(procData)
        if(TOGGLESIGMA == 'Yes'):
           mean, median, std = sigma_clipped_stats(procData)
           print(mean, median, std)
        os.chdir(targetdir)
        procHDU=fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Reduced and Processed Image')
        procHDU.writeto(dirlist1_greenproc[i],overwrite=True)
        os.chdir(path)
    
    for i in range(numHafiles):
        os.chdir(path)
        HDUList=fits.open(dirlist1_Haimages[i])
        primaryHeader=HDUList[0].header
        HDUList.close()
        file=reduced_Haimages_out[i]
        print(file)
        print(path)
    #   The 32768 fudge is to remove an unknown error adding 32768 to all values
        procData=reduced_Haimdata_out[file]#-32768
        procData = np.uint16(procData)
        print(procData)
        if(TOGGLESIGMA == 'Yes'):
           mean, median, std = sigma_clipped_stats(procData)
           print(mean, median, std)
        os.chdir(targetdir)
        procHDU=fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Reduced and Processed Image')
        procHDU.writeto(dirlist1_Haproc[i],overwrite=True)
        os.chdir(path)
    
    for i in range(numlumfiles):
        os.chdir(path)
        HDUList=fits.open(dirlist1_lumimages[i])
        primaryHeader=HDUList[0].header
        HDUList.close()
        file=reduced_lumimages_out[i]
        print(file)
        print(path)
    #   The 32768 fudge is to remove an unknown error adding 32768 to all values
        procData=reduced_lumimdata_out[file]#-32768
        procData = np.uint16(procData)
        print(procData)
        if(TOGGLESIGMA == 'Yes'):
           mean, median, std = sigma_clipped_stats(procData)
           print(mean, median, std)
        os.chdir(targetdir)
        procHDU=fits.PrimaryHDU(procData)
        procHDU.header=primaryHeader
        procHDU.header.add_history('Reduced and Processed Image')
        procHDU.writeto(dirlist1_lumproc[i],overwrite=True)
        os.chdir(path)
    
    ##  Note that prior Alignment Code is available in file Reduction_Alignment Module Only_works.py
    
    ## First, we create Path to find Processed Files
    # path = '/python/processed' + '/' + target
    # os.chdir(path)
    # print(path)
    
    ## WORKS THROUGH HERE!! $$$ Need to Test from Here forward with updated code
    
    
    ##  These lines confirm that Processed Files are of correct shape
    ## We check to be sure that the exposures are the correct size.
        
    ##  THIS IS INITIALIZATION FOR TEST PURPOSES
    
    # This section creates the directory of "writeto" files for the Target Directory
    '''
    dirlist_rproc=glob.glob(targetdir + '/' + '*-r *.fit')+glob.glob(targetdir + '/' + '*-r *.fts')+glob.glob(targetdir + '/' + '*-r-*.fit')+glob.glob(targetdir + '/' + '*-r-*.fts')
    dirlist1_rproc=dirlist_rproc+glob.glob(targetdir + '/' + '*-r *.fits')+glob.glob(targetdir + '/' + '*-r-*.fits')
    # print(dirlist1_rproc)
    
    dirlist_gproc=glob.glob(targetdir + '/' + '*-g *.fit')+glob.glob(targetdir + '/' + '*-g *.fts')+glob.glob(targetdir + '/' + '*-g-*.fit')+glob.glob(targetdir + '/' + '*-g-*.fts')
    dirlist1_gproc=dirlist_gproc+glob.glob(targetdir + '/' + '*-g *.fits')+glob.glob(targetdir + '/' + '*-g-*.fits')
    # print(dirlist1_gproc)
    
    dirlist_uproc=glob.glob(targetdir + '/' + '*-u*.fit')+glob.glob(targetdir + '/' + '*-u*.fts')+glob.glob(targetdir + '/' + '*-u-*.fit')+glob.glob(targetdir + '/' + '*-u-*.fts')
    dirlist1_uproc=dirlist_uproc+glob.glob(targetdir + '/' + '*-u*.fits')+glob.glob(targetdir + '/' + '*-u-*.fits')
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
    
    ##  This is to Display Frames Before Alignment and Average or Summation
    try:
        rimagecube = np.stack([newraw_image_data[science_frame] for science_frame in dirlist1_rproc],axis=0)
        if(TOGGLESIGMA == 'Yes'):
            plot_grid(rimagecube,dirlist1_rproc)
    except:
        pass
    
    try:
        gimagecube = np.stack([newraw_image_data[science_frame] for science_frame in dirlist1_gproc],axis=0)
        if(TOGGLESIGMA == 'Yes'):
            plot_grid(gimagecube,dirlist1_gproc)
    except:
        pass
    
    try:
        uimagecube = np.stack([newraw_image_data[science_frame] for science_frame in dirlist1_uproc],axis=0)
        if(TOGGLESIGMA == 'Yes'):
            plot_grid(uimagecube,dirlist1_uproc)
    except:
        pass
    
    try:
        iimagecube = np.stack([newraw_image_data[science_frame] for science_frame in dirlist1_iproc],axis=0)
        if(TOGGLESIGMA == 'Yes'):
            plot_grid(iimagecube,dirlist1_iproc)
    except:
        pass
    
    try:
        zimagecube = np.stack([newraw_image_data[science_frame] for science_frame in dirlist1_zproc],axis=0)
        if(TOGGLESIGMA == 'Yes'):
            plot_grid(zimagecube,dirlist1_zproc)
    except:
        pass
    
    try:
        redimagecube = np.stack([newraw_image_data[science_frame] for science_frame in dirlist1_redproc],axis=0)
        if(TOGGLESIGMA == 'Yes'):
            plot_grid(redimagecube,dirlist1_redproc)
    except:
        pass
    
    try:
        greenimagecube = np.stack([newraw_image_data[science_frame] for science_frame in dirlist1_greenproc],axis=0)
        if(TOGGLESIGMA == 'Yes'):
            plot_grid(greenimagecube,dirlist1_greenproc)
    except:
        pass
    
    try:
        blueimagecube = np.stack([newraw_image_data[science_frame] for science_frame in dirlist1_blueproc],axis=0)
        if(TOGGLESIGMA == 'Yes'):
            plot_grid(blueimagecube,dirlist1_blueproc)
    except:
        pass
    
    try:
        Haimagecube = np.stack([newraw_image_data[science_frame] for science_frame in dirlist1_Haproc],axis=0)
        if(TOGGLESIGMA == 'Yes'):
            plot_grid(Haimagecube,dirlist1_Haproc)
    except:
        pass
    
    try:
        lumimagecube = np.stack([newraw_image_data[science_frame] for science_frame in dirlist1_lumproc],axis=0)
        if(TOGGLESIGMA == 'Yes'):
            plot_grid(lumimagecube,dirlist1_lumproc)
    except:
        pass
    
    try:
        allimagecube = np.stack([newraw_image_data[science_frame] for science_frame in dirlist1_proc],axis=0)
        plot_grid(allimagecube,dirlist1_proc)
    except:
        pass
    
    plt.show()
    
    print("Image reduction program has concluded.  Your images have been successfully processed.")
    gc.disable()
if __name__ == '__main__':
    main(date,objid)