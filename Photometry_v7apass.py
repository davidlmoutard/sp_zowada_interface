

def main(date, objid, predmag, magrange, fwhmmultiplier, snr, apradius, anninner, annouter, directories):
    # -*- coding: utf-8 -*-
    """
    # # Photometry using astropy photutils and AAVSO comparison stars
    # 
    # ## 1 Import dependencies and setup matplotlib
    
    # All file names MUST end in the <filter>.fits or <filter>.fts to ensure proper sort
    
    # This is current Production Photometry Code 8 - 9 - 20
    """
    import os
    import requests, math, glob
    import pandas as pd
    import numpy as np
    from photutils import DAOStarFinder
    from astropy.stats import mad_std
    from astropy.stats import sigma_clipped_stats
    from astropy.io import fits
    from astropy.coordinates import SkyCoord
    from astropy import coordinates as coord
    from astropy.wcs import WCS
    import astropy.units as u
    import matplotlib.pyplot as plt
    from photutils import aperture_photometry, CircularAperture, CircularAnnulus
    from astroquery.simbad import Simbad
    from astroquery.vizier import Vizier
    from datetime import datetime
    import warnings
    from scipy import stats
    warnings.filterwarnings('ignore')
    Vizier.ROW_LIMIT=-1
    import time
    import json
    '''
    print("It is important to determine whether you will be running Photometry on a single series of")
    print("images or performing a detailed analysis of alternative apertures.")
    trigger = input('Do you want to perform a simple Photometry analysis?  Please answer "Yes" or "No"')
    time.sleep(3)
    '''
    #targetsub = 'J005705.02+120635.5'
    
    ###
    targetsub = objid
    Date = date
    # TARGET_RAN ='93.3246275' #AT2020znl
    # TARGET_DECN = '-15.8919653'
    # TARGET_RAN ='128.4261105' #AT2020adow
    # TARGET_DECN = '27.7122696'
    
    estflux = predmag
    fluxrange =  magrange
    ###
    
    with open(directories) as f:
        dirs = json.load(f)
    snr = snr
    fwhmmultiplier = fwhmmultiplier
    
    trigger = 'Yes'
    if(trigger == 'Yes'):
        trigger = 'off'
    else:
        trigger = 'on'
    
    ANNULUSTOGGLE = 'Median'
    
    # trigger is "off" to ask questions and collect values interactively; trigger is "on" if values are provided by the PhotoWrapper script
    # "Typical" settings 2 - 1.5 - 2     
    if(trigger == 'off'):
        aperture_radius_factor = apradius
        annular_radius_factor_inner = anninner
        annular_radius_factor_outer = annouter
        fieldofview=0.40
    else:
        from PhotoWrapper2 import aperture_radius_factor 
        from PhotoWrapper2 import annular_radius_factor_inner 
        from PhotoWrapper2 import annular_radius_factor_outer 
        from PhotoWrapper2 import fieldofview
    
    print('TestLoop end = ',aperture_radius_factor, annular_radius_factor_inner, annular_radius_factor_outer)
    
    print(trigger)
    
    ##### Temporarily commented out for Troubleshooting
    if(trigger=='off'):
        targetdir = targetsub
    #    targetdir = input('Please type in the subdirectory of your "Processed" files directory which holds your processed images:  ') 
    #    time.sleep(3)
    else:
        from PhotoWrapper2 import targetdir 
    
    calibrationpath = dirs['Calibration']
    flatspath = dirs['Flats']
    processedpath = dirs['Templates']
    #targetpath= processedpath + '/' + targetdir
    
    if(trigger=='off'):
    #    TOGGLESTACK = 'Averaged'
        TOGGLESTACK = 'Median'
    #    TOGGLESTACK=input("What stacking method did you perform during alignment?  Answer: 'Summed', 'Averaged', 'Mean' or 'Median'?  You have 3 seconds.  ")
    #    time.sleep(3)
    else:
        from PhotoWrapper2 import TOGGLESTACK 
    
    #TOGGLESTACK = 'Averaged'
    if (TOGGLESTACK == 'Summed'):
        targetpath= processedpath + '/' + targetdir + '/' + 'summed'
    elif (TOGGLESTACK == 'Averaged'):
        targetpath= processedpath + '/' + targetdir + '/' + 'averaged'
    elif (TOGGLESTACK == 'Mean'):
        targetpath= processedpath + '/' + targetdir + '/' + 'nanmean'
    elif (TOGGLESTACK == 'Median'):
        targetpath= processedpath + '/' + targetdir + '/' + 'nanmedian'
    
    targetpath= processedpath + targetsub+'/' + Date
    findloc = targetpath+"/"+os.listdir(targetpath)[0]
    HDUList = fits.open(findloc) 
    lightheader = HDUList[0].header
    HDUList.close()
    sc = SkyCoord(str(str(lightheader['RA']) + str(lightheader['DEC'])) , unit=(u.hourangle, u.deg))
    TARGET_RAN = str(sc.ra.deg)
    TARGET_DECN = str(sc.dec.deg)
    print(targetpath)
    
    targetlistAr=glob.glob(targetpath + '/' + '*-r.fit')+glob.glob(targetpath + '/' + '*-r.fts')+glob.glob(targetpath + '/' + '*-r.fits')
    print(targetlistAr)
    
    targetlistAu=glob.glob(targetpath + '/' + '*-u.fit')+glob.glob(targetpath + '/' + '*-u.fts')+glob.glob(targetpath + '/' + '*-u.fits')
    print(targetlistAu)
    
    targetlistAg=glob.glob(targetpath + '/' + '*-g.fit')+glob.glob(targetpath + '/' + '*-g.fts')+glob.glob(targetpath + '/' + '*-g.fits')
    print(targetlistAg)
    
    targetlistAi=glob.glob(targetpath + '/' + '*-i.fit')+glob.glob(targetpath + '/' + '*-i.fts')+glob.glob(targetpath + '/' + '*-i.fits')
    print(targetlistAi)
    
    targetlistAz=glob.glob(targetpath + '/' + '*-z.fit')+glob.glob(targetpath + '/' + '*-z.fts')+glob.glob(targetpath + '/' + '*-z.fits')
    print(targetlistAz)
    
    targetlistAred=glob.glob(targetpath + '/' + '*-red.fit')+glob.glob(targetpath + '/' + '*-red.fts')+glob.glob(targetpath + '/' + '*-red.fits')
    print(targetlistAred)
    
    targetlistAblue=glob.glob(targetpath + '/' + '*-blue.fit')+glob.glob(targetpath + '/' + '*-blue.fts')+glob.glob(targetpath + '/' + '*-blue.fits')
    print(targetlistAblue)
    
    targetlistAgreen=glob.glob(targetpath + '/' + '*-green.fit')+glob.glob(targetpath + '/' + '*-green.fts')+glob.glob(targetpath + '/' + '*-green.fits')
    print(targetlistAgreen)
    
    targetlistAHa=glob.glob(targetpath + '/' + '*-Ha.fit')+glob.glob(targetpath + '/' + '*-Ha.fts')+glob.glob(targetpath + '/' + '*-Ha.fits')
    print(targetlistAHa)
    
    targetlistAlum=glob.glob(targetpath + '/' + '*-lum.fit')+glob.glob(targetpath + '/' + '*-lum.fts')+glob.glob(targetpath + '/' + '*-lum.fits')
    print(targetlistAlum)
    
    targetlistA=targetlistAblue+targetlistAred+targetlistAgreen+targetlistAr+targetlistAg \
    +targetlistAi+targetlistAu+targetlistAz+targetlistAlum+targetlistAHa
    print(targetlistA)
    
    targetlistB=targetlistAblue+targetlistAred+targetlistAgreen+targetlistAr+targetlistAg \
    +targetlistAi+targetlistAu+targetlistAz+targetlistAlum+targetlistAHa
    
    
    targetdata = {}
    for targetname in targetlistA: targetdata[targetname] = fits.getdata(targetname)
    
    
    for targetname in targetlistA: print (targetname, targetdata[targetname].shape)
    
    print('Echo', aperture_radius_factor, annular_radius_factor_inner, annular_radius_factor_outer)
    
    # This section builds Menu with listing of Files, Objects and Filters
    
    numtargetlistA=len(targetlistA)
    
    menuA = []
    menuB = []
    menu = []
    menu1 = []
    menu2 = []
    menu3 = []
    menu4 = []
    menu5 = []
    
    for i in range(numtargetlistA): #targetlistA:
        name = targetlistA[i]
        FITS_FILE = targetlistA[i]
        hdulist = fits.open(FITS_FILE)
        data = hdulist[0].data.astype(float)
        header = hdulist[0].header
        wcs = WCS(header)
        filter=hdulist[0].header['Filter']
        target = hdulist[0].header['Object']
        exposure = hdulist[0].header['Exposure']
        try:
            stackcount = hdulist[0].header['Stack_Count']
        except:
            stackcount = 1
        try:
            FWHMtest = hdulist[0].header['FWHM']
        except:
            FWHMtest = 'nan'
    #
        FWHMA = FWHMtest
        menuA.append({'frame_no': i, 'Reduced_Image': name, 'target': target, 'filter': filter, 'Stack Count': stackcount, 'Exposure': exposure, 'SNR': snr,'FWHMA': FWHMA})
        hdulist.close()        
    
    menuB = menuA
    
    menuA = pd.DataFrame(menuA)
    targetlistC = pd.DataFrame(targetlistA)
    
    for image in menuB:    
        filename = image['Reduced_Image']
        filename1 = filename.replace(targetpath,"")
        filename2 = filename1.replace(image['target'],"")
        filename3 = filename2[7:10]
        filename4 = filename3.strip()
        filename4 = filename4.rsplit('-',1)[0]
        image['Science_Image'] = filename4
    
    menuy = pd.DataFrame(menuB)
    menuy
    
    ## THIS SECTION ELIMINATE CRAPPY FILES AND NOTIFIES THE USER
    
    menuC = menuB
    targetlistD = targetlistA
    
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    pd.set_option('max_colwidth', -1)
    
    print("We will now check for files with inadequate quality for analysis based upon FWHM and will not be processed:")
    time.sleep(3)
    iterationA=0
    countnan=0
    iterationB=0
    #print("targetlistA at start")
    #print(targetlistA)
    #print()
    numtargetlistA = len(targetlistA)
    
    for file in targetlistA:
    #    print("current file in loop", file)
    #    print("targetlistB at beginning of loop",targetlistB)
    #    print("at beginning of loop","iterationA", iterationA, "countnan", countnan, "iterationB", iterationB)
    #    print(menuA['FWHMA'][iterationA])
        if(menuA['FWHMA'][iterationA] == 'nan'):
            print("This", file, "has inadequate quality for processing.")
    #        print("delete this file",file)
            targetlistB.pop(iterationB)
            menuB.pop(iterationB)
            countnan=countnan+1
    #    print("targetlistB after loop",targetlistB)
        iterationA = iterationA + 1
        iterationB = iterationA - countnan
    #    print("new iterationA", iterationA, "new iterationB", iterationB)
        print()
    
    if(countnan==0):
        print("No worries.  All of your files had adequate FWHM")
    
    # menu = menu.query('FWHM1 != "nan"')
    # menu = menu.reset_index(drop=True)
    
    targetlist1 = targetlistB
    numtargetlist1=len(targetlist1)
    
    if (numtargetlist1 == 0):
        print("I am sorry: All of your files have inadequate FWHM")
        print("I regret to inform you that further analysis will not be possible.")
        quit()
    
    for i in range(numtargetlist1): 
        name = targetlist1[i]
        FITS_FILE = targetlist1[i]
        hdulist = fits.open(FITS_FILE)
        data = hdulist[0].data.astype(float)
        header = hdulist[0].header
        wcs = WCS(header)
        filter=hdulist[0].header['Filter']
        target = hdulist[0].header['Object']
        exposure = hdulist[0].header['Exposure']
        try:
            stackcount = hdulist[0].header['Stack_Count']
        except:
            stackcount = 1
        FWHM1 = hdulist[0].header['FWHM']
        date=hdulist[0].header['Date']
        menu.append({'frame_no': i, 'Reduced_Image': name, 'FITS_target': target, 'filter': filter, 'Stack Count': stackcount, 'Exposure': exposure,'FWHM1': FWHM1,'SNR':snr})
    
    for image in menu:    
        filename = image['Reduced_Image']
        filename1 = filename.replace(targetpath,"")
        filename2 = filename1.replace(image['FITS_target'],"")
        filename3 = filename2[7:10]
        filename4 = filename3.strip()
        filename4 = filename4.rsplit('-',1)[0]
        image['Science_Image'] = filename4
    
    
    menu1 = menu
    menu = pd.DataFrame(menu)
    
    #print(menu)
    
    ## THIS SECTION IDENTIFIES THE TARGET AND ITS RA AND DEC. ##  MAY NEED TO CORRECT LOGIC 
    
    ##### Temporarily commented out for Troubleshooting
    target = menu['FITS_target'][0]
    Flag = 1
    
    if(trigger=='off'):
    #    TOGGLEACCEPT = 'Yes'
        TOGGLEACCEPT = 'No'
    #    TOGGLEACCEPT=input("Do you wish to accept the FITS File Title Object as your Target?  Answer 'Yes' or 'No'. You have 3 seconds.  ")
    #    time.sleep(3)
    else:
        from PhotoWrapper2 import TOGGLEACCEPT 
    
    #TOGGLEACCEPT='No'
    # time.sleep(3)
    if(TOGGLEACCEPT=='No'):
        if(trigger=='off'):
    #        TOGGLENAME=input("Will you be providing an accepted SDSS or SIMBAD designation?  Answer 'Yes' or 'No'. You have 3 seconds. ")
    #        time.sleep(3)
            TOGGLENAME = 'No'
            if(TOGGLENAME=='Yes'):
                target = input("Please provide the SDSS or SIMBAD designation for your object.  You have 5 seconds.  ")
                time.sleep(5)
            else:
    #            TARGET_RAN=input("Please provide the RA for your object as precisely as possible in decimal form. (Example: 348.16308100) You have 5 seconds.  ")
    #           time.sleep(5)
    #            TRA=TARGET_RAN
    #            TARGET_DECN=input("Please provide the DEC for your object as precisely as possible in decimal form. (Example: 11.04738500) You have 5 seconds.  ")
    #            time.sleep(5)
    #            TDEC=TARGET_DECN
    #            TARGET_FLUX_R=input("Please provide the estimated magnitude of your object. (Example: 16) You have 3 seconds.  ")
    #            TARGET_FLUX_R=float(TARGET_FLUX_R)
    #            time.sleep(3)
    #            flux_range=input("Please provide your uncertainty of magnitude for your object. (1.5 is the Program default) You have 3 seconds.  ")        
    #            flux_range=float(flux_range)
    #            time.sleep(3)
                TARGET_RAN = TARGET_RAN
                TARGET_DECN = TARGET_DECN
                TARGET_FLUX_R = estflux
                flux_range = fluxrange
                Flag=2
                TRA = TARGET_RAN
                TDEC = TARGET_DECN
        else:
            from PhotoWrapper2 import TOGGLENAME
            if(TOGGLENAME == 'Yes'):
                from PhotoWrapper2 import target
            else:
                from PhotoWrapper2 import TARGET_RAN
                TRA = TARGET_RAN
                from PhotoWrapper2 import TARGET_DECN
                TDEC = TARGET_DECN
                from PhotoWrapper2 import TARGET_FLUX_R
                TARGET_FLUX_R=float(TARGET_FLUX_R)
                from PhotoWrapper2 import flux_range
                flux_range=float(flux_range)
                from PhotoWrapper2 import Flag
    else:
        flux_range = fluxrange
    
    
    if(trigger=='off'):
        print("Do you wish to see 'All', 'Selected' or 'No' Maps and Regressions?")
        print("If you select 'Selected', you will see Maps and Regressions for R and U Filters Osnly")
    #    TOGGLEOPTIONS=input("Answer 'All', 'Selected' or 'No'. You have 3 seconds.  ")
        TOGGLEOPTIONS = 'No'
    #    time.sleep(3)
    else:
        from PhotoWrapper2 import TOGGLEOPTIONS
    
    print("Do you wish to print 'Results' for Work-in-Process filtering of the data to Excel files for debugging and/or follow-on analysis?")
    #ResultsPrint = input("Answer 'Yes' or 'No'. You have 3 seconds.  ")
    ResultsPrint = 'No'
    #time.sleep(3)
    
    STAR_NAME = target
    
    for image in menu1:    
        image['Photometry_Object'] = STAR_NAME
    
    #Toggle='No' #Placeholder for SDSS Status
    
    # For SIMBAD Stars
    # Use astroquery to locate the RA/DEC of our target star. # Disconnect from actual centroid location fixed later in program
    # Astroquery will not work for SDSS Object
    
    if (Flag == 1):
        Simbad.add_votable_fields('flux(R)')
        astroquery_results = Simbad.query_object(STAR_NAME)
        if(astroquery_results is not None): #  for SIMBAD STARS Results will be Present
            print('One')
            Toggle='No' #Not SDSS Object (but SIMBAD)
            print(target, " is an Astroquery object")
    #        print(astroquery_results.colnames)
            TARGET_FLUX_R=str(astroquery_results[0]['FLUX_R'])
            if(TARGET_FLUX_R=='--'):
    #            TARGET_FLUX_R=18.5
                if(trigger=='off'):
                    TARGET_FLUX_R=input("Please provide the estimated magnitude of your object. (Example: 16) You have 3 seconds.  ")
    #                TARGET_FLUX_R = float(14.5)
                    TARGET_FLUX_R=float(TARGET_FLUX_R)
                    time.sleep(3)
                else:
                    from PhotoWrapper2 import TARGET_FLUX_R
                    TARGET_FLUX_R=float(TARGET_FLUX_R)
    #            astroquery_results['FLUX_R'][0]=TARGET_FLUX_R
            TARGET_RA = str(astroquery_results[0]['RA'])
            TARGET_DEC = str(astroquery_results[0]['DEC']).replace(".*"," *0").replace('+','')
            star_coord = SkyCoord(TARGET_RA,TARGET_DEC, unit=(u.hourangle, u.deg))
            new0=star_coord.to_string()
            new01=new0.split()
            TRA=new01[0]
            TDEC=new01[1]
            TARGET_RAN=TRA
            TARGET_DECN=TDEC
            TARGET_FLUX_R=float(TARGET_FLUX_R)
    #
            Flag=1 # for Function selection of SIMBAD/Other Results (1=SIMBAD)
    #   
    #        print(astroquery_results)
    #        print(target,TARGET_RAN, TARGET_DECN, TARGET_FLUX_R)
    #
        else:
            print('two')
            if(trigger=='off'):
                Toggle = 'Yes'
    #            Toggle = input("Is this an SDSS Object?  Please answer Yes or No. (You have 3 sec)   ")
    #            time.sleep(3)
            else:
                from PhotoWrapper2 import TOGGLE 
             # Toggle='Yes'
            print('three') 
            print(target, " is an apass9 Catalog object")
            field_of_view =0.001
            Vizier.ROW_LIMIT = -1
            ToggleVal = Vizier.query_region(STAR_NAME, radius=field_of_view * u.deg, catalog = 'apass9') 
            if(ToggleVal is not None):
                targetobjresults=[]
                targetobj={}
                targetobj = Vizier.query_region(STAR_NAME, radius=field_of_view * u.deg, catalog = 'apass9') 
    #            print ('apass9', targetobj[0]['apass9'][0], '; RA = ', targetobj[0]['RA_ICRS'][0],'; DEC = ', targetobj[0]['DE_ICRS'][0],'; rmag = ', targetobj[0]['rmag'][0],'; umag = ', targetobj[0]['umag'][0],'; gmag = ', targetobj[0]['gmag'][0],\
    #                               ';imag = ', targetobj[0]['imag'][0],'; zmag = ', targetobj[0]['zmag'][0])
                targetobjresults.append({'auid':'apass9'+str(targetobj[0]['apass9'][0]),'ra': targetobj[0]['RA_ICRS'][0],'dec':targetobj[0]['DE_ICRS'][0], 'rmag': targetobj[0]['rmag'][0], 'umag': targetobj[0]['umag'][0], 'gmag': targetobj[0]['gmag'][0],\
                                                'imag': targetobj[0]['imag'][0], 'zmag': targetobj[0]['zmag'][0]})
    #
                targetobjresults = pd.DataFrame(targetobjresults)
    #  
                TRA = targetobjresults.ra[0]
                TDEC = targetobjresults.dec[0]
                TARGET_FLUX_R=targetobjresults.rmag[0]
                TARGET_RAN=TRA
                TARGET_DECN=TDEC
                Flag=2  # for function selection of SIMBAD/Other Results (2=Other)
        #
    #            print('target: ',target)
    #            print('RA: ',TARGET_RAN)
    #            print ('DEC: ',TARGET_DECN)
    #            print('Target Flux in R Band', TARGET_FLUX_R)
    #This is new 
    else:
        Simbad.add_votable_fields('flux(R)')
        if(trigger=='off'):
            TARGET_FLUX_R=float(TARGET_FLUX_R)
        else:
            from PhotoWrapper2 import TARGET_FLUX_R 
            TARGET_FLUX_R=float(TARGET_FLUX_R)
    
    '''
    def get_results(TRA, TDEC, TARGET_FLUX_R, Flag,flux_range,fieldofview): 
    #
        field_of_view=fieldofview
    #
        # Run for SIMBAD Objects Only -- Need logic to choose
        if(Flag==1):
            print('alpha')
            calibrators={}
            Vizier.ROW_LIMIT = -1
            calibrators = Vizier(columns=['*','Q','class','flags','q_mode','mode'],row_limit = -1).query_region(STAR_NAME, radius=field_of_view * u.deg, catalog = 'SDSS12') 
            if len(calibrators)>1:
                for i in range(0, len(calibrators[0]['SDSS12'])):
                    if (calibrators[0]['rmag'][i] > astroquery_results['FLUX_R'][0] - flux_range) and (calibrators[0]['rmag'][i] < astroquery_results['FLUX_R'][0] + flux_range): # If the r mag is roughly within a range of +/- 2 mag from target
                        results.append({'auid':'apass9'+str(calibrators[0]['SDSS12'][i]),'ra': calibrators[0]['RA_ICRS'][i],'dec':calibrators[0]['DE_ICRS'][i], 'rmag': calibrators[0]['rmag'][i], 'umag': calibrators[0]['umag'][i], 'gmag': calibrators[0]['gmag'][i],\
                                        'imag': calibrators[0]['imag'][i], 'zmag': calibrators[0]['zmag'][i],'objecttype': calibrators[0]['class'][i],'quality':calibrators[0]['Q'][i],'flags':hex(calibrators[0]['flags'][i]),'q_mode': calibrators[0]['q_mode'][i],'mode': calibrators[0]['mode'][i]})
            else:
                print('beta')
                calibrators={}
                calibrators = Vizier(columns=['*','Q','class','flags','q_mode','mode'],row_limit = -1).query_region(coord.SkyCoord(ra=TRA, dec=TDEC, unit=(u.deg, u.deg),frame='icrs'), radius=field_of_view * u.deg, catalog = 'SDSS12') 
                for i in range(0, len(calibrators[0]['SDSS12'])):
                    if (calibrators[0]['rmag'][i] > TARGET_FLUX_R - flux_range) and (calibrators[0]['rmag'][i] < TARGET_FLUX_R + flux_range): # If the r mag is roughly within a range of +/- 2 mag from target
                        results.append({'auid':'SDSS12'+str(calibrators[0]['SDSS12'][i]),'ra': calibrators[0]['RA_ICRS'][i],'dec':calibrators[0]['DE_ICRS'][i], 'rmag': calibrators[0]['rmag'][i], 'umag': calibrators[0]['umag'][i], 'gmag': calibrators[0]['gmag'][i],\
                                        'imag': calibrators[0]['imag'][i], 'zmag': calibrators[0]['zmag'][i],'objecttype': calibrators[0]['class'][i],'quality':calibrators[0]['Q'][i],'flags':hex(calibrators[0]['flags'][i]),'q_mode': calibrators[0]['q_mode'][i],'mode': calibrators[0]['mode'][i]})
        else:
            print('gamma')
            calibrators={}
            calibrators = Vizier(columns=['*','Q','class','flags','q_mode','mode'],row_limit = -1).query_region(coord.SkyCoord(ra=TRA, dec=TDEC, unit=(u.deg, u.deg),frame='icrs'), radius=field_of_view * u.deg, catalog = 'SDSS12') 
            a = len(calibrators[0]['SDSS12'])
            print(a)
            if len(calibrators)>0:
                for i in range(0, len(calibrators[0]['SDSS12'])):
                    if (calibrators[0]['rmag'][i] > TARGET_FLUX_R - flux_range) and (calibrators[0]['rmag'][i] < TARGET_FLUX_R + flux_range): # If the r mag is roughly within a range of +/- 2 mag from target
                        results.append({'auid':'SDSS12'+str(calibrators[0]['SDSS12'][i]),'ra': calibrators[0]['RA_ICRS'][i],'dec':calibrators[0]['DE_ICRS'][i], 'rmag': calibrators[0]['rmag'][i], 'umag': calibrators[0]['umag'][i], 'gmag': calibrators[0]['gmag'][i],\
                                        'imag': calibrators[0]['imag'][i], 'zmag': calibrators[0]['zmag'][i],'objecttype': calibrators[0]['class'][i],'quality':calibrators[0]['Q'][i],'flags':hex(calibrators[0]['flags'][i]),'q_mode': calibrators[0]['q_mode'][i],'mode': calibrators[0]['mode'][i]})
        return(results)    
    '''
    def get_results(TRA, TDEC, TARGET_FLUX_R, Flag,flux_range,fieldofview): 
    #
        field_of_view=fieldofview
        if(Flag==1):
            print('alpha')
            calibrators={}
            Vizier.ROW_LIMIT = -1
            calibrators = Vizier.query_region(STAR_NAME, radius=field_of_view * u.deg, catalog = 'apass9') 
            if len(calibrators)>1:
                for i in range(0, len(calibrators[0])):
                    if (calibrators[0]['r_mag'][i] > TARGET_FLUX_R - flux_range) and (calibrators[0]['r_mag'][i] < TARGET_FLUX_R + flux_range): # If the r mag is roughly within a range of +/- 2 mag from target
                        results.append({'auid':'apass9'+'_' + str(calibrators[0]['recno'][i]),'ra': calibrators[0]['RAJ2000'][i],'dec':calibrators[0]['DEJ2000'][i], 'rmag': calibrators[0]['r_mag'][i], 'umag': calibrators[0]['Bmag'][i], 'gmag': calibrators[0]['g_mag'][i],\
                                        'imag': calibrators[0]['i_mag'][i], 'zmag': calibrators[0]['Vmag'][i],'objecttype': '6','quality':'3','flags':'0x100001010000200','q_mode': '+','mode': '1'})                    
            else:
                print('beta')
                calibrators={}
                calibrators = Vizier.query_region(coord.SkyCoord(ra=TRA, dec=TDEC, unit=(u.deg, u.deg),frame='icrs'), radius=field_of_view * u.deg, catalog = 'apass9') 
                for i in range(0, len(calibrators[0])):
                    if (calibrators[0]['r_mag'][i] > TARGET_FLUX_R - flux_range) and (calibrators[0]['r_mag'][i] < TARGET_FLUX_R + flux_range): # If the r mag is roughly within a range of +/- 2 mag from target
                        results.append({'auid':'apass9'+'_' +str(calibrators[0]['recno'][i]),'ra': calibrators[0]['RAJ2000'][i],'dec':calibrators[0]['DEJ2000'][i], 'rmag': calibrators[0]['r_mag'][i], 'umag': calibrators[0]['Bmag'][i], 'gmag': calibrators[0]['g_mag'][i],\
                                            'imag': calibrators[0]['i_mag'][i], 'zmag': calibrators[0]['Vmag'][i],'objecttype': '6','quality':'3','flags':'0x100001010000200','q_mode': '+','mode': '1'})
        else:
            print('gamma')
            calibrators={}
            calibrators = Vizier.query_region(coord.SkyCoord(ra=TRA, dec=TDEC, unit=(u.deg, u.deg),frame='icrs'), radius=field_of_view * u.deg, catalog = 'II/336/apass9') 
            a = len(calibrators[0])
            print(a)
            if len(calibrators)>0:
                for i in range(0, len(calibrators[0])):
                    if (calibrators[0]['r_mag'][i] > TARGET_FLUX_R - flux_range) and (calibrators[0]['r_mag'][i] < TARGET_FLUX_R + flux_range): # If the r mag is roughly within a range of +/- 2 mag from target
                        results.append({'auid':'apass9'+'_' +str(calibrators[0]['recno'][i]),'ra': calibrators[0]['RAJ2000'][i],'dec':calibrators[0]['DEJ2000'][i], 'rmag': calibrators[0]['r_mag'][i], 'umag': calibrators[0]['Bmag'][i], 'gmag': calibrators[0]['g_mag'][i],\
                                    'imag': calibrators[0]['i_mag'][i], 'zmag': calibrators[0]['Vmag'][i],'objecttype': '6','quality':'3','flags':'0x100001010000200','q_mode': '+','mode': '1'})                    
    #
        return(results)    
    
    
    results=[]
    results.append({'auid': 'target', 'ra': TRA, 'dec': TDEC})
    
    results=get_results(TARGET_RAN, TARGET_DECN,TARGET_FLUX_R, Flag,flux_range,fieldofview)
    
    resultz = pd.DataFrame(results)
        
    if (float(TARGET_FLUX_R)>13):
        BRIGHTEST_COMPARISON_STAR_MAG = float(TARGET_FLUX_R)-flux_range
        DIMMEST_COMPARISON_STAR_MAG = float(TARGET_FLUX_R)+flux_range
    else:
        BRIGHTEST_COMPARISON_STAR_MAG = float(TARGET_FLUX_R)-2*flux_range
        DIMMEST_COMPARISON_STAR_MAG = float(TARGET_FLUX_R)+2*flux_range  
    
        
    print(BRIGHTEST_COMPARISON_STAR_MAG, DIMMEST_COMPARISON_STAR_MAG)
    
    # Creates 13 Character Names for Results for later elimination of duplicates
    i=0
    for star in results:    
        name = star['auid']
        name0 = name[0:16]
        star['auidshort'] = name0          
        if (name == 'target'):
            star['q_mode'] = '+'
            star['objecttype'] = 6
            
    #print(results)
    
    #  This will be where we loop through images
    
    def get_map(i,stars,results):
        FITS_FILE = targetlist1[i]
        hdulist = fits.open(FITS_FILE)
        data = hdulist[0].data.astype(float)
        plt.figure(figsize=(18, 16), dpi= 80, edgecolor='k')
    #    fig=plt.figure(figsize=(18, 16), dpi= 80, edgecolor='k')
        #fig=plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
        count=0
        #plt.imshow(data, cmap='gray', origin='lower', vmin=25000, vmax=50000)
        #plt.imshow(data, cmap='gray_r', origin='lower', vmin=0, vmax=50000)    
        #plt.show()
        while count<stars:
            if (results['x'][count]=='nan'):
                count=count+1
    #            print(count)
    #            print(x)
            else:
                positions=(results['x'][count],results['y'][count])
                count=count+1
                position_apps = CircularAperture(positions, r=12.)    
        # target_app = CircularAperture(target_xy, r=20.)    
                position_apps.plot(color='red', lw=1.5, alpha=0.5)
                plt.imshow(data, cmap='gray', origin='lower', vmin=25000, vmax=50000)
        #        plt.imshow(data, cmap='gray_r', origin='lower', vmin=0, vmax=2500)
        #plt.show()
        # target_app.plot(color='blue', lw=1.5, alpha=0.5)
        while count<stars:
            if (results['x'][count]=='nan'):
                count=count+1
    #            print(count)
    #            print(x)
            else:
                positions=(results['x'][count],results['y'][count])
                count=count+1
                position_apps = CircularAperture(positions, r=12.)    
        # target_app = CircularAperture(target_xy, r=20.)    
                position_apps.plot(color='red', lw=1.5, alpha=0.5)
                plt.imshow(data, cmap='gray', origin='lower', vmin=0, vmax=50000)
        #        plt.imshow(data, cmap='gray_r', origin='lower', vmin=0, vmax=2500)
        #plt.show()
        # target_app.plot(color='blue', lw=1.5, alpha=0.5)
    #    to_plot = results.query('peak < 50000 and rmag < 20 and rmag > 0') 
        for to_annotate in results.iterrows():
            plt.annotate('{}'.format(to_annotate[1]['auid']),
                xy=(to_annotate[1]['x'], to_annotate[1]['y']), xycoords='data',
                xytext=(-20, 20), color='white', textcoords='offset points', size='8',
        #        xytext=(-20, 20), textcoords='offset points', size='8',
                arrowprops=dict(color='white', arrowstyle="->"))
        #        arrowprops=dict(arrowstyle="->"))
        plt.imshow(data, cmap='gray', origin='lower', vmin=25000, vmax=50000)
        #plt.imshow(data, cmap='gray_r', origin='lower', vmin=0, vmax=50000)
        plt.show()
        hdulist.close
    
    def get_regression(results,filter2):
        print('Commence Regression - Commence Regression - Commence Regression - Commence Regression')
        # now perform ensemble photometry by linear regression of the comparison stars' instrumental mags
        if(filter2=='r'):
            to_linear_fit = results.query('auid != "target" and rmag > {} and rmag < {}'.format(BRIGHTEST_COMPARISON_STAR_MAG, DIMMEST_COMPARISON_STAR_MAG)) 
            x = to_linear_fit['instrumental_mag'].values
            y = to_linear_fit['rmag'].values
        elif(filter2=='g'):
    ##> NOTE:  These are adjustable to maintain range of 'g' filter stars (because results is function only of r-mag)
            G_BRIGHTEST_COMPARISON_STAR_MAG = BRIGHTEST_COMPARISON_STAR_MAG + 0.0
            G_DIMMEST_COMPARISON_STAR_MAG = DIMMEST_COMPARISON_STAR_MAG + 0.5
            to_linear_fit = results.query('auid != "target" and gmag > {} and gmag < {}'.format(G_BRIGHTEST_COMPARISON_STAR_MAG - 1, G_DIMMEST_COMPARISON_STAR_MAG +1 )) 
            x = to_linear_fit['instrumental_mag'].values
            y = to_linear_fit['gmag'].values
        elif(filter2=='u'):
    ##> NOTE:  These are adjustable to maintain range of 'u' filter stars (because results is function only of r-mag)
            U_BRIGHTEST_COMPARISON_STAR_MAG = BRIGHTEST_COMPARISON_STAR_MAG + 0.5
            U_DIMMEST_COMPARISON_STAR_MAG = DIMMEST_COMPARISON_STAR_MAG + 2.0
            to_linear_fit = results.query('auid != "target" and umag > {} and umag < {}'.format(U_BRIGHTEST_COMPARISON_STAR_MAG+1.5, U_DIMMEST_COMPARISON_STAR_MAG+2.0)) 
            x = to_linear_fit['instrumental_mag'].values
            y = to_linear_fit['umag'].values
    #        print(x,y)
        elif(filter2=='i'):
            to_linear_fit = results.query('auid != "target" and imag > {} and imag < {}'.format(BRIGHTEST_COMPARISON_STAR_MAG, DIMMEST_COMPARISON_STAR_MAG)) 
            x = to_linear_fit['instrumental_mag'].values
            y = to_linear_fit['imag'].values
        elif(filter2=='z'):
            to_linear_fit = results.query('auid != "target" and zmag > {} and zmag < {}'.format(BRIGHTEST_COMPARISON_STAR_MAG, DIMMEST_COMPARISON_STAR_MAG))         
            x = to_linear_fit['instrumental_mag'].values
            y = to_linear_fit['zmag'].values
    #### Need = debug track
        step1 = filter2
        print('step1=', step1)    
    #    debug1.append({'step1': step1},ignore_index=True)
        fit, residuals, rank, singular_values, rcond = np.polyfit(x, y, 1, full=True)
        fit_fn = np.poly1d(fit) 
        # fit_fn is now a function which takes in x and returns an estimate for y, 
        # we use this later to calculate the target mag
        #generate a plot to see the quality of the data
        step2 = filter2
        print('step2:',step2)
     #   debug1.append({'step2': step2},ignore_index=True)
        #
        # now perform ensemble photometry by linear regression of the comparison stars' instrumental mags
        target_instrumental_magnitude = results[results.auid=='target']['instrumental_mag'].values[0]
        if len(results)<3:
            target_instrumental_magnitude = results[results.auid=='target']['instrumental_mag'].values[0]
            target_magnitude = results.vmag[0]-2.5*math.log10(results[results.auid=='target']['aperture_sum'].values[0]/results.aperture_sum[0])
        else:
            target_instrumental_magnitude = results[results.auid=='target']['instrumental_mag'].values[0]
            target_magnitude = fit_fn(target_instrumental_magnitude)
    #        print('Magnitude estimate: {} error from residuals {}'.format(target_magnitude, residuals))
        if(TOGGLEOPTIONS == 'All'):
            x = np.append(x,target_instrumental_magnitude)
            y = np.append(y,target_magnitude)
            plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
            plt.plot(x,y, 'yo', x, fit_fn(x), '-k')    
            plt.xlim(np.min(x)-0.2, np.max(x)+0.2)
            plt.ylim(np.min(y)-0.2, np.max(y)+0.2)
            plt.annotate('Target: MAG {} Error {}'.format(target_magnitude, residuals),
                        xy=(target_instrumental_magnitude, target_magnitude), xycoords='data',
                        xytext=(-20, -85), textcoords='offset points', size='16',
                        arrowprops=dict(arrowstyle="->"))
            if(filter2=='r'):
                for to_annotate in to_linear_fit.iterrows():
                    plt.annotate('{}'.format(to_annotate[1]['auid']),
                        xy=(to_annotate[1]['instrumental_mag'], to_annotate[1]['rmag']), xycoords='data',
                        xytext=(-150, 130), textcoords='offset points', size='11',
                        arrowprops=dict(arrowstyle="->"))
            elif(filter2=='g'):
                for to_annotate in to_linear_fit.iterrows():
                    plt.annotate('{}'.format(to_annotate[1]['auid']),
                        xy=(to_annotate[1]['instrumental_mag'], to_annotate[1]['gmag']), xycoords='data',
                        xytext=(-150, 130), textcoords='offset points', size='11',
                        arrowprops=dict(arrowstyle="->"))
            elif(filter2=='u'):
                for to_annotate in to_linear_fit.iterrows():
                    plt.annotate('{}'.format(to_annotate[1]['auid']),
                        xy=(to_annotate[1]['instrumental_mag'], to_annotate[1]['umag']), xycoords='data',
                        xytext=(-150, 130), textcoords='offset points', size='11',
                        arrowprops=dict(arrowstyle="->"))
            elif(filter2=='i'):
                for to_annotate in to_linear_fit.iterrows():
                    plt.annotate('{}'.format(to_annotate[1]['auid']),
                        xy=(to_annotate[1]['instrumental_mag'], to_annotate[1]['imag']), xycoords='data',
                        xytext=(-150, 130), textcoords='offset points', size='11',
                        arrowprops=dict(arrowstyle="->"))
            elif(filter2=='z'):
                for to_annotate in to_linear_fit.iterrows():
                    plt.annotate('{}'.format(to_annotate[1]['auid']),
                        xy=(to_annotate[1]['instrumental_mag'], to_annotate[1]['zmag']), xycoords='data',
                        xytext=(-150, 130), textcoords='offset points', size='11',
                        arrowprops=dict(arrowstyle="->"))
        elif((TOGGLEOPTIONS == 'Selected' and filter2 == 'r') or (TOGGLEOPTIONS == 'Selected' and filter2 == 'u')):
            x = np.append(x,target_instrumental_magnitude)
            y = np.append(y,target_magnitude)
            plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
            plt.plot(x,y, 'yo', x, fit_fn(x), '-k')    
            plt.xlim(np.min(x)-0.2, np.max(x)+0.2)
            plt.ylim(np.min(y)-0.2, np.max(y)+0.2)
            plt.annotate('Target: MAG {} Error {}'.format(target_magnitude, residuals),
                        xy=(target_instrumental_magnitude, target_magnitude), xycoords='data',
                        xytext=(-20, -85), textcoords='offset points', size='16',
                        arrowprops=dict(arrowstyle="->"))
            if(filter2=='r'):
                for to_annotate in to_linear_fit.iterrows():
                    plt.annotate('{}'.format(to_annotate[1]['auid']),
                        xy=(to_annotate[1]['instrumental_mag'], to_annotate[1]['rmag']), xycoords='data',
                        xytext=(-150, 130), textcoords='offset points', size='11',
                        arrowprops=dict(arrowstyle="->"))
            elif(filter2=='u'):
                for to_annotate in to_linear_fit.iterrows():
                    plt.annotate('{}'.format(to_annotate[1]['auid']),
                        xy=(to_annotate[1]['instrumental_mag'], to_annotate[1]['umag']), xycoords='data',
                        xytext=(-150, 130), textcoords='offset points', size='11',
                        arrowprops=dict(arrowstyle="->"))
        step3 = filter2
        print('step3:',step3)
    #    debug1.append({'step3': step3},ignore_index=True)
        plt.show()
        coefficients = fit_fn.c
        slope = coefficients[0].item(0)
        yintercept= coefficients[1].item(0)
        print(target_magnitude,residuals,slope,yintercept)
        return(target_magnitude, residuals, slope, yintercept)
    
    # Following code creates duplicate AUID values for "near" valued objects to support later deletion of duplicates
    #    results(i) = pd.DataFrame(results)
    #    print(results(i))
    
    # Pairs stars from database to stars in science frame
    
    def get_sources(i,snr,fwhmmultiplier):
        FITS_FILE = targetlist1[i]
        hdulist = fits.open(FITS_FILE)
        data = hdulist[0].data.astype(float)
    #    header = hdulist[0].header
    #    fwhm = hdulist[0].header['FWHM']*2.5
        fwhm = hdulist[0].header['FWHM']*fwhmmultiplier
        source_snr = snr #1.5 #2 #5 # 10
        bkg_sigma = mad_std(data)    
        daofind = DAOStarFinder(fwhm=fwhm, threshold=source_snr*bkg_sigma)    
        sources = daofind(data)
        #pd.set_option('display.max_columns', None4
        #pd.set_option('display.max_rows', None)
        #pd.set_option('display.max_colwidth', -1)
        print(sources)
        hdulist.close
        return sources
    
    for i in range(numtargetlist1): #targetlist1:
        name = 'Sources' + str(i)
        sourcename = 'sources' + str(i)
        print(i,sourcename)
        vars()[sourcename] = get_sources(i,snr,fwhmmultiplier)
    
    
    #pd.set_option('display.max_columns', None)
    #pd.set_option('display.max_rows', None)
    #pd.set_option('max_colwidth', -1)
    
    # FOLLOWING IS INSANELY USEFUL FOR SOURCE - CATALOG RESULTS DEBUGGING!!!!
    '''
    loopa = 0
    spots=[]
    spots1=[]
    spots2=[]
    for source in vars()[sourcename]:         
        loopa = loopa + 1
        namenum = loopa
        xcenter = source['xcentroid']
        ycenter = source['ycentroid']
        peakk = source['peak']
        roundness1 = source['roundness1']
        roundness2 = source['roundness2']
        fluxx = source['flux']
        rax = 'na'
        decx = 'na'
        spots.append({'namenum' : namenum, 'xcentroid': xcenter, 'ycentroid': ycenter, 'ra':rax, 'dec' : decx, 'peak': peakk, 'flux':fluxx, 'roundness1' : roundness1, 'roundness2' : roundness2})
    
    spots1 = spots
    spotz = pd.DataFrame(spots)
    spotz = spotz[['namenum','xcentroid','ycentroid','ra','dec','peak','flux','roundness1','roundness2']]
    
    clock = datetime.today()
    SID = str(clock.year) + '-' + str(clock.month) + '-' + str(clock.day) + '-' + 'h' + str(clock.hour) + 'm' + str(clock.minute) + 's' + str(clock.second)
    
    name = STAR_NAME + '_Frame_' + 'Sources' + '_' + SID + '.xlsx'
    
    os.chdir(processedpath)
    name=name.replace(" ","_")
    
    name=name.replace("/","-")
    name=name.replace("*","star")
    spotz.to_excel(name)
    
    loopb = loopa
    results1 = results
    for star1 in results1:
        loopb = loopb +1
        namenum = star1['auid']
        star1['x'] = 0
        star1['y'] = 0
        star_coord = SkyCoord(star1['ra'],star1['dec'], unit=(u.deg, u.deg))
        xy1 = SkyCoord.to_pixel(star_coord, wcs=wcs, origin=1)
        x1 = xy1[0].item(0)
        y1 = xy1[1].item(0)
        star1['x'] = x1
        star1['y'] = y1
        xcenter = star1['x']
        ycenter = star1['y']
        peakk = 'na'
        roundness1 = 'na'
        roundness2 = 'na'
        fluxx = 'na'
        rax = star1['ra']
        decx = star1['dec']
        spots1.append({'namenum' : namenum, 'xcentroid': xcenter, 'ycentroid': ycenter, 'ra': rax, 'dec': decx,'peak': peakk, 'flux':fluxx, 'roundness1' : roundness1, 'roundness2' : roundness2})
    
    spots2 = spots1
    spotz1 = pd.DataFrame(spots1)
    spotz1 = spotz1[['namenum','xcentroid','ycentroid','ra','dec','peak','flux','roundness1','roundness2']]
    
    clock = datetime.today()
    SID = str(clock.year) + '-' + str(clock.month) + '-' + str(clock.day) + '-' + 'h' + str(clock.hour) + 'm' + str(clock.minute) + 's' + str(clock.second)
    
    name = STAR_NAME + '_Frame_' + 'Sources and Results' + '_' + SID + '.xlsx'
    
    os.chdir(processedpath)
    name=name.replace(" ","_")
    
    name=name.replace("/","-")
    name=name.replace("*","star")
    spotz1.to_excel(name)
    
    ''' 
    
    
    # This section creates multiple copies of results, because these get sequantially turned into datafiles
    #for i in range(numtargetlist1): #targetlist1:
    #    results[i] = results
    #    print(i, results[i])
    #    print()    
    
    resultsduplicate=results
    
    debug = []
    #debug1 = []
    
    for i in range(numtargetlist1):
        print(i,numtargetlist1)    
        FITS_FILE = targetlist1[i]
        print(FITS_FILE)
        hdulist = fits.open(FITS_FILE)
        data = hdulist[0].data.astype(float)
        header = hdulist[0].header
        wcs = WCS(header)
        ObsDate = hdulist[0].header['DATE-OBS']
        sourcename = 'sources' + str(i)
    #    print(sourcename)
        results = resultsduplicate
        if (menu['filter'][i] == 'Red'):
            filter2 = 'r'
        else:
            filter2 = menu['filter'][i]
        print(FITS_FILE,filter2)
    #    print(wcs)
    #    print(filter2)
    #    print(results)
        for star in results:
            star['x'] = 0
            star['y'] = 0
            star['peak'] = 0
            star_coord = SkyCoord(star['ra'],star['dec'], unit=(u.deg, u.deg))
            xy = SkyCoord.to_pixel(star_coord, wcs=wcs, origin=1)
            x = xy[0].item(0)
            y = xy[1].item(0)
    #        print('star in results', star, star_coord, x,y) #Following 5 Lines are the issue!
        #   star['xact']=x
        #   star['yact']=y
            loop = 0
            for source in vars()[sourcename]:         
                loop = loop + 1 
    #            print(source, sourcename)        
                if (source['xcentroid']-2 < x < source['xcentroid']+2) and (source['ycentroid']-2 < y < source['ycentroid']+2):
    #                star['x'] = x
    #                star['y'] = y
                    print('xcentroid = ', source['xcentroid'], 'ycentroid = ', source['ycentroid'], 'source peak = ',source['peak'])
                    star['x'] = source['xcentroid']
                    star['y'] = source['ycentroid']
                    star['peak'] = source['peak']
                    print(star['auid'])
                    if (star['auid'] == 'target'):
    #                    print('xcentroid = ', source['xcentroid'], 'ycentroid = ', source['ycentroid'], 'source peak = ',source['peak'])
                        targetdec = star['dec']
                        targetra = star['ra']
    # checking logic for star in results:
            #            print(star['auid'], star['dec'], star['ra'])  
                    if (star['auid'] != 'target'):
    #                    print('Not', star['auid'], star['dec'], star['ra'], targetdec, targetra)        
            #                if (targetdec - 0.0001 < star['dec'] < targetdec + 0.0001) and (targetra - 0.0001 <  star['ra'] < targetra + 0.0001):
                        if (((float(targetdec) - 0.0001) < star['dec'] < (float(targetdec) + 0.0001)) and ((float(targetra) - 0.005) < star['ra'] < (float(targetra) + 0.005))):
                            star['auid']='target'
                            break
        print(i, sourcename, loop)
        results = pd.DataFrame(results)
    #    print(results)
    #    print(results)
    #
    # Section drops duplicate identifications of auid 
        results.drop_duplicates(subset='auid', keep='first', inplace=True)
        results = results.reset_index(drop=True)
        results.dropna(subset=['auid'])
        results.dropna(subset=['peak'])
        results = results.query('x > 0') 
        results.drop_duplicates(subset=['x','y','peak'], keep='first', inplace=True)
        results = results.reset_index(drop=True)
        print(results)
    #
        stars=len(results)
        print("stars = ",stars)
            #
        if(TOGGLEOPTIONS == 'All'):
            get_map(i,stars,results)
        elif((TOGGLEOPTIONS == 'Selected' and filter2 == 'r') or (TOGGLEOPTIONS == 'Selected' and filter2 == 'u')):
            get_map(i,stars,results)     
    #
    # Was 2.5 for brighter stars 
        aperture_radius = hdulist[0].header['FWHM']*aperture_radius_factor
        count=0
        positions=()
        while count<stars:
            if (math.isnan(results['x'][count])):
                print(count, results['x'][count], positions)
                count=count+1
            else:
                position=(results['x'][count],results['y'][count])
                positions=position+positions
                print(count, results['x'][count], positions)
                count=count+1
    #
        it=iter(positions)
        pos=list(zip(it,it))
        pos=list(reversed(pos))            
        results = results.query('x > 1 and y > 1')
    #    print(results)
    #
    #    annular_radius_factor_inner = 1.5
    #    annular_radius_factor_outer = 2.0
        apertures= CircularAperture(pos, r=aperture_radius)    
        #Note these were previously 11/8/19 1.5 / 2
        aperture_in=aperture_radius*annular_radius_factor_inner
        aperture_out=aperture_radius*annular_radius_factor_outer
        print('Echo2', aperture_radius_factor, annular_radius_factor_inner, annular_radius_factor_outer)
        menu1[i]['ap_radius'] = aperture_radius
        menu1[i]['ann_radius_inner']=aperture_in
        menu1[i]['ann_radius_outer']=aperture_out
        annulus_apertures = CircularAnnulus(pos, r_in=aperture_in, r_out=aperture_out)
        apers=[apertures, annulus_apertures]
    #        
        phot_table = aperture_photometry(data, apers)  
        for col in phot_table.colnames:
            phot_table[col].info.format = '%.8g'  # for consistent table output
    #    print(phot_table)
    #
        bkg_mean = phot_table['aperture_sum_1'] / annulus_apertures.area
        bkg_sum = bkg_mean * apertures.area
        final_sum = phot_table['aperture_sum_0'] - bkg_sum
        phot_table['residual_aperture_sum'] = final_sum
        phot_table['residual_aperture_sum'].info.format = '%.8g'  # for consistent table output
    #    print(phot_table['residual_aperture_sum'])  
        results['residual_aperture_sum'] = phot_table['residual_aperture_sum']
    #
    # This section calculates median rather than mean annulus flux    
        annulus_masks = annulus_apertures.to_mask(method='center')
    #
        bkg_median = []
        for mask in annulus_masks:
            annulus_data = mask.multiply(data)
            annulus_data_1d = annulus_data[mask.data > 0]
            _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_sigclip)
    #    
        bkg_median = np.array(bkg_median)
    #   
        bkg_mediansum = bkg_median * apertures.area
    #   
        final_mediansum = phot_table['aperture_sum_0'] - bkg_mediansum
    #    
        phot_table1 = phot_table#aperture_photometry(data, apertures)
    #   
        phot_table1['residual_aperture_median_sum'] = final_mediansum
        phot_table1['residual_aperture_median_sum'].info.format = '%.8g'  # for consistent table output
        print(phot_table1['residual_aperture_median_sum'])  
    #   
        #phot_table1['annulus_median'] = bkg_median
        #phot_table1['aper_bkg'] = bkg_median * annulus_apertures.area()
        #phot_table1['aper_sum_bkgsub'] = phot['aperture_sum'] - phot['aper_bkg']
        for col in phot_table1.colnames:
            phot_table1[col].info.format = '%.8g'  # for consistent table output
    #    
        print(phot_table1)
    #    
        results['residual_aperture_mean_sum'] = phot_table1['residual_aperture_sum']
        results['residual_aperture_median_sum'] = phot_table1['residual_aperture_median_sum']
        print(results['residual_aperture_median_sum'] - results['residual_aperture_mean_sum'])
    #
        if(ANNULUSTOGGLE == 'Median'):
            results['residual_aperture_sum'] = phot_table1['residual_aperture_median_sum']
        else:
            results['residual_aperture_sum'] = results['residual_aperture_mean_sum']
    #
        # NOTE THAT THIS LINE IS MANUALLY SET TO REMOVE NEGATIVE VALUES (11/8/19 remark worth revisiting)
    #    results = results.query('residual_aperture_sum > 1000') 
        results = results.query('residual_aperture_sum > 10') 
        results = results.reset_index(drop=True)
        # add a col with calculation for instrumental mag
        results['instrumental_mag'] = results.apply(lambda x: -2.5 * math.log10(x['residual_aperture_sum']), axis = 1)
    #    print(results)
        debug.append({'file_no': i, 'filter2': filter2})
    #
    # Set Filter Switch for "CLEAN" objects
    #
        results['flagint'] = results['flags']
        results['Clean'] = results['flags']
        #
        test1 = '0x8100000c00a4' # This is the Hexidecimal flag for "Clean"
        test1 = int(test1, 16) # This converts the Hexidecimal flag to a base 10 integer
        #
        for j in range(len(results)):
            place = results['flags'][j]
        #        print(place)
            try:
                if (math.isnan(results['flags'][j])):
                    results['flagint'][j] = 'False'
                    results['Clean'][j] = 'True'
            except:
                results['flagint'][j] = int(place,16)
                flagz = results['flagint'][j]
                if (bool(flagz & test1) == True):  # This Boolean test compares the binary values of the interpreted base 10 flags and the "Clean" test
                    results['Clean'][j] = 'False'
                else:
                    results['Clean'][j] = 'True'            
            print(results['flags'][j], results['flagint'][j], results['Clean'][j])
    #
        resultz1 = pd.DataFrame(results)
        #
    #
        if(TARGET_FLUX_R> 9):   
            #
        # First Filter of Results:  Provides Raw Outcomes, before any editing of data to remove Outliers or image misidentifications    
            #
            try:
                target_magnitude, residuals, slope, yintercept = get_regression(results,filter2)
            #
                menu1[i]['Observation Date'] = ObsDate
                menu1[i]['Num Stars NoFltr'] = len(results)
                menu1[i]['Mag_Estimate_NoFltr'] = target_magnitude
                residuals = residuals[0].item(0)
                menu1[i]['Residuals_NoFltr'] = residuals
                stddev = math.sqrt(residuals/(len(results)-1))
                menu1[i]['NoFltr_Std_Dev'] = stddev
            except:
                target_magnitude = 'nan'
                menu1[i]['Observation Date'] = ObsDate
                menu1[i]['Num Stars NoFltr'] = len(results)
                menu1[i]['Mag_Estimate_NoFltr'] = target_magnitude
                residuals = 'nan'
                menu1[i]['Residuals_NoFltr'] = residuals         
                menu1[i]['NoFltr_Std_Dev'] = 'nan'
        #    
            if(ResultsPrint == 'Yes'):
                resultz1 = pd.DataFrame(results)
        #
                clock = datetime.today()
                ID1 = str(clock.year) + '-' + str(clock.month) + '-' + str(clock.day) + '-' + 'h' + str(clock.hour) + 'm' + str(clock.minute) + 's' + str(clock.second)
                name1 = 'Results1_' + STAR_NAME + '_' + ID1 + '-' + filter2 + '.xlsx'
                resultz1 = resultz1[['auid','auidshort','rmag','gmag','imag','zmag','umag','quality','q_mode','objecttype','flags','flagint','Clean','instrumental_mag','residual_aperture_sum','peak','ra','dec','x','y']]
                print(resultz1)
                os.chdir(processedpath)
                name1=name1.replace(" ","_")
                name1=name1.replace("/","-")
                name1=name1.replace("*","star")
                resultz1.to_excel(name1)
        #
        #
            menu2 = menu1 
    #
            target_magnitude1 = 'nan'
            menu2[i]['Num Stars_QMode'] = 'nan'
            menu2[i]['Mag_Estimate_QMode'] = target_magnitude1
            residuals1 = 'nan'
            menu2[i]['Residuals_QMode'] = residuals1     
            menu2[i]['QMode_Std_Dev'] = 'nan'
        #
            target_magnitude2 = 'nan'
            menu2[i]['Num Stars_Clean'] = 'nan'
            menu2[i]['Mag_Estimate_Clean'] = target_magnitude2
            residuals2 = 'nan'
            menu2[i]['Residuals_Clean'] = residuals2     
            menu2[i]['Clean_Std_Dev'] = 'nan'
            #    
            target_magnitude3 = 'nan'
            menu2[i]['Num Stars_Q & Clean'] = 'nan'
            menu2[i]['Mag_Estimate_Q & Clean'] = target_magnitude3
            residuals3 = 'nan'
            menu2[i]['Residuals_Q & Clean'] = residuals3     
            menu2[i]['Q & Clean_Std_Dev'] = 'nan'
            #    
            menu3 = menu2
            resultfinal = []
            resultfinal = results
            resultfinal['dev_from_regr'] = resultfinal['instrumental_mag']
            #
            try:
                if(filter2 == 'r'): 
                    rmagstdev = math.sqrt(residuals/(len(resultfinal)-1))
                    resultfinal['dev_from_regr'] = abs((resultfinal['rmag'] - (resultfinal['instrumental_mag']*slope + yintercept))/rmagstdev)
                    standard = rmagstdev
                elif(filter2 == 'g'):
                    gmagstdev = math.sqrt(residuals/(len(resultfinal)-1))
                    resultfinal['dev_from_regr'] = abs((resultfinal['gmag'] - (resultfinal['instrumental_mag']*slope + yintercept))/gmagstdev)
                    standard = gmagstdev
                elif(filter2 == 'u'):
                    umagstdev = math.sqrt(residuals/(len(resultfinal)-1))
                    resultfinal['dev_from_regr'] = abs((resultfinal['umag'] - (resultfinal['instrumental_mag']*slope + yintercept))/umagstdev)
                    standard = umagstdev
                elif(filter2 == 'i'):
                    imagstdev = math.sqrt(residuals/(len(resultfinal)-1))
                    resultfinal['dev_from_regr'] = abs((resultfinal['imag'] - (resultfinal['instrumental_mag']*slope + yintercept))/imagstdev)
                    standard = imagstdev
                elif(filter2=='z'):
                    zmagstdev = math.sqrt(residuals/(len(resultfinal)-1))
                    resultfinal['dev_from_regr'] = abs((resultfinal['zmag'] - (resultfinal['instrumental_mag']*slope + yintercept))/zmagstdev)
                    standard = zmagstdev
                #
                resultfinal = resultfinal.query('dev_from_regr == "NaN" or dev_from_regr<3.0') 
                resultfinal = resultfinal.reset_index(drop=True)
                #
                target_magnitude4, residuals4, slope4, yintercept4 = get_regression(resultfinal,filter2)
                print("Got Back Even One More Time!")
                menu3[i]['Num Stars_Final'] = len(resultfinal)
                menu3[i]['Mag_Estimate_Final'] = target_magnitude4
                residuals4 = residuals4[0].item(0)
                menu3[i]['Residuals_Final'] = residuals4
                stddev = math.sqrt(residuals4/(len(resultfinal)-1))
                menu3[i]['Final_Std_Dev'] = stddev
        ##
            except:
                target_magnitude2 = 'nan'
                menu3[i]['Num Stars_Final'] = len(resultfinal)
                menu3[i]['Mag_Estimate_Final'] = target_magnitude4
                residuals4 = 'nan'
                menu3[i]['Residuals_Final'] = residuals4     
                menu3[i]['Final_Std_Dev'] = 'nan'      
        # 
            if(ResultsPrint == 'Yes'):
                resultz5 = pd.DataFrame(resultfinal)
            #
                clock = datetime.today()
                ID5 = str(clock.year) + '-' + str(clock.month) + '-' + str(clock.day) + '-' + 'h' + str(clock.hour) + 'm' + str(clock.minute) + 's' + str(clock.second)
                name5 = 'Results5_' + STAR_NAME + '_' + ID5 + '-' + filter2 + '.xlsx'
                resultz5 = resultz5[['auid','auidshort','rmag','gmag','imag','zmag','umag','quality','q_mode','objecttype','flags','flagint','Clean','instrumental_mag','residual_aperture_sum','peak','ra','dec','x','y']]
                print(resultz5)
                os.chdir(processedpath)
                name5=name5.replace(" ","_")
                name5=name5.replace("/","-")
                name5=name5.replace("*","star")
                resultz5.to_excel(name5)
            menu6 = menu3
    #        
        else:
            print('There is insufficient sample size to perform regression')
            menu6 = menu1 
            
    menu7 = menu6
    pd.set_option('display.max_columns', None)
    menu7 = pd.DataFrame(menu6)
    
    clock = datetime.today()
    
    #ID = str(clock.year) + '-' + str(clock.month) + '-' + str(clock.day) + '-' + 'h' + str(clock.hour) + 'm' + str(clock.minute) + 's' + str(clock.second)
    
    name = STAR_NAME + '_' + Date + '.xlsx'
    
    menu7 = menu7[['frame_no','Reduced_Image', 'Observation Date','FITS_target', 'Photometry_Object' , 'filter' , 'Science_Image' ,'Stack Count', 'Exposure' , 'FWHM1', 'SNR', 'ap_radius', 'ann_radius_inner', 'ann_radius_outer',\
        'Num Stars NoFltr', 'Mag_Estimate_NoFltr', 'Residuals_NoFltr','NoFltr_Std_Dev',\
        'Num Stars_QMode', 'Mag_Estimate_QMode' , 'Residuals_QMode', 'QMode_Std_Dev','Num Stars_Clean', 'Mag_Estimate_Clean', 'Residuals_Clean','Clean_Std_Dev',\
        'Num Stars_Q & Clean', 'Mag_Estimate_Q & Clean', 'Residuals_Q & Clean','Q & Clean_Std_Dev','Num Stars_Final', 'Mag_Estimate_Final', 'Residuals_Final','Final_Std_Dev']]
    
    print(menu7)
    
    os.chdir(dirs['Photometry']+targetsub)
    #os.chdir('/home/david/ZowadaTransients/Photometry/teststar')
    name=name.replace(" ","_")
    
    name=name.replace("/","-")
    name=name.replace("*","star")
    menu7.to_excel(name)
    
    print("Photometry program has concluded.  Your measurements have been successfully completed.  Congratulations!")

if __name__ == "__main__":
    main(date, objid, predmag, magrange, fwhmmultiplier, snr, apradius, anninner, annouter,directories)