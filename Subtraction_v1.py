date = '20210506'
objid = 'ZOW2021aab'

def main(date, objid):
    from astropy.io import fits
    import os
    import numpy as np
    from astropy.wcs import WCS
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from astropy.stats import sigma_clipped_stats, SigmaClip
    import matplotlib.pyplot as plt
    import subprocess
    import gc
    from scipy.optimize import minimize
    
    gc.enable()
    templateDir = "/home/david/ZowadaTransients/Templates/"
    sourcelocs =  objid
    date = date
    
    lightdir = templateDir + sourcelocs + "/"+ date
    os.chdir(lightdir)
    
    #Measure PSF for matching
    filtlist = []
    for file in os.listdir(lightdir): #gets the filter from the name of the files.
        if file[0:12] == "lightaligned":
            filtlist.append(file[-6])
    for filt in filtlist:
        os.system('cp /usr/share/sextractor/prepsf.sex prepsf.sex')
        os.system('cp /usr/share/psfex/default.psfex default.psfex')
        os.system('sex lightaligned-'+filt+'.fits -c prepsf.sex ')
        os.system('psfex prepsf.cat')
        
        
        HDUList = fits.open('prepsf.psf') 
        lightheader = HDUList[1].header
        HDUList.close()
        
        lightfwhm = lightheader['PSF_FWHM']
        lightsigma = lightfwhm/2.355
        
        os.system('rm prepsf.psf')
        os.system('rm prepsf.cat')
        os.system('sex templatealigned-'+filt+'.fits -c prepsf.sex')
        os.system('psfex prepsf.cat')
        
        
        
        HDUList = fits.open('prepsf.psf') 
        tmplheader = HDUList[1].header
        HDUList.close()
        
        tmplfwhm = tmplheader['PSF_FWHM']
        tmplsigma = tmplfwhm/2.355
        
        os.system('rm prepsf.psf')
        os.system('rm prepsf.cat')
        
        HDUList = fits.open('lightaligned-'+filt+'.fits') 
        sciProcdata = HDUList[0].data
        HDUList.close()
        # ilactual = str(int(np.nanmin(sciProcdata)))
        iu = str(int(np.nanmax(sciProcdata))) 
        clipper = SigmaClip(sigma = 3)
        datatest = clipper(sciProcdata)
        il = str(np.nanmin(datatest))
        # plt.hist(datatest)
        # plt.show()
    
    
     
        
        #From PSFex, Gaussians are measured. Used to generate kernel and blur the template 
        diffsig = np.sqrt(lightsigma**2 - tmplsigma**2)
        
        HDUList = fits.open('lightaligned-'+filt+'.fits') 
        sciProcdata = HDUList[0].data
        HDUList.close()
        # ilactual = str(int(np.nanmin(sciProcdata)))
        iu = int(np.nanmax(sciProcdata))
        if iu > 60000:
            iu = 60000
        clipper = SigmaClip(sigma = 3)
        datatest = clipper(sciProcdata)
        il = str(int(np.nanmin(datatest)))
        iu = str(iu)
        os.system('/home/david/hotpants/hotpants -inim lightaligned-'+filt+'.fits \
                  -tmplim templatealigned-'+filt+'.fits -outim sub-'+filt+'.fits  \
                        -tu 150000  -tl -300 -iu '+iu+' -il '+il+' -n i -ng 3 6 ' +str(0.5*diffsig)+\
                          ' 4 '+ str(diffsig)+ ' 2 ' +str(2*diffsig) + \
                              ' -tmi tmplmask-'+filt+'.fits -imi scimask-'+filt+'.fits -v 0')
        
        
        sub = "sub-"+filt+".fits"
        HDUList = fits.open(sub) 
        subheader = HDUList[0].header
        HDUList.close()    
        if filt == 'r':
            X2r = float(subheader['X2NRM00'])
            ilr = il
        elif filt == 'g':
            X2g = float(subheader['X2NRM00'])
            ilg = il
    
    
        
        
            
        sub = "sub-"+filt+".fits"
        HDUList = fits.open(sub) 
        subprocData = HDUList[0].data
        subheader = HDUList[0].header
        HDUList.close()    
    
            
        light = "lightaligned-"+filt+".fits"
        HDUList = fits.open(light) 
        lightprocData = HDUList[0].data
        lightheader = HDUList[0].header
        HDUList.close()
        
    
        
        lightprocData = lightprocData.astype('int16')
        w = WCS(lightheader)
        sc = SkyCoord(str(str(lightheader['RA']) + str(lightheader['DEC'])) , unit=(u.hourangle, u.deg))
        objy, objx = w.world_to_pixel(sc)
        
        if np.mean(subprocData[int(objx-80):int(objx+80), int(objy-80):int(objy+80)]) != (np.nan or 0):
          
            lightprocData[int(objx-80):int(objx+80), int(objy-80):int(objy+80)] = \
                subprocData[int(objx-80):int(objx+80), int(objy-80):int(objy+80)] 
            
            
            lightprocData[lightprocData < -1000] = 0 
            
            lighthdu = fits.PrimaryHDU(lightprocData, uint = False)
            lighthdul = fits.HDUList([lighthdu])
            header = lightheader
            header.set('SUB','T')
            fits.writeto(templateDir + sourcelocs+"/"+ date+"/subcutout-" +filt+ ".fits", data = lighthdul[0].data, header = header, overwrite = True)
            lighthdul.close()
            os.system('rm *_prepsf.fits')
            
            HDUList = fits.open(templateDir + sourcelocs+ "/"+date+"/subcutout-" +filt+ ".fits") 
            HDUList.close()
        else:
            print("Subtraction failed- no object found!")
    #os.system('rm lightaligned* templatealigned* sub-* tmplmask-* scipsf* tmplpsf* scimask* pre* psfex* default*')
    gc.disable()  
    
if __name__ == "__main__":
    main(date, objid)