date = '20210506'
objid = 'ZOW2021aab'

def main(date, objid):
    from astropy.io import fits
    from reproject import reproject_interp
    import os
    from numpy import exp, linspace, nan, zeros, nanmax, nanmin
    from photutils import Background2D, MedianBackground
    import re
    from astropy.stats import sigma_clipped_stats, SigmaClip
    import gc
    
    # for toppath in os.listdir('/home/david/ZowadaTransients/Processed/SN2020uyn'):
    
    gc.enable()
    date = date #date of lights to align to
    targetname =  objid
    
    
    templateDir = "/home/david/ZowadaTransients/Templates/"+targetname #file location of templates
    
    path = os.path.join(templateDir,date)
    try:
        os.mkdir(templateDir)
        print("Top Directory made")
    except FileExistsError:
        print("Top Directory already exists")
    try:
        os.mkdir(path)
        print("Date Directory made")
    except FileExistsError:
        print("Date Directory already exists")
    os.chdir(templateDir)
    
    for filename in os.listdir("/home/david/ZowadaTransients/Processed/" + targetname + "/"+ date + "/nanmedian/"):
        lightFile = "/home/david/ZowadaTransients/Processed/" + targetname + "/"+ date + "/nanmedian/" +filename
        filt = lightFile[-5]
        os.system('cp '+ lightFile +  ' /home/david/ZowadaTransients/Templates/' +targetname + "/"+date +"/lightaligned-"+filt+".fits")
    #Pulls up and manipulates template for alignement, creates a template if one does not exist
        HDUList = fits.open('/home/david/ZowadaTransients/Templates/' +targetname + "/"+date +"/lightaligned-"+filt+".fits") 
        lightprocData = HDUList[0].data
        lightheader = HDUList[0].header
        HDUList.close()
        
        
        sigma_clip = SigmaClip(sigma=3.)
        bkg_estimator = MedianBackground()
        bkg = Background2D(lightprocData, (50, 50), filter_size=(3, 3),\
                      sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    
        lightprocData = lightprocData - bkg.background
        # sizex = len(lightprocData[0])
        # sizey = len(lightprocData[1])
        # nstamps = 6
    
        
        # stampsx = linspace(0,sizex, nstamps)
        # stampsy = linspace(0,sizey,nstamps)
        # stampsx = [int(i) for i in stampsx]
        # stampsy = [int(i) for i in stampsy] 
        
        # for i in range(len(stampsx)-1):
        #     for j in range(len(stampsy)-1):
        #         lightprocData[stampsx[i]:stampsx[i+1], stampsy[j]:stampsy[j+1]] = \
        #             lightprocData[stampsx[i]:stampsx[i+1], stampsy[j]:stampsy[j+1]]\
        #                 - sigma_clipped_stats(lightprocData[stampsx[i]:stampsx[i+1], \
        #                                                     stampsy[j]:stampsy[j+1]], \
        #                                       sigma = 3)[1]
        
        # lightprocData = lightprocData.astype('int16')
        # lightprocData[lightprocData < -1000] = 0
    
        lighthdu = fits.PrimaryHDU(lightprocData)
        lighthdul = fits.HDUList([lighthdu])
        header = lightheader
        fits.writeto('/home/david/ZowadaTransients/Templates/' +targetname + "/"+date +"/lightaligned-"+filt+".fits"\
                     , data = lighthdul[0].data, header = header, overwrite = True)
        
        HDUList = fits.open(lightFile) 
        lightprocData = HDUList[0].data
        lightheader = HDUList[0].header
        HDUList.close()
        
        
        RA = str(lightheader['RA'])
        DEC = str(lightheader['DEC'])
        DEC =re.sub('\s+', ':', DEC)
        RA = re.sub('\s+', ':', RA)
        testlist = []
        for directory in os.listdir(templateDir):
            if directory[0:len(targetname)] == targetname and directory[-6] == filt:
                testlist.append(1)
        
        
        if len(testlist) == 0:
                    
            print('No ',filt,' template found, Creating template, this may take a moment.') 
            os.system('panstamps --width=90 --filters='+filt+ \
                      ' stack ' + RA + ' ' + DEC)
            for file in os.listdir(templateDir):
                if os.path.isdir(file):
                    for item in os.listdir(templateDir + '/'+file):
                        if  item[0:5] == 'stack':
                            itemname = templateDir+'/'+file+'/'+item
                            movedname = templateDir + '/' + targetname+'_unaligned-'\
                                +filt+'.fits'
                            os.system('mv ' +itemname+ ' '+movedname)
                            os.system('rmdir ' + file)
                        
    
    
        Templatename = targetname+"_unaligned-"+filt+".fits"
        HDUtmpl = fits.open(templateDir+"/"+Templatename)[0] 
        HDUsci = fits.open(templateDir+"/" + date + "/lightaligned-"+filt+".fits")[0] 
    
        try:
            array, footprint = reproject_interp(HDUtmpl, HDUsci.header)
            
            fits.writeto(templateDir+"/" + date + "/templatealigned-"+filt+".fits", data = array, header = HDUsci.header, overwrite = True)
            print(filt, "template aligned.")
        except ValueError:
            os.chdir(templateDir)
            os.system('rm -r ' + date)
        
        objname = targetname
        fileloc = '/home/david/ZowadaTransients/Templates/'+objname+'/'+date+'/'+'templatealigned-'+filt+'.fits'
        HDUList = fits.open(fileloc)
        tmpldata = HDUList[0].data
        HDUList.close()
        print('Making '+filt+' Template Mask')
        maskarray = zeros((len(tmpldata),len(tmpldata[1])))
        tmpldata = tmpldata.astype('int32')
        for i in range(len(tmpldata)):
            for j in range(len(tmpldata[i])):
                if tmpldata[i][j] > 175000:
                    maskarray[i][j] = 1
                elif tmpldata[i][j] < -300:
                    maskarray[i][j] = 1
                else:
                    continue
                    
        fits.writeto('/home/david/ZowadaTransients/Templates/' +objname+'/'+date+'/' +"/tmplmask-"+filt+".fits"\
                         , data = maskarray, overwrite = True)
        
        
        fileloc = '/home/david/ZowadaTransients/Templates/'+objname+'/'+date+'/'+'lightaligned-'+filt+'.fits'
        HDUList = fits.open(fileloc)
        scidata = HDUList[0].data
        clipper = SigmaClip(sigma = 2.5)
        datatest = clipper(scidata)
        il = int(nanmin(datatest))
        
        HDUList.close()
        print('Making '+filt+' Science Mask')
        maskarray = zeros((len(scidata),len(scidata[1])))
        scidata = scidata.astype('int32')
        for i in range(len(scidata)):
            for j in range(len(scidata[i])):
                if scidata[i][j] > 60000:
                    maskarray[i][j] = 1
                elif scidata[i][j] < -70:
                    maskarray[i][j] = 1
                else:
                    continue
                    
        fits.writeto('/home/david/ZowadaTransients/Templates/' +objname+'/'+date+'/' +"/scimask-"+filt+".fits"\
                         , data = maskarray, overwrite = True)    
       
    gc.disable()
            
    
if __name__ == "__main__":
    main(date, objid)

            



        
    
