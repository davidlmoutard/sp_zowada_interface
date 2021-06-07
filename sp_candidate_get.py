import requests
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from datetime import date
from astropy.io import ascii
from astropy.table import Table
import datetime
import os

token = 'your-token-here'

groupData = requests.get(
    'http://desi2.lbl.gov:5000/api/groups',
    headers={'Authorization': f'token {token}'})

today = datetime.date.today()

startDate = today - datetime.timedelta(days=1)
dateStr = str(startDate)
#startDate = "2021-05-20"
numPerPage = 500
groupIDs = [] 
pageNum = 1


groupList = ['DESIDIFF Color', 'DESIDIFF TDE', 'DESIDIFF Hydrogen','DESITRIP']
for groupName in groupList:
    for group in groupData.json()['data']['user_groups']:
        if groupName in group.values():
            groupIDs.append(group['id'])


params =  {"startDate":startDate, "numPerPage":numPerPage, "redshiftRange":"(0,0.4)",\
                        "includePhotometry":True, 'groupIDs':groupIDs,'pageNumber':pageNum}


response = requests.get(
    'http://desi2.lbl.gov:5000/api/candidates', params = params,
    headers={'Authorization': f'token {token}'}
)

print(f'HTTP code: {response.status_code}, {response.reason}')

outputDir = '/home/david/ZowadaTransients/SkyPortal/' #Change to desired Directory
def istransient(candidate):
    if candidate['transient'] == False:
        return False
    else:
        return True
    
def isvisible(candidate, airmassThresh = 2.00):
    obsLoc = EarthLocation(lat=31.836944*u.deg, lon=-109.031667*u.deg, height=1257*u.m) #Zowada
    utcoffset = -6*u.hour
    RA = candidate['ra']
    DEC = candidate['dec']
    obj = SkyCoord(RA,DEC, unit = u.deg)      
    midnight = Time(str(date.today())+" 00:00:00") - utcoffset
    delta_midnight = np.linspace(-3, 6, 100)*u.hour
    objaltaz = obj.transform_to(AltAz(obstime=midnight+delta_midnight, location=obsLoc))         
    minairMass = min([i for i in objaltaz.secz if i > 0])
    if  minairMass < airmassThresh:
        return True
    else:
        return False

def isbright(candidate, magThresh = 20):
    if len(candidate['photometry']) == 0:
        return True #assume if no photometry that it's bright enough
    else:
        flux = candidate['photometry'][-1]['flux'] # Grabs most recent photometry
        filt = candidate['photometry'][-1]['filter']
        mjd = candidate['photometry'][-1]['mjd']
        mag = 22.5 - 2.5*np.log10(flux)
        if mag < magThresh:
            return True
        else:
            return False

        
data = response.json()["data"]["candidates"]
totalCand = len(data)
print("Page Number: " + str(pageNum))
print("Candidates: "+ str(totalCand))

invalid = []
for i in data:
    transient = istransient(i)
    visible = isvisible(i)
    bright = isbright(i)
    if not transient or not visible or not bright:
        invalid.append(i)

for k in invalid:
    data.remove(k)
  
print("Valid Transients: "+ str(len(data)))
RAs = []
DECs = []
mags = []
filts = []
mjds = []
redshifts = []
IDs  = []
for i in data:
    IDs.append(i['id'])
    redshifts.append(i['redshift'])
    RAs.append(i['ra'])
    DECs.append(i['dec'])
    if len(i['photometry']) > 0:
        mags.append(22.5 - 2.5*np.log10(i['photometry'][-1]['flux']))
        mjds.append(i['photometry'][-1]['mjd'])
        filts.append(i['photometry'][-1]['filter'])
    else:
        mags.append(np.nan)
        mjds.append(np.nan)
        filts.append(np.nan)
        

transientOutput = Table([IDs,RAs,DECs,redshifts, mjds, mags, filts], \
                        names = ('ID','RA', 'DEC', 'Redshift','MJD','Mag','Filter'))   

print(transientOutput)
groupStr = ''
for ID in groupIDs:
    groupStr += '_'+str(ID)
ascii.write(transientOutput, outputDir+'SkyPortal_'+dateStr+groupStr+'.csv', format = 'csv', overwrite =True)       

pageNum = 2
while totalCand == 500:
    params =  {"startDate":startDate, "numPerPage":numPerPage, "redshiftRange":"(0,0.4)",\
                            "includePhotometry":True, 'groupIDs':groupIDs,'pageNumber':pageNum}
    
    
    response = requests.get(
        'http://desi2.lbl.gov:5000/api/candidates', params = params,
        headers={'Authorization': f'token {token}'}
    )
    
    print(f'HTTP code: {response.status_code}, {response.reason}')
    
    outputDir = '/home/david/ZowadaTransients/SkyPortal/' #Change to desired Directory
    def istransient(candidate):
        if candidate['transient'] == False:
            return False
        else:
            return True
        
    def isvisible(candidate, airmassThresh = 2.00):
        obsLoc = EarthLocation(lat=31.836944*u.deg, lon=-109.031667*u.deg, height=1257*u.m) #Zowada
        utcoffset = -6*u.hour
        RA = candidate['ra']
        DEC = candidate['dec']
        obj = SkyCoord(RA,DEC, unit = u.deg)      
        midnight = Time(str(date.today())+" 00:00:00") - utcoffset
        delta_midnight = np.linspace(-3, 6, 100)*u.hour
        objaltaz = obj.transform_to(AltAz(obstime=midnight+delta_midnight, location=obsLoc))         
        minairMass = min([i for i in objaltaz.secz if i > 0])
        if  minairMass < airmassThresh:
            return True
        else:
            return False
    
    def isbright(candidate, magThresh = 20):
        if len(candidate['photometry']) == 0:
            return True #assume if no photometry that it's bright enough
        else:
            flux = candidate['photometry'][-1]['flux'] # Grabs most recent photometry
            filt = candidate['photometry'][-1]['filter']
            mjd = candidate['photometry'][-1]['mjd']
            mag = 22.5 - 2.5*np.log10(flux)
            if mag < magThresh:
                return True
            else:
                return False
    
            
    data = response.json()["data"]["candidates"]
    totalCand = len(data)
    print("Page Number: " + str(pageNum))
    print("Candidates: "+ str(totalCand))
    
    invalid = []
    for i in data:
        transient = istransient(i)
        visible = isvisible(i)
        bright = isbright(i)
        if not transient or not visible or not bright:
            invalid.append(i)
    
    for k in invalid:
        data.remove(k)
      
    print("Valid Transients: "+ str(len(data)))
    RAs = []
    DECs = []
    mags = []
    filts = []
    mjds = []
    redshifts = []
    IDs  = []
    for i in data:
        IDs.append(i['id'])
        redshifts.append(i['redshift'])
        RAs.append(i['ra'])
        DECs.append(i['dec'])
        if len(i['photometry']) > 0:
            mags.append(22.5 - 2.5*np.log10(i['photometry'][-1]['flux']))
            mjds.append(i['photometry'][-1]['mjd'])
            filts.append(i['photometry'][-1]['filter'])
        else:
            mags.append(np.nan)
            mjds.append(np.nan)
            filts.append(np.nan)
            
    
    transientOutput = Table([IDs,RAs,DECs,redshifts, mjds, mags, filts], \
                            names = ('ID','RA', 'DEC', 'Redshift','MJD','Mag','Filter'))       

    pageNum += 1
    
    with open(outputDir+'SkyPortal_'+dateStr+groupStr+'.csv', mode='a') as f:
        f.seek(0, os.SEEK_END)  # Some platforms don't automatically seek to end when files opened in append mode
        transientOutput.write(f, format='ascii.no_header', delimiter = ',')
