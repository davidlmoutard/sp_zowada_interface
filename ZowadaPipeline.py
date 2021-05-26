import argparse

parser = argparse.ArgumentParser(description='Pipeline for processing Zowada images')

parser.add_argument("-d","--date", help="Date of images to process", required = True)
parser.add_argument("-o","--objid", help="Name of object you are processing", required = True)
parser.add_argument("-m","--mag", help="Predicted magnitude of object", required = True)
parser.add_argument("-r","--magrange", help="Range of magnitudes to fit", required = True)
parser.add_argument("-s","--snr", help="SNR for input into photometry portion", default = 2, type = float)
parser.add_argument("-f","--fwhmmultiplier", help="Multiplier for full width at half max", default = 2, type = float)
parser.add_argument("--apradius", help="Aperture Radius Factor", default = 2, type = float)
parser.add_argument("--anninner", help="Inner ring of annulus", default = 1.5, type = float)
parser.add_argument("--annouter", help="Outer ring of annulus", default = 2, type = float)
parser.add_argument("-p","--photonly", help = "If you only need photometry and do not wish to do processing, set to True", default = False, type = bool)

args = parser.parse_args()
##to run in command line
date = args.date
objid = args.objid
predmag = args.mag
magrange = args.magrange
snr = args.snr
fwhmmultiplier = args.fwhmmultiplier
apradius = args.apradius
anninner = args.anninner
annouter = args.annouter
photonly = args.photonly
## to run in spyder/IDE
# date = 
# objid = 
# predmag =
# magrange = 
# snr = 2
# fwhmmultiplier = 2
# apradius = 2
# anninner = 1.5
# annouter = 2


import Reduction_v2 as reduction
import Alignment_v8 as alignment
import TemplateAlignment_v3 as templatealignment
import Subtraction_v1 as subtraction
import Photometry_v8e as SDSSphotometry
import Photometry_v7apass as APASSphotometry

def main(date, objid, predmag, magrange, fwhmmultiplier, snr, \
                       apradius, anninner, annouter, photonly):
    
    if photonly == False:
        reduction.main(date,objid)
        alignment.main(date, objid)
        templatealignment.main(date, objid)
        subtraction.main(date, objid)
        try:
            SDSSphotometry.main(date, objid, predmag, magrange, fwhmmultiplier, snr, \
                           apradius, anninner, annouter)
        except ValueError:
            APASSphotometry.main(date, objid, predmag, magrange, fwhmmultiplier, snr, \
                           apradius, anninner, annouter)
    else:
        try:
            SDSSphotometry.main(date, objid, predmag, magrange, fwhmmultiplier, snr, \
                           apradius, anninner, annouter)
        except ValueError:
            APASSphotometry.main(date, objid, predmag, magrange, fwhmmultiplier, snr, \
                           apradius, anninner, annouter) 

if __name__ == "__main__":
    main(date, objid, predmag, magrange, fwhmmultiplier, snr, \
                       apradius, anninner, annouter, photonly)
        
    

    
