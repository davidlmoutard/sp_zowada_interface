# sp_zowada_interface

ZowadaPipeline.py is a script that contains all steps to produce photometry, and references the other scripts in the repository. 
It is executable from the command line. This will be amended in future versions.  This pipeline uses the software HOPTANTS for subtraction, as well as popular
pieces of astronomy software like SExtractor and PSFex.It currently requires that you change the path to the location of hotpants in the Subtraction file.  Much of the pipeline is adapted from scripts written by Robert Carr. 

The required inputs are:

-D, --directories: file containing the directories used in photometry. A different directory in this scheme is used to store Calibration, Flats, Templates, etc.
-d, --date: date of images, in the form mmddyyyy. Each night of images is stored in a directory like: "topdir/ObjectName/Date"
-o, --objid: Name of the object
-m, --mag: predicted magnitude of object
-r, --magrange: range of magnitudes to use in Photometry for comparison.

Optional inputs
-f,--fwhmmultiplier: Multiplier for full width at half max
--apradius: Aperture Radius Factor
--anninner: Inner ring of annulus default = 1.5
--annouter: Outer ring of annulus, default = 2
-p,--photonly: If you only need photometry and do not wish to do processing, set to True, default = False

sp_candidate_get.py is a script which allows you to query SkyPortal, and will return a csv file to specified location of any valid transients for use in other 
follow-up algorithms. This script also requires that you change the outputDir variable to your desireed directory for the csv file, and you must also change the 
obsLoc variable to match your observatory's latitude,longitude, and altitude. 
