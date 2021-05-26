# sp_zowada_interface

ZowadaPipeline.py is a script that contains all steps to produce photometry, and references the other scripts in the repository. 
It is executable from the command line. It currently requires that you change the paths in each step of the script to the paths where you keep your calibration 
files, your flats, your raw images, etc. This will be amended in future versions.  This pipeline uses the software HOPTANTS for subtraction, as well as popular
pieces of astronomy software like SExtractor and PSFex. Much of the pipeline is adapted from scripts written by Robert Carr.


sp_candidate_get.py is a script which allows you to query SkyPortal, and will return a csv file to specified location of any valid transients for use in other 
follow-up algorithms. This script also requires that you change the outputDir variable to your desireed directory for the csv file, and you must also change the 
obsLoc variable to match your observatory's latitude,longitude, and altitude. 
