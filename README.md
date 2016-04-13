# astrometric
Code used to analyze astrometric errors in LSST images.

TO TEST USING IDL:  

idl

.compile astrometric.pro

astrometric

TO TEST USING PYTHON:
python astrometric.py


Both programs will produce a plot of astrometric error as a function of separation.

These programs are designed to take in catalogs from 10 successive images of a star field generated by PhoSim (https://bitbucket.org/phosim/phosim_release/wiki/Home).

To generate the images using PhoSim, run the shell file in /examples, called phosim.sh.  It assumes you are running phosim from the /examples directory.  

Astrometric runs on SExtractor catalogs of the original FITS images.  I include SExtractor catalogs for two different PhoSim atmospheres so you can see the difference.  
I include 10 FITS images for one of those atmospheres.  Note that the images were generated using 10 sequential 15 second snaps.  

If you have questions, e-mail matthewwiesner AT aol.com.


