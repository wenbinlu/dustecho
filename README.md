# dustecho
The goal of this project is to compute the dust echo from a jetted transient source. The resulting lightcurves, as well as astrometric shift of the flux centroid, for a number of different dust column densities and observer's viewing angles are in files named like 'nH3_lam4.0um_theobs0.52_Ldnu_xcentr.txt', where the first number means hydrogen number density n_H = 3/cm^3, the second number is for observing wavelength of lambda = 4.0 micro-meter, and the third number is for observer's viewing angle of 0.52 radian with respect to the jet axis. A few of the files have 'on_axis' in their names, and they are for the case where the observer's line of sight coincident with the jet axis.

Inside each file, there are four rows (separated by '\n'). The first row shows the info for the observer's time grid, which is taken to be logarithmic grid between 'tobsmin' and 'tobsmax' with 'Ntobs' number of grid points and the right boundary is not included. For Python users, the time grid is given by

tobs = numpy.logspace(tobsmin, tobsmax, Ntobs)

The second line shows the meaning of the later two rows of numbers: Ldnu is the specific luminosity [in erg/s/Hz] at wavelength lambda specificed in the filename; xcentr is the projected distance [in pc] between the flux centroid and the center of explosion on the sky. The third row shows Ldnu and the last row shows xcentr. When Ldnu = 0 (dust echo hasn't arrived yet), we have xcentr = numpy.nan (does not exist).

A final note is that the dust echo from the counter jet (which propagates away from the observer) is not included. The authors plan to publish the computing code (written in Python 3.7) after the paper is accepted for publication.
